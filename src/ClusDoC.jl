module ClusDoC

using Clustering
using Contour
using CSV
using DataFrames
using ImageFiltering
using Interpolations
using InvertedIndices
using LocalizationMicroscopy
using NearestNeighbors
using PolygonOps
using StatsBase
using XLSX
using Gtk.ShortNames, GtkObservables, NativeFileDialog, Plots, LocalizationMicroscopy, ImageCore, ImageIO, ColorTypes # find out which subpackages of Images I need
using JLD2
using Statistics

export clusdoc, load_raw_results

include("types/ChannelResult.jl")
include("types/ClusDoCParameters.jl")
include("doc.jl")
include("dbscan.jl")
include("smooth.jl")
include("output.jl")
const sourcedir = @__DIR__

const defaultdocparameters = DoCParameters(20, 500, 10, 0.4)
const defaultclusterparameters = ClusterParameters(20, 3, true, 15, 10)
const defaultoutputcsvclusters = false;
const defaultclustercombinechannels = false;
const gladepath = joinpath(sourcedir, "../gui/ClusDoC.glade")


"""
    clusdoc()

Initialize GUI.

    clusdoc(channelnames, localizations, roiarea)

Run ClusDoC on a single set of localizations.

    clusdoc(inputfiles, rois, localizations, outputfolder, update_callback)
"""
clusdoc() = (include(joinpath(sourcedir, "gui.jl")); nothing)

# profiling - 2300 frame in nearest neighbors `inrange`, 3000 frames in `imfilter`, 200 spent in `rankcorr` (out of 6000)
# should check out rankcorr, it's the one I haven't tried to optimize at all - checked and it's about minimal
function clusdoc(channelnames, localizations, roiarea, docparameters::DoCParameters, clusterparameters::Vector{ClusterParameters}, combine_channels_for_clustering = false)
    result = doc(channelnames, localizations, docparameters.localradius, docparameters.radiusmax, docparameters.radiusstep, roiarea)
    dbscan!(result, clusterparameters, combine_channels_for_clustering)
    smooth!(result, clusterparameters, combine_channels_for_clustering)
    calculate_colocalization_data!(result, docparameters, clusterparameters, combine_channels_for_clustering)
    return result
end

function clusdoc(inputfiles, rois, localizations, outputfolder, colors = defaultcolors, docparameters = defaultdocparameters, clusterparameters = nothing, output_clusters_to_csv = false, combine_channels_for_clustering = false, update_callback = () -> nothing)
    isempty(rois) && return
    println("Starting ClusDoC")

    starttime = time_ns()

    for inputfile ∈ inputfiles
        println("Working on $inputfile")
        filename = basename(inputfile)
        locs = localizations[filename]
        chnames = sort(unique(keys(locs)))

        if clusterparameters === nothing
            fileclusterparameters = fill(defaultclusterparameters, length(chnames))
        else
            fileclusterparameters = clusterparameters
        end

        results = ROIResult[]
        
        if haskey(rois, filename) && !isempty(rois[filename])
            filerois = copy(rois[filename])
        else
            filerois = []
        end

        scalefactor = get_scale_factor(locs)

        for (i, roi) ∈ enumerate(filerois)
            println("    ROI $i")
            roi_starttime = time_ns()
            roi = [(x, 1 - y) for (x, y) ∈ roi] # invert y to match localization coordinates - but actually I might need to invert the original image instead
            roilocalizations = get_roi_localizations(locs, chnames, roi)
            result = clusdoc(chnames, roilocalizations, abs(PolygonOps.area(roi) * scalefactor ^ 2), docparameters, fileclusterparameters, combine_channels_for_clustering)
            push!(results, result)
            generate_roi_output(result, outputfolder, filename, i, chnames, colors)
            roi_endtime = time_ns()
            roi_elapsed_s = (roi_endtime - roi_starttime) / 1_000_000_000
            println("         finished in $roi_elapsed_s s")
            if output_clusters_to_csv
                memberscolumns = [Symbol("$(cn)members") for cn ∈ chnames]
                clusterdatacopy = copy(result.clusterdata)[!, Not([:cluster, :contour, memberscolumns...])]
                CSV.write(joinpath(outputfolder, "$filename cluster data.csv"), clusterdatacopy)
            end
            update_callback()
        end

        resulttablepath = joinpath(outputfolder, "$(filename) ClusDoC Results.xlsx")
        # try
            writeresultstables(results, docparameters, fileclusterparameters, resulttablepath, combine_channels_for_clustering)
        # catch
        #     println("Failed to save results tables. Please ensure the path ($resulttablepath) is valid and the file is not open in another program.")
        # end
        save(joinpath(outputfolder, "$(filename) raw data.jld2"), "results", results)
        # should save: localization map with ROIs shown
        # should have rois automatically save and load (if possible)
        # consider padding that rois can be drawn in
    end

    endtime = time_ns()
    elapsed_s = (endtime - starttime) / 1_000_000_000

    println("Done in $elapsed_s s")
end

function get_bounds(locs)
    xmin, xmax, ymin, ymax = (Inf, -Inf, Inf, -Inf)
    for (i, chname) ∈ enumerate(sort(collect(keys(locs))))
        chpoints = extractcoordinates(locs[chname])
        xmin = min(minimum(chpoints[1, :]), xmin)
        xmax = max(maximum(chpoints[1, :]), xmax)
        ymin = min(minimum(chpoints[2, :]), ymin)
        ymax = max(maximum(chpoints[2, :]), ymax)
    end
    return (xmin, xmax), (ymin, ymax)
end

function get_scale_factor(locs)
    # scale factor to map the relative line positions and the image coordinates
    (xmin, xmax), (ymin, ymax) = get_bounds(locs)
    return max(xmax - xmin, ymax - ymin)
end

function get_roi_localizations(locs, chnames, roi)
    roilocalizations = Vector{Vector{Localization}}()
    scalefactor = get_scale_factor(locs)

    for chname ∈ chnames
        coords = extractcoordinates(locs[chname])
        roichlocalizationsmask = inpolygon.(eachcol(coords ./ scalefactor), Ref(roi)) .!= 0
        push!(roilocalizations, locs[chname][roichlocalizationsmask])
    end
    return roilocalizations
end

function generate_roi_output(cr, outputfolder, filename, i, chnames, colors)
    generate_localization_maps(cr, outputfolder, filename, i, chnames, colors)
    generate_doc_maps(cr, outputfolder, filename, i, chnames)
    generate_cluster_maps(cr, outputfolder, filename, i, chnames, colors)
    generate_doc_histograms(cr, outputfolder, filename, i, chnames, colors)
end

function calculate_colocalization_data!(result::ROIResult, docparameters, clusterparameters, combine_channels_for_clustering)
    nchannels = result.nchannels
    for i ∈ 1:nchannels
        for j ∈ 1:nchannels
            docscores12 = filter(x -> x.channel == i, result.pointdata)[!, Symbol(:docscore, j)]
            result.pointschannelresults[i].fraction_colocalized[j] = count(docscores12 .> docparameters.colocalized_threshold) / result.pointschannelresults[i].nlocalizations_abovethreshold
        end
    end

    if combine_channels_for_clustering
        summarize_interaction_data!(result.clusterdata, result.pointdata, result.clusterresults[1], result.sigclusterresults[1], result.pointschannelresults, docparameters, result.channelnames, 0)
        summarize_cocluster_data!(result.clusterdata, result.pointdata, result, result.pointschannelresults, docparameters, clusterparameters[1], result.channelnames)
    else
        expandedclusterdata = DataFrame()
        for i ∈ 1:result.nchannels
            channelclusterdata = filter(x -> x.channel == i, result.clusterdata)
            summarize_interaction_data!(channelclusterdata, result.pointdata, result.clusterresults[i], result.sigclusterresults[i], result.pointschannelresults, docparameters, result.channelnames, i)
            summarize_cocluster_data!(channelclusterdata, result.pointdata, result, result.pointschannelresults, docparameters, clusterparameters[i], result.channelnames)
            append!(expandedclusterdata, channelclusterdata)
        end
        result.clusterdata = expandedclusterdata
    end
    return nothing

    ### In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
    ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
    ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
end

function summarize_cluster_data!(clusterresult::ClustersResult, clusterdata::DataFrame, pointdata::DataFrame, pointschannelresults::Vector{PointsChannelResult}, docparameters::DoCParameters, channelnames::Vector{String}, j::Int)
    clusterresult.meanclusterarea = mean(clusterdata.area)
    clusterresult.meanclustercircularity = mean(clusterdata.circularity)

    nchannels = length(channelnames)

    # statistics for each channel within the coclusters
    for k ∈ 1:nchannels
        countcolumn = Symbol("$(channelnames[k])count")
        absolutedensitycolumn = Symbol("$(channelnames[k])absolutedensity")
        relativedensitycolumn = Symbol("$(channelnames[k])density")

        meanclustersize = mean(clusterdata[!, countcolumn])
        meanclusterabsolutedensity = mean(clusterdata[!, absolutedensitycolumn])
        meanclusterdensity = mean(clusterdata[!, relativedensitycolumn])
        
        if j > 0
            memberscolumn = Symbol("$(channelnames[j])members")
            allclustered = union(clusterdata[!, memberscolumn]...)
            allclusteredcount = length(allclustered)
            fraction_clustered = allclusteredcount / pointschannelresults[j].nlocalizations_abovethreshold
            docscores12 = filter(x -> x.channel == j, pointdata)[!, Symbol(:docscore, k)]
            fraction_interactions_clustered12 = count(>(docparameters.colocalized_threshold), docscores12[allclustered]) / pointschannelresults[j].nlocalizations_abovethreshold
            ch2memberscolumn = Symbol("$(channelnames[k])members")
            ch2allclustered = union(clusterdata[!, ch2memberscolumn]...)
            docscores21 = filter(x -> x.channel == k, pointdata)[!, Symbol(:docscore, j)]
            fraction_interactions_clustered21 = count(>(docparameters.colocalized_threshold), docscores21[ch2allclustered]) / pointschannelresults[j].nlocalizations_abovethreshold
        else
            fraction_clustered = NaN
            fraction_interactions_clustered12 = NaN
            fraction_interactions_clustered21 = NaN
        end

        clusterchannelresult = ClustersChannelResult(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, fraction_interactions_clustered12, fraction_interactions_clustered21)
        push!(clusterresult.channelresults, clusterchannelresult)
    end
end

function summarize_interaction_data!(clusterdata::DataFrame, pointdata::DataFrame, clustersresult::ClustersResult, sigclustersresult::ClustersResult, pointschannelresults::Vector{PointsChannelResult}, docparameters::DoCParameters, channelnames::Vector{String}, i::Int)
    nclusters = clustersresult.nclusters
    nchannels = length(channelnames)
    for j ∈ 1:nchannels
        # count members of each cluster
        ch2pointdata = filter(x -> x.channel == j, pointdata)
        ch2countcolumn = Symbol("$(channelnames[j])count")
        ch2interactingcountcolumns = [Symbol("$(channelnames[j])-$(cn)interactingcount") for cn ∈ channelnames]
        ch2memberscolumn = Symbol("$(channelnames[j])members")
        clusterdata[!, ch2countcolumn] = Vector{Int}(undef, nclusters)
        for k ∈ 1:nchannels
            k == j && continue
            clusterdata[!, ch2interactingcountcolumns[k]] = Vector{Int}(undef, nclusters)
        end
        clusterdata[!, ch2memberscolumn] = Vector{Vector{Int}}(undef, nclusters)
        for (ii, cluster) ∈ enumerate(eachrow(clusterdata))
            inmembers = inpolygon.(eachcol(pointschannelresults[j].coordinates), Ref(cluster.contour)) .!= 0
            members = findall(inmembers)
            clusterdata[ii, ch2memberscolumn] = members
            ninteracting = [count(ch2pointdata[!, Symbol(:docscore, k)][members] .> docparameters.colocalized_threshold) for k ∈ eachindex(channelnames)]
            for k ∈ 1:nchannels
                k == j && continue
                clusterdata[ii, ch2interactingcountcolumns[k]] = ninteracting[k]
            end
            clusterdata[ii, ch2countcolumn] = length(members)
        end

        ch2absolutedensitycolumn = Symbol("$(channelnames[j])absolutedensity")
        clusterdata[!, ch2absolutedensitycolumn] = clusterdata[!, ch2countcolumn] ./ (clusterdata.area ./ 1_000_000)

        ch2relativedensitycolumn = Symbol("$(channelnames[j])density")
        clusterdata[!, ch2relativedensitycolumn] = [calculate_relative_density(pointdata.density[clusterdata[i, ch2memberscolumn]], pointschannelresults[j].roidensity) for i ∈ eachindex(eachrow(clusterdata))]
    end

    for j ∈ 1:nchannels
        summarize_cluster_data!(clustersresult, clusterdata, pointdata, pointschannelresults, docparameters, channelnames, j)
        sigclusterdata = filter(x -> x.issignificant, clusterdata)
        summarize_cluster_data!(sigclustersresult, sigclusterdata, pointdata, pointschannelresults, docparameters, channelnames, j)
    end
end

function calculate_relative_density(pointdensities, roidensity)
    x = mean(pointdensities) / roidensity
    return isnan(x) ? 1.0 : x # If there are no points, we can say the points are evenly distributed.
end

function summarize_cocluster_data!(clusterdata::DataFrame, pointdata::DataFrame, result::ROIResult, pointschannelresults::Vector{PointsChannelResult}, docparameters::DoCParameters, clusterparameters::ClusterParameters, channelnames::Vector{String})
    nchannels = length(channelnames)
    for i ∈ 1:nchannels
        all_colocalized_indexes = []
        clusterdata_abovecutoff = filter(x -> x.cluster.size > clusterparameters.minsigclusterpoints, clusterdata)
    
        interactingcountcolumns = [Symbol("$(channelnames[i])-$(channelnames[j])interactingcount") for j ∈ 1:nchannels]

        nchannels = length(channelnames)

        # coclusters
        for j ∈ 1:nchannels
            i == j && continue

            colocalized_indexes = findall([c[interactingcountcolumns[j]] ≥ clusterparameters.minsigclusterpoints for c ∈ eachrow(clusterdata_abovecutoff)])
            union!(all_colocalized_indexes, colocalized_indexes)
            ncoclusters = length(colocalized_indexes)
            clusterresult = ClustersResult(ncoclusters, ncoclusters / result.roiarea)
            push!(result.coclusterresults, clusterresult)
            length(colocalized_indexes) > 0 || continue

            coclusterdata = clusterdata_abovecutoff[colocalized_indexes, :]
            summarize_cluster_data!(clusterresult, coclusterdata, pointdata, pointschannelresults, docparameters, channelnames, j)
        end

        # intermediate coclusters
        for j ∈ 1:nchannels
            i == j && continue

            colocalized_indexes = findall([0 < c[interactingcountcolumns[j]] < clusterparameters.minsigclusterpoints for c ∈ eachrow(clusterdata_abovecutoff)])
            union!(all_colocalized_indexes, colocalized_indexes)
            nintcoclusters = length(colocalized_indexes)
            clusterresult = ClustersResult(nintcoclusters, nintcoclusters / result.roiarea)
            push!(result.intermediatecoclusterresults, clusterresult)
            length(colocalized_indexes) > 0 || continue

            intcoclusterdata = clusterdata_abovecutoff[colocalized_indexes, :]
            summarize_cluster_data!(clusterresult, intcoclusterdata, pointdata, pointschannelresults, docparameters, channelnames, j)
        end

        # noncolocalized
        noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ clusterparameters.minsigclusterpoints for c ∈ clusterdata_abovecutoff.cluster]), all_colocalized_indexes)
        nnoninteractingclusters = length(noncolocalized_indexes)
        clusterresult = ClustersResult(nnoninteractingclusters, nnoninteractingclusters / result.roiarea)
        result.noncolocalizedclusterresults = clusterresult
        length(all_colocalized_indexes) < length(clusterdata_abovecutoff.cluster) || continue

        nonclusterdata = clusterdata_abovecutoff[noncolocalized_indexes, :]
        summarize_cluster_data!(clusterresult, nonclusterdata, pointdata, pointschannelresults, docparameters, channelnames, 0)
    end
end


load_raw_results(path) = load(path)["results"]

end # module
