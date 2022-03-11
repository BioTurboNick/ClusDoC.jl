module ClusDoC

using Clustering
using Contour
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

const defaultdocparameters = DoCParameters(20, 500, 10)
const defaultclusterparameters = ClusterParameters(20, 3, true, 15, 10)

"""
    clusdoc()

Open the ClusDoC GUI

    clusdoc(channelnames, localizations, roiarea)

Run ClusDoC on a single set of localizations.

    clusdoc(inputfiles, rois, localizations, outputfolder, update_callback)
"""
clusdoc() = (include("src/gui.jl"); nothing)

# profiling - 2300 frame in nearest neighbors `inrange`, 3000 frames in `imfilter`, 200 spent in `rankcorr` (out of 6000)
# should check out rankcorr, it's the one I haven't tried to optimize at all - checked and it's about minimal
function clusdoc(channelnames, localizations, roiarea, docparameters::DoCParameters, clusterparameters::Vector{ClusterParameters})
    cr = doc(channelnames, localizations, docparameters.localradius, docparameters.radiusmax, docparameters.radiusstep, roiarea)
    dbscan!(cr, clusterparameters, docparameters.localradius)
    smooth!(cr, clusterparameters)
    calculate_colocalized_cluster_data!(cr, clusterparameters)
    return cr
end

function clusdoc(inputfiles, rois, localizations, outputfolder, colors = defaultcolors, docparameters = defaultdocparameters, clusterparameters = nothing, update_callback = () -> nothing)
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

        results = Vector{ClusDoC.ChannelResult}[]
        
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
            cr = clusdoc(chnames, roilocalizations, abs(PolygonOps.area(roi) * scalefactor ^ 2), docparameters, fileclusterparameters)
            push!(results, cr)
            generate_roi_output(cr, outputfolder, filename, i, chnames, colors)
            roi_endtime = time_ns()
            roi_elapsed_s = (roi_endtime - roi_starttime) / 1_000_000_000
            println("         finished in $roi_elapsed_s s")
            update_callback()
            calculate_pooled_cluster_statistics(cr)
        end

        writeresultstables(results, docparameters, fileclusterparameters, joinpath(outputfolder, "$(filename) ClusDoC Results.xlsx"))
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

const colocalized_threshold = 0.4
const minclusterpoints = 10 # clusters with fewer points than this are ignored, and clusers must have at least this many colocalized points to be considered colocalized

function calculate_colocalized_cluster_data!(cr::Vector{ChannelResult}, clusterparameters)
    for (i, channel) ∈ enumerate(cr)
        all_colocalized_indexes = []
        clusterdata_abovecutoff = filter(x -> x.cluster.size > clusterparameters[i].pointscutoff, channel.clusterdata)
        for (j, channel2) ∈ enumerate(cr)
            i == j && continue
            clusterdata2_abovecutoff = filter(x -> x.cluster.size > clusterparameters[j].pointscutoff, channel2.clusterdata)
            colocalized_indexes = findall([c.ninteracting[j] ≥ minclusterpoints for c ∈ eachrow(clusterdata2_abovecutoff)])
            union!(all_colocalized_indexes, colocalized_indexes)
            channel.ncoclusters[j] = length(colocalized_indexes)
            length(colocalized_indexes) > 0 || continue

            channel.meancoclustersize[j] = mean(clusterdata2_abovecutoff.size[colocalized_indexes])
            channel.meancoclusterarea[j] = mean(clusterdata2_abovecutoff.area[colocalized_indexes])
            channel.meancoclustercircularity[j] = mean(clusterdata2_abovecutoff.circularity[colocalized_indexes])
            channel.meancoclusterdensity[j] = mean([mean(channel2.pointdata.density[c.core_indices]) / (channel2.roidensity / 1_000_000) for (i, c) ∈ enumerate(clusterdata2_abovecutoff.cluster[colocalized_indexes])])

            incluster = falses(channel2.nlocalizations)
            for cluster ∈ eachrow(clusterdata_abovecutoff)
                incluster .|= inpolygon.(eachcol(channel2.coordinates), Ref(cluster.contour)) .!= 0
            end
            docscoresj = channel.pointdata[!, Symbol(:docscore, j)]
            channel.fraction_colocalized[j] = count(docscoresj .> colocalized_threshold) / count(channel.pointdata.abovethreshold)

            docscoresi = channel2.pointdata[!, Symbol(:docscore, i)]
            channel.fraction_interactions_clustered[j] = count(>(colocalized_threshold), docscoresi[incluster]) / count(>(colocalized_threshold), docscoresi)
        end

        noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ minclusterpoints for c ∈ clusterdata_abovecutoff.cluster]), all_colocalized_indexes)
        channel.ncoclusters[i] = length(noncolocalized_indexes)
        length(all_colocalized_indexes) < length(clusterdata_abovecutoff.cluster) || continue
        channel.meancoclustersize[i] = mean(clusterdata_abovecutoff.size[noncolocalized_indexes])
        channel.meancoclusterarea[i] = mean(clusterdata_abovecutoff.area[noncolocalized_indexes])
        channel.meancoclustercircularity[i] = mean(clusterdata_abovecutoff.circularity[noncolocalized_indexes])
        channel.meancoclusterdensity[i] = mean([mean(channel.pointdata.density[c.core_indices]) / (channel.roidensity / 1_000_000) for (i, c) ∈ enumerate(clusterdata_abovecutoff.cluster[noncolocalized_indexes])])
    end
    return cr

    ### In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
    ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
    ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
end

function calculate_pooled_cluster_statistics(cr::Vector{ChannelResult})
    # Calculate statistics on clusters under the assumption that clusters can be pooled between paired channels

    for (i, channel) ∈ enumerate(cr)
        for (j, channel2) ∈ enumerate(cr)
            no_interaction = filter(x -> x.ninteracting[j] == 0, channel.clusterdata)
            low_interaction = filter(x -> 0 < x.ninteracting[j] ≤ 5, channel.clusterdata)
            high_interaction = filter(x -> x.ninteracting[j] > 5, channel.clusterdata)

            for clustertype ∈ (no_interaction, low_interaction, high_interaction, channel.clusterdata)
                allclusters = []
                for cluster ∈ eachrow(clustertype)
                    incluster = inpolygon.(eachcol(channel2.coordinates), Ref(cluster.contour)) .!= 0
                    inclustercount = count(incluster)
                    points = (1:size(channel2.coordinates, 2))[incluster]

                    push!(allclusters, (; cluster, points, inclustercount))
                end
                if length(allclusters) > 0
                    println("$(length(allclusters)), $(mean(x -> x.cluster.area, allclusters)), $(mean(x -> x.cluster.size, allclusters)), \
                        $(mean(x -> x.inclustercount, allclusters))")
                end
            end
            
        end
    end
end


load_raw_results(path) = load(path)["results"]

end # module

#  TODO: New parameter in UI, Check cocluster results