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
const defaultclusterparameters = ClusterParameters(20, 3, true, 15)

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
    calculate_colocalized_cluster_data!(cr)
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
            filerois = [create_whole_image_roi(locs, chnames)]
        end

        for (i, roi) ∈ enumerate(filerois)
            println("    ROI $i")
            roi_starttime = time_ns()
            roi = [(x, 1 - y) for (x, y) ∈ roi] # invert y to match localization coordinates - but actually I might need to invert the original image instead
            roilocalizations = get_roi_localizations(locs, chnames, roi)
            cr = clusdoc(chnames, roilocalizations, abs(PolygonOps.area(roi) * 40960 * 40960), docparameters, fileclusterparameters)
            push!(results, cr)
            generate_roi_output(cr, outputfolder, filename, i, chnames, colors)
            roi_endtime = time_ns()
            roi_elapsed_s = (roi_endtime - roi_starttime) / 1_000_000_000
            println("         finished in $roi_elapsed_s s")
            update_callback()
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

function create_whole_image_roi(locs, chnames)
    coords = reduce(hcat, [extractcoordinates(locs[chname]) for chname ∈ chnames])
    xmin, xmax = extrema(coords[1, :]) ./ 40960
    ymin, ymax = extrema(coords[2, :]) ./ 40960
    return [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]
end

function get_roi_localizations(locs, chnames, roi)
    roilocalizations = Vector{Vector{Localization}}()
    for chname ∈ chnames
        coords = extractcoordinates(locs[chname])
        roichlocalizationsmask = inpolygon.(eachcol(coords ./ 40960), Ref(roi)) .!= 0
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

function calculate_colocalized_cluster_data!(cr::Vector{ChannelResult})
    for (i, channel) ∈ enumerate(cr)
        channel.fraction_colocalized = [count(channel.pointdata[!, Symbol(:docscore, j)] .> colocalized_threshold) / length(channel.pointdata[!, Symbol(:docscore, j)]) for j ∈ eachindex(cr)]

        all_colocalized_indexes = []
        for (j, channel2) ∈ enumerate(cr)
            i == j && continue
            colocalized_indexes = findall([count(channel2.pointdata[!, Symbol(:docscore, j)][c.core_indices] .> 0.4) ≥ minclusterpoints for c ∈ channel2.clusterdata.cluster])
            union!(all_colocalized_indexes, colocalized_indexes)
            channel.ncoclusters[j] = length(colocalized_indexes)
            length(colocalized_indexes) > 0 || continue

            channel.meancoclustersize[j] = mean(channel2.clusterdata.size[colocalized_indexes])
            channel.meancoclusterarea[j] = mean(channel2.clusterdata.area[colocalized_indexes])
            channel.meancoclustercircularity[j] = mean(channel2.clusterdata.circularity[colocalized_indexes])
            channel.meancoclusterdensity[j] = mean([mean(channel2.pointdata.density[c.core_indices]) / (channel2.roidensity / 1_000_000) for (i, c) ∈ enumerate(channel2.clusterdata.cluster[colocalized_indexes])])
        end

        noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ minclusterpoints for c ∈ channel.clusterdata.cluster]), all_colocalized_indexes)
        channel.ncoclusters[i] = length(noncolocalized_indexes)
        length(all_colocalized_indexes) < length(channel.clusterdata.cluster) || continue
        channel.meancoclustersize[i] = mean(channel.clusterdata.size[noncolocalized_indexes])
        channel.meancoclusterarea[i] = mean(channel.clusterdata.area[noncolocalized_indexes])
        channel.meancoclustercircularity[i] = mean(channel.clusterdata.circularity[noncolocalized_indexes])
        channel.meancoclusterdensity[i] = mean([mean(channel.pointdata.density[c.core_indices]) / (channel.roidensity / 1_000_000) for (i, c) ∈ enumerate(channel.clusterdata.cluster[noncolocalized_indexes])])
    end
    return cr

    ###  In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
    ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
    ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
end

load_raw_results(path) = load(path)["results"]

end # module
