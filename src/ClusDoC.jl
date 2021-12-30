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
using Gtk.ShortNames, GtkObservables, NativeFileDialog, Plots, LocalizationMicroscopy, ImageCore, ImageIO # find out which subpackages of Images I need
using JLD2
using Statistics

export clusdoc, load_raw_results

include("types/ChannelResult.jl")
include("doc.jl")
include("dbscan.jl")
include("smooth.jl")
include("output.jl")

clusdoc() = (include("src/gui.jl"); nothing)

function clusdoc(channelnames, localizations, roiarea)
    cr = doc(channelnames, localizations, 20, 500, 10, roiarea)
    dbscan!(cr, 20, 3, true, 20)
    smooth!(cr, 20, 15)
    calculate_colocalized_cluster_data!(cr)
    return cr
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
            channel.meancoclusterdensity[j] = mean([mean(channel2.pointdata.density[c.core_indices]) / channel2.roidensity for (i, c) ∈ enumerate(channel2.clusterdata.cluster[colocalized_indexes])])
        end

        noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ minclusterpoints for c ∈ channel.clusterdata.cluster]), all_colocalized_indexes)
        channel.ncoclusters[i] = length(noncolocalized_indexes)
        length(all_colocalized_indexes) < length(channel.clusterdata.cluster) || continue
        channel.meancoclustersize[i] = mean(channel.clusterdata.size[noncolocalized_indexes])
        channel.meancoclusterarea[i] = mean(channel.clusterdata.area[noncolocalized_indexes])
        channel.meancoclustercircularity[i] = mean(channel.clusterdata.circularity[noncolocalized_indexes])
        channel.meancoclusterdensity[i] = mean([mean(channel.pointdata.density[c.core_indices]) / channel.roidensity for (i, c) ∈ enumerate(channel.clusterdata.cluster[noncolocalized_indexes])])
    end
    return cr
end

load_raw_results(path) = load(path)["results"]

end # module
