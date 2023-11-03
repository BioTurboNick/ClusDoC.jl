
mutable struct PointsChannelResult
    coordinates::Matrix{Float64}
    nlocalizations::Int
    nlocalizations_abovethreshold::Int
    roidensity::Float64
    fraction_colocalized::Float64
end

# Store data specific to a channel within a cluster.
mutable struct ClustersChannelResult
    meanclustersize::Float64
    meanclusterabsolutedensity::Float64
    meanclusterdensity::Float64
    fraction_clustered::Float64
    fraction_of_interacting_points12::Float64
    fraction_of_interacting_points21::Float64

    ClustersChannelResult(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, fraction_interacting12, fraction_interacting21) =
        new(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, fraction_interacting12, fraction_interacting21)
end

mutable struct ClustersResult
    nclusters::Int
    roiclusterdensity::Float64
    meanclusterarea::Float64
    meanclustercircularity::Float64

    channelresults::Vector{ClustersChannelResult}

    ClustersResult(nclusters, roiclusterdensity) =
        new(nclusters, roiclusterdensity, NaN, NaN, ClustersChannelResult[])
end

mutable struct ROIResult
    roiarea::Float64
    
    clusterdata::Union{Nothing, DataFrame}
    pointdata::Union{Nothing, DataFrame}

    nchannels::Int
    channelnames::Vector{String}

    clusterresults::Vector{ClustersResult}
    sigclusterresults::Vector{ClustersResult}
    coclusterresults::Vector{ClustersResult}
    intermediatecoclusterresults::Vector{ClustersResult}
    noncolocalizedclusterresults::Union{Nothing, ClustersResult}

    pointschannelresults::Vector{PointsChannelResult}

    ROIResult(roiarea, roidensity, channelnames, coordinates, nlocalizations) =
        new(roiarea, nothing, nothing, length(channelnames), channelnames, ClustersResult[], ClustersResult[], ClustersResult[], ClustersResult[], nothing,
        [PointsChannelResult(c, n, -1, d, NaN) for (c, n, d) ∈ zip(coordinates, nlocalizations, roidensity)])
end



# mutable struct ChannelResult
#     channelname::String
#     coordinates::Matrix{Float64}
#     nlocalizations::Int
#     roiarea::Float64
#     roidensity::Float64
    
#     # nclusters::Int
#     # roiclusterdensity::Float64
#     # meanclustersize::Float64
#     # meanclusterarea::Float64
#     # meanclustercircularity::Float64
#     # meanclusterabsolutedensity::Float64
#     # meanclusterdensity::Float64
#     # fraction_clustered::Float64

#     # nsigclusters::Int
#     # roisigclusterdensity::Float64
#     # meansigclustersize::Float64
#     # meansigclusterarea::Float64
#     # meansigclustercircularity::Float64
#     # meansigclusterabsolutedensity::Float64
#     # meansigclusterdensity::Float64
#     # fraction_sig_clustered::Float64
    
#     # pointdata::Union{Nothing, DataFrame}
#     # clusterdata::Union{Nothing, DataFrame} # ["cluster", "size", "area", "circularity", "contour", "ninteracting", "notherchannels"]
#     fraction_colocalized::Vector{Float64}

#     # cocluster (at least minsigclusterpoints interacting) summary data. Entry where index is the same as this channel contains noncolocalized clusters
#     meancoclustersize::Vector{Float64}
#     meancoclusterarea::Vector{Float64}
#     meancoclustercircularity::Vector{Float64}
#     meancoclusterdensity::Vector{Float64}
#     ncoclusters::Vector{Int}
#     fraction_interactions_clustered::Vector{Float64}

#     # intermediate cocluster (1 to minsigclusterpoints interacting) summary data. Entry where index is the same as this channel is not used
#     meancoclustersize_int::Vector{Float64}
#     meancoclusterarea_int::Vector{Float64}
#     meancoclustercircularity_int::Vector{Float64}
#     meancoclusterdensity_int::Vector{Float64}
#     ncoclusters_int::Vector{Int}
#     fraction_interactions_clustered_int::Vector{Float64}

#     ChannelResult(a, b, c, d, e, nchannels) = new(a, b, c, d, e, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, nothing, nothing,
#         fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(0, nchannels), fill(NaN, nchannels),
#                               fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(0, nchannels), fill(NaN, nchannels))
# end
