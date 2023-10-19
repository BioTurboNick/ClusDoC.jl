mutable struct ROIResult
    roiarea::Float64
    roidensity::Float64
    
    clusterdata::Some{DataFrame}
    pointdata::Some{DataFrame}

    nchannels::Int
    channelnames::Vector{String}

    clusterresults::Vector{ClustersResult}
    sigclusterresults::Vector{ClustersResult}
    coclusterresults::Vector{ClustersResult}
    intermediatecoclusterresults::Vector{ClustersResult}

    pointschannelresults::Vector{PointsChannelResult}

    ChannelResult(roiarea, roidensity, channelnames, coordinates, nlocalizations) =
        new(roiarea, roidensity, length(channelnames), channelnames, ClustersResult[], ClustersResult[], ClustersResult[], ClustersResult[],
        [PointsChannelResult(c, n, NaN, NaN) for (c, n) ∈ zip(coordinates, nlocalizations)])
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

# Store data specific to a channel within a cluster.
mutable struct ClustersChannelResult
    meanclustersize::Float64
    meanclusterabsolutedensity::Float64
    meanclusterdensity::Float64
    fraction_of_interacting_points::Float64

    ClusterChannelResult(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered) =
        new(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, NaN, NaN)
end

mutable struct PointsChannelResult
    coordinates::Matrix{Float64}
    nlocalizations::Int
    fraction_colocalized::Float64
    fraction_clustered::Float64
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
