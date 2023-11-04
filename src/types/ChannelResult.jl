
mutable struct PointsChannelResult
    coordinates::Matrix{Float64}
    nlocalizations::Int
    nlocalizations_abovethreshold::Int
    roidensity::Float64
    fraction_colocalized::Vector{Float64}
end

# Store data specific to a channel within a cluster.
mutable struct ClustersChannelResult
    meanclustersize::Float64
    meanclusterabsolutedensity::Float64
    meanclusterdensity::Float64
    fraction_clustered::Float64
    fraction_of_interacting_points::Float64

    ClustersChannelResult(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, fraction_interacting) =
        new(meanclustersize, meanclusterabsolutedensity, meanclusterdensity, fraction_clustered, fraction_interacting)
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
    coclusterresults::Vector{Vector{ClustersResult}}
    intermediatecoclusterresults::Vector{Vector{ClustersResult}}
    noncolocalizedclusterresults::Vector{ClustersResult}

    pointschannelresults::Vector{PointsChannelResult}

    ROIResult(roiarea, roidensity, channelnames, coordinates, nlocalizations) =
        new(roiarea, nothing, nothing, length(channelnames), channelnames, ClustersResult[], ClustersResult[], Vector{ClustersResult}[], Vector{ClustersResult}[], ClustersResult[],
        [PointsChannelResult(c, n, -1, d, fill(NaN, length(channelnames))) for (c, n, d) ∈ zip(coordinates, nlocalizations, roidensity)])
end
