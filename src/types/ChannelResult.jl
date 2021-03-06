mutable struct ChannelResult
    channelname::String
    coordinates::Matrix{Float64}
    nlocalizations::Int
    roiarea::Float64
    roidensity::Float64
    
    nclusters::Int
    roiclusterdensity::Float64
    meanclustersize::Float64
    meanclusterarea::Float64
    meanclustercircularity::Float64
    meanclusterabsolutedensity::Float64
    meanclusterdensity::Float64
    fraction_clustered::Float64

    nsigclusters::Int
    roisigclusterdensity::Float64
    meansigclustersize::Float64
    meansigclusterarea::Float64
    meansigclustercircularity::Float64
    meansigclusterabsolutedensity::Float64
    meansigclusterdensity::Float64
    fraction_sig_clustered::Float64
    
    pointdata::Union{Nothing, DataFrame}
    clusterdata::Union{Nothing, DataFrame} # ["cluster", "size", "area", "circularity", "contour", "ninteracting", "notherchannels"]
    fraction_interactions_clustered::Vector{Float64}

    # cocluster summary data. Entry where index is the same as this channel contains noncolocalized clusters
    meancoclustersize::Vector{Float64}
    meancoclusterarea::Vector{Float64}
    meancoclustercircularity::Vector{Float64}
    meancoclusterdensity::Vector{Float64}
    ncoclusters::Vector{Int}
    fraction_colocalized::Vector{Float64}


    ChannelResult(a, b, c, d, e, nchannels) = new(a, b, c, d, e, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, nothing, nothing,
        fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(0, nchannels), fill(NaN, nchannels))
end
