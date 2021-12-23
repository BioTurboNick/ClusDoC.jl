mutable struct ChannelResult
    channelname::String
    coordinates::Matrix{Float64}
    roidensity::Float64
    pointdata::Union{Nothing, DataFrame}
    clusterdata::Union{Nothing, DataFrame}

    # cocluster summary data. Entry where index is the same as this channel contains noncolocalized clusters
    meancoclustersize::Vector{Float64}
    meancoclusterarea::Vector{Float64}
    meancoclustercircularity::Vector{Float64}
    meancoclusterdensity::Vector{Float64}
    ncoclusters::Vector{Int}
    fraction_colocalized::Vector{Float64}

    ChannelResult(a, b, c, nchannels) = new(a, b, c, nothing, nothing,
        fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(NaN, nchannels), fill(0, nchannels), fill(NaN, nchannels))
end
