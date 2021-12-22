mutable struct ChannelResult
    channelname
    coordinates
    roidensity
    pointdata::Union{Nothing, DataFrame}
    clusterdata::Union{Nothing, DataFrame}

    ChannelResult(a, b, c) = new(a, b, c, nothing, nothing)
end
