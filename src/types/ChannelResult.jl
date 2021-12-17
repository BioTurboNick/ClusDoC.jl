mutable struct ChannelResult
    channelname
    coordinates
    densities
    equivalentradius
    abovethreshold
    docscores
    roidensity
    clusters
    clusternpoints
    clusterareas
    clustercircularities
    clustercontours
end

ChannelResult(a, b, c, d, e, f, g) =
        ChannelResult(a, b, c, d, e, f, g, nothing, nothing, nothing, nothing, nothing)
