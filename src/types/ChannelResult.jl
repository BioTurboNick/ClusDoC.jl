mutable struct ChannelResult
    channelname
    coordinates
    densities
    equivalentradius
    abovethreshold
    docscores
    clusters
    clusternpoints
    clusterimages
    clusterareas
    clustercircularities
    clustercontours
    clusterboxes
    clustercutoffpoints
end

ChannelResult(a, b, c, d, e, f, g) =
        ChannelResult(a, b, c, d, e, f, g, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
