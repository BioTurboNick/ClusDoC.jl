# if channelresult already has abovethreshold set, it'll use that instead of recalculating.

function dbscan!(channels::Vector{ChannelResult}, epsilon, minpoints, uselocalradius_threshold, localradius)
    for c ∈ channels
        if uselocalradius_threshold && isnothing(c.docdata.abovethreshold)
            allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
            allneighbortree = BallTree(allcoordinates)
            nneighbors = inrangecount(allneighbortree, c.coordinates, localradius, true)
            equivalentradii = equivalentradius.(nneighbors, ntotal, roiarea)
            c.docdata.abovethreshold = equivalentradius .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        end
        coordinates = uselocalradius_threshold ? c.coordinates[:, c.docdata.abovethreshold] : c.coordinates
        clusters = Clustering.dbscan(coordinates, epsilon, min_cluster_size = minpoints)
        c.clusterdata = DataFrame(:cluster => clusters, :size => [cluster.size for cluster ∈ clusters])
    end
end