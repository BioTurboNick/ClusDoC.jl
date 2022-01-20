# if channelresult already has abovethreshold set, it'll use that instead of recalculating.

function dbscan!(channels::Vector{ChannelResult}, clusterparameters, localradius)
    for (i, c) ∈ enumerate(channels)
        if clusterparameters[i].uselocalradius_threshold && isnothing(c.pointdata.abovethreshold)
            allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
            allneighbortree = BallTree(allcoordinates)
            nneighbors = inrangecount(allneighbortree, c.coordinates, localradius, true)
            equivalentradii = equivalentradius.(nneighbors, ntotal, roiarea)
            c.pointdata.abovethreshold = equivalentradii .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        end
        coordinates = clusterparameters[i].uselocalradius_threshold ? c.coordinates[:, c.pointdata.abovethreshold] : c.coordinates
        clusters = Clustering.dbscan(coordinates, clusterparameters[i].epsilon, min_cluster_size = clusterparameters[i].minpoints)
        c.clusterdata = DataFrame(:cluster => clusters, :size => [cluster.size for cluster ∈ clusters])
        c.nclusters = length(clusters)
        c.roiclusterdensity = c.nclusters / c.roiarea
        c.meanclustersize = mean(c.clusterdata.size)
        c.fraction_clustered = (c.nlocalizations - sum(c.clusterdata.size)) / c.nlocalizations
    end
end