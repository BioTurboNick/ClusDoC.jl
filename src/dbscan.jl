#=
[~, ClusterCh, ~, classOut, ~, ~, ~, ResultCell{roiIter, cellIter}] = DBSCANHandler(Data_DoC1(:,5:6), ...
                        dbscanParams, cellIter, roiIter, true, false, clusterColor, Data_DoC1(:, NDatacolumns + 2), Data_DoC1(:, NDatacolumns + 6), ...
                        Data_DoC1(:, NDatacolumns + 4));



Algorithm
1. Use Lr on all points to threshold (optionally)


                        =#

# if channelresult already has abovethreshold set, it'll use that instead of recalculating.

function dbscan!(channels::Vector{ChannelResult}, epsilon, minpoints, uselocalradius_threshold, localradius)
    for c ∈ channels
        if uselocalradius_threshold && isnothing(c.abovethreshold)
            allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
            allneighbortree = BallTree(allcoordinates)
            nneighbors = inrangecount(allneighbortree, c.coordinates, localradius, true)
            c.equivalentradius = equivalentradius.(nneighbors, ntotal, roiarea)
            c.abovethreshold = c.equivalentradius .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        end
        
        c.clusters = Clustering.dbscan(c.coordinates, epsilon, min_cluster_size = minpoints)
    end
end