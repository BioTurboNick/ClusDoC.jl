#=
[~, ClusterCh, ~, classOut, ~, ~, ~, ResultCell{roiIter, cellIter}] = DBSCANHandler(Data_DoC1(:,5:6), ...
                        dbscanParams, cellIter, roiIter, true, false, clusterColor, Data_DoC1(:, NDatacolumns + 2), Data_DoC1(:, NDatacolumns + 6), ...
                        Data_DoC1(:, NDatacolumns + 4));



Algorithm
1. Use Lr on all points to threshold (optionally)


                        =#

function dbscan(channels::Vector{Channel}, epsilon, minpoints, uselocalradius_threshold, localradius)

    for c ∈ channels
        if uselocalradius_threshold
            allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
            allneighbortree = BallTree(allcoordinates)
            nneighbors = inrangecount(allneighbortree, c.coordinates, localradius, true)
            c.equivalentradius = equivalentradius.(nneighbors, ntotal, roiarea)
            c.abovethreshold = c.equivalentradius .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        end
        
        c.clusters = dbscan(coordinates, epsilon, minpoints)
    end
end