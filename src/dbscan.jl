#=
[~, ClusterCh, ~, classOut, ~, ~, ~, ResultCell{roiIter, cellIter}] = DBSCANHandler(Data_DoC1(:,5:6), ...
                        dbscanParams, cellIter, roiIter, true, false, clusterColor, Data_DoC1(:, NDatacolumns + 2), Data_DoC1(:, NDatacolumns + 6), ...
                        Data_DoC1(:, NDatacolumns + 4));



Algorithm
1. Use Lr on all points to threshold (optionally)


                        =#

function dbscan(channels::Vector{Channel})
    for c ∈ channels
        c.clusters = dbscan(coordinates, eps, minpts)
    end
end