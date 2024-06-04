
# Overwrite Clustering method to remove error when there are few points.
function Clustering.dbscan(points::AbstractMatrix, radius::Real; metric = Clustering.Euclidean(), min_neighbors::Integer = 1, min_cluster_size::Integer = 1, nntree_kwargs...)
    0 <= radius || throw(ArgumentError("radius $radius must be ≥ 0"))                                                  

    if metric !== nothing
        # points are point coordinates
        dim, num_points = size(points)
        #num_points <= dim && throw(ArgumentError("points has $dim rows and $num_points columns. Must be a D x N matric with D < N"))
        kdtree = Clustering.KDTree(points, metric; nntree_kwargs...)
        data = (kdtree, points)
    else
        # points is a distance matrix
        num_points = size(points, 1)
        size(points, 2) == num_points || throw(ArgumentError("When metric=nothing, points must be a square distance matrix ($(size(points)) given)."))
        num_points >= 2 || throw(ArgumentError("At least two data points are required ($num_points given)."))
        data = points
    end
    clusters = Clustering._dbscan(data, num_points, radius, min_neighbors, min_cluster_size)
    return Clustering.DbscanResult(clusters, num_points)
end

function dbscan!(result::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        coordinates = hcat(collect(pcr.coordinates for pcr ∈ result.pointschannelresults)...)
        coordinates = clusterparameters[1].uselocalradius_threshold ? coordinates[:, result.pointdata.abovethreshold] : coordinates
        clusterdata, clustersresult, sigclustersresult = dbscan!(result, coordinates, clusterparameters[1], 0)
        push!(result.clusterresults, clustersresult)
        push!(result.sigclusterresults, sigclustersresult)
        result.clusterdata = clusterdata
    else
        clustersdata = DataFrame()
        for i ∈ 1:result.nchannels
            coordinates = result.pointschannelresults[i].coordinates
            channelpointdata = filter(x -> x.channel == i, result.pointdata)
            coordinates = clusterparameters[i].uselocalradius_threshold ? coordinates[:, channelpointdata.abovethreshold] : coordinates
            clusterdata, clustersresult, sigclustersresult = dbscan!(result, coordinates, clusterparameters[i], i)
            append!(clustersdata, clusterdata)
            push!(result.clusterresults, clustersresult)
            push!(result.sigclusterresults, sigclustersresult)
        end
        result.clusterdata = clustersdata
    end
end

function dbscan!(result::ROIResult, coordinates::Matrix{Float64}, clusterparameters::ClusterParameters, i::Int)
    length(coordinates) > 0 || return
    res = Clustering.dbscan(coordinates, clusterparameters.epsilon; min_cluster_size = clusterparameters.minpoints)
    clusters = res.clusters
    issignificant = getfield.(clusters, :size) .> clusterparameters.minsigclusterpoints

    clusterdata = DataFrame(
        :channel => i,
        :cluster => clusters,
        :size => [cluster.size for cluster ∈ clusters],
        :issignificant => issignificant)
    
    nclusters = length(clusters)
    roiclusterdensity = nclusters / result.roiarea

    sigclusters = clusters[issignificant]
    nsigclusters = length(sigclusters)
    roisigclusterdensity = nsigclusters / result.roiarea
    return clusterdata, ClustersResult(nclusters, roiclusterdensity), ClustersResult(nsigclusters, roisigclusterdensity)
end