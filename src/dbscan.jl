function dbscan!(result::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        coordinates = hcat(pcr.coordinates for pcr ∈ result.pointschannelresults)
        coordinates = clusterparameters.uselocalradius_threshold ? coordinates[:, pointdata.abovethreshold] : coordinates
        clustersresult, sigclustersresult = dbscan!(coordinates, clusterparameters[1], i)
        push!(result.clusterresults, clustersresult)
        push!(result.sigclusterresults, sigclustersresult)
    else
        for i ∈ 1:result.nchannels
            coordinates = result.pointschannelresults[i].coordinates
            coordinates = clusterparameters.uselocalradius_threshold ? coordinates[:, pointdata.abovethreshold[pointdata.channel .== i]] : coordinates
            clustersresult, sigclustersresult = dbscan!(coordinates, clusterparameters[i], i)
            push!(result.clusterresults, clustersresult)
            push!(result.sigclusterresults, sigclustersresult)
        end
    end
end

function dbscan!(coordinates::Matrix{Float64}, clusterparameters::ClusterParameters, i::Int)
    length(coordinates) > 0 || return
    clusters = Clustering.dbscan(coordinates, clusterparameters.epsilon, min_cluster_size = clusterparameters.minpoints)
    abovethreshold = map(x -> x.size > clusterparameters.minsigclusterpoints, clusters)
    issignificant = getfield.(clusters, :size) .> clusterparameters.minsigclusterpoints

    c.clusterdata = DataFrame(
        :channel => i,
        :cluster => clusters,
        :size => [cluster.size for cluster ∈ clusters],
        :abovethreshold => abovethreshold,
        :issignificant => issignificant)
    
    nclusters = length(clusters)
    roiclusterdensity = c.nclusters / c.roiarea
    #fraction_clustered = sum(c.clusterdata.size) / c.nlocalizations

    sigclusters = clusters[issignificant]
    nsigclusters = length(sigclusters)
    roisigclusterdensity = c.nsigclusters / c.roiarea
    #fraction_sig_clustered = length(sigclusters) == 0 ? 0 : sum(c.size for c ∈ clusters_abovethreshold) / c.nlocalizations
    return ClustersResult(nclusters, roiclusterdensity), ClustersResult(nsigclusters, roisigclusterdensity)
end