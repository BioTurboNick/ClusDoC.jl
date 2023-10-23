function dbscan!(result::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        coordinates = hcat(pcr.coordinates for pcr ∈ result.pointschannelresults)
        coordinates = clusterparameters.uselocalradius_threshold ? coordinates[:, result.pointdata.abovethreshold] : coordinates
        clusterdata, clustersresult, sigclustersresult = dbscan!(result, coordinates, clusterparameters[1], i)
        push!(result.clusterresults, clustersresult)
        push!(result.sigclusterresults, sigclustersresult)
        result.clusterdata = [clusterdata]
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
    clusters = Clustering.dbscan(coordinates, clusterparameters.epsilon, min_cluster_size = clusterparameters.minpoints)
    abovethreshold = map(x -> x.size > clusterparameters.minsigclusterpoints, clusters)
    issignificant = getfield.(clusters, :size) .> clusterparameters.minsigclusterpoints

    clusterdata = DataFrame(
        :channel => i,
        :cluster => clusters,
        :size => [cluster.size for cluster ∈ clusters],
        :abovethreshold => abovethreshold,
        :issignificant => issignificant)
    
    nclusters = length(clusters)
    roiclusterdensity = nclusters / result.roiarea
    #fraction_clustered = sum(c.clusterdata.size) / c.nlocalizations

    sigclusters = clusters[issignificant]
    nsigclusters = length(sigclusters)
    roisigclusterdensity = nsigclusters / result.roiarea
    #fraction_sig_clustered = length(sigclusters) == 0 ? 0 : sum(c.size for c ∈ clusters_abovethreshold) / c.nlocalizations
    return clusterdata, ClustersResult(nclusters, roiclusterdensity), ClustersResult(nsigclusters, roisigclusterdensity)
end