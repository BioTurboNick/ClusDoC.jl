function dbscan!(channels::Vector{ChannelResult}, clusterparameters, localradius)
    for (i, c) ∈ enumerate(channels)
        coordinates = clusterparameters[i].uselocalradius_threshold ? c.coordinates[:, c.pointdata.abovethreshold] : c.coordinates
        length(coordinates) > 0 || continue
        clusters = Clustering.dbscan(coordinates, clusterparameters[i].epsilon, min_cluster_size = clusterparameters[i].minpoints)
        abovethreshold = map(x -> x.size > clusterparameters[i].minsigclusterpoints, clusters)
        ninteracting = [[count(c.pointdata[!, Symbol(:docscore, j)][cluster.core_indices] .> 0.4) for j ∈ eachindex(channels)] for cluster ∈ clusters]
        c.clusterdata = DataFrame(:cluster => clusters, :size => [cluster.size for cluster ∈ clusters], :ninteracting => ninteracting, :abovethreshold => abovethreshold)
        c.nclusters = length(clusters)
        c.roiclusterdensity = c.nclusters / c.roiarea
        c.meanclustersize = mean(c.clusterdata.size)
        c.fraction_clustered = sum(c.clusterdata.size) / c.nlocalizations

        clusters_abovethreshold = filter(x -> x.size > clusterparameters[i].minsigclusterpoints, clusters)
        c.nsigclusters = length(clusters_abovethreshold)
        c.roisigclusterdensity = c.nsigclusters / c.roiarea
        c.meansigclustersize = length(clusters_abovethreshold) == 0 ? NaN : mean(c.size for c ∈ clusters_abovethreshold)
        c.fraction_sig_clustered = length(clusters_abovethreshold) == 0 ? 0 : sum(c.size for c ∈ clusters_abovethreshold) / c.nlocalizations
    end
end