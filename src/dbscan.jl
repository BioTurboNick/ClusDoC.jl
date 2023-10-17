function dbscan!(channels::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        combined = ChannelResult(
            "combinedtemp",
            hcat([c.coordinates for c ∈ channels]...),
            sum(c.nlocalizations for c ∈ channels),
            sum(c.roiarea for c ∈ channels),
            channels[1].roiarea,
            0)
        combined.pointdata = copy(channels[1].pointdata)
        for c ∈ channels[2:end]
            append!(combined.pointdata, c.pointdata)
        end
        dbscan!(combined, clusterparameters[1])

        for c ∈ channels
            c.clusterdata = combined.clusterdata
            c.nclusters = combined.nclusters
            c.roiclusterdensity = combined.roiclusterdensity
            c.meanclustersize = combined.meanclustersize
            c.fraction_clustered = combined.fraction_clustered
            c.nsigclusters = combined.nsigclusters
            c.roisigclusterdensity = combined.roisigclusterdensity
            c.meansigclustersize = combined.meansigclustersize
            c.fraction_sig_clustered = combined.fraction_sig_clustered
        end
    else
        for (i, c) ∈ enumerate(channels)
            dbscan!(c, clusterparameters[i])
        end
    end
end

function dbscan!(channel::ChannelResult, clusterparameters)
    c = channel
    coordinates = clusterparameters.uselocalradius_threshold ? c.coordinates[:, c.pointdata.abovethreshold] : c.coordinates
    length(coordinates) > 0 || return
    clusters = Clustering.dbscan(coordinates, clusterparameters.epsilon, min_cluster_size = clusterparameters.minpoints)
    abovethreshold = map(x -> x.size > clusterparameters.minsigclusterpoints, clusters)
    c.clusterdata = DataFrame(
        :cluster => clusters,
        :size => [cluster.size for cluster ∈ clusters],
        :abovethreshold => abovethreshold)
    c.nclusters = length(clusters)
    c.roiclusterdensity = c.nclusters / c.roiarea
    c.meanclustersize = mean(c.clusterdata.size)
    c.fraction_clustered = sum(c.clusterdata.size) / c.nlocalizations

    clusters_abovethreshold = filter(x -> x.size > clusterparameters.minsigclusterpoints, clusters)
    c.nsigclusters = length(clusters_abovethreshold)
    c.roisigclusterdensity = c.nsigclusters / c.roiarea
    c.meansigclustersize = length(clusters_abovethreshold) == 0 ? NaN : mean(c.size for c ∈ clusters_abovethreshold)
    c.fraction_sig_clustered = length(clusters_abovethreshold) == 0 ? 0 : sum(c.size for c ∈ clusters_abovethreshold) / c.nlocalizations
end
