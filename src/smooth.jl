function smooth!(cr)
    for (i, c) ∈ enumerate(cr)
        Nb = [cluster.size for cluster ∈ c.clusters]
        #sigmas = Smothing Radius parameter; not sure what that is
        clustercoordinates = [@view c.coordinates[:, union(cluster.core_indices, cluster.boundary_indices)] for cluster ∈ c.clusters]
        bounds =  extrema.(clustercoordinates, dims = 2)
        lengths = [last.(b) .- first.(b) for b in bounds]
        boxsizes = 0.5 .* maximum.(lengths) .+ (epsilon + 10)

        centers = [(last.(b) .+ first.(b)) ./ 2 for b ∈ bounds]
        #epsilon = dbscan parameter
        boxmins = [cen .- bs for (cen, bs) ∈ zip(centers, boxsizes)]
        boxmaxes = [cen .+ bs for (cen, bs) ∈ zip(centers, boxsizes)]

        # create the grid
        grids = [UnitRange.(floor.(bmin), ceil.(bmax)) for (bmin, bmax) ∈ zip(boxmins, boxmaxes)]

        
    end
end