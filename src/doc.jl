mutable struct Channel
    coordinates
    Lr
    allrelativedensity
    abovethreshold
    docscores
end


# calculate degree of colocalization (DoC) scores
function doc!(channels, localradius, radiusmax, step, roiarea)
    0 < localradius < radiusmax ||
        throw(ArgumentError("$(:localradius) must be positive and less than $(:radiusmax); got $localradius, $radiusmax"))
    π * radiusmax ^ 2 ≤ roiarea ||
        throw(ArgumentError("$(:radiusmax) must describe a circle with area smaller than $(:roiarea); got $radiusmax"))

    #=
    The algorithm for coordinate-based colocalization (10.1007/s00418-011-0880-5) is:
    1. For each localization, count the number of localizations (other than itself) within a given radius for each channel.
    2. Divide by the area examined to get the density for each radius and each point.
    3. For each point, divide by the density within the maximum radius for each point.
    4. Steps 1-3 simplify to dividing the number observed in a given radius for each point, divided by the number
       observed within the maximum radius for each point, multiplied by the maximum radius divided by the given radius.
    5. Spearman's rank correlation coefficient is calculated between the self-self neighbor distribution and the self-other
       neighbor distribution.
    
    The ClusDoC algorithm has an initial threshold step to decide which points to include in the calculation.
    However, it seems to be broken; either the Lr calculation is producing too-high values, or the threshold is too low.
    Essentially, all points that have any neighbors are above the threshold. I think what it is supposed to be is a calculation
    of the number of points expected to be within Lr_radius of a point assuming all were evenly distributed across the ROI. And the
    Lr for a single point would then just be the number of neighbors within Lr_radius of it. But the Lr calculation appears instead to
    calculate the radius of a circle in which that fraction of points would be expected to appear, with the ROI area doubled. (holdover
    from a previous iteration in which this value was the side of the ROI square?)

    I *think* then that the threshold is really just the Lr_radius. Which is what the original code complained about.
    =#



    println("Segment clustered points from background...")

    # This threshold was apparently simplified from a more convoluted one which
    # a comment said basically equalled Lr_radius
    # Lr_Threshold should be number of points within Lr_r for a random
    # distrubution of the same number of points in the current ROI


    # The threshold is the number of points that would fall within Lr_radius if randomly distributed over the ROI area.
    # originally, this value was treated the same as the Lr function does.
    allcoordinates = [c.coordinates for c ∈ channels]
    
    allneighbortree = BallTree(allcoordinates) # original uses KDTree, should compare
    ctrees = BallTree.(c.coordinates for c ∈ channels)
    radiussteps = (1:ceil(radiusmax / step)) .* step
    for (i, c) ∈ enumerate(channels)
        ineighbors = inrange(allneighbortree, c.coordinates, radius, true)
        nneighbors = length.(ineighbors) .- 1 # remove self
        ntotal = length(allneighbortree) - 1 # remove self
        c.equivalentradius = equivalentradius.(nneighbors, ntotal, roiarea)
        c.density = density.(nneighbors, localradius)
        c.abovethreshold = c.equivalentradius .> localradius # maybe can replace with simple number threshold though?

        distributions = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            k = Int(i == j) # factor to remove self
            dr = [(length.(inrange(ctrees[j], c.coordinates, r)) .- k) ./ r ^ 2 for r ∈ radiussteps]
            distributions[j] = dr ./ dr[end]
        end

        for j ∈ eachindex(channels)
            spearmancoefficient = corspearman.(c.distributions[i], distributions[j])
            _, nearestdistance = nn.(ctrees[j], c.coordinates)
            c.docscore = spearmancoefficient .* exp(-nearestdistance / radiussteps[end])
        end
    end

    # original sets any NaNs in SA1/2 to 0

    return channels
    # original ends a stopwatch here
end

"""
    density(nneighbors::Int, radius)

Calculate the density of points in a circle of a given radius.
"""
function density(nneighbors::Int, radius)
    nneighbors ≥ 0 ||
        throw(ArgumentError("$(:nneighbors) must be greater than equal to zero; got $nneighbors"))
    radius > 0 ||
        throw(ArgumentError("$(:radius) must be positive; got $radius"))
    return nneighbors / (π * radius ^ 2)
end

"""
    equivalentradius(nneighbors::Int, ntotal::Int, roiarea)

Calculate the radius of a circle that would contain the given number of neighbors, assuming the neighbors
were evenly distributed across the ROI.
"""
function equivalentradius(nneighbors::Int, ntotal::Int, roiarea) #Lr
    0 < ntotal ||
        throw(ArgumentError("$(:ntotal) must be positive; got $notal"))
    0 ≤ nneighbors ≤ ntotal ||
        throw(ArgumentError("$(:nneighbors) and $(:ntotal) must be greater than or equal to zero, and $(:ntotal) must be greater \
                             than or equal to $(:nneighbors); got $nneighbors, $ntotal"))
    roiarea > 0 ||
        throw(ArgumentError("$(:roiarea) must be positive; got $roiarea"))
    fneighbors = nneighbors / ntotal
    equivalentarea = fneighbors * roiarea
    equivalentradius = sqrt(equivalentarea / π)
    return equivalentradius
end
