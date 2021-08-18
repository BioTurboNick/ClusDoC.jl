mutable struct Channel
    coordinates
    density
    equivalentradius
    abovethreshold
    docscore
end

channel1pointsX = [28739.6, 29635, 28894.8, 29151.3, 28824.9, 29491.4, 29129.6, 29131.9, 28574.5, 29492.4, 29132.5,
                       28580.6, 28561.3, 28577.5, 28739.6, 29633.5, 29131.6]
channel1pointsY = [37181.4, 37849.4, 37446.7, 37957, 37026.1, 37287.2, 37786.3, 37959, 37558.8, 37309.7, 37771.8,
                    37555, 37571.5, 37568.2, 37194.2, 37851.2, 37761]
coordinates = repeat(permutedims([channel1pointsX channel1pointsY]), inner = (1, 1000))

channels = Vector(undef, 2)
channels[1] = Channel(coordinates, nothing, nothing, nothing, nothing)
channels[2] = Channel(coordinates, nothing, nothing, nothing, nothing)



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
    allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
    
    allneighbortree = BallTree(allcoordinates) # original uses KDTree, should compare
    ctrees = BallTree.(c.coordinates for c ∈ channels)
    radiussteps = (1:ceil(radiusmax / step)) .* step
    for (i, c) ∈ enumerate(channels)
        # determine which localizations have more neighbors than expected by chance
        ineighbors = inrange(allneighbortree, c.coordinates, localradius, true)
        nneighbors = length.(ineighbors) .- 1 # remove self
        ntotal = size(allcoordinates, 2) - 1 # remove self
        c.equivalentradius = equivalentradius.(nneighbors, ntotal, roiarea)
        c.abovethreshold = c.equivalentradius .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        c.density = density.(nneighbors, localradius)

        # calculate density gradient for each point vs. neighbors in each other channel
        distributions = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            k = Int(i == j) # factor to remove self
            # I think this is the major bottleneck*********************************************
            dr = [(length.(inrange(ctrees[j], (@view c.coordinates[:, c.abovethreshold]), r)) .- k) ./ r ^ 2 for r ∈ radiussteps]
            distributions[j] = hcat((drx ./ dr[end] for drx ∈ dr)...)
        end

        # compute degree of colocalization
        c.docscore = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            spearmancoefficient = [corspearman(distributions[i][k, :], distributions[j][k, :]) for k ∈ 1:count(c.abovethreshold)]
            _, nearestdistance = nn(ctrees[j], c.coordinates[:, c.abovethreshold])
            c.docscore[j] = spearmancoefficient .* exp.(-nearestdistance ./ radiussteps[end])
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
        throw(ArgumentError("$(:nneighbors) and $(:ntotal) must be greater than or equal to zero, and $(:ntotal) must be greater than or equal to $(:nneighbors); got $nneighbors, $ntotal"))
    roiarea > 0 ||
        throw(ArgumentError("$(:roiarea) must be positive; got $roiarea"))
    fneighbors = nneighbors / ntotal
    equivalentarea = fneighbors * roiarea
    equivalentradius = sqrt(equivalentarea / π)
    return equivalentradius
end
