"""
    doc(channels::Vector{Channel}, localradius, radiusmax, radiusstep, roiarea)

Calculate the degree of colocalization for all points in each channel against all other channels. Computes a density gradient
within `radiusmax` of each point with steps of size `radiusstep`. `localradius` is used to select foreground points with more
neighbors than expected by chance. Results are added to the `Channel` objects. 
"""
function doc(localizations, localradius, radiusmax, radiusstep, roiarea)
    0 < localradius < radiusmax ||
        throw(ArgumentError("$(:localradius) must be positive and less than $(:radiusmax); got $localradius, $radiusmax"))
    0 < radiusstep < radiusmax ||
        throw(ArgumentError("$(:radiusstep) must be positive and less than $(:radiusmax); got $radiuistep, $radiusmax"))
    π * radiusmax ^ 2 ≤ roiarea ||
        throw(ArgumentError("$(:radiusmax) must describe a circle with area smaller than $(:roiarea); got $radiusmax"))

    channels = ChannelResult.(extractcoordinates.(localizations), nothing, nothing, nothing, nothing, nothing)
    #=
    The algorithm for coordinate-based colocalization (doi: 10.1007/s00418-011-0880-5) is:
    1. For each localization, count the number of localizations (other than itself) within a given radius for each channel.
    2. Divide by the area examined to get the density for each radius and each point.
    3. For each point, divide by the density within the maximum radius for each point.
    4. Steps 1-3 simplify to dividing the number observed in a given radius for each point, divided by the number
       observed within the maximum radius for each point, multiplied by the maximum radius divided by the given radius.
    5. Spearman's rank correlation coefficient is calculated between the self-self neighbor distribution and the self-other
       neighbor distribution.
    6. The coefficient is multipled by an exponential accounting for the distance to the nearest neighbor as a fraction of the
       `radiusmax`.
    
    There is an initial threshold step to ignore points that have fewer neighbors than expected based on a random
    distrubtion across the ROI.
    =#

    allcoordinates = reduce(hcat, c.coordinates for c ∈ channels)
    allneighbortree = BallTree(allcoordinates) # original uses KDTree, should compare
    ctrees = BallTree.(c.coordinates for c ∈ channels)
    radiussteps = (1:ceil(radiusmax / radiusstep)) .* radiusstep
    for (i, c) ∈ enumerate(channels)
        # determine which localizations have more neighbors than expected by chance
        nneighbors = NearestNeighbors.inrangecount(allneighbortree, c.coordinates, localradius, true)
        ntotal = size(allcoordinates, 2) - 1 # remove self
        c.equivalentradius = equivalentradius.(nneighbors, ntotal, roiarea)
        c.abovethreshold = c.equivalentradius .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        c.density = pointdensity.(nneighbors, localradius)

        # calculate density gradient for each point vs. neighbors in each other channel
        distributions = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            k = Int(i == j) # factor to remove self
            dr = [(NearestNeighbors.inrangecount(ctrees[j], (@view c.coordinates[:, c.abovethreshold]), r) .- k) ./ r ^ 2 for r ∈ radiussteps]
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
function pointdensity(nneighbors::Int, radius)
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
