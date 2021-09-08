# channel1pointsX = [28739.6, 29635, 28894.8, 29151.3, 28824.9, 29491.4, 29129.6, 29131.9, 28574.5, 29492.4, 29132.5,
#                        28580.6, 28561.3, 28577.5, 28739.6, 29633.5, 29131.6]
# channel1pointsY = [37181.4, 37849.4, 37446.7, 37957, 37026.1, 37287.2, 37786.3, 37959, 37558.8, 37309.7, 37771.8,
#                     37555, 37571.5, 37568.2, 37194.2, 37851.2, 37761]
# coordinates = repeat(permutedims([channel1pointsX channel1pointsY]), inner = (1, 1000))

# channels = Vector(undef, 2)
# channels[1] = ChannelResult("test1", coordinates, nothing, nothing, nothing, nothing, nothing)
# channels[2] = ChannelResult("test2", coordinates, nothing, nothing, nothing, nothing, nothing)

"""
    doc(channels::Vector{Channel}, localradius, radiusmax, radiusstep, roiarea)

Calculate the degree of colocalization for all points in each channel against all other channels. Computes a density gradient
within `radiusmax` of each point with steps of size `radiusstep`. `localradius` is used to select foreground points with more
neighbors than expected by chance. Results are added to the `Channel` objects. 
"""
function doc(channelnames, localizations, localradius, radiusmax, radiusstep, roiarea)
    0 < localradius < radiusmax ||
        throw(ArgumentError("$(:localradius) must be positive and less than $(:radiusmax); got $localradius, $radiusmax"))
    0 < radiusstep < radiusmax ||
        throw(ArgumentError("$(:radiusstep) must be positive and less than $(:radiusmax); got $radiuistep, $radiusmax"))
    π * radiusmax ^ 2 ≤ roiarea ||
        throw(ArgumentError("$(:radiusmax) must describe a circle with area smaller than $(:roiarea); got $radiusmax"))

    channels = ChannelResult.(channelnames, extractcoordinates.(localizations), nothing, nothing, nothing, nothing, nothing)
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
    allneighbortree = BallTree(allcoordinates) # original uses KDTree; I timed it and it is worse
    ctrees = BallTree.(c.coordinates for c ∈ channels)
    radiussteps = (1:ceil(radiusmax / radiusstep)) .* radiusstep
    for (i, c) ∈ enumerate(channels)
        # determine which localizations have more neighbors than expected by chance
        nneighbors = inrangecount(allneighbortree, c.coordinates, localradius)
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
        c.docscores = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            spearmancoefficient = [corspearman(distributions[i][k, :], distributions[j][k, :]) for k ∈ 1:count(c.abovethreshold)]
            _, nearestdistance = nn(ctrees[j], c.coordinates[:, c.abovethreshold])
            c.docscores[j] = spearmancoefficient .* exp.(-nearestdistance ./ radiussteps[end])
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
