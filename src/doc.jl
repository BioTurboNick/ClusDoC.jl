using GeometryBasics: Point

"""
    doc(channels::Vector{Channel}, localradius, radiusmax, radiusstep, roiarea)

Calculate the degree of colocalization for all points in each channel against all other channels. Computes a density gradient
within `radiusmax` of each point with steps of size `radiusstep`. `localradius` is used to select foreground points with more
neighbors than expected by chance. Results are added to the `Channel` objects. Units are nanometers.
"""
function doc(channelnames, localizations, localradius, radiusmax, radiusstep, roiarea)
    length(channelnames) == length(localizations) ||
        throw(ArgumentError("$(:channelnames) must be the same length as $(:localizations)"))
    0 < localradius < radiusmax ||
        throw(ArgumentError("$(:localradius) must be positive and less than $(:radiusmax); got $localradius, $radiusmax"))
    0 < radiusstep < radiusmax ||
        throw(ArgumentError("$(:radiusstep) must be positive and less than $(:radiusmax); got $radiuistep, $radiusmax"))

    coordinates = map(localizations) do l
        if length(l) == 0
            fill(0.0, 3, 1) # downstream fails if there isn't at least one point.
        else
            extractcoordinates(l)
        end
    end
    channels = ChannelResult.(channelnames, coordinates, length.(localizations), roiarea / 1_000_000,
        length.(localizations) ./ (roiarea / 1_000_000), length(channelnames))
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
        nneighbors = count_inrange(allneighbortree, c.coordinates, localradius)
        ntotal = size(allcoordinates, 2) - 1 # remove self
        equivalentradii = equivalentradius.(nneighbors .- 1, ntotal, roiarea * 1_000_000)
        abovethreshold = equivalentradii .> localradius # maybe can replace with simple number threshold though, if don't need to compare across channels
        densities = pointdensity.(nneighbors, localradius)
        # calculate density gradient for each point vs. neighbors in each other channel
        distributions = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            k = Int(i == j) # factor to remove self
            dr = [(NearestNeighbors.count_inrange(ctrees[j], (@view c.coordinates[:, abovethreshold]), r) .- k) ./ r ^ 2 for r ∈ radiussteps]
            distributions[j] = hcat((drx ./ dr[end] for drx ∈ dr)...)
        end

        # compute degree of colocalization
        # original sets any NaNs to 0. I'm setting them to -1 because "nothing anywhere nearby" is as anticorrelated as it can get.
        docscores = Vector(undef, length(channels))
        for j ∈ eachindex(channels)
            length(channels[j].coordinates) > 0 || continue
            docscore = fill(NaN, length(abovethreshold))
            spearmancoefficient = [corspearman(distributions[i][k, :], distributions[j][k, :]) for k ∈ 1:count(abovethreshold)]
            _, nearestdistance = nn(ctrees[j], c.coordinates[:, abovethreshold])
            docscore[abovethreshold] = spearmancoefficient .* exp.(-nearestdistance ./ radiussteps[end])
            docscore[abovethreshold][isnan.(docscore[abovethreshold])] .= -1
            docscores[j] = docscore
        end

        channels[i].pointdata = DataFrame(:abovethreshold => abovethreshold, :density => densities, [Symbol(:docscore, j) => docscores[j] for j ∈ eachindex(channels)]...)
    end

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
    equivalentradii = sqrt(equivalentarea / π)
    return equivalentradii
end
