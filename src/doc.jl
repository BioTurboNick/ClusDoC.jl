# calculate degree of colocalization (DoC) scores
function doc!(channels, Lr_radius, radiusmax, step, roiarea)
    println("Segment clustered points from background...")
    
    # This threshold was apparently simplified from a more convoluted one which
    # a comment said basically equalled Lr_radius
    # Lr_Threshold should be number of points within Lr_r for a random
    # distrubution of the same number of points in the current ROI
    Lr_threshold = length(allcoordinates) * π * (Lr_radius ^ 2) / roiarea; # original divided by ROI size
    
    allcoordinates = [c.coordinates for c ∈ channels]
    allneighbortree = BallTree(allcoordinates) # original uses KDTree, should compare
    ctrees = BallTree.(c.coordinates for c ∈ channels)
    for (i, c) ∈ enumerate(channels)
        c.Lr, c.allrelativedensity = lr(allneighbortree, c.coordinates, lr_radius)
        _, c.relativedensity = lr(ctrees[i], c.coordinates, lr_radius)
        c.abovethreshold = Lr .> Lr_threshold
    end
    
    # assumes two channels
    # not sure that these denominators are all correct, but they are what the original code used...
    D1max = count(channels[1].abovethreshold) / roiarea ^ 2 # should it be squared?
    D2max = count(channels[2].abovethreshold) / roiarea ^ 2
    D2maxCh2Ch1 = length.(first(lr(channels[2].coordinates, channels[1].coordinates, radiusmax))) / radiusmax ^ 2 # check for Ch2 points that are within radiusmax of each Ch1 point

    # original starts a stopwatch index here

    println("Calculating DoC scores...")

    # original does this in a parallel loop
    nsteps = ceil(radiusmax / step)
    N11 = Vector(undef, nsteps)
    N12 = Vector(undef, nsteps)
    N22 = Vector(undef, nsteps)
    N21 = Vector(undef, nsteps)
    for (i, r) ∈ enumerate((1:nsteps) .* step)
        r2 = r ^ 2
        N11[i] = inrange(ctrees[1], c[1].coordinates, r) .- 1 ./ (D1max * r2)
        N12[i] = inrange(ctrees[2], c[1].coordinates, r) ./ (D2max * r2)
        N22[i] = inrange(ctrees[2], c[2].coordinates, r) .- 1 ./ (D1max * r2)
        N21[i] = inrange(ctrees[1], c[2].coordinates, r) ./ (D2maxCh2Ch1 * r2)
    end

    println("Correlating coefficients...")

    SA1 = corspearman.(N11, N12)
    SA2 = corspearman.(N22, N21)

    # original sets any NaNs in SA1/2 to 0

    # DoC Assignment
    # finds distance to the nearest neighbor between each channel
    _, nearestdistance1 = nn(ctrees[2], channels[1].coordinates)
    _, nearestdistance2 = nn(ctrees[2], channels[2].coordinates)

    # DoC is the correlation score multiplied by an exponential that scales with the distance to the nearest neighbor as a fraction of radiusmax
    channels[1].doc = SA1 .* exp(-nearestdistance1 / radiusmax)
    channels[2].doc = SA2 .* exp(-nearestdistance2 / radiusmax)

    return channels
    # original ends a stopwatch here
end

lr(coordinates1::AbstractVector, coordinates2::AbstractVector, radius) = lr(BallTree(coordinates1), coordinates2, radius)
function lr(neighbortree::NNTree, coordinates2::AbstractVector, radius, roiarea)
    ineighbors = inrange(neighbortree, coordinates2, radius, true)
    kfuncans = length.(ineighbors)
    density = length.(ineighbors) ./ (π * radius ^ 2)

    Lr = ((roiarea ^ 2) * kfuncans / (length(coordinates2) - 1) / π) .^ 0.5 # numerator multiplied by the area of the ROI

    return Lr, density # original returns Lr, indexes, displacement, and density
end