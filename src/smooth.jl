function smooth!(result::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        coordinates = hcat(collect(pcr.coordinates for pcr ∈ result.pointschannelresults)...)
        coordinates = clusterparameters[1].uselocalradius_threshold ? coordinates[:, result.pointdata.abovethreshold] : coordinates
        area, circularity, contour = smooth!(result.clusterdata, coordinates, result.clusterresults[1], result.sigclusterresults[1], clusterparameters[1])
        result.clusterdata.area = area
        result.clusterdata.circularity = circularity
        result.clusterdata.contour = contour
    else
        areas = Float64[]
        circularities = Float64[]
        contours = []
        for i ∈ 1:result.nchannels
            coordinates = result.pointschannelresults[i].coordinates
            channelpointdata = filter(x -> x.channel == i, result.pointdata)
            coordinates = clusterparameters[i].uselocalradius_threshold ? coordinates[:, channelpointdata.abovethreshold] : coordinates
            channelclusterdata = filter(x -> x.channel == i, result.clusterdata)
            area, circularity, contour = smooth!(channelclusterdata, coordinates, result.clusterresults[i], result.sigclusterresults[i], clusterparameters[i])
            append!(areas, area)
            append!(circularities, circularity)
            append!(contours, contour)
        end
        result.clusterdata.area = areas
        result.clusterdata.circularity = circularities
        result.clusterdata.contour = contours
    end
end

function smooth!(clusterdata::DataFrame, coordinates::Matrix{Float64}, clusterresult::ClustersResult, sigclusterresult::ClustersResult, clusterparameters::ClusterParameters)
    sigmas = clusterparameters.smoothingradius
    clustercoordinates = [@view coordinates[:, union(cluster.core_indices, cluster.boundary_indices)] for cluster ∈ clusterdata.cluster]
    is2d = !any(size(clustercoordinates, 1) == 3 && (@view clustercoordinates[3, :]) .!= 0)
    if is2d
        clustercoordinates = [@view cluster[1:2, :] for cluster ∈ clustercoordinates]
    end
    # NOTE: Currently, 3d may result in OOM errors
    bounds1 = extrema.(clustercoordinates, dims = 2)
    lengths = [last.(b) .- first.(b) for b in bounds1]
    boxsizes = 0.5 .* maximum.(lengths) .+ (clusterparameters.epsilon + 10)

    centers = [(last.(b) .+ first.(b)) ./ 2 for b ∈ bounds1]
    boxmins = [cen .- bs for (cen, bs) ∈ zip(centers, boxsizes)]
    boxmaxes = [cen .+ bs for (cen, bs) ∈ zip(centers, boxsizes)]

    # create the grid
    boxes = [UnitRange.(floor.(bmin), ceil.(bmax)) for (bmin, bmax) ∈ zip(boxmins, boxmaxes)]

    areas = Vector{Float64}(undef, length(boxes))
    circularities = Vector{Float64}(undef, length(boxes))
    contours = Vector{Any}(undef, length(boxes))
    for (ii, box) ∈ enumerate(boxes)
        coords1 = clustercoordinates[ii]
        bincounts = create_histogram(coords1, box)
        clusimage = convolve_histogram(bincounts, length(box), sigmas)
        intensities = interpolate_intensities(clusimage, coords1, box)

        # My attempt to keep this dimension-agnostic breaks here; no generic way to do contour/surface? Look into MDBM.jl

        contour, contourarea = find_contour(clusimage, box, intensities)
        areas[ii] = contourarea
        circularities[ii] = calculate_circularity(contourarea, contour)
        contours[ii] = contour
    end
    
    clusterresult.meanclusterarea = mean(areas)
    clusterresult.meanclustercircularity = mean(circularities)

    sigclusters = findall(clusterdata.issignificant)
    sigclusterresult.meanclusterarea = mean(areas[sigclusters])
    sigclusterresult.meanclustercircularity = mean(circularities[sigclusters])
    return areas, circularities, contours
end

function create_histogram(coords, box)
    # create histogram of the cluster with a resolution of 1 unit
    bincounts = zeros(Int, length.(box)...) # opportunity for sparse matrix? But current base implementation only 2d
    for k ∈ axes(coords, 2)
        binindex = ()
        for j ∈ eachindex(box)
            binindex = (binindex..., findfirst(t -> t > coords[j, k], box[j]))
        end
        bincounts[binindex...] += 1
    end
    return bincounts
end

function convolve_histogram(bincounts, ndims, sigmas)
    # convolve histogram with gaussian kernel using its separability
    q = 3 # half-size of the kernel in terms of number of standard deviations (6 sigmas= 99.8% of the Normal distribution)
    nbinpoints = ceil(q * sigmas) * 2 + 1
    aux1 = LinRange(-q * sigmas, q * sigmas, nbinpoints) # x ^ 2
    Ik = exp.(-aux1 .^ 2 / (2 * sigmas ^ 2)) # exp( - (x ^ 2 + y ^ 2) / (2 * sigma ^ 2) )
    clusimage = imfilter(bincounts, reshape(Ik, length(Ik), ntuple(x -> 1, Val(ndims - 1))...))

    for d ∈ 2:ndims
        dims = (ntuple(x -> 1, Val(d - 1))..., length(Ik), ntuple(x -> 1, Val(ndims - d))...)
        imfilter!(clusimage, clusimage, reshape(Ik, dims))
    end
    return clusimage
end

function interpolate_intensities(clusimage, coords, box)
    # create interpolation of the grid to assess value at the points (real positions)
    itp = interpolate((box...,), clusimage, Gridded(Linear()))
    intensities = [itp(p...) for p ∈ eachcol(coords)]
    return intensities
end


#=
# Potential 3d code (start, anyway)
points, _ = isosurface(clusimage, MarchingCubes(iso = cutoff))
unique!(points)
offset = size(clusimage) ./ 2
midxypoints = filter(x -> x[3] == 0.0, points)
# note: points are unsorted

midxypoints = reduce(hcat, midxypoints)

midxypoints .*= offset
midxypoints .+= offset
=#
function find_contour(clusimage, box, intensities)
    # Choose the smallest contour taking all the points            
    cutoff = minimum(intensities)
    # take the biggest cluster
    clusimagesubset = @view clusimage[:, :, size(clusimage, 3) ÷ 2 + 1]
    contour = Contour.contour(box[1], box[2], clusimagesubset, cutoff)
    points = [Point.(zip(coordinates(line)...)) for line ∈ Contour.lines(contour)]
    contourarea, maxline = findmax([abs(PolygonOps.area(pts)) for pts ∈ points]) # original implementation chooses largest perimeter
    return points[maxline], contourarea
end

function calculate_circularity(contourarea, maxpoints)
    # perimeter circularity calculation
    dx = diff([first(p) for p ∈ maxpoints])
    dy = diff([last(p) for p ∈ maxpoints])
    perimeter = sum(sqrt.(dx .^ 2 + dy .^ 2))
    circularity = 4 * π * contourarea / (perimeter ^ 2)
    return circularity
end