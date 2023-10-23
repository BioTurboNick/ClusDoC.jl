function smooth!(result::ROIResult, clusterparameters, combinechannels)
    if combinechannels
        coordinates = hcat(pcr.coordinates for pcr ∈ result.pointschannelresults)
        coordinates = clusterparameters.uselocalradius_threshold ? coordinates[:, pointdata.abovethreshold] : coordinates
        smooth!(result.clusterdata, coordinates, result.clusterresults[1], result.sigclusterresults[1], clusterparameters[1])
    else
        areas = Float64[]
        circularities = Float64[]
        contours = []
        for i ∈ 1:result.nchannels
            coordinates = result.pointschannelresults[i].coordinates
            coordinates = clusterparameters.uselocalradius_threshold ? coordinates[:, pointdata.abovethreshold[pointdata.channel .== i]] : coordinates
            channelclusterdata = filter(x -> x.channel == i, result.clusterdata)
            area, circularity, contour = smooth!(channelclusterdata, coordinates, result.clusterresults[i], result.sigclusterresults[i], clusterparameters[i])
            append!(areas, area)
            append!(circularities, circularity)
            append!(contours, contour)
        end
        clusterdata.area = areas
        clusterdata.circularity = circularities
        clusterdata.contour = contours
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

    area = Vector{Float64}(undef, length(boxes))
    circularity = Vector{Float64}(undef, length(boxes))
    contour = Vector{Any}(undef, length(boxes))
    for (ii, box) ∈ enumerate(boxes)
        coords1 = clustercoordinates[ii]
        bincounts = create_histogram(coords1, box)
        clusimage = convolve_histogram(bincounts, length(box), sigmas)
        intensities = interpolate_intensities(clusimage, coords1, box)

        # My attempt to keep this dimension-agnostic breaks here; no generic way to do contour/surface? Look into MDBM.jl

        contour, contourarea = find_contour(clusimage, box, intensities)
        area[ii] = contourarea
        circularity[ii] = calculate_circularity(contourarea, contour)
        contour[ii] = contour
    end
    
    clusterresult.meanclusterarea = mean(area)
    clusterresult.meanclustercircularity = mean(circularity)

    sigclusters = findall(x -> x.issignificant, clusterdata)
    sigclusterresult.meansigclusterarea = mean(area[sigclusters])
    sigclusterresult.meansigclustercircularity = mean(circularity[sigclusters])
    return area, circularity, contour
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