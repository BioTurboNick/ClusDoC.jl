
# huge memory usage here, and taking a long time.
function smooth!(cr::Vector{ChannelResult}, epsilon, smoothingradius)
    for (i, cc) ∈ enumerate(cr)
        sigmas = smoothingradius
        cccoordinates = @view cc.coordinates[:, cc.pointdata.abovethreshold]
        clustercoordinates = [@view cccoordinates[:, union(cluster.core_indices, cluster.boundary_indices)] for cluster ∈ cc.clusterdata.cluster]
        is2d = any((@view clustercoordinates[3, :]) .!= 0)
        if is2d
            clustercoordinates = [@view cluster[1:2, :] for cluster ∈ clustercoordinates]
        end
        # NOTE: Currently, 3d may result in OOM errors
        bounds1 =  extrema.(clustercoordinates, dims = 2)
        lengths = [last.(b) .- first.(b) for b in bounds1]
        boxsizes = 0.5 .* maximum.(lengths) .+ (epsilon + 10)

        centers = [(last.(b) .+ first.(b)) ./ 2 for b ∈ bounds1]
        boxmins = [cen .- bs for (cen, bs) ∈ zip(centers, boxsizes)]
        boxmaxes = [cen .+ bs for (cen, bs) ∈ zip(centers, boxsizes)]

        # create the grid
        boxes = [UnitRange.(floor.(bmin), ceil.(bmax)) for (bmin, bmax) ∈ zip(boxmins, boxmaxes)]

        cc.clusterdata.area = Vector{Float64}(undef, length(boxes))
        cc.clusterdata.circularity = Vector{Float64}(undef, length(boxes))
        cc.clusterdata.contour = Vector{Any}(undef, length(boxes))
        for (ii, box) ∈ enumerate(boxes)
            # create histogram of the cluster with a resolution of 1 unit
            bincounts = zeros(Int, length.(box)...) # opportunity for sparse matrix? But current base implementation only 2d
            coords1 = clustercoordinates[ii]
            for k ∈ 1:size(coords1, 2)
                binindex = ()
                for j ∈ 1:length(box)
                    binindex = (binindex..., findfirst(t -> t > coords1[j, k], box[j]))
                end
                bincounts[binindex...] += 1
            end

            # convolve histogram with gaussian kernel using its separability
            q = 3 # half-size of the kernel in terms of number of standard deviations (6 sigmas= 99.8% of the Normal distribution)
            nbinpoints = ceil(q * sigmas) * 2 + 1
            aux1 = LinRange(-q * sigmas, q * sigmas, nbinpoints) # x ^ 2
            Ik = exp.(-aux1 .^ 2 / (2 * sigmas ^ 2)) # exp( - (x ^ 2 + y ^ 2) / (2 * sigma ^ 2) )
            clusimage = imfilter(bincounts, reshape(Ik, length(Ik), ntuple(x -> 1, Val(length(box) - 1))...))

            for d ∈ 2:length(box)
                dims = (ntuple(x -> 1, Val(d - 1))..., length(Ik), ntuple(x -> 1, Val(length(box) - d))...)
                imfilter!(clusimage, clusimage, reshape(Ik, dims))
            end

            # create interpolation of the grid to assess value at the points (real positions)
            itp = interpolate((box...,), clusimage, Gridded(Linear()))
            intensities = [itp(p...) for p ∈ eachcol(coords1)]
            
            # Choose the smallest contour taking all the points            
            cutoff = minimum(intensities)

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

            # My attempt to keep this dimension-agnostic breaks here; no generic way to do contour/surface? Look into MDBM.jl
            # take the biggest cluster
            clusimage = @view clusimage[:, :, size(clusimage, 3) ÷ 2 + 1]
            contour = Contour.contour(box[1], box[2], clusimage, cutoff)
            points = [Point.(zip(coordinates(line)...)) for line ∈ Contour.lines(contour)]
            contourarea, maxline = findmax([abs(PolygonOps.area(pts)) for pts ∈ points]) # original implementation chooses largest perimeter

            # perimeter circularity calculation
            maxpoints = points[maxline]
            dx = diff([first(p) for p ∈ maxpoints])
            dy = diff([last(p) for p ∈ maxpoints])
            perimeter = sum(sqrt.(dx .^ 2 + dy .^ 2))
            circularity = 4 * π * contourarea / (perimeter ^ 2)

            cc.clusterdata.area[ii] = contourarea
            cc.clusterdata.circularity[ii] = circularity
            cc.clusterdata.contour[ii] = points[maxline]
        end
    end
end
