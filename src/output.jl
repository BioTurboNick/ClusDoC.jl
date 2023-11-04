const defaultcolors = (colorant"blue", colorant"orange", colorant"purple")

function generate_whole_localization_map(locs::Dict{String, Vector{Localization}}, outputpath, filename, colors)
    plot()
    for (i, chname) ∈ enumerate(sort(collect(keys(locs))))
        chpoints = extractcoordinates(locs[chname])
        Plots.scatter!(chpoints[1, :], chpoints[2, :], markercolor = colors[i], markeralpha = 0.5, markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
    end
    Plots.plot!(ticks=:none, legend = :none, axis = false, widen = false, margin=-2(Plots.mm)) # change margin when Plots is updated
    imagepath = joinpath(outputpath, filename * ".png")
    Plots.savefig(imagepath)
end

rect(w, h, x, y) = Shape(x .+ [0, w, w, 0, 0], y .+ [0, 0, h, h, 0])

function add_scalebar!(xmin, xmax, ymin, ymax)
    baselength = max(xmax - xmin, ymax - ymin)
    scalebar_height = baselength / 100
    scalebar_width = round(baselength / 10)
    scalebar_y = ymin - 1.5 * scalebar_height
    annotation_y = scalebar_y - 2 * scalebar_height
    plot!(rect(scalebar_width, scalebar_height, xmin, scalebar_y), c = :black)
    annotate!([(xmin, annotation_y, text(string(scalebar_width), :left, :bottom, 24))])
    ylims!(annotation_y, ymax)
end

function generate_localization_maps(result::ROIResult, outputpath, filename, i, chnames, colors)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[1,:] for c ∈ result.pointschannelresults))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[2,:] for c ∈ result.pointschannelresults))
    for j ∈ 1:result.nchannels
        coordinates = result.pointschannelresults[j].coordinates
        scatter(coordinates[1, :], coordinates[2, :], markercolor = colors[j], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        add_scalebar!(xmin, xmax, ymin, ymax)
        path = joinpath(outputpath, "localization maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_maps(result::ROIResult, outputpath, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[1,:] for c ∈ result.pointschannelresults))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[2,:] for c ∈ result.pointschannelresults))
    for j ∈ 1:result.nchannels
        for k ∈ 1:result.nchannels
            k != j || continue
            channelpointdata = filter(x -> x.channel == j, result.pointdata)
            coordinates = result.pointschannelresults[j].coordinates
            abovethreshold = channelpointdata.abovethreshold
            scatter(coordinates[1, abovethreshold], coordinates[2, abovethreshold], markerz = channelpointdata[abovethreshold, Symbol(:docscore, k)], markersize = 4, markerstrokewidth = 0, alpha = 0.5, seriescolor = :balance, clims = (-1, 1), tickfontsize = 24)
            plot!(size=(2048,2048), margin = 7Plots.mm, legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
            add_scalebar!(xmin, xmax, ymin, ymax)
            path = joinpath(outputpath, "doc maps")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * " to " * chnames[k] * ".png")
            savefig(imagepath)
        end
    end
end

function generate_cluster_maps(result::ROIResult, outputpath, filename, i, chnames, colors)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[1,:] for c ∈ result.pointschannelresults))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ result.pointschannelresults)), maximum(maximum(c.coordinates[2,:] for c ∈ result.pointschannelresults))

    for j ∈ 1:result.nchannels
        cr = result.pointschannelresults[j]
        scatter(cr.coordinates[1, :], cr.coordinates[2, :], color = :gray, markersize = 4, alpha = 0.1)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        [plot!(ai, lw = 5, linecolor = colors[j]) for ai in result.clusterdata.contour[result.clusterdata.channel .== j]]
        add_scalebar!(xmin, xmax, ymin, ymax)
        path = joinpath(outputpath, "cluster maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_histograms(result::ROIResult, outputpath, filename, i, chnames, colors)
    # what to do about huge spike at -1?
    for j ∈ 1:result.nchannels
        for k ∈ 1:result.nchannels
            j != k || continue
            channelpointdata = filter(x -> x.channel == j, result.pointdata)
            histogram(channelpointdata[!, Symbol(:docscore, k)], fillcolor = colors[j], bins = 100, xlims = (-1, 1), size = (1024, 256), legend = :none,
                xlabel = "Degree of Colocalization", ylabel = "Frequency", margin = 6Plots.mm, widen = false)
            path = joinpath(outputpath, "doc histograms")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
            savefig(imagepath)
        end
    end
end

# XLSX creates a fixable error in the output with NaN values
replacenan(data) = isnan(data) ? "" : data

function writeresultstables(roiresults::Vector{ROIResult}, docparameters, clusterparameters, path, combine_channels)
    XLSX.openxlsx(path, mode = "w") do xf
        writeresultstables_colocalization(xf, roiresults)
        
        for (r, result) ∈ enumerate(roiresults)
            if combine_channels
                writeresultstables_clustering(xf, r, 0, result)
            else
                for i ∈ 1:result.nchannels
                    writeresultstables_clustering(xf, r, i, result)
                end
            end
        end
        for (r, result) ∈ enumerate(roiresults)
            if combine_channels
                writeresultstables_clusdoc(xf, r, 0, result)
            else
                for i ∈ 1:result.nchannels
                    writeresultstables_clusdoc(xf, r, i, result)
                end
            end
        end

        writeresultstables_parameters(xf, roiresults[1].channelnames, docparameters, clusterparameters)
    end
end

function writeresultstables_colocalization(xf, roiresults)
    sheet = xf[1]
    XLSX.rename!(sheet, "DoC Results")
    sheet["A1"] = "Percentage of colocalized molecules"
    sheet["A2"] = ["$x -> $y" for x ∈ roiresults[1].channelnames for y ∈ roiresults[1].channelnames if x != y]
    for (k, roiresult) ∈ enumerate(roiresults)
        for i ∈ 1:roiresult.nchannels
            for j ∈ 1:roiresult.nchannels
                i == j && continue
                sheet[2 + k, 1 + (j - 1) * i] = roiresult.pointschannelresults[i].fraction_colocalized[j]
            end
        end
    end
end

function writeresultstables_clustering(xf, r, i, result::ROIResult)
    # Clustering results
    if r == 1
        XLSX.addsheet!(xf)
        if i == 0 # combined
            sheet = xf[2]
            XLSX.rename!(sheet, "Clustering Results")
        else
            sheet = xf[i + 1]
            XLSX.rename!(sheet, "Clustering Results $(result.channelnames[i])")
        end
        sheet["A1"] = "ROI area (μm^2)"
        sheet["B1"] = "Number of localizations in ROI"
        sheet["C1"] = "Number of clusters"
        sheet["D1"] = "Density of clusters (clusters / μm^2)"
        sheet["E1"] = "Cluster area (nm^2)"
        sheet["F1"] = "Cluster circularity"

        col = 7
        for j ∈ 1:result.nchannels
            chname = result.channelnames[j]
            sheet[1, col] = "$chname Fraction of localizations in clusters"
            sheet[1, col + 1] = "$chname Number of localizations per cluster"
            sheet[1, col + 2] = "$chname Absolute density in clusters (molecules / μm^2)"
            sheet[1, col + 3] = "$chname Relative density in clusters"
            col += 4
        end
        
        sheet[1, col] = "Number of significant clusters"
        sheet[1, col + 1] = "Density of significant clusters (clusters / μm^2)"
        sheet[1, col + 2] = "Significant cluster area (nm^2)"
        sheet[1, col + 3] = "Significant cluster circularity"
        col += 4

        for j ∈ 1:result.nchannels
            chname = result.channelnames[j]
            sheet[1, col] = "$chname Fraction of localizations in significant clusters"
            sheet[1, col + 1] = "$chname Number of localizations per significant cluster"
            sheet[1, col + 2] = "$chname Absolute density in significant clusters (molecules / μm^2)"
            sheet[1, col + 3] = "$chname Relative density in significant clusters"
            col += 4
        end
    else
        if i == 0 # combined
            sheet = xf[2]
        else
            sheet = xf[i + 1]
        end
    end

    i = i == 0 ? 1 : i # combined

    sheet[1 + r, 1] = result.roiarea
    sheet[1 + r, 2] = result.pointschannelresults[i].nlocalizations_abovethreshold

    clusterresult = result.clusterresults[i]
    sheet[1 + r, 3] = clusterresult.nclusters
    sheet[1 + r, 4] = clusterresult.roiclusterdensity
    sheet[1 + r, 5] = clusterresult.meanclusterarea
    sheet[1 + r, 6] = clusterresult.meanclustercircularity    

    col = 7
    for j ∈ 1:result.nchannels
        chresult = clusterresult.channelresults[j]
        sheet[1 + r, col] = chresult.fraction_clustered
        sheet[1 + r, col + 1] = chresult.meanclustersize
        sheet[1 + r, col + 2] = chresult.meanclusterabsolutedensity
        sheet[1 + r, col + 3] = chresult.meanclusterdensity
        col += 4
    end

    sigclusterresult = result.sigclusterresults[i]
    sheet[1 + r, col] = sigclusterresult.nclusters
    sheet[1 + r, col + 1] = sigclusterresult.roiclusterdensity
    sheet[1 + r, col + 2] = sigclusterresult.meanclusterarea
    sheet[1 + r, col + 3] = sigclusterresult.meanclustercircularity
    col += 4

    for j ∈ 1:result.nchannels
        chresult = sigclusterresult.channelresults[j]
        sheet[1 + r, col] = chresult.fraction_clustered
        sheet[1 + r, col + 1] = chresult.meanclustersize
        sheet[1 + r, col + 2] = chresult.meanclusterabsolutedensity
        sheet[1 + r, col + 3] = chresult.meanclusterdensity
        col += 4
    end
end

function writeresultstables_clusdoc(xf, r, i, result::ROIResult)
    clusterinfo_rowlength = 6 + 2 * 5 * result.nchannels
    sheetoffset = result.nchannels
    # Clustering-colocalization results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1 + sheetoffset]
        XLSX.rename!(sheet, "Clus-DoC Results $(result.channelnames[i])")
        sheet["A1"] = "Properties of significant clusters by type"
        sheet["A2"] = "Noncolocalized"
        sheet["A3"] = "Number of clusters"
        sheet["B3"] = "Area (nm^2)"
        sheet["C3"] = "Circularity"
        col = 4

        for j ∈ 1:result.nchannels
            chname = result.channelnames[j]
            sheet[3, col] = "$chname Number of localizations per cluster"
            sheet[3, col + 1] = "$chname Absolute density (molecules / μm^2)"
            sheet[3, col + 2] = "$chname Relative density"
            sheet[3, col + 3] = "$chname Fraction of localizations in clusters"
            sheet[3, col + 4] = "Fraction of $chname interacting inside clusters"
            col += 5
        end
    else
        sheet = xf[i + 1 + sheetoffset]
    end

    k = 1
    for j ∈ 1:result.nchannels
        i == j && continue
        baseoffset = clusterinfo_rowlength * k + 1
        if r == 1
            offset = baseoffset
            sheet[2, offset] = "Colocalization with $(result.channelnames[j])"
            sheet[3, offset] = "Number of clusters"
            sheet[3, offset + 1] = "Area (nm^2)"
            sheet[3, offset + 2] = "Circularity"
            offset += 3

            for jj ∈ 1:result.nchannels
                chname = result.channelnames[jj]
                sheet[3, offset] = "$chname Number of localizations per cluster"
                sheet[3, offset + 1] = "$chname Absolute density (molecules / μm^2)"
                sheet[3, offset + 2] = "$chname Relative density"
                sheet[3, offset + 3] = "$chname Fraction of localizations in clusters"
                sheet[3, offset + 4] = "Fraction of $chname interacting inside clusters"
                offset += 5
            end

            sheet[2, offset] = "Intermediate Colocalization with $(result.channelnames[j])"
            sheet[3, offset] = "Number of clusters"
            sheet[3, offset + 1] = "Area (nm^2)"
            sheet[3, offset + 2] = "Circularity"
            offset += 3

            for jj ∈ 1:result.nchannels
                chname = result.channelnames[jj]
                sheet[3, offset] = "$chname Number of localizations per cluster"
                sheet[3, offset + 1] = "$chname Absolute density (molecules / μm^2)"
                sheet[3, offset + 2] = "$chname Relative density"
                sheet[3, offset + 3] = "$chname Fraction of localizations in clusters"
                sheet[3, offset + 4] = "Fraction of $chname interacting inside clusters"
                offset += 5
            end
        end

        offset = baseoffset
        coclusterresult = result.coclusterresults[i][j]
        sheet[3 + r, offset] = coclusterresult.nclusters
        sheet[3 + r, offset + 1] = coclusterresult.meanclusterarea
        sheet[3 + r, offset + 2] = coclusterresult.meanclustercircularity
        offset += 3

        for jj ∈ 1:result.nchannels
            chresult = coclusterresult.channelresults[jj]
            sheet[3 + r, offset] = chresult.meanclustersize
            sheet[3 + r, offset + 1] = chresult.meanclusterabsolutedensity
            sheet[3 + r, offset + 2] = chresult.meanclusterdensity
            sheet[3 + r, offset + 3] = chresult.fraction_clustered
            sheet[3 + r, offset + 4] = chresult.fraction_of_interacting_points
            offset += 5
        end

        intcoclusterresult = result.intermediatecoclusterresults[i][j]
        sheet[3 + r, offset] = intcoclusterresult.nclusters
        sheet[3 + r, offset + 1] = intcoclusterresult.meanclusterarea
        sheet[3 + r, offset + 2] = intcoclusterresult.meanclustercircularity
        offset += 3

        for jj ∈ 1:result.nchannels
            chresult = intcoclusterresult.channelresults[jj]
            sheet[3 + r, offset] = chresult.meanclustersize
            sheet[3 + r, offset + 1] = chresult.meanclusterabsolutedensity
            sheet[3 + r, offset + 2] = chresult.meanclusterdensity
            sheet[3 + r, offset + 3] = chresult.fraction_clustered
            sheet[3 + r, offset + 4] = chresult.fraction_of_interacting_points
            offset += 5
        end

        k += 1
    end

    noncoclusterresult = result.noncolocalizedclusterresults[i]
    sheet[3 + r, 1] = noncoclusterresult.nclusters
    sheet[3 + r, 2] = replacenan(noncoclusterresult.meanclusterarea)
    sheet[3 + r, 3] = replacenan(noncoclusterresult.meanclustercircularity)
    offset = 4

    for jj ∈ 1:result.nchannels
        chresult = noncoclusterresult.channelresults[jj]
        sheet[3 + r, offset] = replacenan(chresult.meanclustersize)
        sheet[3 + r, offset + 1] = chresult.meanclusterabsolutedensity
        sheet[3 + r, offset + 2] = chresult.meanclusterdensity
        sheet[3 + r, offset + 3] = chresult.fraction_clustered
        sheet[3 + r, offset + 4] = chresult.fraction_of_interacting_points
        offset += 5
    end
end

function writeresultstables_parameters(xf, chnames, docparameters, clusterparameters)
    sheet = XLSX.addsheet!(xf)
    XLSX.rename!(sheet, "Algorithm parameters")
    sheet["A1"] = "DoC parameters"
    sheet["A2"] = "Local radius (nm)"
    sheet["B2"] = "Radius max (nm)"
    sheet["C2"] = "Radius step (nm)"
    sheet["D2"] = "Colocalized DoC threshold"

    sheet["A3"] = docparameters.localradius
    sheet["B3"] = docparameters.radiusmax
    sheet["C3"] = docparameters.radiusstep
    sheet["D3"] = docparameters.colocalized_threshold

    sheet["A5"] = "Clustering parameters"
    sheet["A6"] = "Channel"
    sheet["B6"] = "Epsilon"
    sheet["C6"] = "Min points"
    sheet["D6"] = "Use local radius threshold"
    sheet["E6"] = "Smoothing radius (nm)"
    sheet["F6"] = "Min significant cluster points"

    for (i, c) ∈ enumerate(clusterparameters)
        i > length(chnames) && break
        sheet[6 + i, 1] = chnames[i]
        sheet[6 + i, 2] = c.epsilon
        sheet[6 + i, 3] = c.minpoints
        sheet[6 + i, 4] = c.uselocalradius_threshold
        sheet[6 + i, 5] = c.smoothingradius
        sheet[6 + i, 6] = c.minsigclusterpoints
    end
end
