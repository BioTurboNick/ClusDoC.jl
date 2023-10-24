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
        #writeresultstables_colocalization(xf, roiresults)
        
        for (r, result) ∈ enumerate(roiresults)
            # if combine_channels
            #     writeresultstables_clustering_combined(xf, r, i, channel[1])
            # else
                for i ∈ 1:result.nchannels
                    writeresultstables_clustering(xf, r, i, result)
                end
            # end
        end
        for (r, result) ∈ enumerate(roiresults)
            for i ∈ 1:result.nchannels
                #writeresultstables_clusdoc(xf, r, i, result)
            end
        end

        writeresultstables_parameters(xf, [c.channelname for c ∈ roiresults[1]], docparameters, clusterparameters)
    end
end

function writeresultstables_colocalization(xf, roiresults)
    # sheet = xf[1]
    # XLSX.rename!(sheet, "DoC Results")
    # sheet["A1"] = "Percentage of colocalized molecules"
    # sheet["A2"] = ["$x -> $y" for x ∈ [cr.channelname for cr ∈ roiresults[1]] for y ∈ [cr.channelname for cr ∈ roiresults[1]] if x != y]
    # for k ∈ eachindex(roiresults)
    #     for (i, cr) ∈ enumerate(roiresults[k])
    #         for j ∈ eachindex(roiresults[k])
    #             i == j && continue
    #             sheet[2 + k, 1 + (j - 1) * i] = cr.fraction_colocalized[j]
    #         end
    #     end
    # end
end

function writeresultstables_clustering(xf, r, i, result::ROIResult)
    # Clustering results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1]
        XLSX.rename!(sheet, "Clustering Results $(result.channelnames[i])")
        sheet["A1"] = "ROI area (μm^2)"
        sheet["B1"] = "Number of clusters"
        sheet["C1"] = "Density of clusters (clusters / μm^2)"
        sheet["D1"] = "Cluster area (nm^2)"
        sheet["E1"] = "Cluster circularity"
        sheet["F1"] = "Number of significant clusters"
        sheet["G1"] = "Density of significant clusters (clusters / μm^2)"
        sheet["H1"] = "Significant cluster area (nm^2)"
        sheet["I1"] = "Significant cluster circularity"
        sheet["J1"] = "Number of localizations in ROI"
        sheet["K1"] = "Fraction of localizations in clusters"
        sheet["L1"] = "Number of localizations per cluster"
        # sheet["M1"] = "Absolute density in clusters (molecules / μm^2)"
        # sheet["N1"] = "Relative density in clusters"
        sheet["O1"] = "Fraction of localizations in significant clusters"
        sheet["P1"] = "Number of localizations per significant cluster"
        # sheet["Q1"] = "Absolute density in significant clusters (molecules / μm^2)"
        # sheet["R1"] = "Relative density in significant clusters"
    else
        sheet = xf[i + 1]
    end

    sheet[1 + r, 1] = result.roiarea
    sheet[1 + r, 2] = result.clusterresults[i].nclusters
    sheet[1 + r, 3] = result.clusterresults[i].roiclusterdensity
    sheet[1 + r, 4] = result.clusterresults[i].meanclusterarea
    sheet[1 + r, 5] = result.clusterresults[i].meanclustercircularity
    sheet[1 + r, 6] = result.sigclusterresults[i].nclusters
    sheet[1 + r, 7] = result.sigclusterresults[i].roiclusterdensity
    sheet[1 + r, 8] = result.sigclusterresults[i].meanclusterarea
    sheet[1 + r, 9] = result.sigclusterresults[i].meanclustercircularity
    sheet[1 + r, 10] = result.pointschannelresults[i].nlocalizations
    sheet[1 + r, 11] = result.pointschannelresults[i].fraction_clustered
    # sheet[1 + r, 12] = channel.meanclustersize
    # sheet[1 + r, 13] = channel.meanclusterabsolutedensity
    #sheet[1 + r, 14] = channel.meanclusterdensity
    #sheet[1 + r, 15] = channel.fraction_sig_clustered
    #sheet[1 + r, 16] = channel.meansigclustersize
    # sheet[1 + r, 17] = channel.meansigclusterabsolutedensity
    # sheet[1 + r, 18] = channel.meansigclusterdensity
end

function writeresultstables_clusdoc(xf, r, i, result::ROIResult)
    clusterinfo_rowlength = 12
    sheetoffset = result.nchannels
    # Clustering-colocalization results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1 + sheetoffset]
        XLSX.rename!(sheet, "Clus-DoC Results $(result.channelnames[i])")
        sheet["A1"] = "Properties of significant clusters by type"
        sheet["A2"] = "Noncolocalized"
        sheet["A3"] = "Number of clusters"
        sheet["B3"] = "Number of localizations per cluster"
        sheet["C3"] = "Area (nm^2)"
        sheet["D3"] = "Circularity"
        sheet["E3"] = "Relative density"
    else
        sheet = xf[i + 1 + sheetoffset]
    end

    k = 1
    for j ∈ 1:result.nchannels
        i == j && continue
        offset = clusterinfo_rowlength * k + 1
        if r == 1
            sheet[2, offset] = "Colocalization with $(result.channelnames[i])"
            sheet[3, offset] = "Number of clusters"
            sheet[3, offset + 1] = "Number of localizations per cluster"
            sheet[3, offset + 2] = "Area (nm^2)"
            sheet[3, offset + 3] = "Circularity"
            sheet[3, offset + 4] = "Relative density"
            sheet[3, offset + 5] = "Fraction of $(result.channelnames[i]) interacting inside clusters"

            sheet[2, offset + 6] = "Intermediate colocalization with $(result.channelnames[j])"
            sheet[3, offset + 6] = "Number of clusters"
            sheet[3, offset + 7] = "Number of localizations per cluster"
            sheet[3, offset + 8] = "Area (nm^2)"
            sheet[3, offset + 9] = "Circularity"
            sheet[3, offset + 10] = "Relative density"
            sheet[3, offset + 11] = "Fraction of $(result.channelnames[i]) interacting inside clusters"
        end

        sheet[3 + r, offset] = result.coclusterresults[j].nclusters
        #sheet[3 + r, offset + 1] = replacenan(result.coclusterresults[j].meanclustersize)
        sheet[3 + r, offset + 2] = replacenan(result.coclusterresults[j].meanclusterarea)
        sheet[3 + r, offset + 3] = replacenan(result.coclusterresults[j].meanclustercircularity)
        #sheet[3 + r, offset + 4] = replacenan(result.coclusterresults[j].meancoclusterdensity[j])
        #sheet[3 + r, offset + 5] = replacenan(result.coclusterresults[j].fraction_interactions_clustered[j])

        sheet[3 + r, offset + 6] = result.intermediatecoclusterresults[j].nclusters
        #sheet[3 + r, offset + 7] = replacenan(result.intermediatecoclusterresults[j].meanclustersize)
        sheet[3 + r, offset + 8] = replacenan(result.intermediatecoclusterresults[j].meanclusterarea)
        sheet[3 + r, offset + 9] = replacenan(result.intermediatecoclusterresults[j].meanclustercircularity)
        #sheet[3 + r, offset + 10] = replacenan(result.intermediatecoclusterresults[j].meanclusterdensity)
        #sheet[3 + r, offset + 11] = replacenan(result.intermediatecoclusterresults[j].fraction_interactions_clustered)
        
        k += 1
    end

    sheet[3 + r, 1] = result.noncolocalizedclusterresults.nclusters
    #sheet[3 + r, 2] = replacenan(result.noncolocalizedclusterresults[i].meanclustersize)
    sheet[3 + r, 3] = replacenan(result.noncolocalizedclusterresults.meanclusterarea)
    sheet[3 + r, 4] = replacenan(result.noncolocalizedclusterresults.meanclustercircularity)
    #sheet[3 + r, 5] = replacenan(result.noncolocalizedclusterresults.meanclusterdensity)
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
