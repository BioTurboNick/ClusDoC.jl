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

function generate_localization_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames, colors)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], markercolor = colors[j], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        path = joinpath(outputpath, "localization maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        for k ∈ eachindex(cr)
            k != j || continue
            scatter(cr[j].coordinates[1, cr[j].pointdata.abovethreshold], cr[j].coordinates[2, cr[j].pointdata.abovethreshold], markerz = cr[j].pointdata[!, Symbol(:docscore, k)], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
            plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
            path = joinpath(outputpath, "doc maps")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * " to " * chnames[k] * ".png")
            savefig(imagepath)
        end
    end
end

function generate_cluster_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames, colors)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], color = :gray, markersize = 4, alpha = 0.1)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        [plot!(ai, lw = 5, linecolor = colors[j]) for ai in cr[j].clusterdata.contour]
        path = joinpath(outputpath, "cluster maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_histograms(cr::Vector{ChannelResult}, outputpath, filename, i, chnames, colors)
    # what to do about huge spike at -1?
    for j ∈ eachindex(cr)
        for k ∈ eachindex(cr)
            j != k || continue
            histogram(cr[j].pointdata[!, Symbol(:docscore, k)], fillcolor = colors[j], bins = 100, size = (1024, 256), legend = :none,
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

function writeresultstables(roiresults::Vector{Vector{ChannelResult}}, docparameters, clusterparameters, path)
    XLSX.openxlsx(path, mode = "w") do xf
        writeresultstables_colocalization(xf, roiresults)
        
        for (r, roichannels) ∈ enumerate(roiresults)
            for (i, channel) ∈ enumerate(roichannels)
                writeresultstables_clustering(xf, r, i, channel)
            end
        end
        for (r, roichannels) ∈ enumerate(roiresults)
            for (i, channel) ∈ enumerate(roichannels)
                writeresultstables_clusdoc(xf, r, roichannels, i, channel)
            end
        end

        writeresultstables_parameters(xf, docparameters, clusterparameters)
    end
end

function writeresultstables_colocalization(xf, roiresults)
    sheet = xf[1]
    XLSX.rename!(sheet, "DoC Results")
    sheet["A1"] = "Percentage of colocalized molecules"
    sheet["A2"] = ["$x -> $y" for x ∈ [cr.channelname for cr ∈ roiresults[1]] for y ∈ [cr.channelname for cr ∈ roiresults[1]] if x != y]
    for k ∈ eachindex(roiresults)
        for (i, cr) ∈ enumerate(roiresults[k])
            for j ∈ eachindex(roiresults[k])
                i == j && continue
                sheet[2 + k, 1 + (j - 1) * i] = cr.fraction_colocalized[j]
            end
        end
    end
end

function writeresultstables_clustering(xf, r, i, channel)
    # Clustering results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1]
        XLSX.rename!(sheet, "Clustering Results $(channel.channelname)")
        sheet["A1"] = "ROI area (μm^2)"
        sheet["B1"] = "Number of clusters"
        sheet["C1"] = "Density of clusters (clusters / μm^2)"
        sheet["D1"] = "Cluster area (nm^2)"
        sheet["E1"] = "Cluster circularity"
        sheet["F1"] = "Number of localizations in ROI"
        sheet["G1"] = "Fraction of localizations in clusters"
        sheet["H1"] = "Number of localizations per cluster"
        sheet["I1"] = "Absolute density in clusters (molecules / μm^2)"
        sheet["J1"] = "Relative density in clusters"
    else
        sheet = xf[i + 1]
    end

    sheet[1 + r, 1] = channel.roiarea
    sheet[1 + r, 2] = channel.nclusters
    sheet[1 + r, 3] = channel.roiclusterdensity
    sheet[1 + r, 4] = channel.meanclusterarea
    sheet[1 + r, 5] = channel.meanclustercircularity
    sheet[1 + r, 6] = channel.nlocalizations
    sheet[1 + r, 7] = channel.fraction_clustered
    sheet[1 + r, 8] = channel.meanclustersize
    sheet[1 + r, 9] = channel.meanclusterabsolutedensity
    sheet[1 + r, 10] = channel.meanclusterdensity
end

function writeresultstables_clusdoc(xf, r, roichannels, i, channel)
    clusterinfo_rowlength = 5
    sheetoffset = length(roichannels)
    # Clustering-colocalization results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1 + sheetoffset]
        XLSX.rename!(sheet, "Clus-DoC Results $(channel.channelname)")
        sheet["A1"] = "Properties of clusters by type"
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
    for (j, channel2) ∈ enumerate(roichannels)
        i == j && continue
        offset = clusterinfo_rowlength * k + 1
        if r == 1
            sheet[2, offset] = "Colocalized with $(channel2.channelname)"
            sheet[3, offset] = "Number of clusters"
            sheet[3, offset + 1] = "Number of localizations per cluster"
            sheet[3, offset + 2] = "Area (nm^2)"
            sheet[3, offset + 3] = "Circularity"
            sheet[3, offset + 4] = "Relative density"
        end

        sheet[3 + r, offset] = channel.ncoclusters[j]
        sheet[3 + r, offset + 1] = replacenan(channel.meancoclustersize[j])
        sheet[3 + r, offset + 2] = replacenan(channel.meancoclusterarea[j])
        sheet[3 + r, offset + 3] = replacenan(channel.meancoclustercircularity[j])
        sheet[3 + r, offset + 4] = replacenan(channel.meancoclusterdensity[j])
        
        k += 1
    end

    sheet[3 + r, 1] = channel.ncoclusters[i]
    sheet[3 + r, 2] = replacenan(channel.meancoclustersize[i])
    sheet[3 + r, 3] = replacenan(channel.meancoclusterarea[i])
    sheet[3 + r, 4] = replacenan(channel.meancoclustercircularity[i])
    sheet[3 + r, 5] = replacenan(channel.meancoclusterdensity[i])

    ###  In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
    ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
    ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
end

function writeresultstables_parameters(xf, docparameters, clusterparameters)
    # # Clustering results
    # if r == 1
    #     XLSX.addsheet!(xf)
    #     sheet = xf[end]
    #     XLSX.rename!(sheet, "Parameters $(channel.channelname)")
    #     sheet["A1"] = "ROI area (μm^2)"
    #     sheet["A2"] = "Number of clusters"
    #     sheet["A3"] = "Density of clusters (clusters / μm^2)"
    #     sheet["A4"] = "Cluster area (nm^2)"
    #     sheet["E1"] = "Cluster circularity"
    #     sheet["F1"] = "Number of localizations in ROI"
    #     sheet["G1"] = "Fraction of localizations in clusters"
    #     sheet["H1"] = "Number of localizations per cluster"
    #     sheet["I1"] = "Absolute density in clusters (molecules / μm^2)"
    #     sheet["J1"] = "Relative density in clusters"
    # else
    #     sheet = xf[i + 1]
    # end

    # sheet[1 + r, 1] = channel.roiarea
    # sheet[1 + r, 2] = channel.nclusters
    # sheet[1 + r, 3] = channel.roiclusterdensity
    # sheet[1 + r, 4] = channel.meanclusterarea
    # sheet[1 + r, 5] = channel.meanclustercircularity
    # sheet[1 + r, 6] = channel.nlocalizations
    # sheet[1 + r, 7] = channel.fraction_clustered
    # sheet[1 + r, 8] = channel.meanclustersize
    # sheet[1 + r, 9] = channel.meanclusterabsolutedensity
    # sheet[1 + r, 10] = channel.meanclusterdensity
end
