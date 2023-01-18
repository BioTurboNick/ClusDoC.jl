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

function generate_localization_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames, colors)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], markercolor = colors[j], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        add_scalebar!(xmin, xmax, ymin, ymax)
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
            abovethreshold = cr[j].pointdata.abovethreshold
            scatter(cr[j].coordinates[1, abovethreshold], cr[j].coordinates[2, abovethreshold], markerz = cr[j].pointdata[abovethreshold, Symbol(:docscore, k)], markersize = 4, markerstrokewidth = 0, alpha = 0.5, seriescolor = :balance, clims = (-1, 1), tickfontsize = 24)
            plot!(size=(2048,2048), margin = 7Plots.mm, legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
            add_scalebar!(xmin, xmax, ymin, ymax)
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
        cr[j].clusterdata !== nothing || continue
        [plot!(ai, lw = 5, linecolor = colors[j]) for ai in cr[j].clusterdata.contour]
        add_scalebar!(xmin, xmax, ymin, ymax)
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
            histogram(cr[j].pointdata[!, Symbol(:docscore, k)], fillcolor = colors[j], bins = 100, xlims = (-1, 1), size = (1024, 256), legend = :none,
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

        writeresultstables_parameters(xf, [c.channelname for c ∈ roiresults[1]], docparameters, clusterparameters)
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
        sheet["F1"] = "Number of significant clusters"
        sheet["G1"] = "Density of significant clusters (clusters / μm^2)"
        sheet["H1"] = "Significant cluster area (nm^2)"
        sheet["I1"] = "Significant cluster circularity"
        sheet["J1"] = "Number of localizations in ROI"
        sheet["K1"] = "Fraction of localizations in clusters"
        sheet["L1"] = "Number of localizations per cluster"
        sheet["M1"] = "Absolute density in clusters (molecules / μm^2)"
        sheet["N1"] = "Relative density in clusters"
        sheet["O1"] = "Fraction of localizations in significant clusters"
        sheet["P1"] = "Number of localizations per significant cluster"
        sheet["Q1"] = "Absolute density in significant clusters (molecules / μm^2)"
        sheet["R1"] = "Relative density in significant clusters"
    else
        sheet = xf[i + 1]
    end

    sheet[1 + r, 1] = channel.roiarea
    sheet[1 + r, 2] = channel.nclusters
    sheet[1 + r, 3] = channel.roiclusterdensity
    sheet[1 + r, 4] = channel.meanclusterarea
    sheet[1 + r, 5] = channel.meanclustercircularity
    sheet[1 + r, 6] = channel.nsigclusters
    sheet[1 + r, 7] = channel.roisigclusterdensity
    sheet[1 + r, 8] = channel.meansigclusterarea
    sheet[1 + r, 9] = channel.meansigclustercircularity
    sheet[1 + r, 10] = channel.nlocalizations
    sheet[1 + r, 11] = channel.fraction_clustered
    sheet[1 + r, 12] = channel.meanclustersize
    sheet[1 + r, 13] = channel.meanclusterabsolutedensity
    sheet[1 + r, 14] = channel.meanclusterdensity
    sheet[1 + r, 15] = channel.fraction_sig_clustered
    sheet[1 + r, 16] = channel.meansigclustersize
    sheet[1 + r, 17] = channel.meansigclusterabsolutedensity
    sheet[1 + r, 18] = channel.meansigclusterdensity
end

function writeresultstables_clusdoc(xf, r, roichannels, i, channel)
    clusterinfo_rowlength = 12
    sheetoffset = length(roichannels)
    # Clustering-colocalization results
    if r == 1
        XLSX.addsheet!(xf)
        sheet = xf[i + 1 + sheetoffset]
        XLSX.rename!(sheet, "Clus-DoC Results $(channel.channelname)")
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
    for (j, channel2) ∈ enumerate(roichannels)
        i == j && continue
        offset = clusterinfo_rowlength * k + 1
        if r == 1
            sheet[2, offset] = "Colocalization with $(channel2.channelname)"
            sheet[3, offset] = "Number of clusters"
            sheet[3, offset + 1] = "Number of localizations per cluster"
            sheet[3, offset + 2] = "Area (nm^2)"
            sheet[3, offset + 3] = "Circularity"
            sheet[3, offset + 4] = "Relative density"
            sheet[3, offset + 5] = "Fraction of $(channel.channelname) interacting inside clusters"

            sheet[2, offset + 6] = "Intermediate colocalization with $(channel2.channelname)"
            sheet[3, offset + 6] = "Number of clusters"
            sheet[3, offset + 7] = "Number of localizations per cluster"
            sheet[3, offset + 8] = "Area (nm^2)"
            sheet[3, offset + 9] = "Circularity"
            sheet[3, offset + 10] = "Relative density"
            sheet[3, offset + 11] = "Fraction of $(channel.channelname) interacting inside clusters"
        end

        sheet[3 + r, offset] = channel.ncoclusters[j]
        sheet[3 + r, offset + 1] = replacenan(channel.meancoclustersize[j])
        sheet[3 + r, offset + 2] = replacenan(channel.meancoclusterarea[j])
        sheet[3 + r, offset + 3] = replacenan(channel.meancoclustercircularity[j])
        sheet[3 + r, offset + 4] = replacenan(channel.meancoclusterdensity[j])
        sheet[3 + r, offset + 5] = replacenan(channel.fraction_interactions_clustered[j])

        sheet[3 + r, offset + 6] = channel.ncoclusters_int[j]
        sheet[3 + r, offset + 7] = replacenan(channel.meancoclustersize_int[j])
        sheet[3 + r, offset + 8] = replacenan(channel.meancoclusterarea_int[j])
        sheet[3 + r, offset + 9] = replacenan(channel.meancoclustercircularity_int[j])
        sheet[3 + r, offset + 10] = replacenan(channel.meancoclusterdensity_int[j])
        sheet[3 + r, offset + 11] = replacenan(channel.fraction_interactions_clustered_int[j])
        
        k += 1
    end

    sheet[3 + r, 1] = channel.ncoclusters[i]
    sheet[3 + r, 2] = replacenan(channel.meancoclustersize[i])
    sheet[3 + r, 3] = replacenan(channel.meancoclusterarea[i])
    sheet[3 + r, 4] = replacenan(channel.meancoclustercircularity[i])
    sheet[3 + r, 5] = replacenan(channel.meancoclusterdensity[i])
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
