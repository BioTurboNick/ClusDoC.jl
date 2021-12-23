function generate_whole_localization_map(locs::Dict{String, Vector{Localization}}, outputpath, filename)
    plot()
    for (i, chname) ∈ enumerate(sort(collect(keys(locs))))
        chpoints = extractcoordinates(locs[chname])
        Plots.scatter!(chpoints[1, :], chpoints[2, :], markercolor = colors[i], markeralpha = 0.5, markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
    end
    Plots.plot!(ticks=:none, legend = :none, axis = false, widen = false, margin=-2(Plots.mm)) # change margin when Plots is updated
    imagepath = joinpath(outputpath, filename * ".png")
    Plots.savefig(imagepath)
end

function generate_localization_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames)
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
            scatter(cr[j].coordinates[1, cr[j].abovethreshold], cr[j].coordinates[2, cr[j].abovethreshold], markerz = cr[j].pointdata[!, Symbol(:docscore, k)], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
            plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
            path = joinpath(outputpath, "doc maps")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * " to " * chnames[k] * ".png")
            savefig(imagepath)
        end
    end
end

function generate_cluster_maps(cr::Vector{ChannelResult}, outputpath, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], color = :gray, markersize = 4, alpha = 0.1)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        [plot!(ai, lw = 5, linecolor = colors[j]) for ai in cr[j].clustercontours]
        path = joinpath(outputpath, "cluster maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_histograms(cr::Vector{ChannelResult}, outputpath, filename, i, chnames)
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

function writeresultstables(roiresults::Vector{Vector{ChannelResult}}, path)
    XLSX.openxlsx(path, mode = "w") do xf
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

        clusterinfo_rowlength = 5

        for (r, roichannels) ∈ enumerate(roiresults)
            for (i, channel) ∈ enumerate(roichannels)
                if r == 1
                    XLSX.addsheet!(xf)
                    sheet = xf[i + 1]
                    XLSX.rename!(sheet, "Clus-DoC Results $(channel.channelname)")
                    sheet["A1"] = "Properties of clusters by type"
                    sheet["A2"] = "Noncolocalized"
                    sheet["A3"] = "Number of clusters"
                    sheet["B3"] = "Number of localizations per cluster"
                    sheet["C3"] = "Area"
                    sheet["D3"] = "Circularity"
                    sheet["E3"] = "Relative density / granularity"
                else
                    sheet = xf[i + 1]
                end

                k = 1
                for (j, channel2) ∈ enumerate(roichannels)
                    i == j && continue
                    offset = clusterinfo_rowlength * k + 1
                    if r == 1
                        sheet[2, offset] = "Colocalized with $(channel2.channelname)"
                        sheet[3, offset] = "Number of clusters"
                        sheet[3, offset + 1] = "Number of localizations per cluster"
                        sheet[3, offset + 2] = "Area"
                        sheet[3, offset + 3] = "Circularity"
                        sheet[3, offset + 4] = "Relative density / granularity"
                    end

                    sheet[3 + r, offset + 1] = channel.ncoclusters[j]
                    sheet[3 + r, offset + 2] = replacenan(channel.meancoclustersize[j])
                    sheet[3 + r, offset + 3] = replacenan(channel.meancoclusterarea[j])
                    sheet[3 + r, offset + 4] = replacenan(channel.meancoclustercircularity[j])
                    sheet[3 + r, offset + 5] = replacenan(channel.meancoclusterdensity[j])
                    
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
        end
    end
end

