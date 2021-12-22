function generate_localization_maps(cr::Vector{ClusDoC.ChannelResult}, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], markercolor = colors[j], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        path = joinpath(outputfolder[], "localization maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_maps(cr::Vector{ClusDoC.ChannelResult}, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        for k ∈ eachindex(cr)
            k != j || continue
            scatter(cr[j].coordinates[1, cr[j].abovethreshold], cr[j].coordinates[2, cr[j].abovethreshold], markerz = cr[j].pointdata.docscore[k], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
            plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
            path = joinpath(outputfolder[], "doc maps")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * " to " * chnames[k] * ".png")
            savefig(imagepath)
        end
    end
end

function generate_cluster_maps(cr::Vector{ClusDoC.ChannelResult}, filename, i, chnames)
    xmin, xmax = minimum(minimum(c.coordinates[1,:] for c ∈ cr)), maximum(maximum(c.coordinates[1,:] for c ∈ cr))
    ymin, ymax = minimum(minimum(c.coordinates[2,:] for c ∈ cr)), maximum(maximum(c.coordinates[2,:] for c ∈ cr))
    for j ∈ eachindex(cr)
        scatter(cr[j].coordinates[1, :], cr[j].coordinates[2, :], color = :gray, markersize = 4, alpha = 0.1)
        plot!(size=(2048,2048), legend = :none, aspectratio = :equal, axis = false, ticks = false, xlims = (xmin, xmax), ylims = (ymin, ymax))
        [plot!(ai, lw = 5, linecolor = colors[j]) for ai in cr[j].clustercontours]
        path = joinpath(outputfolder[], "cluster maps")
        mkpath(path)
        imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
        savefig(imagepath)
    end
end

function generate_doc_histograms(cr::Vector{ClusDoC.ChannelResult}, filename, i, chnames)
    # what to do about huge spike at -1?
    for j ∈ eachindex(cr)
        for k ∈ eachindex(cr)
            j != k || continue
            histogram(cr[j].docscores[k], fillcolor = colors[j], bins = 100, size = (1024, 256), legend = :none,
                xlabel = "Degree of Colocalization", ylabel = "Frequency", margin = 6Plots.mm, widen = false)
            path = joinpath(outputfolder[], "doc histograms")
            mkpath(path)
            imagepath = joinpath(path, filename * " region $i " * chnames[j] * ".png")
            savefig(imagepath)
        end
    end
end


function writeresultstables(roiresults::Vector{Vector{ClusDoC.ChannelResult}}, path)
    XLSX.openxlsx(path, mode = "w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "DoC Results")
        sheet["A1"] = "Percentage of colocalized molecules"
        sheet["A2"] = ["$x -> $y" for x ∈ [cr.channelname for cr ∈ roiresults[1]] for y ∈ [cr.channelname for cr ∈ roiresults[1]] if x != y]
        for k ∈ eachindex(roiresults)
            for (i, cr) ∈ enumerate(roiresults[k])
                for j ∈ eachindex(roiresults[k])
                    i == j && continue
                    sheet[2 + k, 1 + (j - 1) * i] = count(cr.pointdata.docscore[j] .> 0.4) / length(cr.pointdata.docscore[j])
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

                minclusterpoints = 10 # clusters with fewer points than this are ignored, and clusers must have at least this many colocalized points to be considered colocalized
                all_colocalized_indexes = []
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

                    #clusterpoints = union(c.core_indices, c.po)
                    colocalized_indexes = findall([count(channel.pointdata.docscore[j][c.core_indices] .> 0.4) ≥ minclusterpoints for c ∈ channel.clusterdata.cluster])
                    union!(all_colocalized_indexes, colocalized_indexes)
                    sheet[3 + r, offset] = length(colocalized_indexes)
                    length(colocalized_indexes) > 0 || continue

                    meansize = mean(channel.clusterdata.size[colocalized_indexes])
                    sheet[3 + r, offset + 1] = isnan(meansize) ? "" : meansize # XLSX creates a fixable error in the output with NaN values
                    meanarea = mean(channel.clusterdata.area[colocalized_indexes])
                    sheet[3 + r, offset + 2] = isnan(meanarea) ? "" : meanarea
                    meancircularity = mean(channel.clusterdata.circularity[colocalized_indexes])
                    sheet[3 + r, offset + 3] = isnan(meancircularity) ? "" : meancircularity
                    meandensity = mean([mean(channel.pointdata.density[c.core_indices]) / channel.roidensity for (i, c) ∈ enumerate(channel.clusterdata.cluster[colocalized_indexes])])
                    sheet[3 + r, offset + 4] = isnan(meandensity) ? "" : meandensity

                    k += 1
                end

                noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ minclusterpoints for c ∈ channel.clusterdata.cluster]), all_colocalized_indexes)
                sheet[3 + r, 1] = length(noncolocalized_indexes)
                length(all_colocalized_indexes) < length(channel.clusters) || continue
                meansize = mean(chanel.clusterdata.size[noncolocalized_indexes])
                sheet[3 + r, 2] = isnan(meansize) ? "" : meansize # XLSX creates a fixable error in the output with NaN values
                meanarea = mean(channel.clusterdata.area[noncolocalized_indexes])
                sheet[3 + r, 3] = isnan(meanarea) ? "" : meanarea
                meancircularity = mean(channel.clusterdata.circularity[noncolocalized_indexes])
                sheet[3 + r, 4] = isnan(meancircularity) ? "" : meancircularity
                meandensity = mean([mean(channel.pointdata.density[c.core_indices]) / channel.roidensity for (i, c) ∈ enumerate(channel.clusterdata.cluster[noncolocalized_indexes])])
                sheet[3 + r, 5] = isnan(meandensity) ? "" : meandensity

                ###  In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
                ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
                ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
            end
        end
    end
end

