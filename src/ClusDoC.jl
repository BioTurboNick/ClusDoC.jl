module ClusDoC

using Clustering
using Contour
using ImageFiltering
using Interpolations
using InvertedIndices
using LocalizationMicroscopy
using Makie
using WGLMakie
using NearestNeighbors
#using Plots
using PolygonOps
using StatsBase
using XLSX

include("src/types/ChannelResult.jl") # remove src
include("src/doc.jl")
include("src/dbscan.jl")
include("src/smooth.jl")

path = "test/realtest.bin.txt"
outputpath = "output"

function clusdoc(path, outputpath)
    # load file, allow selection of ROI, run algorithm with default parameters
    locs = loadlocalizations(path, LocalizationMicroscopy.nikonelementstext)
    ch1locs = getlocalizations(locs, "488", 1, 11000, 100, 10)
    ch2locs = getlocalizations(locs, "647", 11001, 11000, 100, 10)

    p = Makie.scatter(extractcoordinates(ch1locs)[1, :], extractcoordinates(ch1locs)[2, :], color = :green, markersize = 2)
    Makie.scatter!(extractcoordinates(ch2locs)[1, :], extractcoordinates(ch2locs)[2, :], color = :red, markersize = 2)
    clicks = Node(Point2f0[(0, 0)])
    on(events(p).mousebutton, priority = 0) do event
        if event.action == Mouse.press && event.button == Mouse.left
            pos = mouseposition(event)
            push!(clicks, push!(clicks[], pos))
        end
    end
    Makie.scatter!(clicks, color = :black, marker = '+', markersize = 20, show_axis = false)

    # on(events(p).mousebutton, priority = 0) do buttons
    #     if ispressed(p, Mouse.left)
    #         pos = AbstractPlotting.mouseposition(p)
    #         push!(clicks, push!(clicks[], pos))
    #     end
    #     return
    # end
    # scatter!(p, clicks, color = :red, marker = '+', markersize = 20, show_axis = false)
    # RecordEvents(p, "roi_drawing")
    # p

    cr = doc(["488", "647"], [ch1locs, ch2locs], 20, 500, 10, 40960*40960)
    dbscan!(cr, 20, 3, true, 20)
    smooth!(cr, 20, 15)
    mkpath(outputpath)
    writeresultstables(cr, joinpath(outputpath, "ClusDoC Results.xlsx"))
    histogram(cr[1].docscores[2]) # contains NaNs, histogram ignores
    savefig(joinpath(outputpath, "ch1dochist.png"))
    histogram(cr[2].docscores[1])
    savefig(joinpath(outputpath, "ch2dochist.png"))
    Plots.scatter(cr[1].coordinates[1, :], cr[1].coordinates[2, :], markercolor = cgrad(:thermal)[cr[1].docscores[2]], markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
    savefig(joinpath(outputpath, "ch1doc.png"))
    Plots.scatter(cr[2].coordinates[1, :], cr[2].coordinates[2, :], markercolor = cgrad(:thermal)[cr[2].docscores[1]], markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
    savefig(joinpath(outputpath, "ch2doc.png"))
    return cr
end

#=
using WGLMakie
fig = Figure()

ax = Axis(fig[1, 1])
fig[2, 1] = buttongrid = GridLayout(tellwidth = false)

counts = Node([1, 4, 3, 7, 2])

buttonlabels = [@lift("Count: $($counts[i])") for i in 1:5]

buttons = buttongrid[1, 1:5] = [Button(fig, label = l) for l in buttonlabels]

for i in 1:5
    on(buttons[i].clicks) do n
        counts[][i] += 1
        notify(counts)
    end
end

barplot!(counts, color = cgrad(:Spectral)[LinRange(0, 1, 5)])
ylims!(ax, 0, 20)

fig
=#

#=
julia> @time clusdoc(path)
134.263914 seconds (15.54 M allocations: 6.483 GiB, 1.32% gc time)
=#

# plot


function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1 # activate these limits

    localizations = filter(l -> l.channel == channelname, alllocalizations)
end

function writeresultstables(channels::Vector{ChannelResult}, path)
    XLSX.openxlsx(path, mode = "w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "DoC Results")
        sheet["A1"] = "Percentage of colocalized molecules"
        sheet["A2"] = "from\\to"
        sheet["A3", dim = 1] = [cr.channelname for cr ∈ channels]
        sheet["B2"] = [cr.channelname for cr ∈ channels]
        for (i, cr) ∈ enumerate(channels)
            for j ∈ eachindex(channels)
                i == j && continue
                sheet[2 + j, 1 + i] = count(cr.docscores[j] .> 0.4) / length(cr.docscores[j])
            end
        end

        clusterinfo_rowlength = 5

        for (i, channel) ∈ enumerate(channels)
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

            all_colocalized_indexes = []
            k = 1
            for (j, channel2) ∈ enumerate(channels)
                i == j && continue
                offset = clusterinfo_rowlength * k + 1
                sheet[2, offset] = "Colocalized with $(channel2.channelname)"
                sheet[3, offset] = "Number of clusters"
                sheet[3, offset + 1] = "Number of localizations per cluster"
                sheet[3, offset + 2] = "Area"
                sheet[3, offset + 3] = "Circularity"
                sheet[3, offset + 4] = "Relative density / granularity"

                #clusterpoints = union(c.core_indices, c.po)
                colocalized_indexes = findall([count(channel.docscores[j][c.core_indices] .> 0.4) > 5 for c ∈ channel.clusters])
                union!(all_colocalized_indexes, colocalized_indexes)
                sheet[4, offset] = length(colocalized_indexes)
                sizes = [c.size for c ∈ channel.clusters[colocalized_indexes]]
                meansize = mean(sizes)
                sheet[4, offset + 1] = isnan(meansize) ? "" : meansize
                sheet[4, offset + 2] = mean(channel.clusterareas[colocalized_indexes])
                sheet[4, offset + 3] = mean(channel.clustercircularities[colocalized_indexes])
                #meandensities = [chann.densities[] for c ∈ channel.clusters[colocalized_indexes]]
                #sheet[3, offset + 4] = next step is how to get relative densities back...

                k += 1
            end

            noncolocalized_indexes = Not(all_colocalized_indexes)            
            sheet["A4"] = length(channel.clusters) - length(all_colocalized_indexes)
            sizes = [c.size for c ∈ channel.clusters[noncolocalized_indexes]]
            meansize = mean(sizes)
            sheet["B4"] = isnan(meansize) ? "" : meansize
            sheet[4, 3] = mean(channel.clusterareas[noncolocalized_indexes])
            sheet[4, 4] = mean(channel.clustercircularities[noncolocalized_indexes])
            # # sheet[3, offset + 4] = 
        end
    end
end

end # module
