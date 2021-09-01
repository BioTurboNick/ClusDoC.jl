module ClustDoC

using Clustering
using Makie
using LocalizationMicroscopy
using NearestNeighbors
using StatsBase
using XLSX

include("src/types/ChannelResult.jl") # remove src
include("src/doc.jl")
include("src/dbscan.jl")

path = "test/realtest.bin.txt"

function clusdoc(path)
    # load file, allow selection of ROI, run algorithm with default parameters
    locs = loadlocalizations(path, LocalizationMicroscopy.nikonelementstext)
    ch1locs = getlocalizations(locs, "488", 1, 11000, 100, 10)
    ch2locs = getlocalizations(locs, "647", 11001, 11000, 100, 10)
    cr = doc([ch1locs, ch2locs], 20, 500, 10, 40960*40960)
    dbscan!(cr, 0.5, )
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
histogram(cr[1].docscore[2]) # contains NaNs, histogram ignores
histogram(cr[2].docscore[1])
Plots.scatter(cr[1].coordinates[1, :], cr[1].coordinates[2, :], markercolor = cgrad(:thermal)[cr[1].docscore[2]], markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
savefig("ch1doc.png")
Plots.scatter(cr[2].coordinates[1, :], cr[2].coordinates[2, :], markercolor = cgrad(:thermal)[cr[2].docscore[1]], markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
savefig("ch2doc.png")

function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1

    localizations = filter(l -> l.channel == channelname, alllocalizations)
end

end # module
