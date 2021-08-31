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
    docscores = doc([ch1locs, ch2locs], 20, 500, 10, 40960*40960)
end

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
