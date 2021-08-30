module ClustDoC

using Clustering
using Makie
using LocalizationMicroscopy
using NearestNeighbors
using StatsBase
using XLSX

include("types/ChannelResult.jl")
include("doc.jl")
include("dbscan.jl")

path = "test/nikontestdata.bin.txt"

function clusdoc(path)
    # load file, allow selection of ROI, run algorithm with default parameters
    locs = loadlocalizations(path, LocalizationMicroscopy.nikonelementstext)
    ch1locs = getlocalizations(locs, "488", 1, 11000, 100, 10)
    ch2locs = getlocalizations(locs, "647", 11001, 11000, 100, 10)
    doc = doc([ch1locs, ch2locs])
end

function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1

    localizations = filter(l -> l.channel == channelname, alllocalizations)
end

end # module
