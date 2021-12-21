using Gtk.ShortNames, GtkObservables, NativeFileDialog, Printf, Plots, LocalizationMicroscopy, Images, ImageIO # find out which subpackages of Images I need
using PolygonOps
using ClusDoC
using XLSX
using Statistics
using InvertedIndices

# independent observables
inputfiles = Observable([""])
outputfolder = Observable("")
localizations = Observable(Dict{String, Dict{String, Vector{Localization}}}()) # vector for each channel in each image
selectedimg = Observable{Union{Nothing, Matrix{RGB{N0f8}}}}(nothing)
rois = Observable(Dict{String, Any}())
activedrawing = Observable(false)
polyroibuilder = Observable([])
nextlineposition = Observable{Union{Nothing, NTuple{2, Float64}}}(nothing)
selectedroi = Observable{Union{Nothing, Int}}(nothing)

# initialize UI elements
b = Gtk.GtkBuilder(filename="gui/clusdoc.glade")
imgcanvas = canvas(UserUnit, 500, 500)
push!(b["canvasbox"], imgcanvas)
win = b["mainwin"]

inputbtn = button(widget = b["inputbtn"])
inputtxt = textbox(String; widget = b["inputtxt"])
outputbtn = button(widget = b["outputbtn"])
textbox(String; widget = b["outputtxt"], observable = outputfolder)
fileselector = dropdown([], widget = b["fileselector"])
addroibtn = button(widget = b["roiadd"])
deleteroibtn = button(widget = b["roidelete"])
savebtn = button(widget = b["roisave"])
loadbtn = button(widget = b["roiload"])
runbtn = button(widget = b["runbtn"])

Gtk.showall(win)


# define event handlers
function loadfiles(_)
    files = pick_multi_file()
    if length(files) > 0
        inputfiles[] = files
    end
end

function vec_to_text(obs)
    if length(obs) > 0
        txt = "\"" * join(obs, "\" \"") * "\""
    else
        txt = ""
    end
    if txt != inputtxt[]
        inputtxt[] = txt
    end
end

function text_to_vec(obs)
    files = filter!(x -> !contains(x, r"^\s*$"), split(obs, "\"", keepempty = false))
    if files != inputfiles[]
        inputfiles[] = files
    end
end

function set_outputfolder(_)
    path = pick_folder()
    if path != ""
        outputfolder[] = joinpath(path, "ClusDoC Results")
    end
end

function populate_fileselector(obs)
    @idle_add begin
        empty!(fileselector)
        append!(fileselector, basename.(obs))
        fileselector[] = basename(inputfiles[][1])
    end
end

function load_data(obs)
    obs === nothing && return
    empty!(localizations[])
    for f ∈ inputfiles[]
        locs = loadlocalizations(f, LocalizationMicroscopy.nikonelementstext)
        fname = basename(f)
        localizations[][fname] = Dict{String, Vector{Localization}}()
        for chname ∈ unique(l.channel for l ∈ locs)
            localizations[][fname][chname] = filter(l -> l.channel == chname, locs)
        end
    end
    empty!(rois[])
    notify!(localizations)
    notify!(rois)
end

const colors = (:blue, :orange, :purple)

function drawplots(_)
    outputfolder[] == "" && return
    for f ∈ inputfiles[]
        filename = basename(f)
        haskey(localizations[], filename) || continue
        locs = localizations[][filename]
        plot()
        for (i, chname) ∈ enumerate(sort(collect(keys(locs))))
            chpoints = extractcoordinates(locs[chname])
            Plots.scatter!(chpoints[1, :], chpoints[2, :], markercolor = colors[i], markeralpha = 0.5, markersize = 4, aspectratio = :equal, size=(2048, 2048), markerstrokewidth = 0)
        end
        Plots.plot!(ticks=:none, legend = :none, axis = false, widen = false, margin=-2(Plots.mm)) # change margin when Plots is updated
        imagepath = joinpath(outputfolder[], filename * ".png")
        Plots.savefig(imagepath)
        # could probably generate plots, but delay saving until an output folder selected.
    end
    selectedimg[] = load_image(fileselector[])
end

function load_image(obs)
    obs === nothing && return
    try
        # may fire before images created
        selectedimg[] = load(joinpath(outputfolder[], obs * ".png"))
    catch ex
        if !(ex isa ArgumentError)
            rethrow()
        end
    end
end

function change_file(obs)
    selectedroi[] = nothing
end

function start_roi_drawing(obs)
    activedrawing[] = true
    selectedroi[] = nothing
end

function delete_roi(obs)
    deleteat!(rois[][fileselector[]], selectedroi[])
    selectedroi[] = nothing
end

function save_rois(obs)
    outputfolder !== "" || return
    roifile = joinpath(outputfolder[], "roicoordinates.txt")

    open(roifile, "w") do f
        for (filename, filerois) ∈ rois[]
            write(f, "#", filename, "\n")
            for roi ∈ filerois
                write(f, "\t")
                for p ∈ @view roi[1:end-1]
                    write(f, string(p), ",")
                end
                write(f, string(roi[end]), "\n")
            end
            write(f, "\n")
        end
    end
end

function load_rois(obs)
    # if file doesn't exist, ask for location
    outputfolder !== "" || return
    roifile = joinpath(outputfolder[], "roicoordinates.txt")

    if !isfile(roifile)
        roifile = pick_file()
    end
    isfile(roifile) || return
    roidict = Dict{String, Any}()
    open(roifile, "r") do f
        lines = readlines(f)
        filename = ""
        filerois = []
        for line ∈ lines
            if startswith(line, "#")
                if filename != ""
                    roidict[filename] = filerois
                    filerois = []
                end
                filename = only(split(line, "#", keepempty = false))
            elseif startswith(line, "\t") && filename != ""
                points = split.(
                    split(only(split(line, ('\t'), keepempty = false)), "),("),
                    Ref(('(', ',', ')')), keepempty = false)
                roi = [tuple(parse.(Float64, y)...) for y ∈ points]
                push!(filerois, roi)
            end
        end
        filename != "" && (roidict[filename] = filerois)
    end
    rois[] = roidict
end

resultsvar = nothing

function run_clusdoc(obs)
    try
        for inputfile ∈ inputfiles[]
            filename = basename(inputfile)
            locs = localizations[][filename]
            chnames = sort(unique(keys(locs)))

            results = Vector{ClusDoC.ChannelResult}[]
            for (i, roi) ∈ enumerate(rois[][filename])
                roi = [(x, 1 - y) for (x, y) ∈ roi] # invert y to match localization coordinates - but actually I might need to invert the original image instead
                roilocalizations = Vector{Vector{Localization}}()
                for chname ∈ chnames
                    coords = extractcoordinates(locs[chname])
                    roichlocalizationsmask = inpolygon.(eachcol(coords ./ 40960), Ref(roi)) .!= 0
                    push!(roilocalizations, locs[chname][roichlocalizationsmask])
                end
                cr = clusdoc(chnames, roilocalizations, abs(PolygonOps.area(roi) * 40960 * 40960))
                push!(results, cr)
                global resultsvar = results
                generate_localization_maps(cr, filename, i, chnames)
                generate_doc_maps(cr, filename, i, chnames)
                generate_cluster_maps(cr, filename, i, chnames)
                generate_doc_histograms(cr, filename, i, chnames)
            end
            
            writeresultstables(results, joinpath(outputfolder[], "$(filename) ClusDoC Results.xlsx"))
            # should save: localization map with ROIs shown
            # next: generate plots of DoC scores etc.
        end
    catch ex
        println(ex) # why do errors get swallowed?
    end
end

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
            scatter(cr[j].coordinates[1, cr[j].abovethreshold], cr[j].coordinates[2, cr[j].abovethreshold], markerz = cr[j].docscores[k], markersize = 4, alpha = 0.5, markerstrokewidth = 0)
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
            histogram(cr[j].docscores[k], fillcolor = colors[j], bins = 50, size = (1024, 256), legend = :none,
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
                    sheet[2 + k, 1 + (j - 1) * i] = count(cr.docscores[j] .> 0.4) / length(cr.docscores[j])
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
                    colocalized_indexes = findall([count(channel.docscores[j][c.core_indices] .> 0.4) ≥ minclusterpoints for c ∈ channel.clusters])
                    union!(all_colocalized_indexes, colocalized_indexes)
                    sheet[3 + r, offset] = length(colocalized_indexes)
                    length(colocalized_indexes) > 0 || continue

                    sizes = [c.size for c ∈ channel.clusters[colocalized_indexes]]
                    meansize = mean(sizes)
                    sheet[3 + r, offset + 1] = isnan(meansize) ? "" : meansize # XLSX creates a fixable error in the output with NaN values
                    meanarea = mean(channel.clusterareas[colocalized_indexes])
                    sheet[3 + r, offset + 2] = isnan(meanarea) ? "" : meanarea
                    meancircularity = mean(channel.clustercircularities[colocalized_indexes])
                    sheet[3 + r, offset + 3] = isnan(meancircularity) ? "" : meancircularity
                    meandensity = mean([mean(channel.densities[c.core_indices]) / channel.roidensity for (i, c) ∈ enumerate(channel.clusters[colocalized_indexes])])
                    sheet[3 + r, offset + 4] = isnan(meandensity) ? "" : meandensity

                    k += 1
                end

                noncolocalized_indexes = setdiff!(findall([length(c.core_indices) ≥ minclusterpoints for c ∈ channel.clusters]), all_colocalized_indexes)
                sheet[3 + r, 1] = length(noncolocalized_indexes)
                sizes = [c.size for c ∈ channel.clusters[noncolocalized_indexes]]
                length(all_colocalized_indexes) < length(channel.clusters) || continue
                meansize = mean(sizes)
                sheet[3 + r, 2] = isnan(meansize) ? "" : meansize # XLSX creates a fixable error in the output with NaN values
                meanarea = mean(channel.clusterareas[noncolocalized_indexes])
                sheet[3 + r, 3] = isnan(meanarea) ? "" : meanarea
                meancircularity = mean(channel.clustercircularities[noncolocalized_indexes])
                sheet[3 + r, 4] = isnan(meancircularity) ? "" : meancircularity
                meandensity = mean([mean(channel.densities[c.core_indices]) / channel.roidensity for (i, c) ∈ enumerate(channel.clusters[noncolocalized_indexes])])
                sheet[3 + r, 5] = isnan(meandensity) ? "" : meandensity

                ###  In the original, the average density was determined by the rectangular range between min and max point coordinates for the area
                ### A benefit of that is that if your ROI overshoots, you will end up with the same area regardless. But it also means that non-square ROIs
                ### will underestimate average density by a lot. I should just stick with the actual area, which I'm doing now.
                
                # maybe I should add a bit of padding to the image to allow ROIs to catch the edge
            end
        end
    end
end

function draw_canvas(c, img, rois, newroi, nextlineposition, selectedroi)
    img !== nothing || return
    copy!(c, img)
    set_coordinates(c, BoundingBox(0, 1, 0, 1))
    ctx = Gtk.getgc(c)
    if haskey(rois, fileselector[])
        for (i, roi) ∈ enumerate(rois[fileselector[]])
            if i == selectedroi
                drawroi(ctx, roi, colorant"red")
            else
                drawroi(ctx, roi, colorant"gray")
            end
        end
    end
    nextpoly = nextlineposition === nothing ? newroi : [newroi; nextlineposition]
    drawroi(ctx, nextpoly, colorant"red", false)
end

function drawroi(ctx, roi, color, close = true)
    isempty(roi) && return
    p1 = first(roi)
    p = p1
    Gtk.move_to(ctx, p[1], p[2])
    Gtk.set_source(ctx, color)
    for i ∈ 2:length(roi)
        p = roi[i]
        Gtk.line_to(ctx, p[1], p[2])
    end
    close && Gtk.line_to(ctx, p[1], p[2])
    Gtk.stroke(ctx)
end

function onmouseclick(btn)
    if activedrawing[]
        if btn.button == 1 && btn.modifiers == 0
            push!(polyroibuilder[], (btn.position.x.val, btn.position.y.val))
            if btn.clicktype == DOUBLE_BUTTON_PRESS
                # two single-click events occur when a double-click event occurs
                pop!(polyroibuilder[])
                pop!(polyroibuilder[])
                push!(polyroibuilder[], polyroibuilder[][1])
                if haskey(rois[], fileselector[])
                    push!(rois[][fileselector[]], polyroibuilder[])
                else
                    rois[][fileselector[]] = [polyroibuilder[]]
                end
                polyroibuilder[] = []
                activedrawing[] = false
                nextlineposition[] = nothing
                notify!(rois)
            else
                notify!(polyroibuilder)
            end
        end
    else
        if haskey(rois[], fileselector[])
            filerois = collect(rois[][fileselector[]])
            pos = (btn.position.x.val, btn.position.y.val)
            matches = inpolygon.(Ref(pos), filerois) .!= 0
            if selectedroi[] === nothing
                selectedroi[] = findfirst(matches)
            else
                selectedroi[] = findnext(matches, selectedroi[] + 1)
                if selectedroi[] === nothing
                    selectedroi[] = findfirst(matches)
                end
            end
        end
    end
end

function onmousemove(btn)
    if activedrawing[]
        nextlineposition[] = (btn.position.x.val, btn.position.y.val)
    end
end


# hook up event handlers
on(loadfiles, inputbtn)
on(vec_to_text, inputfiles)
on(load_data, inputfiles)
on(populate_fileselector, inputfiles)
on(text_to_vec, inputtxt)
on(set_outputfolder, outputbtn)
on(drawplots, outputfolder)
on(drawplots, localizations)
on(load_image, fileselector)
on(change_file, fileselector)
on(start_roi_drawing, addroibtn)
on(delete_roi, deleteroibtn)
on(save_rois, savebtn)
on(load_rois, loadbtn)
on(run_clusdoc, runbtn)

draw(draw_canvas, imgcanvas, selectedimg, rois, polyroibuilder, nextlineposition, selectedroi)

on(onmouseclick, imgcanvas.mouse.buttonpress)
on(onmousemove, imgcanvas.mouse.motion)



# should have a colorset selector?

