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
colors = Observable{NTuple{3, Colorant}}(defaultcolors)
docparameters = Observable{DoCParameters}(defaultdocparameters)
clusterparameters = Observable{Vector{ClusterParameters}}(fill(defaultclusterparameters, 3))
isrunning = Observable(false)

# initialize UI elements
b = Gtk.GtkBuilder(filename=gladepath)
imgcanvas = canvas(UserUnit, 500, 500)
push!(b["canvasbox"], imgcanvas)
win = b["mainwin"]
Gtk.showall(win)

inputbtn = button(widget = b["inputbtn"])
inputtxt = textbox(String; widget = b["inputtxt"], gtksignal = "focus-out-event")
outputbtn = button(widget = b["outputbtn"])
textbox(String; widget = b["outputtxt"], observable = outputfolder, gtksignal = "focus-out-event")
fileselector = dropdown([], widget = b["fileselector"])
addroibtn = button(widget = b["roiadd"])
deleteroibtn = button(widget = b["roidelete"])
savebtn = button(widget = b["roisave"])
loadbtn = button(widget = b["roiload"])
settingsbtn = button(widget = b["opensettings"])
runbtn = button(widget = b["runbtn"])
ch1label = label("Ch 1", widget = b["ch1label"])
ch2label = label("Ch 2", widget = b["ch2label"])
ch3label = label("Ch 3", widget = b["ch3label"])

ch1colorbtn = colorbutton(defaultcolors[1]; widget = b["ch1colorbutton"])
ch2colorbtn = colorbutton(defaultcolors[2]; widget = b["ch2colorbutton"])
ch3colorbtn = colorbutton(defaultcolors[3]; widget = b["ch3colorbutton"])

localradiustxt = textbox(defaultdocparameters.localradius; widget = b["localradiusinput"], gtksignal = "focus-out-event")
radiusmaxtxt = textbox(defaultdocparameters.radiusmax; widget = b["radiusmaxinput"], gtksignal = "focus-out-event")
radiussteptxt = textbox(defaultdocparameters.radiusstep; widget = b["radiusstepinput"], gtksignal = "focus-out-event")
colocthreshtxt = textbox(defaultdocparameters.colocalized_threshold; widget = b["colocalizedthreshinput"], gtksignal = "focus-out-event")
epsilontxt1 = textbox(defaultclusterparameters.epsilon; widget = b["epsiloninput"], gtksignal = "focus-out-event")
minpointstxt1 = textbox(defaultclusterparameters.minpoints; widget = b["minpointsinput"], gtksignal = "focus-out-event")
usethreshholdcheckbox1 = checkbox(defaultclusterparameters.uselocalradius_threshold, widget = b["usethresholdinput"])
smoothingradiustxt1 = textbox(defaultclusterparameters.smoothingradius; widget = b["smoothingradiusinput"], gtksignal = "focus-out-event")
minsigclusterpointstxt1 = textbox(defaultclusterparameters.minsigclusterpoints; widget = b["minsigclusterpointsinput"], gtksignal = "focus-out-event")
epsilontxt2 = textbox(defaultclusterparameters.epsilon; widget = b["epsiloninput1"], gtksignal = "focus-out-event")
minpointstxt2 = textbox(defaultclusterparameters.minpoints; widget = b["minpointsinput1"], gtksignal = "focus-out-event")
usethreshholdcheckbox2 = checkbox(defaultclusterparameters.uselocalradius_threshold, widget = b["usethresholdinput1"])
smoothingradiustxt2 = textbox(defaultclusterparameters.smoothingradius; widget = b["smoothingradiusinput1"], gtksignal = "focus-out-event")
minsigclusterpointstxt2 = textbox(defaultclusterparameters.minsigclusterpoints; widget = b["minsigclusterpointsinput1"], gtksignal = "focus-out-event")
epsilontxt3 = textbox(defaultclusterparameters.epsilon; widget = b["epsiloninput2"], gtksignal = "focus-out-event")
minpointstxt3 = textbox(defaultclusterparameters.minpoints; widget = b["minpointsinput2"], gtksignal = "focus-out-event")
usethreshholdcheckbox3 = checkbox(defaultclusterparameters.uselocalradius_threshold, widget = b["usethresholdinput2"])
smoothingradiustxt3 = textbox(defaultclusterparameters.smoothingradius; widget = b["smoothingradiusinput2"], gtksignal = "focus-out-event")
minsigclusterpointstxt3 = textbox(defaultclusterparameters.minsigclusterpoints; widget = b["minsigclusterpointsinput2"], gtksignal = "focus-out-event")
settingsok = button(widget = b["settingsok"])
settingsreset = button(widget = b["settingsreset"])

outputcsvclustercheckbox = checkbox(defaultoutputcsvclusters; widget = b["outputcsvcluster"])


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
        outputpath = joinpath(path, "Clus-DoC Results")
        mkpath(outputpath)
        outputfolder[] = outputpath
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
    notify(localizations)
    notify(rois)
end

function drawplots(_)
    outputfolder[] == "" && return
    @idle_add Gtk.start(b["runspinner"])
    draw_task = Threads.@spawn begin
        for f ∈ inputfiles[]
            filename = basename(f)
            haskey(localizations[], filename) || continue
            locs = localizations[][filename]
            generate_whole_localization_map(locs, outputfolder[], filename, colors[])
            # could probably generate plots, but delay saving until an output folder selected.
        end
    end
    wait(draw_task)
    load_image(fileselector[])
    @idle_add Gtk.stop(b["runspinner"])
end

function load_image(obs)
    if obs === nothing
        visible(b["ch1box"], false)
        visible(b["ch2box"], false)
        visible(b["ch3box"], false)
        return
    end
    try
        # may fire before images created
        selectedimg[] = load(joinpath(outputfolder[], obs * ".png"))

        chnames = sort(unique(keys(localizations[][fileselector[]])))
        ch1label[] = chnames[1]
        visible(b["ch1box"], true)
        if length(chnames) > 1
            ch2label[] = chnames[2]
            visible(b["ch2box"], true)
        else
            visible(b["ch2box"], false)
        end
        if length(chnames) > 2
            ch3label[] = chnames[3]
            visible(b["ch3box"], true)
        else
            visible(b["ch3box"], false)
        end
    catch ex
        if !(ex isa ArgumentError)
            rethrow()
        end
    end
end

function change_file(_)
    selectedroi[] = nothing
end

function start_roi_drawing(_)
    activedrawing[] = true
    selectedroi[] = nothing
end

function delete_roi(_)
    fileselector[] !== nothing && selectedroi[] !== nothing || return
    deleteat!(rois[][fileselector[]], selectedroi[])
    selectedroi[] = nothing
end

function save_rois(_)
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

function load_rois(_)
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

function edit_settings(_)
    Gtk.showall(b["settingsdialog"])
end

function save_settings(_)
    try
        docparameters[] = DoCParameters(localradiustxt[], radiusmaxtxt[], radiussteptxt[], colocthreshtxt[])
        clusterparameters[] = [ClusterParameters(epsilontxt1[], minpointstxt1[], usethreshholdcheckbox1[], smoothingradiustxt1[], minsigclusterpointstxt1[]),
                               ClusterParameters(epsilontxt2[], minpointstxt2[], usethreshholdcheckbox2[], smoothingradiustxt2[], minsigclusterpointstxt2[]),
                               ClusterParameters(epsilontxt3[], minpointstxt3[], usethreshholdcheckbox3[], smoothingradiustxt3[], minsigclusterpointstxt3[])]
        Gtk.hide(b["settingsdialog"])
    catch
        cancel_settings(nothing)
    end
end

function reset_settings(_)
    docparameters[] = defaultparameters
    clusterparameters[] = fill(defaultclusterparameters, 3)
    cancel_settings(nothing)
end

function cancel_settings(_)
    localradiustxt[] = docparameters[].localradius
    radiusmaxtxt[] = docparameters[].radiusmax
    radiussteptxt[] = docparameters[].radiusstep
    colocthreshtxt[] = docparameters[].colocalized_threshold
    epsilontxt1[] = clusterparameters[][1].epsilon
    minpointstxt1[] = clusterparameters[][1].minpoints
    usethreshholdcheckbox1[] = clusterparameters[][1].uselocalradius_threshold
    smoothingradiustxt1[] = clusterparameters[][1].smoothingradius
    minsigclusterpointstxt1[] = clusterparameters[][1].minsigclusterpoints
    epsilontxt2[] = clusterparameters[][2].epsilon
    minpointstxt2[] = clusterparameters[][2].minpoints
    usethreshholdcheckbox2[] = clusterparameters[][2].uselocalradius_threshold
    smoothingradiustxt2[] = clusterparameters[][2].smoothingradius
    minsigclusterpointstxt2[] = clusterparameters[][2].minsigclusterpoints
    epsilontxt3[] = clusterparameters[][3].epsilon
    minpointstxt3[] = clusterparameters[][3].minpoints
    usethreshholdcheckbox3[] = clusterparameters[][3].uselocalradius_threshold
    smoothingradiustxt3[] = clusterparameters[][3].smoothingradius
    minsigclusterpointstxt3[] = clusterparameters[][3].minsigclusterpoints
end

function on_settingsdialog_delete(d, _)
    cancel_settings(nothing)
    Gtk.hide(d)
    return true # prevent dialog from being destroyed
end

function run_clusdoc(_)
    isrunning[] && return
    isrunning[] = true
    save_rois(nothing)
    roicount = sum((n = length(get(rois[], basename(filename), [])); n == 0 ? 1 : n) for filename ∈ inputfiles[])
    progress = progressbar(0:roicount; widget = b["statusprogress"]) # won't update live unless user enables more than one thread
    @idle_add Gtk.start(b["runspinner"])
    clusdoc_task = Threads.@spawn clusdoc(inputfiles[], rois[], localizations[], outputfolder[], colors[], docparameters[], clusterparameters[], outputcsvclustercheckbox[], () -> @idle_add progress[] += 1)
    wait(clusdoc_task)
    @idle_add Gtk.stop(b["runspinner"])
    isrunning[] = false
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
    # should activate spinner here when possible
    isempty(roi) && return
    p1 = first(roi)
    p = p1
    Gtk.move_to(ctx, p[1], p[2])
    Gtk.set_source(ctx, color)
    Gtk.set_line_width(ctx, 3)
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
                notify(rois)
            else
                notify(polyroibuilder)
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
on(edit_settings, settingsbtn)
on(run_clusdoc, runbtn)
on(c -> colors[] = (c, colors[][2:3]...), ch1colorbtn)
on(c -> colors[] = (colors[][1], c, colors[][3]), ch2colorbtn)
on(c -> colors[] = (colors[][1:2]..., c), ch3colorbtn)
on(drawplots, colors)
on(save_settings, settingsok)
on(reset_settings, settingsreset)
signal_connect(on_settingsdialog_delete, b["settingsdialog"], "delete-event")

draw(draw_canvas, imgcanvas, selectedimg, rois, polyroibuilder, nextlineposition, selectedroi)

on(onmouseclick, imgcanvas.mouse.buttonpress)
on(onmousemove, imgcanvas.mouse.motion)
