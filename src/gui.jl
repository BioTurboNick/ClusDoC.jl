using Gtk.ShortNames, GtkObservables, NativeFileDialog, Printf, Plots, LocalizationMicroscopy, Images, ImageIO # find out which subpackages of Images I need

#utility functions
function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1 # activate these limits

    return filter(l -> l.channel == channelname, alllocalizations)
end


# independent observables
inputfiles = Observable([""])
outputfolder = Observable("")
localizations = Observable(Vector{Vector{Localization}}[]) # vector for each channel in each image
selectedimg = Observable{Union{Nothing, Matrix{RGB{N0f8}}}}(nothing)


# initialize UI elements
b = Gtk.GtkBuilder(filename="gui/clusdoc.glade")
imgcanvas = canvas(DeviceUnit, 500, 500)
push!(b["canvasbox"], imgcanvas)
win = b["mainwin"]

inputbtn = button(widget = b["inputbtn"])
inputtxt = textbox(String; widget = b["inputtxt"])
outputbtn = button(widget = b["outputbtn"])
textbox(String; widget = b["outputtxt"], observable = outputfolder)
fileselector = dropdown([], widget = b["fileselector"])

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
    outputfolder[] = path == "" ? "" : joinpath(path, "ClusDoC Results")
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
        ch1 = getlocalizations(locs, "488", 1, 11000, 100, 10)
        ch2 = getlocalizations(locs, "647", 11001, 11000, 100, 10) # need to generalize, obviously
        push!(localizations[], [ch1, ch2])
    end
    notify(localizations)
end

function drawplots(_)
    outputfolder[] == "" && return
    for (i, locs) ∈ enumerate(localizations[])
        ch1 = extractcoordinates(locs[1])
        ch2 = extractcoordinates(locs[2]) # generalize
        Plots.scatter(ch1[1, :], ch1[2, :], markercolor = RGBA(1.0, 0.0, 0.0, 0.5), markersize = 2, aspectratio = :equal, size=(1024, 1024), markerstrokewidth = 0)
        Plots.scatter!(ch2[1, :], ch2[2, :], markercolor = RGBA(0.0, 1.0, 0.0, 0.5), markersize = 2, aspectratio = :equal, size=(1024, 1024), markerstrokewidth = 0)
        Plots.plot!(ticks=:none, legend = :none, axis = false, margin=0 * Plots.mm)
        path = joinpath(outputfolder[], "localizationmaps")
        mkpath(path)
        imagepath = joinpath(path, basename(inputfiles[][i]) * ".png")
        Plots.savefig(imagepath)
        # could probably generate plots, but delay saving until an output folder selected.
    end
    selectedimg[] = load_image(fileselector[])
end

function load_image(obs)
    try
        # may fire before images created
        selectedimg[] = load(joinpath(outputfolder[], "localizationmaps", obs * ".png"))
    catch
    end
end

function draw_canvas(_)
    selectedimg[] !== nothing || return
    copy!(imgcanvas, selectedimg[])
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
on(draw_canvas, selectedimg)
draw(draw_canvas, imgcanvas)

# should have a colorset selector?


function cdraw(widget)
    @idle_add begin
        ctx = Gtk.getgc(widget)
        h = Gtk.height(widget)
        w = Gtk.width(widget)
        ENV["GKS_WSTYPE"] = "142"
        ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))
        plt = gcf()
        plt[:size] = (w, h)
        if fileselector[] !== nothing
            
            GR.imshow(joinpath(outputfolder[], "localizationmaps", fileselector[] * ".svg"))
        end
        delete!(ENV, "GKS_WSTYPE")
        delete!(ENV, "GKSconid")
    end
end

draw(cdraw, canvas) # provide function to run when canvas redraws

showall(win)

filebox2 = Box(:h)
push!(fileboxes, filebox2)
inputrois = Signal("")
add_button!(filebox2, "Input ROIs", () -> push!(inputrois, pick_file()))



