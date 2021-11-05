using Gtk, Gtk.ShortNames, GtkObservables, NativeFileDialog, GR, Printf, LocalizationMicroscopy

b = GtkBuilder(filename="gui/clusdoc.glade")
canvas = Canvas(500, 500)
push!(b["canvasbox"], canvas)

win = b["mainwin"]

inputfolder = Observable("")
inputfiles = Observable([""])
on(button(widget = b["inputbtn"])) do _
    files = pick_multi_file()
    if length(files) > 0
        inputfiles[] = files
    end
end
inputtxt = textbox(String; widget = b["inputtxt"])
on(inputfiles) do obs
    if length(obs) > 0
        txt = "\"" * join(obs, "\" \"") * "\""
    else
        txt = ""
    end
    if txt != inputtxt[]
        inputtxt[] = txt
    end
end
on(inputtxt) do obs
    files = filter!(x -> !contains(x, r"^\s*$"), split(obs, "\"", keepempty = false))
    if files != inputfiles[]
        inputfiles[] = files
    end
end

outputfolder = Observable("")
on(button(widget = b["outputbtn"])) do _
    outputfolder[] = joinpath(pick_folder(), "ClusDoC Results")
end
textbox(String; widget = b["outputtxt"], observable = outputfolder)

selectedfile = Observable{Union{Nothing, String}}(nothing)
fileselector = dropdown([], widget = b["fileselector"])
on(inputfiles) do obs
    @idle_add begin
        empty!(fileselector)
        append!(fileselector, basename.(obs))
        fileselector[] = basename(inputfiles[][1])
        selectedfile[] = inputfiles[][1]
    end
end
on(fileselector) do obs
    selectedfile[] = obs !== nothing ? inputfiles[][findfirst(x -> basename(x) == obs, inputfiles[])] : nothing
end

selectedlocs = Observable((Matrix{Float64}(undef, 2, 0), Matrix{Float64}(undef, 2, 0)))

function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1 # activate these limits

    localizations = filter(l -> l.channel == channelname, alllocalizations)
end

on(selectedfile) do obs
    obs === nothing && return
    locs = loadlocalizations(obs, LocalizationMicroscopy.nikonelementstext)
    ch1 = extractcoordinates(getlocalizations(locs, "488", 1, 11000, 100, 10))
    ch2 = extractcoordinates(getlocalizations(locs, "647", 11001, 11000, 100, 10))
    selectedlocs[] = (ch1, ch2)
end

function cdraw(widget)
    @idle_add begin
    ctx = Gtk.getgc(canvas)
    h = Gtk.height(widget)
    w = Gtk.width(widget)
    ENV["GKS_WSTYPE"] = "142"
    ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))
    plt = gcf()
    plt[:size] = (w, h)
    #scatter([1, 2, 3], [4, 5, 6])
    scatter(selectedlocs[][1][1,:], selectedlocs[][1][2,:])
    # this draws, but it's so much data it's probably better to generate image files once and then display the image
    end
end
draw(cdraw, canvas) # provide function to run when canvas redraws

on(selectedlocs) do _
    cdraw(canvas)
end
showall(win)

filebox2 = Box(:h)
push!(fileboxes, filebox2)
inputrois = Signal("")
add_button!(filebox2, "Input ROIs", () -> push!(inputrois, pick_file()))



