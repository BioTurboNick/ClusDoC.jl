using Gtk, Gtk.ShortNames, GtkObservables, NativeFileDialog

b = GtkBuilder(filename="gui/clusdoc.glade")

win = b["mainwin"]
showall(win)

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

on(selectedfile) do obs
    obs === nothing && return
    locs = loadlocalizations(obs, LocalizationMicroscopy.nikonelementstext)
    ch1locs = getlocalizations(locs, "488", 1, 11000, 100, 10)
    ch2locs = getlocalizations(locs, "647", 11001, 11000, 100, 10)
    
end


filebox2 = Box(:h)
push!(fileboxes, filebox2)
inputrois = Signal("")
add_button!(filebox2, "Input ROIs", () -> push!(inputrois, pick_file()))

is_initializing = false
Gtk.showall(win)



win = Window("Gtk") |> (bx = Box(:v))
lb = Label("(-, -)")
canvas = Canvas(600, 450)
push!(bx, lb, canvas)
@guarded function draw(widget)
    ctx = Gtk.getgc(widget)
    w = Gtk.width(widget)
    h = Gtk.height(widget)
    #ENV["GKS_WSTYPE"] = "142"
    #ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))
    Gtk.rectangle(ctx, 0, 0, w, h)
    Gtk.set_source_rgb(ctx, 1, 0.5, 1)
    Gtk.fill(ctx)
    plot([1, 2, 3])
end

canvas.draw = draw

Gtk.showall(win)


win = Window("Gtk") |> (bx = Box(:v))
c = Canvas(600, 450)
push!(bx, c)
@guarded draw(c) do widget
    ctx = getgc(c)
    # ENV["GKS_WSTYPE"] = "142"
    # ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))
    h = height(c)
    w = width(c)
    Gtk.rectangle(ctx, 0, 0, w, h)
    Gtk.set_source_rgb(ctx, 1, 0.5, 1)
    Gtk.fill(ctx)
    #plot([1, 2, 3], size=(h, w))
end
c.draw = draw
showall(win)
