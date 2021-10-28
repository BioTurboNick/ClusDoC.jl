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
    outputfolder[] = pick_folder()
end
textbox(String; widget = b["outputtxt"], observable = outputfolder)

selectedfile = Observable("")
on(inputfiles) do obs
    dropdown(basename.(obs), widget = b["fileselector"], observable = selectedfile)
end


filebox2 = Box(:h)
push!(fileboxes, filebox2)
inputrois = Signal("")
add_button!(filebox2, "Input ROIs", () -> push!(inputrois, pick_file()))

is_initializing = false
Gtk.showall(win)
