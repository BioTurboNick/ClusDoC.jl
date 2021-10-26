using Gtk, Gtk.ShortNames, GtkReactive, NativeFileDialog

function connect_button(btnwidget, handler)
    btn = button(widget = btnwidget)
    preserve(map(btn) do _
        is_initializing && return
        handler()
    end)
end

function connect_entry(txtwidget, signal::Signal{String})
    textbox(String; widget = txtwidget, signal)
    return nothing
end

function connect_entry(txtwidget, signal::Signal{Vector{String}})
    txtsig = preserve(map(signal) do sig
        # format paths
        is_initializing && return ""
        return "\"" * join(sig, "\" \"") * "\""
    end)
    tb = textbox(String; widget = txtwidget, signal = txtsig)
    preserve(map(tb) do sig
        # parse input file list
        is_initializing && return String[]
        println("parsing $sig")
        return String[]
    end)
end

is_initializing = true

b = GtkBuilder(filename="gui/clusdoc.glade")

win = b["mainwin"]
showall(win)

inputfiles = Signal([""])
connect_button(b["inputbtn"], () -> push!(inputfiles, pick_multi_file()))
connect_entry(b["inputtxt"], inputfiles)

outputfolder = Signal("")
connect_button(b["outputbtn"], () -> push!(outputfolder, pick_folder()))
connect_entry(b["outputtxt"], outputfolder)



filebox2 = Box(:h)
push!(fileboxes, filebox2)
inputrois = Signal("")
add_button!(filebox2, "Input ROIs", () -> push!(inputrois, pick_file()))

is_initializing = false
Gtk.showall(win)
