# ClusDoC.jl
Julia port and extension of the MATLAB-based ClusDoC project: (https://github.com/PRNicovich/ClusDoC)

Installation
============

1. Install Julia (https://julialang.org/downloads/)
2. Open Julia terminal
3. Type `]` to access the package manager prompt
4. Type `add https://github.com/BioTurboNick/ClusDoC.jl` and hit enter to install this package
5. Backspace to return to the main prompt

Usage
=====

1. Open Julia terminal
2. Type `using ClusDoC` and hit enter to load this package
3. Type `clusdoc()` and hit enter to load the GUI
4. Click "Input Files" to select one or more localization list text files to load (currently supported: Nikon Elements output)
5. Click "Output Folder" to choose the folder that the `ClusDoC Results` output folder will be created in.
6. Once both Input Files and Output Folder are selected, images are generated.
7. Optional: Use the controls to add one or more ROIs by clicking to add a point and double-clicking to complete the polygon. Click inside an ROI to select it for deletion. ROIs may be saved in a text file to the chosen output directory or loaded from elsewhere. ROIs are keyed to the image file names.
8. Optional: Change the colors for each channel using the color controls. Each change will regenerate the images.
9. Optional: Use the Settings button to adjust the DoC and DBSCAN clustering settings.
10. Start the analysis with the Run button. Progress will be printed in the terminal window.

Output
======
   - `roicoordinates.txt` containing the ROI shapes
   - Excel (*.xlsx) files for each source image; these contain the main calculated values.
   - Julia data (*.jld2) files for each source image; these contain the raw cluster and localization data which may be used for custom downstream analysis. The `load_raw_results(path)` function may be used to load them.
   - Image (*.png) files for each source image showing the whole field of localizations; this is the image shown in the ROI picker.
   - `cluster maps` contains images of localizations for each ROI and channel showing all the identified clusters and excluded points
   - `doc histograms` contains frequency distributions of the DoC scores for the localizations in each ROI.
   - `doc maps` contains images of localizations for each ROI and channel showing the DoC scores by color.
   - `localization maps` contains images of localizations for each ROI and channel.