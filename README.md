# ClusDoC.jl
Julia port and extension of the MATLAB-based ClusDoC project: (https://github.com/PRNicovich/ClusDoC)

Installation
============

1. Install Julia (https://julialang.org/downloads/)
2. Open Julia terminal
3. Type `]` to access the package manager prompt
4. Type `add https://github.com/BioTurboNick/ClusDoC.jl` and hit enter to install this package
5. Backspace to return to the main prompt

Note: Highly recommend setting the default number of threads Julia uses to at least 2.

Usage
=====

1. Open Julia terminal
2. Type `using ClusDoC` and hit enter to load this package and the GUI.
3. Click "Input Files" to select one or more localization list text files to load (currently supported: Nikon Elements output)
4. Click "Output Folder" to choose the folder that the `ClusDoC Results` output folder will be created in.
5. Once both Input Files and Output Folder are selected, images are generated.
6. Optional: Use the controls to add one or more ROIs by clicking to add a point and double-clicking to complete the polygon. Click inside an ROI to select it for deletion. ROIs may be saved in a text file to the chosen output directory or loaded from elsewhere. ROIs are keyed to the image file names.
7. Optional: Change the colors for each channel using the color controls. Each change will regenerate the images.
8. Optional: Use the Settings button to adjust the DoC and DBSCAN clustering settings.
9. Start the analysis with the Run button. Progress will be printed in the terminal window.

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

*NOTE*: There are some variations in calculations between this version and the MATLAB version. For example, the MATLAB version would use a square ROI area, rather than the computed polygon area. For a more impactful example, the MATLAB version used a separate number of points per cluster cutoff when reporting DBSCAN cluster
results, but not when reporting the colocalized/noncolocalized cluster subsets. In this version, both apply the cutoff for consistency.

Troubleshooting
================
   - If you encounter an error, you may need to restart the Julia session.
