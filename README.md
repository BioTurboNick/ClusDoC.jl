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
2. Type `using ClusDoC` and hit enter to load this package
3. Type `clusdoc()` and hit enter to load the GUI
4. Click "Input Files" to select one or more localization list text files to load (currently supported: Nikon Elements output)
5. Click "Output Folder" to choose the folder that the `ClusDoC Results` output folder will be created in.
6. Once both Input Files and Output Folder are selected, images are generated.
7. Optional: Use the controls to add one or more ROIs by clicking to add a point and double-clicking to complete the polygon. Click inside an ROI to select it for deletion. ROIs may be saved in a text file to the chosen output directory or loaded from elsewhere. ROIs are keyed to the image file names. NOTE: Avoid including too much dead space in each ROI.
8. Optional: Change the colors for each channel using the color controls. Each change will regenerate the images.
9. Optional: Use the Settings button to adjust the DoC and DBSCAN clustering settings (see below).
10. Start the analysis with the Run button. Progress will be printed in the terminal window.

Settings
========
DoC:
* Local Radius (nm): Used to examine the local neighborhood of each localization to determine if there are more nearby than expected by chance within the ROI. Localizations may be filtered 
* Radius Max (nm): When calculating the DoC score for each localization, all points within this radius will be considered.
* Radius Step (nm): The DoC score examines the localizations located within each band extending out from a localization.
* Colocalized DoC Limit: Classify localizations as colocalized if DoC exceeds this value.

DBSCAN Clustering:
* Epsilon (nm): DBSCAN radius, to determine which localizations are neighbors.
* Min Points: The minimum number of localizations needed to form a cluster.
* Use Threshold: Exclude localizations with fewer neighbors than expected by chance, from Local Radius (nm).
* Smoothing Radius (nm): Controls the smoothing of the cluster perimeter.
* Min Points/Sig. Cluster: The lower bound for the number of localizations that comprise a "significant" cluster, and for considering a cluster to be interacting.

Output
======
   - `roicoordinates.txt` containing the ROI shapes
   - Excel (*.xlsx) files for each source image; these contain the main calculated values, one row summarizing each ROI.
   - Julia data (*.jld2) files for each source image; these contain the raw cluster and localization data which may be used for custom downstream analysis. The `load_raw_results(path)` function may be used to load them.
   - Image (*.png) files for each source image showing the whole field of localizations; this is the image shown in the ROI picker.
   - `cluster maps` contains images of localizations for each ROI and channel showing all the identified clusters and excluded points
   - `doc histograms` contains frequency distributions of the DoC scores for the localizations in each ROI.
   - `doc maps` contains images of localizations for each ROI and channel showing the DoC scores by color.
   - `localization maps` contains images of localizations for each ROI and channel.

Data Reference
=================
   - DoC Results: Percentage of colocalized molecules for each channel pair
   - Clustering Results X: DBSCAN output for clustering channel X, for all clusters and significant clusters (defined by `Min Points/Sig. Cluster` setting)
   - Clus-DoC Results X: Colocalization data for signficant clusters in channel X to each other channel. Colocalized clusters are defined by having at least `Min Points/Sig. Cluster` points with DoC scores exceeding `Colocalized DoC Limit`. Noncolocalized clusters have no colocalization, and Intermediate clusters are those with minimal colocalization (below `Min Points/Sig. Cluster`).
   - Algorithm parameters: Stores the settings used to create the data.
   - Columns:
       * ROI area (μm)
       * Number of clusters
       * Density of clusters (clusters / μm^2)
       * Cluster area (nm^2)
       * Cluster circularity
       * Number of localizations in ROI
       * Fraction of localizations in clusters
       * Absolute density in clusters (localizations / μm^2)
       * Relative density in clusters: A measure of granularity; high indicates the cluster is clumpy, low indicates evenly distributed.


*NOTE*: There are some variations in calculations between this version and the MATLAB version. For example, the MATLAB version would use a square ROI area, rather than the computed polygon area. For a more impactful example, the MATLAB version used a separate number of points per cluster cutoff when reporting DBSCAN cluster
results, but not when reporting the colocalized/noncolocalized cluster subsets. In this version, both apply the cutoff for consistency.

Troubleshooting
================
   - If you encounter an error, you may need to restart the Julia session.
   - There is currently an issue on some Mac computers where the file selection dialog is broken. As a workaround, you could collect a list of file paths separated by quotation marks and spaces and paste it into the box. E.g. `"/path/to/file1.txt" "/path/to/file2.txt"` Alternatively, you can directly provide an array of filenames to the observable tracking them as follows: `ClusDoC.inputfiles[] = ["/path/to/file1.txt", "/path/to/file2.txt"]`
