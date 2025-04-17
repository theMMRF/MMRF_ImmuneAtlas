This folder contains 'helper functions' which were used internally across multiple scripts, primarily for progression-based analyses. This includes tasks like loading and incorporating the up-to-date metadata, various differential abundance tests, differential expresssion tests, and so on.

*__DotPlot_custom.R - A modified version of Seurat's dotplot to render as a grid of patchworked plots, incorporate cluster colors, and allow for grouping on the celltypes as well as genes. Does not alter the actual calculations for Seurat's DotPlot
*COMMON_diff_abundance_functions.R - Implementations of Wilcoxon, Beta Regression, and Dirichlet Regression. An implementation to compute a 95% CI for dirichlet means are also included
*DEG_functions.R - Code for computing differentially expressed genes, or for plotting various markers. Used to render DotPlots, FeaturePlots, VlnPlots, and to compute markers with Seurat's FindMarkers (See relevnat figures for limma code)
*EMORY_standard_load_data.R - Loads the current version of the object metadata, incoporates cell annotations, sets colors and sizes for figure rendering, and primary output directory
*general.R - General functions which depend on the above marker and DEG code. Includes code to 'subcluster' a given object into various subclusters and compute top markers, along with code to run Trajectory Analysis
*pdf_and_png.R - Function to conveniently save a .png and .pdf version of a given ggplot object with similar scaling and overall visual appearance
