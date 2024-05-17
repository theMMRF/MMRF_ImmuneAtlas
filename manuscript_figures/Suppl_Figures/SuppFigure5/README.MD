# Longitudinal Analysis of MMRF Data

This  contains an R script for analyzing and visualizing longitudinal data from the Multiple Myeloma Research Foundation (MMRF) study. The script performs various data processing steps, including reading data, merging datasets, calculating proportions, and plotting results.

## Requirements

The following R packages are required to run the script:

- reticulate
- ggplot2
- dplyr
- magrittr
- patchwork
- viridis
- qs
- BiocParallel
- data.table
- IRdisplay
- scales
- readxl
- tidyr
- stringr
- ggbeeswarm
- broom

Additionally, the script uses the `reticulate` package to interface with Python. Ensure that the necessary Python environment and packages are installed.

## Setup

1. Clone this repository to your local machine.
2. Set up the required R and Python packages.
3. Adjust the file paths in the script to point to your data directories.

## Usage

Run the `longitudinal_analysis.R` script to perform the analysis and generate plots. The script performs the following steps:

1. Sets up library paths and loads necessary libraries.
2. Reads data from various sources.
3. Merges datasets and processes data.
4. Calculates proportions of cell types.
5. Generates plots to visualize the results.
6. Saves the plots to the specified directory.

## Output

The script generates and saves various plots, including histograms and trend plots, to the specified output directory. The plots provide insights into the distribution and trends of cell proportions over time.

