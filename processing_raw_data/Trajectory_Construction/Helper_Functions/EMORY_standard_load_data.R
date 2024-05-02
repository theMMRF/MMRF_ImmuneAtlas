# Figure 3, Setup and Data Loading

library(kableExtra)
library(Seurat)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 600 * 1024^3)

setwd("/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/")

source("source/Helper_Functions/general.R", chdir = T)
MEGASEQ <- T
if (dir.exists("/opt/megaseq-data")) {
    dataPath <- "/opt/megaseq-data/wcpilch/MMRF-ImmuneAtlas-Analysis-Sync/data/"
    output <- file.path("/opt/megaseq-data/wcpilch/MMRF-ImmuneAtlas-Analysis-Sync/output", "PAPER_FIGURES")
} else {
    dataPath <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/"
    output <- file.path("/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output", "PAPER_FIGURES")
}
dir.create(output, recursive = T)
# objectPath <- "INTEGRATED_OBJECTS_MMRF/Integration/SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
objectPath <- "objects_stripped/Full/FULL_NO_MD_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
newMD <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/per_cell_md/per_cell_md.rds"
allDimReducs <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata/subumap_dimreduc.rds"

clinical_group <- "progression_group"

reg_auc <- "/opt/megaseq-data/wcpilch/MMRF-ImmuneAtlas-Analysis-Sync/data/SCENIC_FILES/auc_mtx.csv"

output_F2 <- file.path(output, "FIG2_OLD")
dir.create(output_F2, recursive = T)

output_F3 <- file.path(output, "FIG3")
dir.create(output_F3, recursive = T, showWarnings = F)

cluster_naming <- read.csv("/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/custom_MD/annotations_plotting_info_010824.csv", header = T)
lineage_label <- data.frame(lineage_group = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"), lineage_name = c("CD4+ T Cell", "CD8+ T Cell", "B Cell", "Myeloid", "Nk Cell", "Plasma", "Erythroid", "Other", "LQ/Doublet"))
cluster_naming <- dplyr::left_join(cluster_naming, lineage_label, by = "lineage_group")

merged <- readRDS(file.path(dataPath, objectPath))
# merged <- merged %>% fix_public_id()
merged[["umap.sub"]] <- readRDS(allDimReducs) # Adds umap embeddings for all subclusters for display
newMD <- readRDS(newMD)
merged <- AddMetaData(merged, newMD)

merged@meta.data <- merged@meta.data |>
    dplyr::left_join(cluster_naming, by = "subcluster_V03072023", suffix = c("", ".y")) |>
    dplyr::select(-ends_with(".y")) # remove duplicates, keep existing

rownames(merged@meta.data) <- merged$cellname
merged$display_name <- merged$cellID_short

cluster_group <- "subcluster_V03072023"
merged$seurat_cluster_compartment <- merged$subcluster_V03072023_compartment

merged@meta.data$seurat_cluster_subcluster <- merged@meta.data[, cluster_group]
merged@meta.data <- merged@meta.data %>% mutate(treatment_simple = dplyr::case_when(
    d_tx_induction_cat %in% c("imid_pi_steroid", "chemo_imid_pi_steroid") ~ "Triplet",
    .default = "Doublet"
))

merged@meta.data <- merged@meta.data %>% mutate(isDoublet = dplyr::case_when(
    .data[["doublet_pred"]] %in% c("dblet_cluster", "poss_dblet_cluster") ~ "Doublet",
    .default = "Singlet"
))

merged@meta.data$lineage_group <- factor(merged@meta.data$lineage_group, levels = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"))

H_FULL <- 11 / 4.5
H_HALF <- 11 / 9

W_FULL <- 8.5 / 3
W_DBL <- 8.5 / 3 * 2

library(paletteer)

df <- data.frame(
    palette_name = c("lisa::OskarSchlemmer", "NineteenEightyR::miami1", "LaCroixColoR::Lemon", "LaCroixColoR::Orange", "LaCroixColoR::CranRaspberry", "colorBlindness::Blue2Orange8Steps", "khroma::BuRd"),
    fill_NP = c(2, 2, 5, 5, 5, 2, 2),
    fill_RP = c(4, 4, 2, 2, 2, 7, 7),
    color_NP = c(1, 1, 6, 6, 6, 1, 1),
    color_RP = c(5, 5, 1, 1, 1, 8, 8),
    center = c(3, 3, NA, NA, NA, NA, NA)
)

RP_color_pallette <- "lisa::OskarSchlemmer"
NP_color_palette <- "LaCroixColoR::Lemon"
fill_NP <- paletteer_d(NP_color_palette)[[5]] |>
    as.character() |>
    substr(1, 7)
color_NP <- paletteer_d(NP_color_palette)[[6]] |>
    as.character() |>
    substr(1, 7)

fill_RP <- paletteer_d(RP_color_pallette)[[4]] |>
    as.character() |>
    substr(1, 7)
color_RP <- paletteer_d(RP_color_pallette)[[5]] |>
    as.character() |>
    substr(1, 7)

center <- "white"

fill_theme_binary <- c(
    "NP" = fill_NP,
    "RP" = fill_RP
)

color_theme_binary <- c(
    "NP" = color_NP,
    "RP" = color_RP
)

fill_theme_emphasis <- c(
    "NP" = fill_NP,
    "NP*" = color_NP,
    "RP" = fill_RP,
    "RP*" = color_RP,
    "No Difference" = "white"
)

gradientn_colors <- c(fill_theme_binary[["RP"]], "white", fill_theme_binary[["NP"]])

fill_theme_continuous <- list(low = fill_NP, mid = "white", high = color_RP)
