# Loads integrated discovery + validation object, along with extended metadata fields. When subset to the 'discovery' cohort, contents of the count matrix and clinical covariates are identical to the previous object. Updated analyses using new metadata, such as the extended risk analysis, will use this version of the object.
library(kableExtra)
library(Seurat)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 600 * 1024^3)

setwd("MMRF-ImmuneAtlas-Annotation/")

source("source/Helper_Functions/general.R", chdir = T)
MEGASEQ <- T
if (dir.exists("/opt/megaseq-data")) {
    dataPath <- "MMRF-ImmuneAtlas-Analysis-Sync/data/"
    output <- file.path("/opt/megaseq-data/wcpilch/MMRF-ImmuneAtlas-Analysis-Sync/output", "PAPER_FIGURES")
} else {
    dataPath <- "MMRF-ImmuneAtlas-Annotation/data/"
    output <- file.path("MMRF-ImmuneAtlas-Annotation/output", "PAPER_FIGURES")
}
dir.create(output, recursive = T)
# objectPath <- "INTEGRATED_OBJECTS_MMRF/Integration/SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
objectPath <- "objects_stripped/Full/FULL_NO_MD_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds"
newMD <- "MMRF-ImmuneAtlas-Annotation/data/metadata_092424/per_cell_metadata/092424_per_cell_md.rds" # "MMRF-ImmuneAtlas-Annotation/data//per_cell_md.rds"
# allDimReducs <- "MMRF-ImmuneAtlas-Annotation/data/metadata/subumap_dimreduc.rds"

objectPath <- "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/VALIDATION_COHORT/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds"


clinical_group <- "progression_group"

reg_auc <- "/opt/megaseq-data/wcpilch/MMRF-ImmuneAtlas-Analysis-Sync/data/SCENIC_FILES/auc_mtx.csv"

output_F2 <- file.path(output, "FIG2_OLD")
dir.create(output_F2, recursive = T)

output_F3 <- file.path(output, "FIG3")
dir.create(output_F3, recursive = T, showWarnings = F)

# Maps cellID_short -> subcluster_V03072023
cluster_naming <- read.csv("MMRF-ImmuneAtlas-Annotation/data/custom_MD/annotations_plotting_info_010824.csv", header = T)
lineage_label <- data.frame(lineage_group = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"), lineage_name = c("CD4+ T Cell", "CD8+ T Cell", "B Cell", "Myeloid", "Nk Cell", "Plasma", "Erythroid", "Other", "LQ/Doublet"))
cluster_naming <- dplyr::left_join(cluster_naming, lineage_label, by = "lineage_group")

merged <- readRDS(file.path(objectPath))
# merged <- merged %>% fix_public_id()
# merged[["umap.sub"]] <- readRDS(allDimReducs) # Adds umap embeddings for all subclusters for display
merged[["umap.sub"]] <- merged[["umap"]]
newMD <- readRDS(newMD)
merged <- AddMetaData(merged, newMD)

merged@meta.data$LABEL_TRANSFER_subcluster_V03072023 <- merged@meta.data$seurat_subclusters_label_transferring_Yizhe_v1

cluster_naming_mod <- cluster_naming
colnames(cluster_naming_mod) <- paste0("LABEL_TRANSFER_", colnames(cluster_naming_mod))
merged@meta.data <- merged@meta.data |>
    dplyr::left_join(cluster_naming_mod, by = "LABEL_TRANSFER_subcluster_V03072023", suffix = c("", ".y")) |>
    dplyr::select(-ends_with(".y")) # remove duplicates, keep existing

merged@meta.data$INTEGRATION_subcluster_V03072023 <- merged@meta.data$seurat_subclusters_integration_Yizhe_v1
cluster_naming_mod <- cluster_naming
colnames(cluster_naming_mod) <- paste0("INTEGRATION_", colnames(cluster_naming_mod))
merged@meta.data <- merged@meta.data |>
    dplyr::left_join(cluster_naming_mod, by = "INTEGRATION_subcluster_V03072023", suffix = c("", ".y")) |>
    dplyr::select(-ends_with(".y")) # remove duplicates, keep existing


# Formatting the label-mapped data.
# Integrated label data is stored in the default columns. This value is identical to the previous labels for samples in the 'discovery' cohort
merged@meta.data$subcluster_V03072023 <- merged@meta.data$INTEGRATION_subcluster_V03072023
merged@meta.data$lineage_group <- merged@meta.data$INTEGRATION_lineage_group
merged@meta.data$lineage_order <- merged@meta.data$INTEGRATION_lineage_order
merged@meta.data$color <- merged@meta.data$INTEGRATION_color
merged@meta.data$color_sub <- merged@meta.data$INTEGRATION_color_sub
merged@meta.data$primary_id_color <- merged@meta.data$INTEGRATION_primary_id_color

merged@meta.data$subcluster_V03072023_compartment <- merged@meta.data$compartment

merged@meta.data$cellID_short <- merged@meta.data$INTEGRATION_cellID_short
merged@meta.data$cellID_long <- merged@meta.data$INTEGRATION_cellID_long

rownames(merged@meta.data) <- merged$cellname
# merged$display_name <- merged$cellID_short

# cluster_group <- "subcluster_V03072023"
# merged$seurat_cluster_compartment <- merged$subcluster_V03072023_compartment

# merged@meta.data$seurat_cluster_subcluster <- merged@meta.data[, cluster_group]
merged@meta.data <- merged@meta.data %>% mutate(treatment_simple = dplyr::case_when(
    d_tx_induction_cat %in% c("imid_pi_steroid", "chemo_imid_pi_steroid") ~ "Triplet",
    .default = "Doublet"
))

merged@meta.data <- merged@meta.data %>% mutate(isDoublet = dplyr::case_when(
    .data[["INTEGRATION_doublet_pred"]] %in% c("dblet_cluster", "poss_dblet_cluster") ~ "Doublet",
    .default = "Singlet"
))

# merged@meta.data$lineage_group <- factor(merged@meta.data$lineage_group, levels = c("CD4", "CD8", "B", "M", "Nk", "P", "E", "Other", "LQ"))

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
