library(Seurat)
library(harmony)
library(tibble)
library(dplyr)

# merged seurat object using merge function in Seurat
sobj <- readRDS("SeuratObj_rm38_44_in_361_samples_Human_Ref_CB.rds")

# LogNormalize 
raw.count <- GetAssayData(object = sobj[["RNA"]], slot = "counts")

combined <- CreateSeuratObject(counts = raw.count, project = "MMRF", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = sobj@var.genes, npcs = 25, verbose = FALSE)

###############################################
### Harmony batch correction
##############################################
# harmony batch correct by both site and batch
combined@meta.data$siteXbatch <- paste0(combined@meta.data$site,"_",combined@meta.data$batch_num)

out_path <- "path_of_output_folder"
combined <- combined %>% RunHarmony(c("siteXbatch"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(combined, 'harmony')

combined <- combined %>%
    RunUMAP(reduction = "harmony", dims = 1:25) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25, force.recalc = T) %>%
    FindClusters(resolution = 0.5) %>%
    identity()


saveRDS(combined,paste0(out_path,"/xx_Harmony.rds"))
