library(Seurat)
library(SingleR)
library(Matrix)
library(Seurat)
library(rhdf5)
library(reticulate)

### function +++
ExportH5AD <- function(obj, # Seurat object
                       file # Path to save the .h5ad file
){
  sc <- import("scanpy")
  ad <- import("anndata")
  
  expression_matrix <- t(GetAssayData(obj, layer = 'counts'))
  metadata <- obj@meta.data
  feature_names <- data.frame(gene = colnames(expression_matrix))  
  
  # Create an AnnData object with X, obs, and var
  adata <- sc$AnnData(
    X = expression_matrix, 
    obs = metadata,
    var = feature_names  # Add gene names to the 'var' slot
  )
  
  # Write the AnnData object to an .h5ad file
  adata$write_h5ad(file)
}
### +++

setwd('/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas')

### updated meta
#updated_meta_per_cell <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/updateMMRF_yered_william/091524_per_cell_md.rds')
### single cell data
imm_atlas_expr_483 <- readRDS(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/new_data_for_article/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')

### clean metadata
#updated_meta_per_cell <- updated_meta_per_cell[updated_meta_per_cell$cohort %in% 'discovery',]
#rm(updated_meta_per_cell)

### 
rm(list=setdiff(ls(), c('imm_atlas_expr_483','updated_meta_per_cell')))

### store raw data
py_config()

use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)

imm_atlas_expr_483@meta.data <- imm_atlas_expr_483@meta.data[,c('sample_id','seurat_subclusters_label_transferring_Yizhe_v1',
                                                                'celltype_subclusters_label_transferring_Yizhe_v1','compartment','siteXbatch')]
mdata <- imm_atlas_expr_483@meta.data
data <- imm_atlas_expr_483[['RNA']]$counts
rm(imm_atlas_expr_483)

Matrix::writeMM(data, "/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/imm_atlas_expr_discovery_locked_dgCMatrix.mtx")
write.csv(data, "/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/imm_atlas_expr_discovery_locked_metadata.csv")

### Store Everything (VERY large for local)
#ExportH5AD(obj=imm_atlas_expr_483, file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/imm_atlas_expr_discovery_locked.h5ad')

#imm_atlas_expr_483@meta.data <- imm_atlas_expr_483@meta.data[,]
#table(imm_atlas_expr_483$seurat_subclusters_label_transferring_Yizhe_v1)
#table(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1)
#imm_atlas_expr_483$

### subset

#sort(table(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1))
x <- as.data.frame(sort(table(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1)))
x$filter <-  x$Freq > 1000
x$Var1 <- as.character(x$Var1)

my_cells <- list()
for(ii in 1:nrow(x)){
  if (x$filter[ii] == TRUE) {
    my_cells[[ii]] <- sample(which(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1 %in% x$Var1[ii]), 1000)
  } else {
    my_cells[[ii]] <- which(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1 %in% x$Var1[ii])
  }
}

my_subset <- imm_atlas_expr_483[,unlist(my_cells)]

### store subset data
ExportH5AD(obj=my_subset, file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/imm_atlas_expr_discovery_locked_subset.h5ad')
### test (maybe later)
rm(my_subset)




#############################################################################################################################################
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','immunarch',
        'limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','variancePartition','ggbeeswarm')
lapply(libs, require, character.only = TRUE) ; rm(libs)
#############################################################################################################################################
### Work space
### CHANGE AS NECESSARY IN YOUR COMPUTER
setwd("/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/")
#############################################################################################################################################


#############################################################################################################################################
load(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/temp_re_analysis_sep2024.RData')
#############################################################################################################################################

dim(subcluster_all_fulldf)
dim(updated_meta_per_cell_validation)
dim(validation_dataset_w_meta)

plot(density(updated_meta_per_cell_validation$d_dx_amm_age))
plot(density(updated_meta_per_cell_validation$d_dx_amm_age))
#dev.off()
plot(density(subcluster_comp_fulldf$d_dx_amm_age))


#############################################################################################################################################

# 
# Age -> Cell Proportions
# Age -> Gene Expression
# Variables: Prior Treatment / State / OS / PFS

#############################################################################################################################################
