
library(ggplot2)
library(Seurat)

imm_atlas_expr_483 <- readRDS(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/new_data_for_article/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')

Idents(imm_atlas_expr_483) <- imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1
my_levels <- unique(imm_atlas_expr_483$celltype_subclusters_label_transferring_Yizhe_v1)

### looop and merge
irc_pseudo_bulk_pops <- list() ; nlist <- list()
for(i in 1:length(my_levels)) {irc_pseudo_bulk_pops[[i]] <- rowSums(imm_atlas_expr_483@assays$RNA$counts[, which(imm_atlas_expr_483$temp %in% my_levels[i]) ])}
combined_mat <- do.call(cbind,irc_pseudo_bulk_pops)
colnames(combined_mat) <- my_levels
irc_pseudo_bulk_pops <- combined_mat

c('IL12A','IL12B')

ggplot(data=x) + aes(x=var1,y=value) + geom_bar()
