#############################################################################################################################################
### Edgar Gonzalez-Kozlova
#############################################################################################################################################
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','monocle3',
        'limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2')
lapply(libs, require, character.only = TRUE) ; rm(libs)
#############################################################################################################################################
### Work space
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/")
#############################################################################################################################################


### load and prepare single cell data
#############################################################################################################################################
### everything
immune_atlas <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds')
DimPlot(object = immune_atlas, reduction = "umap")
dev.off()
imm_atlas_meta <- immune_atlas@meta.data
imm_atlas_meta$cell_ids <- rownames(imm_atlas_meta)
rownames(imm_atlas_meta) <- NULL
#+
FeaturePlot(immune_atlas, features = 'PDCD1')
#+
VlnPlot(immune_atlas, features = 'PDCD1')
#+
grep( "PDCD1" , rownames(immune_atlas),value = TRUE )
#+
#+
#+
# cells: 1149344

### plasma
tmp_var <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_Plasma_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Plasma_LogNorm_PC25_Harmony_correct_sampleid_plasma_ann.rds')
Idents(tmp_var) <- tmp_var$seurat_clusters
DimPlot(object = tmp_var, reduction = "umap")
dev.off()
colnames(tmp_var@meta.data)[which(colnames(tmp_var@meta.data) %in% 'seurat_clusters')] <- 'SPlasma_Subclusters'
tmp_var@meta.data$cell_ids <- rownames(tmp_var@meta.data)
tmp_var_meta <- tmp_var@meta.data[,which(colnames(tmp_var@meta.data) %in% c('cell_ids','SPlasma_Subclusters'))]
rownames(tmp_var_meta) <- NULL
### merge
imm_atlas_meta <- merge(imm_atlas_meta,tmp_var_meta,by='cell_ids',all=TRUE)

### b_ery
tmp_var <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_B_ery_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_B_ery_PC25_Harmony.rds')
Idents(tmp_var) <- tmp_var$seurat_clusters
DimPlot(object = tmp_var, reduction = "umap")
dev.off()
colnames(tmp_var@meta.data)[which(colnames(tmp_var@meta.data) %in% 'seurat_clusters')] <- 'Sb_ery_Subclusters'
tmp_var@meta.data$cell_ids <- rownames(tmp_var@meta.data)
tmp_var_meta <- tmp_var@meta.data[,which(colnames(tmp_var@meta.data) %in% c('cell_ids','Sb_ery_Subclusters'))]
rownames(tmp_var_meta) <- NULL
### merge
imm_atlas_meta <- merge(imm_atlas_meta,tmp_var_meta,by='cell_ids',all=TRUE)

### NK_T
tmp_var <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_NK_T_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_NK_T_PC30_Harmony.rds')
Idents(tmp_var) <- tmp_var$seurat_clusters
DimPlot(object = tmp_var, reduction = "umap")
dev.off()
colnames(tmp_var@meta.data)[which(colnames(tmp_var@meta.data) %in% 'seurat_clusters')] <- 'SNKT_Subclusters'
tmp_var@meta.data$cell_ids <- rownames(tmp_var@meta.data)
tmp_var_meta <- tmp_var@meta.data[,which(colnames(tmp_var@meta.data) %in% c('cell_ids','SNKT_Subclusters'))]
rownames(tmp_var_meta) <- NULL
### merge
imm_atlas_meta <- merge(imm_atlas_meta,tmp_var_meta,by='cell_ids',all=TRUE)

### myeloid
tmp_var <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_Myeloid_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_MacroMono_PC30_Harmony.rds')
Idents(tmp_var) <- tmp_var$seurat_clusters
DimPlot(object = tmp_var, reduction = "umap")
dev.off()
colnames(tmp_var@meta.data)[which(colnames(tmp_var@meta.data) %in% 'seurat_clusters')] <- 'Smyeloid_Subclusters'
tmp_var@meta.data$cell_ids <- rownames(tmp_var@meta.data)
tmp_var_meta <- tmp_var@meta.data[,which(colnames(tmp_var@meta.data) %in% c('cell_ids','Smyeloid_Subclusters'))]
rownames(tmp_var_meta) <- NULL
### merge
imm_atlas_meta <- merge(imm_atlas_meta,tmp_var_meta,by='cell_ids',all=TRUE)

### eryth
tmp_var <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_Erythrocyte_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Ery_PC15_Harmony.rds')
Idents(tmp_var) <- tmp_var$seurat_clusters
DimPlot(object = tmp_var, reduction = "umap")
dev.off()
colnames(tmp_var@meta.data)[which(colnames(tmp_var@meta.data) %in% 'seurat_clusters')] <- 'Seryth_Subclusters'
tmp_var@meta.data$cell_ids <- rownames(tmp_var@meta.data)
tmp_var_meta <- tmp_var@meta.data[,which(colnames(tmp_var@meta.data) %in% c('cell_ids','Seryth_Subclusters'))]
rownames(tmp_var_meta) <- NULL
### merge
imm_atlas_meta <- merge(imm_atlas_meta,tmp_var_meta,by='cell_ids',all=TRUE)
### 
rm(tmp_var_meta)

###
head(imm_atlas_meta)
imm_atlas_meta$complete_subclusters <- NA
imm_atlas_meta$complete_subclusters[!is.na(imm_atlas_meta$SPlasma_Subclusters)] <- paste('Plasma',imm_atlas_meta$SPlasma_Subclusters[!is.na(imm_atlas_meta$SPlasma_Subclusters)],sep='_')
imm_atlas_meta$complete_subclusters[!is.na(imm_atlas_meta$Sb_ery_Subclusters)] <- paste('BERY',imm_atlas_meta$Sb_ery_Subclusters[!is.na(imm_atlas_meta$Sb_ery_Subclusters)],sep='_')
imm_atlas_meta$complete_subclusters[!is.na(imm_atlas_meta$SNKT_Subclusters)] <- paste('NKT',imm_atlas_meta$SNKT_Subclusters[!is.na(imm_atlas_meta$SNKT_Subclusters)],sep='_')
imm_atlas_meta$complete_subclusters[!is.na(imm_atlas_meta$Smyeloid_Subclusters)] <- paste('Myeloid',imm_atlas_meta$Smyeloid_Subclusters[!is.na(imm_atlas_meta$Smyeloid_Subclusters)],sep='_')
imm_atlas_meta$complete_subclusters[!is.na(imm_atlas_meta$Seryth_Subclusters)] <- paste('Eryth',imm_atlas_meta$Seryth_Subclusters[!is.na(imm_atlas_meta$Seryth_Subclusters)],sep='_')
imm_atlas_meta$complete_subclusters[is.na(imm_atlas_meta$complete_subclusters)] <- c('C23')

### imm_atlas_meta$complete_subclusters
table(is.na(imm_atlas_meta$complete_subclusters))
imm_atlas_meta$seurat_clusters[is.na(imm_atlas_meta$complete_subclusters)]

rm(immune_atlas,tmp_var)
#############################################################################################################################################


### prepare and clean clinical data
#############################################################################################################################################
### pt_metadata <- as.data.frame(readxl::read_excel(path = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/CD138-depleted BMA Cell Samples_Immune Atlas IA-001_Clinical Charactersitics_07302021v1.xlsx',sheet = 1))
pt_metadata <- read.csv('CD138-depleted BMA Cell Samples_Immune Atlas IA-001_Clinical Charactersitics_07302021v1.csv')
pt_metadata$sample_id <- paste(substr(pt_metadata$Public_ID,6,9),substr(pt_metadata$Aliquot.ID,7,12),sep="_")

### For Seurat Clusters
abundances <- unclass(table(imm_atlas_meta$sample_id, imm_atlas_meta$seurat_clusters))

### remove the SM
rownames(abundances) <- substr(rownames(abundances),1,11)
table(rownames(abundances) %in% pt_metadata$sample_id) ## 2 samples unlocated

### filter
abundances <- abundances[which(rownames(abundances) %in% pt_metadata$sample_id),]
pt_metadata <- pt_metadata[which(pt_metadata$sample_id %in% rownames(abundances)),]
### match
pt_metadata <- pt_metadata[match(rownames(abundances),pt_metadata$sample_id),]
identical(pt_metadata$sample_id,rownames(abundances)) ### pass check

colnames(pt_metadata)
colnames(pt_metadata)[1] <- "Batch"
colnames(pt_metadata)[11] <- "first.line.therapy"

### labels
rownames(pt_metadata) <- pt_metadata$sample_id
pt_metadata$Ethnicity[pt_metadata$Ethnicity %in% ''] <- "NA"

pt_metadata$Genetics <- as.character(pt_metadata$Genetics)
pt_metadata$Genetics[pt_metadata$Genetics %in% ''] <- "NA"
pt_metadata$Genetics[pt_metadata$Genetics %in% 'ND'] <- "NA"
pt_metadata$Genetics[is.na(pt_metadata$Genetics)] <- "NA"
pt_metadata$Genetics <- as.factor(pt_metadata$Genetics)

table(pt_metadata$Study.Site)
ix <- which(pt_metadata$Study.Site %in% '')

abundances <- abundances[-ix,]
pt_metadata <- pt_metadata[-ix,]



### For Seurat Clusters
abundances_subclustering <- unclass(table(imm_atlas_meta$sample_id, imm_atlas_meta$complete_subclusters))
rownames(abundances_subclustering) <- substr(rownames(abundances_subclustering),1,11)
abundances_subclustering <- abundances_subclustering[which(rownames(abundances_subclustering) %in% pt_metadata$sample_id),]

pt_metadata_subclustering <- pt_metadata[match(rownames(abundances_subclustering),pt_metadata$sample_id),]
identical(pt_metadata_subclustering$sample_id,rownames(abundances_subclustering)) ### pass check
#############################################################################################################################################



### Analysis
#############################################################################################################################################
identical(pt_metadata$sample_id,rownames(abundances)) ### pass check
y.ab <- DGEList(t(abundances), samples=pt_metadata)

###### not req for DAA
## keep <- filterByExpr(y.ab, group=y.ab$samples$tomato)
## y.ab <- y.ab[keep,]
## summary(keep)

### DREAM
libs<-c("variancePartition","doParallel") ;lapply(libs, require, character.only = TRUE)
registerDoParallel(makeCluster(10)) ; rm(libs)

### Model form
form <- ~ (1|Study.Site) + (1|Genetics) + (1|first.line.therapy) + (1|Ethnicity) + (1|Batch) + (1|Public_ID)

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( y.ab$counts, form, pt_metadata)

### Run model
variance_vodjdream_major_clusters <- fitExtractVarPartModel(vobjDream, form, pt_metadata)

### Plot the figure
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/DAA_Major_Clusters_variance_partition.pdf",width = 3, height = 3)
plotVarPart( sortCols( variance_vodjdream_major_clusters ) ) + ggtitle( '') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
# plotPercentBars( vp[1:10,] )
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### subclusters
identical(pt_metadata_subclustering$sample_id,rownames(abundances_subclustering)) ### pass check

### Model form
form <- ~ (1|Study.Site) + (1|Genetics) + (1|first.line.therapy) + (1|Ethnicity) + (1|Batch) + (1|Public_ID)

# estimate weights using linear mixed model of dream
vobjDream_subclustering = voomWithDreamWeights( t(abundances_subclustering), form, pt_metadata_subclustering)

### Run model
variance_vodjdream_subclustering <- fitExtractVarPartModel(vobjDream_subclustering, form, pt_metadata_subclustering)

### Plot the figure
### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/DAA_SUB_Clusters_variance_partition.pdf",width = 3, height = 3)
plotVarPart( sortCols( variance_vodjdream_subclustering ) ) + ggtitle( '') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
# plotPercentBars( vp[1:10,] )

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### subclusters
identical(pt_metadata_subclustering$sample_id,rownames(abundances_subclustering)) ### pass check

form <- ~ 0 + Study.Site + (1|Genetics) + (1|first.line.therapy) + (1|Ethnicity) + (1|Batch) + (1|Public_ID)

L =  makeContrastsDream( form, pt_metadata_subclustering, 
                         contrasts = c('Study.SiteEMORY - Study.SiteMAYO', 
                                       'Study.SiteEMORY - Study.SiteMSSM' ,  
                                       'Study.SiteEMORY - Study.SiteWUSTL',
                                       'Study.SiteMAYO - Study.SiteMSSM' ,
                                       'Study.SiteMAYO - Study.SiteWUSTL',
                                       'Study.SiteMSSM - Study.SiteWUSTL',
                                       
                                       'Study.SiteEMORY - (Study.SiteMAYO+Study.SiteMSSM+Study.SiteWUSTL)/3', 
                                       'Study.SiteMSSM - (Study.SiteMAYO+Study.SiteEMORY+Study.SiteWUSTL)/3',
                                       'Study.SiteWUSTL - (Study.SiteMAYO+Study.SiteMSSM+Study.SiteEMORY)/3',
                                       'Study.SiteMAYO - (Study.SiteEMORY+Study.SiteMSSM+Study.SiteWUSTL)/3'
                                       
                                       ) )

fitmm = dream( vobjDream_subclustering, form, pt_metadata_subclustering, L)
fitmm = eBayes(fitmm)

lmfreq_results <- list() ; for ( i in 1:10 ) { lmfreq_results[[i]] <- topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(fitmm$coefficients)[i]
lmfreq_results[[i]]$Marker <- rownames( topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") ) }

dream_lmfreq_results <- do.call(rbind,lmfreq_results)
dream_lmfreq_results$nLogFDR <- -log10(dream_lmfreq_results$adj.P.Val)
dream_lmfreq_results$Comparison <- gsub('Study.Site','',dream_lmfreq_results$Comparison)

table(dream_lmfreq_results$Comparison[dream_lmfreq_results$adj.P.Val<0.05 & abs(dream_lmfreq_results$logFC) > 1])
to_test <- unique(dream_lmfreq_results$Marker[dream_lmfreq_results$adj.P.Val<0.05 & abs(dream_lmfreq_results$logFC) > 2])

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ann_rows <- pt_metadata_subclustering[,colnames(pt_metadata_subclustering) %in% c('Batch','Study.Site')]

identical(colnames(abundances_subclustering),rownames(vobjDream_subclustering))

iy <- which(colnames(abundances_subclustering) %in% to_test)

pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/SubClusters_heatmap.pdf",width = 13, height = 4)
pheatmap::pheatmap(mat = t(vobjDream_subclustering$E[,]),
                   scale = 'none',
                   #mat = log2(abundances_subclustering+1)[,iy],
                   angle_col = 90,
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   annotation_row = ann_rows,
                   #annotation_col = ann_rows#,
                   clustering_method = "ward.D2"
                   )
dev.off()
dev.off() ; dev.off() ; dev.off()

pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/SubClusters_barplot.pdf",width = 4, height = 10)
ggplot() + aes(x= as.data.frame(sort(table(imm_atlas_meta$complete_subclusters)))$Freq,
               y= as.data.frame(sort(table(imm_atlas_meta$complete_subclusters)))$Var1) + geom_bar(stat='identity') + 
               labs(x='Cell Number',y='Subcluster')
dev.off()

tmp_melted <- melt(vobjDream_subclustering$E[iy,])
colnames(tmp_melted) <- c('subcluster','sample_id','value')
tmp_melted <- merge(tmp_melted,pt_metadata_subclustering,by='sample_id')
tmp_melted <- tmp_melted[which(tmp_melted$subcluster %in% to_test),]

my_comparison <- list(c('EMORY','MAYO'), 
                      c('EMORY','MSSM'),  
                      c('EMORY','WUSTL'),
                      c('WUSTL','MSSM'),  
                      c('WUSTL','MAYO'))

pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/DAA_SubClusters_boxplots.pdf",width = 5.5, height = 5.5)
#ggplot(data=tmp_melted) + aes(x=Study.Site,y=value,fill=Study.Site) + facet_wrap(~subcluster,ncol=4) + geom_boxplot() + labs(y="Log2(Norm.Cell.Prop)") +
#  stat_compare_means(method='t.test',comparisons = my_comparison)
ggplot(data=tmp_melted) + aes(x=Study.Site,y=value,fill=Study.Site) + facet_wrap(~subcluster,ncol=4) + geom_boxplot() + labs(y="Log2(Norm.Cell.Prop)") +
  stat_compare_means(method='wilcox.test',comparisons = my_comparison) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_bw()
dev.off()

### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#############################################################################################################################################



### Classic DA
#############################################################################################################################################
### abundance
design <- model.matrix(~ factor(Study.Site) + factor(Batch), y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

summary(fit.ab$df.prior)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))

topTags(res)
#############################################################################################################################################


### compositional effects
#############################################################################################################################################
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)
#############################################################################################################################################


###++++
#immune_atlas <- readRDS('RData/WashU_Immune_atlas_prelimiary_analysis_Integration_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet.rds')
#imm_atlas_meta <- imm_atlas_meta[match(rownames(immune_atlas@meta.data),imm_atlas_meta$cell_ids),]
#identical(rownames(immune_atlas@meta.data),imm_atlas_meta$cell_ids)
#immune_atlas@meta.data$subclustering <- imm_atlas_meta$complete_subclusters
#imm_atlas_meta$complete_subclusters
saveRDS(file='RData/metadata_imm_atlas_single_cell_object.rds',imm_atlas_meta) 
rm(immune_atlas)
###++++
dev.off()

#############################################################################################################################################
#save.image('RData/surv_models_and_abbundance.RData')
load('RData/surv_models_and_abbundance.RData')
#############################################################################################################################################

#############################################################################################################################################
### 2025
#############################################################################################################################################
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/DAA_SUB_Clusters_variance_partition_ethnicity_top10.pdf",width = 5, height = 5)
plotPercentBars(variance_vodjdream_subclustering[order(variance_vodjdream_subclustering$Ethnicity,decreasing = T)[1:10],] ) + ggtitle( '') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#############################################################################################################################################
