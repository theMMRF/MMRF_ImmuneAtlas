#############################################################################################################################################
### Edgar Gonzalez-Kozlova
#############################################################################################################################################
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','monocle3',
        'limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2')
lapply(libs, require, character.only = TRUE) ; rm(libs)
#############################################################################################################################################

#############################################################################################################################################
### Work space
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/")
#############################################################################################################################################

#############################################################################################################################################
### FUNCTIONS
#############################################################################################################################################
get_auc  <- function (m,df) {
  predictions <- predict(object = m, df[,-which(colnames(df) %in% 'censpfs')],se.fit = TRUE)
  labels <- df$censpfs
  pred <- prediction(predictions, labels)
  perf <- performance(pred,'tpr','fpr')
  auc_ROCR <- performance(pred, measure = "auc")@y.values[[1]]  
  return(list(perf=perf,auc=auc_ROCR,predictions=predictions,labels=labels,pred=pred))
}
###
get_auc2  <- function (m,df) {
  predictions <- predict(object = m, df[which(!is.na(df$d_dx_amm_iss_stage)),-which(colnames(df) %in% 'censpfs')],se.fit = TRUE,type="fitted.ind")
  labels <- df$censpfs[which(!is.na(df$d_dx_amm_iss_stage))]
  pred <- prediction(predictions, labels)
  perf <- performance(pred,'tpr','fpr')
  auc_ROCR <- performance(pred, measure = "auc")@y.values[[1]]  
  return(list(perf=perf,auc=auc_ROCR,predictions=predictions,labels=labels,pred=pred))
}
get_auc3  <- function (m,df) {
  predictions <- predict(object = m, df[which(!is.na(df$d_dx_amm_iss_stage)),-which(colnames(df) %in% 'censpfs')],se.fit = TRUE,type="fitted.ind",na.action=na.exclude)
  labels <- df$censpfs[which(!is.na(df$d_dx_amm_iss_stage))]
  pred <- prediction(predictions, labels)
  perf <- performance(pred,'tpr','fpr')
  auc_ROCR <- performance(pred, measure = "auc")@y.values[[1]]  
  return(list(perf=perf,auc=auc_ROCR,predictions=predictions,labels=labels,pred=pred))
}
###
get_lrm_validation <- function(m){
  x <- validate(m,B=1000)
  class(x) <- 'matrix'
  x <- as.data.frame(x)
  x$auc <- x[1,1]/2 +0.5
  return(x)
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
#############################################################################################################################################


#############################################################################################################################################
### load data
#############################################################################################################################################
imm_atlas_expr_483 <- readRDS(file = 'RData/new_data_for_article/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')
imm_atlas_meta_483 <- readRDS(file = 'RData/new_data_for_article/Meta_data_SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')
ia_meta_070224 <- readRDS(file = 'RData/new_data_for_article/DATA_FOR_VALIDATION_OBJECT/per_cell_metadata/070224_per_cell_md.rds')
#############################################################################################################################################


#############################################################################################################################################
### REVIEW AND STORE DATA
#############################################################################################################################################
#imm_atlas_expr_483 (validation)
table(imm_atlas_expr_483$cohort)
Idents(imm_atlas_expr_483) <- imm_atlas_expr_483$cohort
imm_atlas_expr_validation <- subset(imm_atlas_expr_483,idents='validation')
rm(imm_atlas_expr_483)

### STORE EXPR DATA
saveRDS(file='RData/new_data_for_article/validation_cohort_only.RDS',imm_atlas_expr_validation)

#table(imm_atlas_meta_483$siteXbatch) ; table(imm_atlas_meta_483$cohort)
x <- table(imm_atlas_meta_483$celltype_subclusters_integration_Yizhe_v1,imm_atlas_meta_483$cohort)
class(x) <- 'matrix'
x <- as.data.frame(x)
head(x)
#
cor.test(x$discovery,x$validation)
cor.test(x$discovery,x$validation,method = 'kendall')
x$pop <- rownames(x)
x$var <- log10(x$discovery-x$validation)/log10(x$discovery+x$validation)

pdf(file ='figures_validation/correlation_pseudobulk_full_cohorts_discovery_validation.pdf',width = 6,height = 5)
ggplot(data=x)+aes(x = log2(discovery) , y = log2(validation), color=var) + geom_point() + #xlim(0,16)+ylim(0,16)+
  geom_smooth(method = 'lm') + theme_classic() +#geom_abline(color='red') +
  geom_text_repel(data=x[log2(x$validation)<3.3,],aes(label=pop),color='black')
dev.off()

pdf(file ='figures_validation/correlation_pseudobulk_full_cohorts_discovery_validation_all_labels.pdf',width = 6,height = 5)
ggplot(data=x)+aes(x = log2(discovery) , y = log2(validation), color=var) + geom_point() + #xlim(0,16)+ylim(0,16)+
  geom_smooth(method = 'lm') + theme_classic() + #geom_abline(color='red') +
  geom_text_repel(data=x,aes(label=pop),color='black')
dev.off()

#imm_atlas_expr_483
table(imm_atlas_meta_483$siteXbatch) ; table(imm_atlas_meta_483$cohort)

#EMORY_B1 EMORY_B2 EMORY_B3 EMORY_B4  MAYO_B1  MAYO_B2  MAYO_B3  MAYO_B4  MSSM_B1  MSSM_B2  MSSM_B3  MSSM_B4 WUSTL_B1 WUSTL_B2 WUSTL_B3 WUSTL_B4 
#135354    96163    69610    87854    58359    38573    39473    52654    68232    40818    61099   107420   149390    82503    90463   219307 
#discovery validation 
#1149344     247928 

table(gsub('_SM','',ia_meta_070224$sample_id) %in% imm_atlas_meta_483$sample_id)
head(ia_meta_070224$sample_id) ; head(imm_atlas_meta_483$sample_id)
tail(ia_meta_070224$sample_id) ; tail(imm_atlas_meta_483$sample_id)
ia_meta_070224$sample_id[9000] ; imm_atlas_meta_483$sample_id[9000]
ia_meta_070224$sample_id[25000] ; imm_atlas_meta_483$sample_id[25000]
ia_meta_070224$sample_id[105000] ; imm_atlas_meta_483$sample_id[105000]
ia_meta_070224$sample_id_old <- ia_meta_070224$sample_id
ia_meta_070224$sample_id <- imm_atlas_meta_483$sample_id

table(imm_atlas_meta_483$sample_id %in% imm_atlas_expr_validation$sample_id)
imm_atlas_meta_validation <- imm_atlas_meta_483[imm_atlas_meta_483$sample_id %in% imm_atlas_expr_validation$sample_id,]
rm(imm_atlas_meta_483)

ia_meta_070224_validation <- ia_meta_070224[ia_meta_070224$sample_id %in% imm_atlas_meta_validation$sample_id,]
rm(ia_meta_070224)

identical(imm_atlas_meta_validation$sample_id,ia_meta_070224_validation$sample_id)
table(colnames(imm_atlas_meta_validation) %in% colnames(ia_meta_070224_validation))
ix <- which(!colnames(imm_atlas_meta_validation) %in% colnames(ia_meta_070224_validation))
test <- cbind(imm_atlas_meta_validation[,ix],ia_meta_070224_validation)

test <- merge(imm_atlas_meta_validation,ia_meta_070224,by='sample_id')
#'censpfs','ttcpfs','censos','ttcos', 'collection_event','sample_id','d_pt_sex','public_id','ecog','d_amm_tx_asct_ever'
# davies_based_risk + d_dx_amm_age + d_amm_tx_asct_ever + d_pt_sex + d_dx_amm_iss_stage
imm_atlas_meta_validation <- test

rm(ia_meta_070224)

### review and clear; then do predictions

### exploring figures

x <- table(imm_atlas_meta_483$celltype_subclusters_integration_Yizhe_v1,paste(imm_atlas_meta_483$cohort,imm_atlas_meta_483$siteXbatch))
class(x) <- 'matrix'
x <- as.data.frame(x)
x<-type.convert(x)

#1
pdf(file ='figures_validation/correlation_subcluster_unscaled_log2_abundance_cohorts_discovery_validation.pdf',width = 15,height = 5)
pheatmap::pheatmap(t(log2(x+1)),angle_col = 90)
dev.off()

#2
pdf(file ='figures_validation/correlation_subcluster_unscaled_log2_abundance_cohorts_discovery_validation.pdf',width = 6,height = 5)
pheatmap::pheatmap(cor(x),angle_col = 90)
dev.off()

#ia_meta_070224
table( colnames(ia_meta_070224) %in% c('davies_based_risk', 'skerget_based_risk', 'nsd2_call', 'ccnd3_call', 
                                'myc_call', 'mafa_call', 'ccnd1_call', 'ccnd2_call', 'maf_call', 'mafb_call', 
                                'seqwgs_cp_hyperdiploid_call', 'seqwgs_cp_1q21_call', 'seqwgs_cp_13q14_call', 
                                'seqwgs_cp_13q34_call', 'seqwgs_cp_17p13_call', 'seqwgs_cp_cdkn2c_call', 'seqwgs_cp_fam46c_call',
                                'seqwgs_cp_rb1_call', 'seqwgs_cp_tp53_call', 'seqwgs_cp_traf3_call') )

grep('SDC1',rownames(imm_atlas_expr_483@assays$RNA))
summary(imm_atlas_expr_483@assays$RNA$counts[3240,])
summary(imm_atlas_expr_483@assays$RNA$counts[647,])
#############################################################################################################################################


#############################################################################################################################################
### validation cohort ALONE
#############################################################################################################################################
identical(as.character(imm_atlas_expr_validation$sample_id),imm_atlas_meta_validation$sample_id)
saveRDS(file='RData/new_data_for_article/validation_meta_cohort_only.RDS',imm_atlas_meta_validation)
#############################################################################################################################################


#############################################################################################################################################
### PREPARE DISCOVERY DATA FOR REPRODUCIBILITY
#############################################################################################################################################

library(rms)
library(ROCR)

### make-model
file_717_per_cell_md <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/1110_per_cell_md.rds')
load(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/clinical_metadata_IA22_Feb03_release.RData')
###
cmia22 <- clinical_metadata_IA22_Feb03_release[which(clinical_metadata_IA22_Feb03_release$public_id %in% file_717_per_cell_md$public_id),]
cmia22 <- cmia22[,colnames(cmia22) %in% c('censpfs','ttcpfs','censos','ttcos', 'collection_event','sample_id','d_pt_sex','public_id','ecog','d_amm_tx_asct_ever')]
cmia22 <- unique(cmia22)
###
my_vars <- c('public_id','davies_based_risk','visit_type','siteXbatch','progression_group','d_dx_amm_age','d_dx_amm_iss_stage','collection_event','sample_id')
###
x <- file_717_per_cell_md[,which(colnames(file_717_per_cell_md) %in% my_vars)]
rownames(x) <- NULL ; x <- unique(x)
x$site <- tidyr::separate(data.frame(x$siteXbatch),1, sep="_", c("a","b",'c','d'))$a
x$batch <- tidyr::separate(data.frame(x$siteXbatch),1, sep="_", c("a","b",'c','d'))$b
x$siteXbatch <- NULL
meta_df <- x
meta_df <- merge(meta_df,cmia22,by='public_id',all=TRUE) ; rm(cmia22)
###
ix <- which(meta_df$davies_based_risk %in% c('standard_risk','high_risk') & meta_df$collection_event %in% 'Baseline')
meta_baseline_df <- meta_df[ix,] ; rm(meta_df)

y <- table(file_717_per_cell_md$subcluster_V03072023,file_717_per_cell_md$sample_id)
class(y) <- 'matrix'
my_doublets <- unique(as.character(file_717_per_cell_md$subcluster_V03072023[file_717_per_cell_md$doublet_pred %in% c('dblet_cluster','poss_dblet_cluster')]))
y <- y[!rownames(y) %in% my_doublets,]
### remove non baseline/relevant
y <- y[,which(colnames(y) %in% meta_baseline_df$sample_id)]
cell_props_all = DR_data(t(y))

identical(meta_baseline_df$sample_id,rownames(cell_props_all))
cell_props_all <- cell_props_all[match(meta_baseline_df$sample_id,rownames(cell_props_all)),]
identical(meta_baseline_df$sample_id,rownames(cell_props_all))

### data for prediction
subcluster_all_fulldf <- cbind(meta_baseline_df,cell_props_all)
subcluster_discovery_fulldf <- subcluster_all_fulldf

### top predictors
predictive_pops=c('NkT.2.0','Myeloid.1','BEry.12','Plasma.4','NkT.2.1','Plasma.2','Myeloid.10','NkT.2.2','Myeloid.8','Myeloid.4','NkT.14')

#############################################################################################################################################
### PREPARE DATA FOR VALIDATION
#############################################################################################################################################

# <-------------- redo and test
#imm_atlas_meta_validation
#imm_atlas_expr_validation

### get the same clusters for proportions
#imm_atlas_expr_validation$predicted_doublet %in% 'True'
ix <- which(imm_atlas_expr_validation$seurat_subclusters_integration_Yizhe_v1 %in% colnames(subcluster_all_fulldf))

y <- table(imm_atlas_expr_validation$seurat_subclusters_integration_Yizhe_v1[ix],imm_atlas_expr_validation$sample_id[ix])
class(y) <- 'matrix'

table(colnames(y) %in% imm_atlas_meta_validation$sample_id)
identical(as.character(colnames(y)),as.character(imm_atlas_meta_validation$sample_id))
imm_atlas_meta_validation <- imm_atlas_meta_validation[match(as.character(colnames(y)),as.character(imm_atlas_meta_validation$sample_id)),]
identical(as.character(colnames(y)),as.character(imm_atlas_meta_validation$sample_id))

table(imm_atlas_meta_validation$collection_event)

### remove non baseline/relevant
ix <- which(imm_atlas_meta_validation$collection_event %in% 'Baseline')
imm_atlas_meta_validation_baseline <- imm_atlas_meta_validation[ix,]
y <- y[,ix]
cell_props_all = DR_data(t(y))
###
identical(as.character(colnames(y)),as.character(imm_atlas_meta_validation_baseline$sample_id))
identical(as.character(rownames(cell_props_all)),as.character(imm_atlas_meta_validation_baseline$sample_id))
rownames(imm_atlas_meta_validation_baseline) <- NULL
imm_atlas_meta_validation_baseline$cellname <- NULL
imm_atlas_meta_validation_baseline$barcode_full <- NULL
class(cell_props_all) <- 'matrix'
subcluster_validation_fulldf <- cbind(imm_atlas_meta_validation_baseline,cell_props_all)
dim(subcluster_validation_fulldf)

subcluster_validation_fulldf$d_amm_tx_asct_ever <- subcluster_validation_fulldf$d_amm_tx_asct_1st

### store base objects
write.csv(file = 'subcluster_discovery_baseline_fulldf.csv',subcluster_discovery_fulldf)
saveRDS(file = 'subcluster_discovery_baseline_fulldf.rds',subcluster_discovery_fulldf)
write.csv(file = 'subcluster_validation_baseline_fulldf.csv',subcluster_validation_fulldf)
saveRDS(file = 'subcluster_validation_baseline_fulldf.rds',subcluster_validation_fulldf)

subcluster_discovery_fulldf <- readRDS(file = 'RData/subcluster_discovery_baseline_fulldf.rds')
subcluster_validation_fulldf <- readRDS(file = 'RData/subcluster_validation_baseline_fulldf.rds')
#############################################################################################################################################
### REVIEW VARIABLE IMPORTANCE
#############################################################################################################################################

dr_all_res_prog <- readRDS(file='RData/dr_all_res_prog.RDS')
k<-readRDS(file='RData/model.test.RDS')

x <- dr_all_res_prog[,] %>% group_by(Marker,log2fc,type) %>% summarize(log2fc=median(log2fc))
test <- x %>% pivot_wider(names_from = "Marker", values_from = "log2fc")
test <- as.data.frame(test) ; rownames(test) <- test$type # test[is.na(test)]<-0 ;
test$type  <- NULL

test <- t(test)
test <- as.data.frame(test)
test$comp[test$comp > 0] <- 1
test$comp[test$comp < 0] <- -1

test$all[test$all > 0] <- 1
test$all[test$all < 0] <- -1
test$marker <- rownames(test)

z <- merge(k,test,by='marker')
z$true <- 'NA'
z$true[z$comp > 0 & z$all > 0] <- 'Worse'
z$true[z$comp < 0 & z$all < 0] <- 'Better'

z <- z[order(z$`Chi-Square`,decreasing = T),]
z
z$d.f.<-NULL
write.csv(file='RData/chisq_top_pops.csv',z)
#############################################################################################################################################
### PACKAGES AND BASICS ON THE DATA
#############################################################################################################################################

library(rms)
library(ROCR)
library(forestmodel)
library(survival)
library(survminer)

### basics before prediction
pdf(file ='figures_validation/km_plot_d_amm_tx_asct_ever_discovery.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ d_amm_tx_asct_ever, data = subcluster_discovery_fulldf), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('#F1C1A1','#FD8A8A')) +labs(title = '')
dev.off()

subcluster_validation_fulldf
pdf(file ='figures_validation/km_plot_d_amm_tx_asct_ever_validation.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ d_amm_tx_asct_1st, data = subcluster_validation_fulldf), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('#F1C1A1','#FD8A8A')) +labs(title = '')
dev.off()

tmp_df_V3 <- subcluster_validation_fulldf

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,d_amm_tx_asct_1st=tmp_df_V3$d_amm_tx_asct_1st,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_extendeddata_fig8k_pfs.csv',dataForRDDPFS)
# fit survival curve
fit <- survfit(Surv(ttcpfs, censpfs == 1) ~ d_amm_tx_asct_1st, data = tmp_df_V3)
# log-rank test (default is unweighted, i.e., Mantel–Cox)
lr <- survdiff(Surv(ttcpfs, censpfs == 1) ~ d_amm_tx_asct_1st, data = tmp_df_V3)
# extract chisq and convert to p-value
chisq <- lr$chisq
df <- length(lr$n) - 1
pval <- 1 - pchisq(chisq, df)
pval
pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval


tmp_df_V3 <- subcluster_discovery_fulldf

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,d_amm_tx_asct_1st=tmp_df_V3$d_amm_tx_asct_ever,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_extendeddata_fig8j_pfs.csv',dataForRDDPFS)
# fit survival curve
fit <- survfit(Surv(ttcpfs, censpfs == 1) ~ d_amm_tx_asct_ever, data = tmp_df_V3)
# log-rank test (default is unweighted, i.e., Mantel–Cox)
lr <- survdiff(Surv(ttcpfs, censpfs == 1) ~ d_amm_tx_asct_ever, data = tmp_df_V3)
# extract chisq and convert to p-value
chisq <- lr$chisq
df <- length(lr$n) - 1
pval <- 1 - pchisq(chisq, df)
pval
pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval



pdf(file ='figures_validation/plot_forrest_discovery.pdf',width = 8,height = 8)
forest_model( coxph(Surv(ttcpfs, censpfs==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + 
                      d_amm_tx_asct_ever + d_dx_amm_iss_stage + siteXbatch, data = subcluster_discovery_fulldf) )

dev.off()

pdf(file ='figures_validation/plot_forrest_validation.pdf',width = 8,height = 8)
forest_model( coxph(Surv(ttcpfs, censpfs==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + 
                      d_amm_tx_asct_ever + d_dx_amm_iss_stage + siteXbatch, data = subcluster_validation_fulldf) )

dev.off()

#############################################################################################################################################
### reproducing DISCOVERY COHORT
#############################################################################################################################################
### form with top pops
# full_form <- "censpfs ~ davies_based_risk + d_dx_amm_age + d_amm_tx_asct_ever + d_pt_sex + d_dx_amm_iss_stage + NkT.2.0 + 
#                        Myeloid.1 + BEry.12 + Plasma.4 + NkT.2.1 + Plasma.2 + Myeloid.10 + NkT.2.2 + Myeloid.8 + Myeloid.4 + NkT.14 +
#                        Plasma.8 + Myeloid.12 + NkT.3.2 + Ery.6 + BEry.6 + BEry.16 + NkT.13 + BEry.15.0 + Ery.8 + NkT.2.3 + NkT.5.0 +
#                        Plasma.12+Myeloid.15+Ery.0+NkT.5.1+Plasma.7+Ery.1+Myeloid.0+NkT.1.5+Myeloid.11.0+Myeloid.11.1 + Myeloid.6 + Ery.3"

#full_form <- "censpfs ~ davies_based_risk + d_dx_amm_age + d_amm_tx_asct_ever + d_pt_sex + d_dx_amm_iss_stage + batch + site + 
#              NkT.2.0 + Myeloid.1 + BEry.12 + Plasma.4 + NkT.2.1 + Plasma.2 + Myeloid.10 + NkT.2.2 + Myeloid.8 + Myeloid.4 + NkT.14"

subcluster_discovery_fulldf$siteXbatch <- paste(subcluster_discovery_fulldf$site,subcluster_discovery_fulldf$batch,sep='_')

full_form <- "censpfs ~ davies_based_risk + d_dx_amm_age + d_amm_tx_asct_ever + d_pt_sex + d_dx_amm_iss_stage + siteXbatch + 
              NkT.2.0 + Myeloid.1 + BEry.12 + Plasma.4 + NkT.2.1 + Plasma.2 + Myeloid.10 + NkT.2.2 + Myeloid.8 + Myeloid.4 + NkT.14"

# CRITICAL FOR ACCURATE MODEL PERFORMANCE
# Can use function m as an argument to Predict or nomogram to
# get predicted means instead of log odds or probabilities
ddist <- datadist(subcluster_discovery_fulldf)
options(datadist='ddist')

mod.bi <- lrm(as.formula(full_form), subcluster_discovery_fulldf, x=TRUE,y=TRUE,na.action = na.delete)
predictions_list <- get_auc2(mod.bi, subcluster_discovery_fulldf)  
pred_validation_list <- get_lrm_validation(mod.bi)
pred_validation_list

pred <- prediction(predictions_list$predictions, predictions_list$labels)
perf <- performance(pred,'tpr','fpr')
predictions_list$auc
#AUC: 0.81

plot(perf,lwd=2,main='ROC curves from 10-fold cross-validation')
dev.off()

x<- perf@x.values
y<- perf@y.values
x <- data.frame(model='discovery',prediction=unlist(x))
x$y <- unlist(y)
x$auc <- predictions_list$auc

subcluster_discovery_V2 <- subcluster_discovery_fulldf[which(!is.na(subcluster_discovery_fulldf$d_dx_amm_iss_stage)),]
subcluster_discovery_V2$p <- as.numeric(predictions_list$predictions) > median(as.numeric(predictions_list$predictions))
subcluster_discovery_V2$p[subcluster_discovery_V2$p == 'TRUE'] <- 'High'
subcluster_discovery_V2$p[subcluster_discovery_V2$p == 'FALSE'] <- 'Low'

subcluster_discovery_V2$p <- factor(subcluster_discovery_V2$p,levels = c('Low','High'))

pdf(file ='figures_validation/km_plot_discovery_cohort_prediction_m11.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_discovery_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','steelblue')) +labs(title = 'M11')
dev.off()

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_discovery_V2) %>% gtsummary::tbl_regression(exp = TRUE)
#HR: 3.84
#CI: 2.64-5.59

dev.off()

#############################################################################################################################################
### VALIDATION COHORT
#############################################################################################################################################
subcluster_validation_fulldf$sample_orig_id<-NULL
subcluster_validation_fulldf$DISCOVERY.doublet_pred<-NULL
subcluster_validation_fulldf$sample_id_old<-NULL

ddist <- datadist(subcluster_validation_fulldf)
options(datadist='ddist')

full_form <- "censpfs ~ davies_based_risk + d_dx_amm_age + d_amm_tx_asct_1st + d_pt_sex + d_dx_amm_iss_stage + siteXbatch + 
              NkT.2.0 + Myeloid.1 + BEry.12 + Plasma.4 + NkT.2.1 + Plasma.2 + Myeloid.10 + NkT.2.2 + Myeloid.8 + Myeloid.4 + NkT.14"


mod.bi <- lrm(as.formula(full_form), subcluster_validation_fulldf, x=TRUE,y=TRUE,na.action = na.delete)
predictions_list <- get_auc2(mod.bi, subcluster_validation_fulldf)  
pred_validation_list <- get_lrm_validation(mod.bi)
pred_validation_list

pred <- prediction(predictions_list$predictions, predictions_list$labels)
perf <- performance(pred,'tpr','fpr')
predictions_list$auc
#AUC:0.95

plot(perf,lwd=2,main='ROC curves from 10-fold cross-validation')
dev.off()

x <- perf@x.values
y <- perf@y.values
x <- data.frame(model='validation',prediction=unlist(x))
x$y <- unlist(y)
x$auc <- predictions_list$auc


subcluster_validation_V2 <- subcluster_validation_fulldf[which(!is.na(subcluster_validation_fulldf$d_dx_amm_iss_stage)),]
subcluster_validation_V2$p <- as.numeric(predictions_list$predictions) > median(as.numeric(predictions_list$predictions))
subcluster_validation_V2$p[subcluster_validation_V2$p == 'TRUE'] <- 'High'
subcluster_validation_V2$p[subcluster_validation_V2$p == 'FALSE'] <- 'Low'

subcluster_validation_V2$p <- factor(subcluster_validation_V2$p,levels = c('Low','High'))

tmp_df_V3 <- subcluster_validation_V2

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=rownames(tmp_df_V3))
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_main_figure_7_j.csv',dataForRDDPFS)

pdf(file ='figures_validation/km_plot_validation_cohort_prediction_m11.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_validation_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','firebrick')) +labs(title = 'M11')
dev.off()

pdf(file ='figures_validation/km_plot_validation_cohort_prediction_m11_fig7_j.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_validation_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','firebrick')) +labs(title = 'M11')
dev.off()

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_validation_V2) %>% gtsummary::tbl_regression(exp = TRUE)
#HR:9.25
#CI:4.17-20.5
dev.off()

#############################################################################################################################################
# fit survival curve
fit <- survfit(Surv(ttcpfs, censpfs == 1) ~ p, data = tmp_df_V3)

# log-rank test (default is unweighted, i.e., Mantel–Cox)
lr <- survdiff(Surv(ttcpfs, censpfs == 1) ~ p, data = tmp_df_V3)

# extract chisq and convert to p-value
chisq <- lr$chisq
df <- length(lr$n) - 1
pval <- 1 - pchisq(chisq, df)

pval

pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval

#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
### CORRECTIONS / NO PRIOR TX
#############################################################################################################################################
subcluster_discovery_fulldf <- readRDS(file='RData/subcluster_discovery_baseline_fulldf.rds')

subcluster_discovery_fulldf$siteXbatch <- paste(subcluster_discovery_fulldf$site,subcluster_discovery_fulldf$batch,sep='_')

full_form <- "censpfs ~ davies_based_risk + d_dx_amm_age + d_pt_sex + d_dx_amm_iss_stage + siteXbatch + 
              NkT.2.0 + Myeloid.1 + BEry.12 + Plasma.4 + NkT.2.1 + Plasma.2 + Myeloid.10 + NkT.2.2 + Myeloid.8 + Myeloid.4 + NkT.14"

# CRITICAL FOR ACCURATE MODEL PERFORMANCE
# Can use function m as an argument to Predict or nomogram to
# get predicted means instead of log odds or probabilities
ddist <- datadist(subcluster_discovery_fulldf)
options(datadist='ddist')

mod.bi <- lrm(as.formula(full_form), subcluster_discovery_fulldf, x=TRUE,y=TRUE,na.action = na.delete)
predictions_list <- get_auc2(mod.bi, subcluster_discovery_fulldf)  
pred_validation_list <- get_lrm_validation(mod.bi)
pred_validation_list

pred <- prediction(predictions_list$predictions, predictions_list$labels)
perf <- performance(pred,'tpr','fpr')
predictions_list$auc

plot(perf,lwd=2,main='ROC curves from 10-fold cross-validation')
dev.off()

x<- perf@x.values
y<- perf@y.values
x <- data.frame(model='discovery',prediction=unlist(x))
x$y <- unlist(y)
x$auc <- predictions_list$auc

subcluster_discovery_V2 <- subcluster_discovery_fulldf[which(!is.na(subcluster_discovery_fulldf$d_dx_amm_iss_stage)),]
subcluster_discovery_V2$p <- as.numeric(predictions_list$predictions) > median(as.numeric(predictions_list$predictions))
subcluster_discovery_V2$p[subcluster_discovery_V2$p == 'TRUE'] <- 'High'
subcluster_discovery_V2$p[subcluster_discovery_V2$p == 'FALSE'] <- 'Low'

subcluster_discovery_V2$p <- factor(subcluster_discovery_V2$p,levels = c('Low','High'))


dataForRDDPFS <- data.frame(ttcpfs=subcluster_discovery_V2$ttcpfs,censpfs=subcluster_discovery_V2$censpfs,
                            prediction=subcluster_discovery_V2$p,sample_id=subcluster_discovery_V2$sample_id)

write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FS4_ex_tx_ASCT.csv',dataForRDDPFS)


pdf(file ='figures_validation/km_plot_discovery_cohort_prediction_m11_CORRECTION.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_discovery_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','steelblue2')) +labs(title = 'M11')
dev.off()

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_discovery_V2) %>% gtsummary::tbl_regression(exp = TRUE)
#HR 4.15
#CI 2.84-6.06
dev.off()

#############################################################################################################################################
# fit survival curve
fit <- survfit(Surv(ttcpfs, censpfs == 1) ~ p, data = subcluster_discovery_V2)

# log-rank test (default is unweighted, i.e., Mantel–Cox)
lr <- survdiff(Surv(ttcpfs, censpfs == 1) ~ p, data = subcluster_discovery_V2)

# extract chisq and convert to p-value
chisq <- lr$chisq
df <- length(lr$n) - 1
pval <- 1 - pchisq(chisq, df)

pval

pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval

#############################################################################################################################################


#############################################################################################################################################
### VALIDATION COHORT CORRECTION
#############################################################################################################################################
ddist <- datadist(subcluster_validation_fulldf)
options(datadist='ddist')

mod.bi <- lrm(as.formula(full_form), subcluster_validation_fulldf, x=TRUE,y=TRUE,na.action = na.delete)
predictions_list <- get_auc2(mod.bi, subcluster_validation_fulldf)  
pred_validation_list <- get_lrm_validation(mod.bi)
pred_validation_list

pred <- prediction(predictions_list$predictions, predictions_list$labels)
perf <- performance(pred,'tpr','fpr')
predictions_list$auc
#AUC: 0.94

plot(perf,lwd=2,main='ROC curves from 10-fold cross-validation')
dev.off()

x <- perf@x.values
y <- perf@y.values
x <- data.frame(model='validation',prediction=unlist(x))
x$y <- unlist(y)
x$auc <- predictions_list$auc

subcluster_validation_V2 <- subcluster_validation_fulldf[which(!is.na(subcluster_validation_fulldf$d_dx_amm_iss_stage)),]
subcluster_validation_V2$p <- as.numeric(predictions_list$predictions) > median(as.numeric(predictions_list$predictions))
subcluster_validation_V2$p[subcluster_validation_V2$p == 'TRUE'] <- 'High'
subcluster_validation_V2$p[subcluster_validation_V2$p == 'FALSE'] <- 'Low'

subcluster_validation_V2$p <- factor(subcluster_validation_V2$p,levels = c('Low','High'))

pdf(file ='figures_validation/km_plot_validation_cohort_prediction_m11_CORRECTION.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_validation_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','firebrick2')) +labs(title = 'M11')
dev.off()

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = subcluster_validation_V2) %>% gtsummary::tbl_regression(exp = TRUE)
#HR:9.79
#CI:4.23-22.7
dev.off()
#############################################################################################################################################

#############################################################################################################################################
write.csv(file='chisq_top_pops.csv',k)
#############################################################################################################################################


#OTHERS
#forest_model( coxph(Surv(ttcpfs, censpfs==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + d_amm_tx_asct_ever + d_dx_amm_iss_stage + site + batch, data = subcluster_all_fulldf_V2) )
#forest_model( coxph(Surv(ttcpfs, censpfs==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + d_amm_tx_asct_ever + site + batch, data = subcluster_all_fulldf_V2) )
#forest_model( coxph(Surv(ttcpfs, censpfs==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + d_amm_tx_asct_ever + site, data = subcluster_all_fulldf_V2) )
#ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ d_amm_tx_asct_ever, data = subcluster_all_fulldf_V2), 
#            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
#            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('#F1C1A1','#FD8A8A')) +labs(title = 'Myeloid.8 [V5]')

#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################