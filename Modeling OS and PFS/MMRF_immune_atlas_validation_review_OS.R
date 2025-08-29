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
#imm_atlas_expr_483 <- readRDS(file = 'RData/new_data_for_article/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')
imm_atlas_meta_483 <- readRDS(file = 'RData/new_data_for_article/Meta_data_SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds')
ia_meta_070224 <- readRDS(file = 'RData/new_data_for_article/DATA_FOR_VALIDATION_OBJECT/per_cell_metadata/070224_per_cell_md.rds')
#############################################################################################################################################


#############################################################################################################################################
### REVIEW AND STORE DATA
#############################################################################################################################################
#imm_atlas_expr_483 (validation)
#table(imm_atlas_expr_483$cohort)
#Idents(imm_atlas_expr_483) <- imm_atlas_expr_483$cohort
#imm_atlas_expr_validation <- subset(imm_atlas_expr_483,idents='validation')
#rm(imm_atlas_expr_483)

### STORE EXPR DATA
#saveRDS(file='RData/new_data_for_article/validation_cohort_only.RDS',imm_atlas_expr_validation)
imm_atlas_expr_validation <- readRDS(file='RData/new_data_for_article/validation_cohort_only.RDS')

#table(imm_atlas_meta_483$siteXbatch) ; table(imm_atlas_meta_483$cohort)
x <- table(imm_atlas_meta_483$celltype_subclusters_integration_Yizhe_v1,imm_atlas_meta_483$cohort)
class(x) <- 'matrix'
x <- as.data.frame(x)
head(x)
#
#cor.test(x$discovery,x$validation)
#cor.test(x$discovery,x$validation,method = 'kendall')
#x$pop <- rownames(x)
#x$var <- log10(x$discovery-x$validation)/log10(x$discovery+x$validation)

#pdf(file ='figures_validation/correlation_pseudobulk_full_cohorts_discovery_validation.pdf',width = 6,height = 5)
#ggplot(data=x)+aes(x = log2(discovery) , y = log2(validation), color=var) + geom_point() + #xlim(0,16)+ylim(0,16)+
#  geom_smooth(method = 'lm') + theme_classic() +#geom_abline(color='red') +
#  geom_text_repel(data=x[log2(x$validation)<3.3,],aes(label=pop),color='black')
#dev.off()

#pdf(file ='figures_validation/correlation_pseudobulk_full_cohorts_discovery_validation_all_labels.pdf',width = 6,height = 5)
#ggplot(data=x)+aes(x = log2(discovery) , y = log2(validation), color=var) + geom_point() + #xlim(0,16)+ylim(0,16)+
#  geom_smooth(method = 'lm') + theme_classic() + #geom_abline(color='red') +
#  geom_text_repel(data=x,aes(label=pop),color='black')
#dev.off()

#imm_atlas_expr_483
#table(imm_atlas_meta_483$siteXbatch) ; table(imm_atlas_meta_483$cohort)

#EMORY_B1 EMORY_B2 EMORY_B3 EMORY_B4  MAYO_B1  MAYO_B2  MAYO_B3  MAYO_B4  MSSM_B1  MSSM_B2  MSSM_B3  MSSM_B4 WUSTL_B1 WUSTL_B2 WUSTL_B3 WUSTL_B4 
#135354    96163    69610    87854    58359    38573    39473    52654    68232    40818    61099   107420   149390    82503    90463   219307 
#discovery validation 
#1149344     247928 

#table(gsub('_SM','',ia_meta_070224$sample_id) %in% imm_atlas_meta_483$sample_id)
#head(ia_meta_070224$sample_id) ; head(imm_atlas_meta_483$sample_id)
#tail(ia_meta_070224$sample_id) ; tail(imm_atlas_meta_483$sample_id)
#ia_meta_070224$sample_id[9000] ; imm_atlas_meta_483$sample_id[9000]
#ia_meta_070224$sample_id[25000] ; imm_atlas_meta_483$sample_id[25000]
#ia_meta_070224$sample_id[105000] ; imm_atlas_meta_483$sample_id[105000]
ia_meta_070224$sample_id_old <- ia_meta_070224$sample_id
ia_meta_070224$sample_id <- imm_atlas_meta_483$sample_id

table(imm_atlas_meta_483$sample_id %in% imm_atlas_expr_validation$sample_id)
imm_atlas_meta_validation <- imm_atlas_meta_483[imm_atlas_meta_483$sample_id %in% imm_atlas_expr_validation$sample_id,]
rm(imm_atlas_meta_483)

ia_meta_070224_validation <- ia_meta_070224[ia_meta_070224$sample_id %in% imm_atlas_meta_validation$sample_id,]

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
#pdf(file ='figures_validation/correlation_subcluster_unscaled_log2_abundance_cohorts_discovery_validation.pdf',width = 15,height = 5)
#pheatmap::pheatmap(t(log2(x+1)),angle_col = 90)
#dev.off()

#2
#pdf(file ='figures_validation/correlation_subcluster_unscaled_log2_abundance_cohorts_discovery_validation.pdf',width = 6,height = 5)
#pheatmap::pheatmap(cor(x),angle_col = 90)
#dev.off()

#ia_meta_070224
#table( colnames(ia_meta_070224) %in% c('davies_based_risk', 'skerget_based_risk', 'nsd2_call', 'ccnd3_call', 
#                                'myc_call', 'mafa_call', 'ccnd1_call', 'ccnd2_call', 'maf_call', 'mafb_call', 
#                                'seqwgs_cp_hyperdiploid_call', 'seqwgs_cp_1q21_call', 'seqwgs_cp_13q14_call', 
#                                'seqwgs_cp_13q34_call', 'seqwgs_cp_17p13_call', 'seqwgs_cp_cdkn2c_call', 'seqwgs_cp_fam46c_call',
#                                'seqwgs_cp_rb1_call', 'seqwgs_cp_tp53_call', 'seqwgs_cp_traf3_call') )

#grep('SDC1',rownames(imm_atlas_expr_483@assays$RNA))
#summary(imm_atlas_expr_483@assays$RNA$counts[3240,])
#summary(imm_atlas_expr_483@assays$RNA$counts[647,])
#############################################################################################################################################


#############################################################################################################################################
### validation cohort ALONE
#############################################################################################################################################
#identical(as.character(imm_atlas_expr_validation$sample_id),imm_atlas_meta_validation$sample_id)
#saveRDS(file='RData/new_data_for_article/validation_meta_cohort_only.RDS',imm_atlas_meta_validation)
imm_atlas_meta_validation <- readRDS(file='RData/new_data_for_article/validation_meta_cohort_only.RDS')
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
#############################################################################################################################################

### store base objects
#write.csv(file = 'subcluster_discovery_baseline_fulldf.csv',subcluster_discovery_fulldf)
subcluster_discovery_fulldf <- readRDS(file = 'RData/subcluster_discovery_baseline_fulldf.rds')
#write.csv(file = 'subcluster_validation_baseline_fulldf.csv',subcluster_validation_fulldf)
subcluster_validation_fulldf <- readRDS(file = 'RData/subcluster_validation_baseline_fulldf.rds')


#############################################################################################################################################
### PACKAGES AND BASICS ON THE DATA
#############################################################################################################################################

library(rms)
library(ROCR)
library(forestmodel)
library(survival)
library(survminer)

### basics before prediction
pdf(file ='figures_validation/km_plot_d_amm_tx_asct_ever_discovery_OS.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcos, censos==1) ~ d_amm_tx_asct_ever, data = subcluster_discovery_fulldf), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "OS Prob.",pval = TRUE,risk.table = TRUE,palette = c('#F1C1A1','#FD8A8A')) +labs(title = '')
dev.off()

pdf(file ='figures_validation/km_plot_d_amm_tx_asct_ever_validation_OS.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcos, censos==1) ~ d_amm_tx_asct_1st, data = subcluster_validation_fulldf), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "OS Prob.",pval = TRUE,risk.table = TRUE,palette = c('#F1C1A1','#FD8A8A')) +labs(title = '')
dev.off()

pdf(file ='figures_validation/plot_forrest_discovery_OS.pdf',width = 8,height = 8)
forest_model( coxph(Surv(ttcos, censos==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + 
                      d_amm_tx_asct_ever + d_dx_amm_iss_stage + siteXbatch, data = subcluster_discovery_fulldf) )

dev.off()

pdf(file ='figures_validation/plot_forrest_validation_OS.pdf',width = 8,height = 8)
forest_model( coxph(Surv(ttcos, censos==1) ~ davies_based_risk + d_pt_sex + d_dx_amm_age + 
                      d_amm_tx_asct_1st + d_dx_amm_iss_stage + siteXbatch, data = subcluster_validation_fulldf) )

dev.off()

#############################################################################################################################################
### VALIDATION COHORT
#############################################################################################################################################
used_features <- readRDS(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/used_features_OS.rds')
###
library(rms)
b <- paste(used_features$R2,'+')
a <- paste("censos ~ davies_based_risk + d_dx_amm_age + d_pt_sex + site + d_dx_amm_iss_stage +")
full_form <- capture.output(cat(a,b))
full_form <- gsub(' \\+$','',full_form)

subcluster_validation_fulldf$sample_orig_id<-NULL
subcluster_validation_fulldf$DISCOVERY.doublet_pred<-NULL
subcluster_validation_fulldf$sample_id_old<-NULL

ddist <- datadist(subcluster_validation_fulldf)
options(datadist='ddist')

mod.bi <- lrm(as.formula(full_form), subcluster_validation_fulldf, x=TRUE,y=TRUE,na.action = na.delete)
predictions_list <- get_auc2(mod.bi, subcluster_validation_fulldf)  
pred_validation_list <- get_lrm_validation(mod.bi)
pred_validation_list

pred <- prediction(predictions_list$predictions, predictions_list$labels)
perf <- performance(pred,'tpr','fpr')
predictions_list$auc
#AUC:0.87

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
tmp_df_V3<- subcluster_validation_V2

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=rownames(tmp_df_V3))
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_extended_sup8t.csv',dataForRDDPFS)


pdf(file ='figures_validation/km_plot_validation_cohort_prediction_m20_OS.pdf',width = 3.5,height = 5)
ggsurvplot( fit = survfit(Surv(ttcos, censos==1) ~ p, data = subcluster_validation_V2), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "OS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','firebrick')) +labs(title = 'M20')
dev.off()

coxph(Surv(ttcos, censos==1) ~ p, data = subcluster_validation_V2) %>% gtsummary::tbl_regression(exp = TRUE)
#HR:34.2
#CI:4.7-60
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
# CRITICAL FOR ACCURATE MODEL PERFORMANCE
# Can use function m as an argument to Predict or nomogram to
# get predicted means instead of log odds or probabilities

#############################################################################################################################################
### VALIDATION COHORT CORRECTION
#############################################################################################################################################

#############################################################################################################################################
write.csv(file='chisq_top_pops.csv',k)
#############################################################################################################################################



#############################################################################################################################################
### Regression per Compartment
library(DirichletReg)

subcluster_all_fulldf <- readRDS(file = 'RData/subcluster_all_fulldf_V2.RDS')

my_res_list<-list()
Compartment <- c('NkT','BEry','Ery','Myeloid','Plasma') #,'Full'

for(iii in Compartment){
  x <- subcluster_all_fulldf
  y <- t(subcluster_all_fulldf)
  ix <- grep(paste('^',iii,sep=''), rownames(y))
  z <- x[,ix]
  
  cell_proportions = DR_data(z)
  x$site_batch <- paste(x$site,x$batch,sep='_')
  x$censos <- factor(x$censos,levels=c('0','1'))
  dr_fit_common <- DirichReg( cell_proportions ~ censos + site_batch, x, model = "common" )
  
  u = summary(dr_fit_common)
  pvals = u$coef.mat[grep('censos1', rownames(u$coef.mat), invert=FALSE), 4]
  v = names(pvals)
  
  prob.ratio = exp( summary(dr_fit_common)$coefficients[paste0(colnames(z),":censos1")] )
  
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('censos1', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  
  dr_comp_prog_res <- data.frame(log2fc = log2(exp( summary(dr_fit_common)$coefficients[paste0(colnames(z),":censos1")] )) , 
                                 pval=colMeans(pvals),
                                 Marker=colnames(z))
  rownames(dr_comp_prog_res) <- gsub(':censos1','',rownames(dr_comp_prog_res))
  
  my_res_list[[iii]] <- dr_comp_prog_res
}
dr_comp_prog_res <- do.call(rbind,my_res_list)

### 
dr_comp_prog_res$fdr <- p.adjust(dr_comp_prog_res$pval,method = 'BH')
dr_comp_prog_res$cell <- rownames(dr_comp_prog_res)
dr_comp_prog_res$type <- 'comp'


ggplot(data=dr_comp_prog_res[dr_comp_prog_res$fdr<0.05,]) + aes(x=Marker,y=type,size=-log10(fdr),color=log2fc) + geom_point()+
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red') + rotate_x_text(angle = 90)

#### STORE
#saveRDS(file='RData/dr_comp_prog_res_OS.RDS',dr_comp_prog_res)

used_features <- readRDS(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/used_features_OS.rds')
###
library(rms)
b <- paste(used_features$R2,'+')
a <- paste("censos ~ davies_based_risk + d_dx_amm_age + d_pt_sex + batch + site + d_dx_amm_iss_stage +")
full_form <- capture.output(cat(a,b))
full_form <- gsub(' \\+$','',full_form)
###
mod.bi <- lrm(as.formula(full_form), subcluster_all_fulldf, x=TRUE,y=TRUE)
###
k <- as.data.frame(as.matrix(anova(mod.bi)))
k$marker <- rownames(k)
k <- k[!k$marker %in% c('TOTAL',k$marker[1:6]),]
k$ap <- p.adjust(k$P,method = 'BH')
k <- k[order(k$P),]
#saveRDS(file='RData/model.test_OS.RDS',k)

#############################################################################################################################################

dr_all_res_prog <- readRDS(file='RData/dr_comp_prog_res_OS.RDS')
k<-readRDS(file='RData/model.test_OS.RDS')
library(tidyverse)

colnames(k)[4] <- 'Marker'
z <- merge(k,dr_all_res_prog,by='Marker',all=TRUE)
z$true <- 'ns'
z$true[z$pval<0.1 & z$log2fc > 0 ] <- 'Worse'
z$true[z$pval<0.1 & z$log2fc < 0 ] <- 'Better'
z <- z[order(z$`Chi-Square`,decreasing = T),]

pdf(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/ROC_subcluster_names_20_variables_DIRECTION_OS.pdf',width =3,height = 3.25)
ggplot(data=z)+aes(x=log10(`Chi-Square`),y=-log10(P))+geom_point(aes(color=true))+ theme_classic()+
  scale_color_manual(values = c('steelblue','grey50','firebrick')) +
  geom_text_repel(data=z[1:20,],aes(label=Marker,color=true),max.overlaps = 1000,force = 50) + geom_hline(yintercept = 1,linetype=2,color='red')+
  labs(x='Log10(Chi-Square)',y='-Log10(P)',color='Outcomes') + theme(legend.position = 'bottom')
dev.off()
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