#############################################################################################################################################
### Edgar Gonzalez-Kozlova
#############################################################################################################################################
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','ggbeeswarm','clusterProfiler',
        'monocle3','limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2','ComplexHeatmap','survRM2',
        "survminer","survival","ggalluvial",'gtsummary','variancePartition','DirichletReg','cowsay','emmeans','forestmodel','ROCR','rms')
lapply(libs, require, character.only = TRUE) ; rm(libs)
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','monocle3',
        'limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2')
lapply(libs, require, character.only = TRUE) ; rm(libs)
#############################################################################################################################################

### data for article

### (1)
#Main Figure 7, Extended Data 8, Extended Data 9 (Old: Supplemental figures 34-38)
#Required: Data frame of patient clinical covariates and abundances which serve an input to the risk stratification model. Needs to have input parameters (Age, Site, ISS, Sex, Batch, Cytogenetics (Davies for Main, IMWG 2024 HR for Supplemental), and pfs/os censoring covariates (ttcpfs and censpfs, or whichever variables were used in the section). Code on the github should indicate how individual panels are derived.
                                                                                                                                                          
#
final_discovery_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_discovery_cohort.RDS')
#write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/Fig7_IA_final_discovery_cohort.csv',final_discovery_cohort)
#
final_validation_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_validation_cohort.RDS')
#write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/Fig7_IA_final_validation_cohort.csv',final_validation_cohort)

### (2)
#ED8a-c, ED8q-r (Supplemental 34a-c, Supplemental 37b-c):
#  Tables for bootstrapped AUC values per model/compartment, for OS and PFS.

### 
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','ggbeeswarm','clusterProfiler',
        'monocle3','limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2','ComplexHeatmap','survRM2',
        "survminer","survival","ggalluvial",'gtsummary','variancePartition','DirichletReg','cowsay','emmeans','forestmodel','ROCR','rms')
lapply(libs, require, character.only = TRUE) ; rm(libs)
### 
setwd("/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/")
###

key_predictions_PFS_list <- readRDS(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/key_predictions_PFS_list.RDS')
load(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/predictions_lists_os_pfs.RData')

pred <- prediction(key_predictions_PFS_list$predictions, key_predictions_PFS_list$labels)
perf <- performance(pred,'tpr','fpr')
#plot(perf,lwd=2,main='ROC curves from 10-fold cross-validation')
x<- perf@x.values
y<- perf@y.values
names(x) <- names(key_predictions_PFS_list$predictions)
names(y) <- names(key_predictions_PFS_list$predictions)
for(i in 1:length(x)){ x[[i]]<- data.frame(model=names(x)[i],x=x[[i]])}
x <- do.call(rbind,x)
x$y <- unlist(y)
rownames(x)<-NULL
### redo aucs
aucs <- list()
for(i in 1:length(predictions_list)) { aucs[[i]] <- data.frame(auc=predictions_list[[i]]$auc,model=names(predictions_list)[i]) }
aucs <- do.call(rbind,aucs)

x <- merge(x,aucs,by='model') ; rm(aucs)

x$type <- x$model
x$type[grep('cyto',x$model)] <- 'cyto'
x$type[grep('NC',x$model)] <- 'NC'
x$type[grep('V1',x$model)] <- 'Single'
x$type[grep('V2',x$model)] <- '+Age+Sex'
x$type[grep('V3',x$model)] <- '+Age+Sex+B*Site'
x$type[grep('V4',x$model)] <- '+Age+Sex+B*Site+ISS+ASTC'
x$type[grep('V5',x$model)] <- '+Cyto+Age+Sex+B*Site+ISS+ASTC'
x$type[grep('V6',x$model)] <- '+NC+Age+Sex+B*Site+ISS+ASTC'

x$type <- factor(x$type,levels = c('cyto','NC','Single','+Age+Sex','+Age+Sex+B*Site','+Age+Sex+B*Site+ISS+ASTC',
                                   '+Cyto+Age+Sex+B*Site+ISS+ASTC','+NC+Age+Sex+B*Site+ISS+ASTC',
                                   "Feat_Sel","Feat_Sel_R1","Feat_Sel_R2","Feat_Sel_R3","Feat_Sel_R4","mod_asct" ))

write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/Fig7_IA_ED8a-c.ED8q-r.AUC.values.per.model.PSF.csv',x)

mypfs <- x


### Simple Figure with key models (PFS)

#### Demographics alone + technical correction (PFS)

#############################################################################################################################################
#tmp_df_V3 <- subcluster_all_fulldf[which(!is.na(subcluster_all_fulldf$d_dx_amm_iss_stage)),]
tmp_df_V3 <- final_discovery_cohort
tmp_df_V3$p <- as.numeric(predictions_list$mod_iss_V5$predictions) > median(as.numeric(predictions_list$mod_iss_V5$predictions))
tmp_df_V3$p[tmp_df_V3$p %in% 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p %in% 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p,levels = c('Low','High'))
#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_demographics_pfs.csv',dataForRDDPFS)

pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_demographics_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#64B200')) +labs(title = 'Demographics') #Age+Sex+ISS+ASCT
dev.off()

#############################################################################################################################################
#### IMWG_2024_HR_pfs
#############################################################################################################################################
tmp_df_V3$p <- tmp_df_V3$IMWG_2024_HR
tmp_df_V3$p[tmp_df_V3$p %in% 'high_risk'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p %in% 'standard_risk'] <- 'Std'
tmp_df_V3$p <- factor(tmp_df_V3$p,levels = c('Std','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,IMWG_2024_HR=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_IMWG_2024_HR_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_IMWG_2024_HR_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#00BD5C')) +labs(title = 'IMWG 2024 HR')
dev.off()
#############################################################################################################################################

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
#############################################################################################################################################

#############################################################################################################################################
#### IA alone
#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Myeloid.8_V6) > median(as.numeric(key_predictions_PFS_list$predictions$Myeloid.8_V6))
tmp_df_V3$p[tmp_df_V3$p %in% 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p %in% 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p,levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_IA_w_covariates_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_IA_w_covariates_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#00C1A7')) +labs(title = 'IA_w_covariates')
dev.off()
#############################################################################################################################################

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
#############################################################################################################################################

#### FSall

#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Feat_Sel) > median(as.numeric(key_predictions_PFS_list$predictions$Feat_Sel))
tmp_df_V3$p[tmp_df_V3$p %in% 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p %in% 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p,levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FSall_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_FSall_covariates_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#EF67EB')) +labs(title = 'FSall')
dev.off()
#############################################################################################################################################

pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval

#### IA+CYTO+CLINICAL (feat selection 1)

#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R1) > median(as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R1))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FS1_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_FS1_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#00C1A7')) +labs(title = 'FS1')
dev.off()
#############################################################################################################################################

#### IA+CYTO+CLINICAL (feat selection 2)

#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R2) > median(as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R2))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FS2_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_FS2_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#00C1A7')) +labs(title = 'FS2')
dev.off()

#############################################################################################################################################

#### IA+CYTO+CLINICAL (feat selection 3)

#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R3) > median(as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R3))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FS3_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_FS3_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#EF67EB')) +labs(title = 'IA+CYTO+PFS')
dev.off()
#############################################################################################################################################

#### IA+CYTO+CLINICAL (feat selection 4)

#############################################################################################################################################
tmp_df_V3$p <- as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R4) > median(as.numeric(key_predictions_PFS_list$predictions$Feat_Sel_R4))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=tmp_df_V3$sample_id)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_discovery_FS4_pfs.csv',dataForRDDPFS)

#############################################################################################################################################
coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)
#############################################################################################################################################
pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/revised_discovery_FS4_pfs.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#EF67EB')) +labs(title = 'IA+CYTO+PFS')
dev.off()
#############################################################################################################################################

#############################################################################################################################################
pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval

# fit survival curve
fit <- survfit(Surv(ttcpfs, censpfs == 1) ~ p, data = tmp_df_V3)

# log-rank test (default is unweighted, i.e., Mantel–Cox)
lr <- survdiff(Surv(ttcpfs, censpfs == 1) ~ p, data = tmp_df_V3)

# extract chisq and convert to p-value
chisq <- lr$chisq
df <- length(lr$n) - 1
pval <- 1 - pchisq(chisq, df)

pval
#############################################################################################################################################

load(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/predictions_lists_os_pfs.RData')
#predictions_list_OS,pred_validation_list_OS,predictions_list,pred_validation_list

x<-list()
for(i in 1:length(predictions_list_OS)){
  x[[i]] <- data.frame(name=names(predictions_list_OS)[i],auc=predictions_list_OS[[i]]$auc) 
}
x<-do.call(rbind,x)
x <- x[order(x$auc),]
tail(x[grep('V1',x$name),])

x$model <- x$name

x$type <- x$model
x$type[grep('cyto',x$model)] <- 'cyto'
x$type[grep('NC',x$model)] <- 'NC'
x$type[grep('V1',x$model)] <- 'Single'
x$type[grep('V2',x$model)] <- '+Age+Sex'
x$type[grep('V3',x$model)] <- '+Age+Sex+B*Site'
x$type[grep('V4',x$model)] <- '+Age+Sex+B*Site+ISS+ASTC'
x$type[grep('V5',x$model)] <- '+Cyto+Age+Sex+B*Site+ISS+ASTC'
x$type[grep('V6',x$model)] <- '+NC+Age+Sex+B*Site+ISS+ASTC'

x$type <- factor(x$type,levels = c('cyto','NC','Single','+Age+Sex','+Age+Sex+B*Site','+Age+Sex+B*Site+ISS+ASTC',
                                   '+Cyto+Age+Sex+B*Site+ISS+ASTC','+NC+Age+Sex+B*Site+ISS+ASTC',
                                   "Feat_Sel","Feat_Sel_R1","Feat_Sel_R2","Feat_Sel_R3","Feat_Sel_R4","mod_asct" ))

head(mypfs)

length(unique(x$model))

y <- x[,c('model','auc','type')]

write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/Fig7_IA_ED8a-c.ED8q-r.AUC.values.per.model.OS.csv',y)

### (3)
#Main Figure 7h, ED8i, ED9h (Main Figure 7h, Supplemental Figure 35f, Supplemental Figure 38h): “Top Predictive Subclusters”
#Underlying data frame for this plot (Chi-Square, P value, whichever metric from the model is used to define ‘Better’ or ‘Worse’ outcomes)

# already exists

### (4)

#Former Supplemental 34D
# from suppl tables.

### For allessandro

### (1)
#Main Figure 7, Extended Data 8, Extended Data 9 (Old: Supplemental figures 34-38)
#Required: Data frame of patient clinical covariates and abundances which serve an input to the risk stratification model. Needs to have input parameters (Age, Site, ISS, Sex, Batch, Cytogenetics (Davies for Main, IMWG 2024 HR for Supplemental), and pfs/os censoring covariates (ttcpfs and censpfs, or whichever variables were used in the section). Code on the github should indicate how individual panels are derived.

#
final_discovery_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_discovery_cohort.RDS')
#
final_validation_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_validation_cohort.RDS')


fdc <- final_discovery_cohort[,16:104]
fvc <- final_validation_cohort[,1:106]

table(colnames(fvc) %in% colnames(fdc))
fvc <- fvc[,which(colnames(fvc) %in% colnames(fdc))]

table( colnames(fdc) %in% colnames(fvc) )

f_all <- rbind(fdc,fvc)

to_correct <- c(paste(final_discovery_cohort$batch,final_discovery_cohort$site),paste(final_validation_cohort$batch,final_validation_cohort$site))


##BiocManager::install("PLSDAbatch")
library(PLSDAbatch)

f_all_corr <- PLSDA_batch(X = f_all, batch = to_correct)
f_all_corr <- PLSDA_batch(X = f_all, 
                         #Y.trt = ad.trt, 
                         Y.bat = to_correct,
                         ncomp.trt = 1, 
                         ncomp.bat = 10)
#summary(rowSums(f_all_corr$X.nobatch))
#crumblr::crumblr(fdcs)

write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_discovery_validation_CorrectedBatchSite.csv',f_all_corr$X.nobatch)
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_discovery_validation_NoCorrection.csv',f_all)

