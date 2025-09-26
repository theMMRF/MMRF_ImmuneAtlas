#############################################################################################################################################
### Edgar Gonzalez-Kozlova
#############################################################################################################################################
libs<-c('scater','loomR','Seurat','patchwork','SeuratDisk','monocle3','ggpubr','cowplot','ggplot2','monocle3',
        'limma','edgeR','fitdistrplus','factoextra','ggrepel','tidyverse','pheatmap','reshape2')
lapply(libs, require, character.only = TRUE) ; rm(libs)
#############################################################################################################################################


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
###
get_lrm_validation <- function(m){
  x <- validate(m,B=1000)
  class(x) <- 'matrix'
  x <- as.data.frame(x)
  x$auc <- x[1,1]/2 +0.5
  return(x)
}
#############################################################################################################################################


#############################################################################################################################################
setwd('/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/')
#############################################################################################################################################


#############################################################################################################################################
final_discovery_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_discovery_cohort.RDS')
final_validation_cohort <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/final_validation_cohort.RDS')
key_predictions_PFS_list <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/key_predictions_PFS_list.RDS')
#############################################################################################################################################


#############################################################################################################################################
my_subcluster_names <- colnames(final_discovery_cohort)[16:104]
x <- list()
for(iii in my_subcluster_names){
  x[[iii]]<- length(table(final_discovery_cohort[[iii]]))
}
my_subcluster_names <- names(sort(unlist(x))[sort(unlist(x)) > 50])
###
ix <- which(colnames(final_discovery_cohort) %in% my_subcluster_names)
###
my_subcluster_names <- names(colSums(final_discovery_cohort[,ix])[colSums(final_discovery_cohort[,ix]) >1])
#############################################################################################################################################


#############################################################################################################################################
clinical_vars <- c("d_dx_amm_age", "d_pt_sex",  "batch" , "site" , "d_dx_amm_iss_stage" , "d_amm_tx_asct_1st")

all_variables <- c(clinical_vars,my_subcluster_names,'censpfs','IMWG_2024_HR','censos','ttcos','ttcpfs')

final_discovery_cohort <- final_discovery_cohort[,which(colnames(final_discovery_cohort) %in% all_variables)]
final_validation_cohort <- final_validation_cohort[,which(colnames(final_validation_cohort) %in% all_variables)]

final_discovery_cohort <- final_discovery_cohort[,order(colnames(final_discovery_cohort))]
final_validation_cohort <- final_validation_cohort[,order(colnames(final_validation_cohort))]
#############################################################################################################################################


#############################################################################################################################################
ddist <- datadist(final_discovery_cohort)
options(datadist='ddist')
#############################################################################################################################################


### formula
#############################################################################################################################################
clinical_vars <- c("d_dx_amm_age", "d_pt_sex",  "batch" , "site" , "d_dx_amm_iss_stage" , "d_amm_tx_asct_1st")
c <- paste( my_subcluster_names[!my_subcluster_names %in% c('')] ,'+') # AI
b <- paste( clinical_vars[!clinical_vars %in% c('')] ,'+') # clinical
a <- paste("censpfs ~ IMWG_2024_HR + ")
full_form <- capture.output(cat(a,b,c))
full_form <- gsub(' \\+$','',full_form)
#############################################################################################################################################


### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################


### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################


### optimized model
#############################################################################################################################################
z <- key_predictions_PFS_list$used_features$R1
z <- paste( z[!z %in% c('')] ,'+') # clinical
full_form <- capture.output(cat(a,b,z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################



### optimized model 0.65
#############################################################################################################################################
z <- key_predictions_PFS_list$used_features$R4
z <- paste( z[!z %in% c('')] ,'+') # clinical
full_form <- capture.output(cat(a,b,z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################



### optimized model 0.65 (no ASCT)
#############################################################################################################################################
z <- key_predictions_PFS_list$used_features$R4
z <- paste( z[!z %in% c('')] ,'+') # clinical
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################


as.data.frame(as.matrix(anova(mod.bi)))


### key variables
#############################################################################################################################################
x<-list()
for(i in 1:length(key_predictions_PFS_list$aucs)){
  x[[i]] <- data.frame(name=names(key_predictions_PFS_list$aucs)[i],auc=as.numeric(key_predictions_PFS_list$aucs[i])) 
}
x<-do.call(rbind,x)
x <- x[order(x$auc),]
#############################################################################################################################################



### optimized model 0.65 (no ASCT) + straglers 
#############################################################################################################################################
z <- key_predictions_PFS_list$used_features$R4
z <- paste( z[!z %in% c('Plasma.2')] ,'+') # clinical
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################



#############################################################################################################################################
library(resample)
ix <- which(colnames(final_validation_cohort) %in% my_subcluster_names)
#colVars(as.matrix(final_validation_cohort[,ix]))
### formula
#############################################################################################################################################
clinical_vars <- c("d_dx_amm_age", "d_pt_sex", "d_dx_amm_iss_stage")
c <- paste( my_subcluster_names[!my_subcluster_names %in% c(colnames(final_validation_cohort)[ix][order(colVars(as.matrix(final_validation_cohort[,ix])),decreasing = F)][1:20])] ,'+') # AI
b <- paste( clinical_vars[!clinical_vars %in% c('')] ,'+') # clinical
a <- paste("censpfs ~ IMWG_2024_HR + ")
full_form <- capture.output(cat(a,b,c))
full_form <- gsub(' \\+$','',full_form)
#############################################################################################################################################
### Feature Selection Modeling
#############################################################################################################################################
#subcluster_all_fulldf <- type_convert(subcluster_all_fulldf)
mod.bi <- lrm(as.formula(full_form), final_validation_cohort, x=TRUE,y=TRUE)
#vif(mod.bi)
names(vif(mod.bi))[vif(mod.bi) > 5]
#get_auc2(mod.bi,final_discovery_cohort)  
#get_lrm_validation(mod.bi)
as.data.frame(as.matrix(anova(mod.bi)))
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################



### elastic net # discovery
#############################################################################################################################################
library(glmnet)

# Prepare the data matrix
x <- model.matrix(full_formula, data = final_discovery_cohort)[, -1]  # Exclude intercept
y <- final_discovery_cohort$censpfs

# Apply LASSO using cross-validation
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)  # Lasso regression

# Best lambda value
best_lambda <- lasso_model$lambda.min

# Fit the final LASSO model
final_lasso <- glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)

# Coefficients of selected variables
#coef(final_lasso)
discovery_glm_vars <- names(final_lasso$beta[,1])[final_lasso$beta[,1] > 0]
#############################################################################################################################################



### elastic net # validation
#############################################################################################################################################
library(glmnet)

### formula
#############################################################################################################################################
clinical_vars <- c("d_dx_amm_age", "d_pt_sex", "d_dx_amm_iss_stage")
c <- paste( my_subcluster_names[!my_subcluster_names %in% c(colnames(final_validation_cohort)[ix][order(colVars(as.matrix(final_validation_cohort[,ix])),decreasing = F)][1:10])] ,'+') # AI
b <- paste( clinical_vars[!clinical_vars %in% c('')] ,'+') # clinical
a <- paste("censpfs ~ IMWG_2024_HR + ")
full_form <- capture.output(cat(a,b,c))
full_form <- gsub(' \\+$','',full_form)

qwrt <- final_validation_cohort[,!colnames(final_validation_cohort) %in% c('site','batch')]
# Prepare the data matrix
x <- model.matrix(as.formula(full_form), data =  qwrt)[, -1]  # Exclude intercept
y <- final_validation_cohort$censpfs

# Apply LASSO using cross-validation
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)  # Lasso regression

# Best lambda value
best_lambda <- lasso_model$lambda.min

# Fit the final LASSO model
final_lasso <- glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)

# Coefficients of selected variables
coef(final_lasso)
validation_glm_vars <- names(final_lasso$beta[,1])[final_lasso$beta[,1] > 0 | final_lasso$beta[,1] < 0]
#############################################################################################################################################



### TEST validation alone
### formula
#############################################################################################################################################
z <- c(key_predictions_PFS_list$used_features$R4,"BEry.2","BEry.10","Myeloid.6","NkT.8","NkT.13")
z <- paste( z[!z %in% c('Plasma.2','BEry.9')] ,'+') # clinical
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_validation_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### Feature Selection Modeling
#############################################################################################################################################
#subcluster_all_fulldf <- type_convert(subcluster_all_fulldf)
mod.bi <- lrm(as.formula(full_form), final_validation_cohort, x=TRUE,y=TRUE)
#vif(mod.bi)
names(vif(mod.bi))[vif(mod.bi) > 5]
#get_auc2(mod.bi,final_discovery_cohort)  
#get_lrm_validation(mod.bi)
as.data.frame(as.matrix(anova(mod.bi)))
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_validation_cohort, x=TRUE,y=TRUE)
### train
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################

#############################################################################################################################################
tmp_df_V3 <- final_validation_cohort
tmp_df_V3$p <- as.numeric(predictions) > median(as.numeric(predictions))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

dataForRDDPFS <- data.frame(ttcpfs=tmp_df_V3$ttcpfs,censpfs=tmp_df_V3$censpfs,prediction=tmp_df_V3$p,sample_id=rownames(tmp_df_V3))
write.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/revised_validation_FS4_pfs.csv',dataForRDDPFS)

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)

pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/km_labels_model_IA+CYTO+CLIN_PFS_FS4_validation_test.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','firebrick2')) +labs(title = 'IA+CYTO+PFS')
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

pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
pval

#############################################################################################################################################


### optimized model 0.65 (no ASCT) + straglers 
#############################################################################################################################################
z <- c(key_predictions_PFS_list$used_features$R4,"BEry.2","BEry.10","Myeloid.6","NkT.8","NkT.13")
z <- paste( z[!z %in% c('Plasma.2','BEry.9')] ,'+') # clinical
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################


#############################################################################################################################################
tmp_df_V3 <- final_validation_cohort
tmp_df_V3$p <- as.numeric(predictions) > median(as.numeric(predictions))
tmp_df_V3$p[tmp_df_V3$p == 'TRUE'] <- 'High'
tmp_df_V3$p[tmp_df_V3$p == 'FALSE'] <- 'Low'
tmp_df_V3$p <- factor(tmp_df_V3$p, levels = c('Low','High'))

coxph(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3) %>% gtsummary::tbl_regression(exp = TRUE)

pdf(file = '/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/km_labels_model_IA+CYTO+CLIN_PFS_FS4_validation.pdf',width = 3.85, height = 4.95)
ggsurvplot( fit = survfit(Surv(ttcpfs, censpfs==1) ~ p, data = tmp_df_V3), 
            type=c("kaplan-meier"),surv.median.line = 'hv',log.rank.weights= "survdiff",pval.method = TRUE,
            xlab = "Days", ylab = "PFS Prob.",pval = TRUE,risk.table = TRUE,palette = c('black','#00C1A7')) +labs(title = 'IA+CYTO+PFS')
dev.off()
#############################################################################################################################################




### TEST
#############################################################################################################################################
z <- c(key_predictions_PFS_list$used_features$R4,"BEry.2","BEry.10","Myeloid.6","NkT.8","NkT.13")
z <- paste( z[!z %in% c('Plasma.2','BEry.9')] ,'+') # clinical
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_validation_cohort, x=TRUE,y=TRUE)
### train
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################

### TEST
#############################################################################################################################################
z <- c(key_predictions_PFS_list$used_features$R4,"BEry.2","BEry.10","Myeloid.6","NkT.8",
       "NkT.13",
       'NkT.2.0','Myeloid.1','NkT.2.1','NkT.2.2','Myeloid.8','Myeloid.4','NkT.14')
z<-unique(z)
z <- paste( z[!z %in% c('Plasma.2','BEry.9')] ,'+') # clinical
#
full_form <- capture.output(cat(a,paste( clinical_vars[!clinical_vars %in% c('d_amm_tx_asct_1st','d_pt_sex','batch','site')] ,'+'),z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
#############################################################################################################################################
### model
#############################################################################################################################################
mod.bi <- lrm(formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)
### train
predictions <- predict(object = mod.bi, 
                       final_discovery_cohort[,-which(colnames(final_discovery_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_discovery_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################
### test
#############################################################################################################################################
predictions <- predict(object = mod.bi, 
                       final_validation_cohort[,-which(colnames(final_validation_cohort) %in% 'censpfs')],
                       se.fit = TRUE,type="fitted.ind")

labels <- final_validation_cohort$censpfs
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')
performance(pred, measure = "auc")@y.values[[1]]  
#############################################################################################################################################



#############################################################################################################################################
updated_meta_per_cell <- readRDS(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/updateMMRF_yered_william/091524_per_cell_md.rds')

updated_meta_per_cell <- updated_meta_per_cell[,c('seurat_subclusters_label_transferring_Yizhe_v1','celltype_subclusters_label_transferring_Yizhe_v1')]
rownames(updated_meta_per_cell) <- NULL
updated_meta_per_cell <- unique(updated_meta_per_cell)
updated_meta_per_cell <- as.data.frame(table(updated_meta_per_cell$seurat_subclusters_label_transferring_Yizhe_v1,updated_meta_per_cell$celltype_subclusters_label_transferring_Yizhe_v1))
updated_meta_per_cell <- updated_meta_per_cell[updated_meta_per_cell$Freq>0,]
updated_meta_per_cell$Freq<-NULL
colnames(updated_meta_per_cell) <- c('cluster','name')

x <- as.data.frame(coef(final_lasso))
x$cluster <- rownames(x)

x <- merge(x,updated_meta_per_cell,by='cluster')
x$rank <- order(x$s0)
colnames(x)[2] <- 'coef'

ggplot(data=x)+aes(x=coef,y=rank)+geom_point()

z <- my_subcluster_names
z <- paste( z[!z %in% c('')] ,'+') # clinical
full_form <- capture.output(cat(a,b,z))
full_form <- gsub(' \\+$','',full_form)
mod.bi <- lrm(as.formula(full_form), final_discovery_cohort, x=TRUE,y=TRUE)

y <- as.data.frame(as.matrix(anova(mod.bi)))
colnames(y) <- c('chi','df','p')
y$cluster <- rownames(y)
y <- merge(y,updated_meta_per_cell,by='cluster')

z <- c(key_predictions_PFS_list$used_features$R4,"BEry.2","BEry.10","Myeloid.6","NkT.8","NkT.13")
z <- z[!z %in% c('Plasma.2','BEry.9')] 

y$key <- NA
y$key[y$name %in% as.character(x$name)[!x$coef == 0]] <- 'x'
y$key[y$cluster %in% z] <- 'x'



y$dir <- 'NA'
y$dir[y$cluster %in% z$marker[z$true %in% 'Better']] <- 'Better'
y$dir[y$cluster %in% z$marker[z$true %in% 'Worse']] <- 'Worse'

ggplot(data=y)+aes(x=chi,y=-log10(p))+
  theme_classic()+
  geom_point() + 
  geom_hline(yintercept = 1.3,color='red',linetype=2)+
  geom_text_repel(data=y[y$key %in% 'x',],aes(label=name))

#y$cluster %in% z$
setwd('')
z <- read.csv(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/RData/chisq_top_pops.csv')
z$cluster <- z$marker

z <- merge(z,x,by='cluster')

z$true[ -log10(z$P) <1 ]   <- 'NA'
z$name <-
z$name[ -log10(z$P) <1 ]

pdf(file='/Users/gonzae34/Documents/projects_gnjatic/MMRF_projects/immune_atlas/figures/ROC_subcluster_names_NEW.pdf',width =3,height = 3.25)
ggplot(data=z)+aes(x=log10(Chi.Square),y=-log10(P))+geom_point(aes(color=true))+ theme_classic()+
  scale_color_manual(values = c('steelblue','firebrick','grey50')) +
  geom_text_repel(data=z[1:11,],aes(label=name,color=true),max.overlaps = 1000,force = 50) + geom_hline(yintercept = 1,linetype=2,color='red')+
  labs(x='Log10(Chi-Square)',y='-Log10(P)',color='Outcomes') + theme(legend.position = 'bottom')
dev.off()
#############################################################################################################################################


