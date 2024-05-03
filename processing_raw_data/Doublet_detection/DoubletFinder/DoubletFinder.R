###############################
###doublet detection from doubletFinder
###############################
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

library(dplyr)
library(DoubletFinder)
library(Seurat)

# create a df to summarize prediction results
DF_total <- as.data.frame(matrix(nrow=1,ncol=4))
colnames(DF_total) <- c("DF_score","DF_doublet","BC","sample_id")
for (sample_id in sample_ids){  
  sobj <- readRDS(file)
  # define the expected number of doublet cellscells.
  nExp <- round(ncol(sobj) * 0.01)  # expect 1% doublets based on scrublet
  data.filt <- doubletFinder_v3(sobj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25,sct=TRUE)
  # name of the DF prediction can change, so extract the correct column name.
  DF.names = colnames(data.filt@meta.data)[grepl("DF.classification|pANN", colnames(data.filt@meta.data))]
  doublet_pred <- data.filt@meta.data[,DF.names]
  doublet_pred$BC <- paste0(rownames(doublet_pred),"_",sample_id)
  doublet_pred$sample_id <- sample_id
  colnames(doublet_pred) <- c("DF_score","DF_doublet","BC","sample_id")
  DF_total <- rbind(DF_total,doublet_pred)
}

DF_total <- DF_total[complete.cases(DF_total),]
out_path <- "path of output folder"
write.table(DF_total,paste0(out_path,"/merged_doublet_finder.txt"),sep="\t",col.names=T,row.names=F,quote=F)

