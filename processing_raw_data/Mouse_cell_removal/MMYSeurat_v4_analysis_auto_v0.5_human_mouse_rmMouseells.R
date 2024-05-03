#!/usr/bin/env Rscript --vanilla
#Ding Lab Seurat Processing Script
#Modified by Lijun Yao, adapted from Reyka Jayasinghe

# load required libraries
library(optparse)
library(Matrix)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(devtools)
library(DropletUtils)
library(stringr)
library(DoubletFinder)
library(gridExtra)
library(viridis)
library(tibble)

# create user options
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("--pre_filter"),
              type="integer",
              default=300,
              help="min number of reads per cell to prefilter",
              metavar="integer"),
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--nfeature_max"),
              type="integer",
              default=50000,
              help="nFeature_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_max"),
              type="integer",
              default=10000,
              help="nCount_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--mito_max"),
              type="double",
              default=20,
              help="maximum allowed mitochondrial percent",
              metavar="double"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character"),
  make_option(c("--pc_num"),
              type="integer",
              default=30,
              help="number of principal components to use",
              metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}

# read in initial arguments
sample_id <- opt$sample_id
sample<-opt$sample_id
out_path <- opt$output
matrix_dir = opt$input

analysis_dir= gsub("raw_feature_bc_matrix","analysis",matrix_dir)
# make output  make output directory
dir.create(out_path) 

out_path<-paste(out_path,"/",sep="")
# get direct paths to data
data <- Read10X(data.dir = paste(matrix_dir))

#EMPTY DROPLET DETECTION USING EMPTY DROPS
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y
# Identify likely cell-containing droplets.
#out <- emptyDrops(data)

#is.cell <- out$FDR <= 0.01 #Might need to mess with this cutoff
#sum(is.cell, na.rm=TRUE)

# Check if p-values are lower-bounded by 'niters'
# (increase 'niters' if any Limited==TRUE and Sig==FALSE)
#table(Sig=is.cell, Limited=out$Limited)

#write("Number of cells/nuclei:",stdout())
#write(sum(is.cell, na.rm=TRUE), stdout())
#trueCells = rownames(subset(out,FDR <= 0.01))
# Filter for good barcodes
#annotateddata <- data[,which(is.cell),drop=FALSE]

# read in matrix
#input <- readMM(file = matrix.path)
#feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE)
#barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
#colnames(input) = barcode.names$V1
#rownames(input) = feature.names$V2

# pre-filter and create initial seurat object
#bc_sums <- colSums(input)
#bg_id <- names(bc_sums[bc_sums >= opt$pre_filter])
#panc = CreateSeuratObject(counts = input[,bg_id],project=opt$sample_id,min.cells = 0)
obj <- CreateSeuratObject(counts = data,min.cells = 0)

### QC
# get percent mitochondrial content
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^GRCh38-MT-")

# plot pre-filter metadata
# plot metadata associations
pdf(paste(out_path,sample,"_prefilter_stats_2.pdf",sep=""))
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

pdf(paste(out_path,sample,"_prefilter_stats_3.pdf",sep=""))
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
# filter step

# Calculate % of cells with Mouse + Human genome information
obj[["human.percent"]] <- PercentageFeatureSet(obj, pattern = "GRCh38-")
obj[["mouse.percent"]] <- PercentageFeatureSet(obj, pattern = "mm10---")

#SUBSET OBJECT
sobj<-subset(x = obj, subset = nFeature_RNA > opt$nfeature_min & nFeature_RNA < opt$nfeature_max & nCount_RNA > opt$ncount_min & nCount_RNA < opt$ncount_max & percent.mt<opt$mito_max & human.percent > 95)

#rm mouse cells from gem classification
gem_class <- read.table(paste0(analysis_dir,"/gem_classification.csv"),sep=",",header=T)
meta_df <- sobj@meta.data %>% rownames_to_column
meta_df$call <- "NA"
meta_df <- meta_df %>% mutate(call= replace(call,rowname %in% subset(gem_class,call=="mm10")$barcode,"mm10"))
meta_df <- meta_df %>% mutate(call= replace(call,rowname %in% subset(gem_class,call=="GRCh38")$barcode,"GRCh38")) %>% column_to_rownames

bc <- sobj@meta.data %>% rownames
sobj@meta.data <- meta_df[bc,]

Idents(sobj) <- "call"
sobj <- subset(sobj,idents="GRCh38")
  
# run sctransform
sobj <- SCTransform(sobj, vars.to.regress = c("percent.mt","nCount_RNA"),return.only.var.genes = F)
# These are now standard steps in the Seurat workflow for visualization and clustering
sobj <- RunPCA(sobj, npcs = opt$pc_num,verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:opt$pc_num, reduction="pca")
sobj <- FindNeighbors(sobj, dims = 1:opt$pc_num, reduction="pca")
sobj <- FindClusters(sobj, resolution = 0.5)
sobj <- RunUMAP(sobj, dims = 1:opt$pc_num)

pdf(paste(out_path,sample,"_filter_dimplot_5.pdf",sep=""))
DimPlot(sobj, label = TRUE) + NoLegend()
dev.off()

#SAVE Filtered Object
saveRDS(sobj,paste(out_path,sample,"_processed_celltype.rds",sep=""))

