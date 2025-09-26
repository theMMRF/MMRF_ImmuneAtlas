#### MMRF Immune Atlas - single cell RNA seq analysis ####

options(stringsAsFactors = F)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(Matrix)
library(openxlsx)
library(glmGamPoi)
library(harmony)
library(patchwork)
library(ComplexHeatmap)
library(unix)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(DescTools)
library(viridis)
library(RColorBrewer)
library(ggpointdensity)
library(ggforce)
library(RColorBrewer)
library(randomcoloR)
library(rstatix)
library(EnhancedVolcano)
library(scatterplot3d)
library(readxl)

#### Data preparation ####

# Load data #

# Full processed object with SCT label transfer
MMRF_full_obj <- readRDS("/data1/datasets/MMRF/immune_atlas/immune_atlas_full_seurat_objects/derivative_data_RDS_object_current_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_SCTNorm_PC25_Harmony_singleR_doublet_subcluster_V03242023.rds")

# Plasma cell reference object (Coffey at al)
plasma_ref_coffey <- readRDS("/data1/users/junia/MMRF/plasmacell/plasmacell_reference_seurat_obj.rds")

# MMRF Immune atlas metadata onject
MMRF <- load("/data1/users/alison/RDS/risk_clinical_data_and_cell_abundances.RData")

# MMRF Seurat object (log norm)
MMRF_full_obj_log <- readRDS("/data1/datasets/MMRF/immune_atlas/immune_atlas_full_seurat_objects/derivative_data_RDS_object_current_SeuratObj_rm38_44_in_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_singleR_doublet_subcluster_V03242023.rds")


load("/data1/users/junia/MMRF/files/MMRF_immuneatlas.RData")



################################################################################
##get some cohort metrics  for immune PSN  cohort #####
load("/data1/users/junia/MMRF/filesSubclusterV03072023_Count_Percentage_PSN-MMRF.RData")
## subset MMRF_OBJ to immune psn samples ###
immune_psn_samples <- unique(PSN_percentage$Public_ID)
Idents(MMRF_full_obj) <- "public_id"
# unique(MMRF_full_obj$public_id)
# MMRF_full_obj$public_id[which(is.na(MMRF_full_obj$public_id))] <- "MMRF_2471"
MMRF_full_obj_immune_Psn <- subset(x = MMRF_full_obj, idents = immune_psn_samples)


## num samples ##
length(unique(MMRF_full_obj_immune_Psn$public_id))
##total cells ##
ncol(MMRF_full_obj_immune_Psn)
### cell specific numbers 
table(MMRF_full_obj_immune_Psn$compartment)

##############################################################################



# Plasma compartment analysis #


# Subset object

Idents(MMRF_full_obj) <- "subcluster_V03072023"

plasma_clusters <- paste0("Plasma.", 0:23)

MMRF_plasma <- subset(MMRF_full_obj, idents = plasma_clusters )

Idents(plasma_ref_coffey) <- "predicted.celltype.l2"

bplasma_clusters <- c("B naive", "B memory","B intermediate", "Plasmablast")

plasma_ref_coffey <- subset(plasma_ref_coffey, idents = bplasma_clusters )

# Prepare objects to label transfer

plasma_ref <- NormalizeData(plasma_ref)
plasma_ref <- FindVariableFeatures(plasma_ref)
plasma_ref <- ScaleData(plasma_ref)
plasma_ref <- RunPCA(plasma_ref)
plasma_ref <- FindNeighbors(plasma_ref, dims = 1:30)
plasma_ref <- FindClusters(plasma_ref)

plasma_ref <- RunUMAP(plasma_ref, dims = 1:30)
DimPlot(plasma_ref, group.by = c("celltytpe", "tech"))


MMRF_plasma <- NormalizeData(MMRF_plasma)

# Label Transfer

plasma_anchors <- FindTransferAnchors(reference = plasma_ref, query = MMRF_plasma, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = plasma_anchors, refdata = plasma_ref$celltype, dims = 1:30)
MMRF_plasma <- AddMetaData(MMRF_plasma, metadata = predictions)

# Check efficiency

MMRF_plasma$prediction.match <- MMRF_plasma$predicted.id == MMRF_plasma$celltype
table(MMRF_plasma$prediction.match)









################################################################################

# MM-PSN analysis #

# Unique patient IDs from lijun's dataframe
unique(lijun_df$Public_ID)
# 263 Unique Patient IDs from batch 1-5 (lijun's dataframe) out of 361 total samples
# "MMRF_1089", "MMRF_1227", "MMRF_1251", "MMRF_1264", "MMRF_1267", 
# "MMRF_1284", "MMRF_1286", "MMRF_1296", "MMRF_1311", "MMRF_1323", 
# "MMRF_1332", "MMRF_1335", "MMRF_1338", "MMRF_1342", "MMRF_1348", 
# "MMRF_1349", "MMRF_1352", "MMRF_1359", "MMRF_1360", "MMRF_1364", 
# "MMRF_1371", "MMRF_1375", "MMRF_1391", "MMRF_1395", "MMRF_1402", 
# "MMRF_1407", "MMRF_1408", "MMRF_1415", "MMRF_1423", "MMRF_1424", 
# "MMRF_1426", "MMRF_1430", "MMRF_1431", "MMRF_1433", "MMRF_1434", 
# "MMRF_1450", "MMRF_1454", "MMRF_1461", "MMRF_1466", "MMRF_1470", 
# "MMRF_1474", "MMRF_1479", "MMRF_1484", "MMRF_1490", "MMRF_1491", 
# "MMRF_1494", "MMRF_1497", "MMRF_1500", "MMRF_1510", "MMRF_1518", 
# "MMRF_1519", "MMRF_1522", "MMRF_1523", "MMRF_1525", "MMRF_1535", 
# "MMRF_1537", "MMRF_1538", "MMRF_1539", "MMRF_1542", "MMRF_1555", 
# "MMRF_1560", "MMRF_1565", "MMRF_1572", "MMRF_1573", "MMRF_1576", 
# "MMRF_1584", "MMRF_1586", "MMRF_1593", "MMRF_1605", "MMRF_1607", 
# "MMRF_1613", "MMRF_1614", "MMRF_1624", "MMRF_1627", "MMRF_1628", 
# "MMRF_1634", "MMRF_1639", "MMRF_1640", "MMRF_1643", "MMRF_1646", 
# "MMRF_1647", "MMRF_1648", "MMRF_1650", "MMRF_1659", "MMRF_1670", 
# "MMRF_1677", "MMRF_1683", "MMRF_1687", "MMRF_1695", "MMRF_1717", 
# "MMRF_1719", "MMRF_1720", "MMRF_1722", "MMRF_1730", "MMRF_1740", 
# "MMRF_1746", "MMRF_1755", "MMRF_1768", "MMRF_1778", "MMRF_1780", 
# "MMRF_1781", "MMRF_1783", "MMRF_1787", "MMRF_1802", "MMRF_1807", 
# "MMRF_1819", "MMRF_1823", "MMRF_1832", "MMRF_1833", "MMRF_1837", 
# "MMRF_1838", "MMRF_1839", "MMRF_1841", "MMRF_1850", "MMRF_1865", 
# "MMRF_1880", "MMRF_1888", "MMRF_1893", "MMRF_1859", "MMRF_1895", 
# "MMRF_1918", "MMRF_1932", "MMRF_1941", "MMRF_1951", "MMRF_1965", 
# "MMRF_1967", "MMRF_1977", "MMRF_1981", "MMRF_1984", "MMRF_1988", 
# "MMRF_1991", "MMRF_2018", "MMRF_2019", "MMRF_2035", "MMRF_2055", 
# "MMRF_2057", "MMRF_2062", "MMRF_2064", "MMRF_2076", "MMRF_2082", 
# "MMRF_2106", "MMRF_2111", "MMRF_2119", "MMRF_2122", "MMRF_2126", 
# "MMRF_2128", "MMRF_2129", "MMRF_2132", "MMRF_2141", "MMRF_2143", 
# "MMRF_2153", "MMRF_2168", "MMRF_2172", "MMRF_2176", "MMRF_2180", 
# "MMRF_2183", "MMRF_2186", "MMRF_2187", "MMRF_2190", "MMRF_2191", 
# "MMRF_2204", "MMRF_2207", "MMRF_2211", "MMRF_2231", "MMRF_2232", 
# "MMRF_2250", "MMRF_2256", "MMRF_2258", "MMRF_2259", "MMRF_2270", 
# "MMRF_2271", "MMRF_2272", "MMRF_2289", "MMRF_2290", "MMRF_2291", 
# "MMRF_2297", "MMRF_2301", "MMRF_2302", "MMRF_2324", "MMRF_2326", 
# "MMRF_2329", "MMRF_2330", "MMRF_2331", "MMRF_2333", "MMRF_2339", 
# "MMRF_2341", "MMRF_2353", "MMRF_2373", "MMRF_2377", "MMRF_2378", 
# "MMRF_2397", "MMRF_2400", "MMRF_2409", "MMRF_2410", "MMRF_2422", 
# "MMRF_2436", "MMRF_2440", "MMRF_2442", "MMRF_2445", "MMRF_2455", 
# "MMRF_2458", "MMRF_2461", "MMRF_2464", "MMRF_2469", "MMRF_2471", 
# "MMRF_2473", "MMRF_2457", "MMRF_2475", "MMRF_2476", "MMRF_2477", 
# "MMRF_2480", "MMRF_2485", "MMRF_2487", "MMRF_2489", "MMRF_2497", 
# "MMRF_2498", "MMRF_2504", "MMRF_2505", "MMRF_2506", "MMRF_2511", 
# "MMRF_2514", "MMRF_2516", "MMRF_2517", "MMRF_2523", "MMRF_2525", 
# "MMRF_2531", "MMRF_2532", "MMRF_2538", "MMRF_2539", "MMRF_2546", 
# "MMRF_2549", "MMRF_2550", "MMRF_2554", "MMRF_2557", "MMRF_2560", 
# "MMRF_2562", "MMRF_2564", "MMRF_2566", "MMRF_2568", "MMRF_2574", 
# "MMRF_2579", "MMRF_2581", "MMRF_2585", "MMRF_2589", "MMRF_2590", 
# "MMRF_2595", "MMRF_2596", "MMRF_2605", "MMRF_2611", "MMRF_2614", 
# "MMRF_2621", "MMRF_2622", "MMRF_1460", "MMRF_1487", "MMRF_1489", 
# "MMRF_1436", "MMRF_1445", "MMRF_1266", "MMRF_1295", "MMRF_1389", 
# "MMRF_1401", "MMRF_1485", "MMRF_1486"

# MSSM samples
unique(edgar_df$Study_Site)
# "EMORY", "MAYO", "WUSTL", "MSSM"

mssm_edgar <- subset(edgar_df, Study_Site == "MSSM")
View(mssm_edgar)
unique(mssm_edgar$Public_ID)
# 45 mssm patients 
## "MMRF_1778", "MMRF_1780", "MMRF_1783", "MMRF_1787", "MMRF_1823", 
## "MMRF_1833", "MMRF_1837", "MMRF_1965", "MMRF_1981", "MMRF_1991", 
## "MMRF_2055", "MMRF_2076", "MMRF_2106", "MMRF_2119", "MMRF_2122", 
## "MMRF_2132", "MMRF_2141", "MMRF_2143", "MMRF_2180", "MMRF_2186", 
## "MMRF_2190", "MMRF_2204", "MMRF_2211", "MMRF_2256", "MMRF_2259", 
## "MMRF_2301", "MMRF_2329", "MMRF_2436", "MMRF_2458", "MMRF_2461", 
## "MMRF_2471", "MMRF_2457", "MMRF_2480", "MMRF_2485", "MMRF_2487", 
## "MMRF_2489", "MMRF_2498", "MMRF_2516", "MMRF_2523", "MMRF_2568",
## "MMRF_2585", "MMRF_2589", "MMRF_2596", "MMRF_2605", "MMRF_2611"

mssm_IDs <- list(c("MMRF_1778", "MMRF_1780", "MMRF_1783", "MMRF_1787", "MMRF_1823", 
                   "MMRF_1833", "MMRF_1837", "MMRF_1965", "MMRF_1981", "MMRF_1991", 
                   "MMRF_2055", "MMRF_2076", "MMRF_2106", "MMRF_2119", "MMRF_2122", 
                   "MMRF_2132", "MMRF_2141", "MMRF_2143", "MMRF_2180", "MMRF_2186", 
                   "MMRF_2190", "MMRF_2204", "MMRF_2211", "MMRF_2256", "MMRF_2259", 
                   "MMRF_2301", "MMRF_2329", "MMRF_2436", "MMRF_2458", "MMRF_2461", 
                   "MMRF_2471", "MMRF_2457", "MMRF_2480", "MMRF_2485", "MMRF_2487", 
                   "MMRF_2489", "MMRF_2498", "MMRF_2516", "MMRF_2523", "MMRF_2568",
                   "MMRF_2585", "MMRF_2589", "MMRF_2596", "MMRF_2605", "MMRF_2611"))

mssm_lijun <- subset(lijun_df, Study_Site == "MSSM")
unique(mssm_lijun$Public_ID)
# 66 pts from mssm (edgar explained this is becasue lijun's is more complete with the batch 5 samples)
mssm_IDs_2 <- list(c("MMRF_1755", "MMRF_1768", "MMRF_1778", "MMRF_1780", "MMRF_1783", 
                     "MMRF_1787", "MMRF_1823", "MMRF_1833", "MMRF_1837", "MMRF_1839", 
                     "MMRF_1841", "MMRF_1888", "MMRF_1918", "MMRF_1965", "MMRF_1981", 
                     "MMRF_1991", "MMRF_2035", "MMRF_2055", "MMRF_2057", "MMRF_2076", 
                     "MMRF_2106", "MMRF_2119", "MMRF_2122", "MMRF_2129", "MMRF_2132", 
                     "MMRF_2141", "MMRF_2143", "MMRF_2168", "MMRF_2180", "MMRF_2186", 
                     "MMRF_2190", "MMRF_2204", "MMRF_2207", "MMRF_2211", "MMRF_2256", 
                     "MMRF_2259", "MMRF_2291", "MMRF_2301", "MMRF_2329", "MMRF_2378", 
                     "MMRF_2400", "MMRF_2436", "MMRF_2458", "MMRF_2461", "MMRF_2464", 
                     "MMRF_2469", "MMRF_2471", "MMRF_2457", "MMRF_2480", "MMRF_2485", 
                     "MMRF_2487", "MMRF_2489", "MMRF_2498", "MMRF_2516", "MMRF_2517", 
                     "MMRF_2523", "MMRF_2532", "MMRF_2539", "MMRF_2562", "MMRF_2564", 
                     "MMRF_2568", "MMRF_2585", "MMRF_2589", "MMRF_2596", "MMRF_2605", 
                     "MMRF_2611"))

# seeing which mssm ids are different between edgar's and lijun's
setdiff(mssm_IDs_2[[1]], mssm_IDs[[1]])
# 21 samples
# "MMRF_1755", "MMRF_1768", "MMRF_1839", "MMRF_1841", "MMRF_1888", 
# "MMRF_1918", "MMRF_2035", "MMRF_2057", "MMRF_2129", "MMRF_2168", 
# "MMRF_2207", "MMRF_2291", "MMRF_2378", "MMRF_2400", "MMRF_2464", 
# "MMRF_2469", "MMRF_2517", "MMRF_2532", "MMRF_2539", "MMRF_2562", 
# "MMRF_2564"

# Loading in PSN IDs and dataframes
# loading PSN_labels
PSN_labels <- read.csv("/data1/users/junia/MMRF/files/PSN_labels.csv")
PSN_labels <- as.data.frame(PSN_labels)
View(PSN_labels)

## creating a main group column
# sub() replaces everything that's not a digit (\D) with an empty string
PSN_labels$MajorGroup <- sub("\\D.*", "", PSN_labels$SubGroup)
head(PSN_labels)

intersect((unique(PSN_labels$Sample)), (unique(lijun_df$Public_ID)))
# PSN sample ids
PSN_IDs <- list(c("MMRF_1408", "MMRF_1342", "MMRF_1395", "MMRF_1364", "MMRF_1360",
                  "MMRF_1371", "MMRF_1089", "MMRF_1323", "MMRF_1284", "MMRF_1466",
                  "MMRF_1445", "MMRF_1436", "MMRF_1450", "MMRF_1460", "MMRF_1489",
                  "MMRF_1485", "MMRF_1519", "MMRF_1510", "MMRF_1538", "MMRF_1522",
                  "MMRF_1573", "MMRF_1593", "MMRF_1614", "MMRF_1624", "MMRF_1639",
                  "MMRF_1640", "MMRF_1643", "MMRF_1695", "MMRF_1730", "MMRF_1648",
                  "MMRF_1659", "MMRF_1720", "MMRF_1719", "MMRF_1768", "MMRF_1802",
                  "MMRF_1423", "MMRF_1859", "MMRF_1833", "MMRF_1839", "MMRF_1850",
                  "MMRF_1895", "MMRF_1537", "MMRF_1967", "MMRF_1981", "MMRF_1880",
                  "MMRF_1932", "MMRF_1977", "MMRF_1888", "MMRF_1991", "MMRF_2018",
                  "MMRF_2035", "MMRF_2057", "MMRF_2126", "MMRF_2055", "MMRF_2076",
                  "MMRF_2082", "MMRF_2111", "MMRF_2119", "MMRF_2141", "MMRF_2143",
                  "MMRF_2180", "MMRF_2204", "MMRF_2183", "MMRF_2187", "MMRF_2211",
                  "MMRF_2232", "MMRF_2168", "MMRF_2258", "MMRF_2329", "MMRF_2291",
                  "MMRF_2290", "MMRF_2302", "MMRF_2301", "MMRF_2330", "MMRF_2353",
                  "MMRF_2377", "MMRF_2409", "MMRF_2378", "MMRF_2400", "MMRF_2440",
                  "MMRF_2455", "MMRF_2422", "MMRF_2464", "MMRF_2436", "MMRF_2445",
                  "MMRF_2461", "MMRF_2475", "MMRF_2523", "MMRF_2489", "MMRF_2498",
                  "MMRF_2506", "MMRF_2514", "MMRF_2516", "MMRF_2525", "MMRF_2532",
                  "MMRF_2549", "MMRF_2557", "MMRF_2562", "MMRF_2581", "MMRF_2590",
                  "MMRF_2554", "MMRF_2595", "MMRF_2614", "MMRF_2596", "MMRF_1988",
                  "MMRF_1251", "MMRF_1389", "MMRF_1391", "MMRF_1359", "MMRF_1267",
                  "MMRF_1433", "MMRF_1424", "MMRF_1454", "MMRF_1461", "MMRF_1349",
                  "MMRF_1491", "MMRF_1497", "MMRF_1560", "MMRF_1572", "MMRF_1586",
                  "MMRF_1605", "MMRF_1613", "MMRF_1634", "MMRF_1627", "MMRF_1683",
                  "MMRF_1778", "MMRF_1823", "MMRF_1832", "MMRF_1783", "MMRF_1780",
                  "MMRF_1965", "MMRF_2062", "MMRF_2106", "MMRF_2122", "MMRF_2172",
                  "MMRF_2176", "MMRF_2231", "MMRF_2250", "MMRF_2271", "MMRF_2272",
                  "MMRF_2297", "MMRF_2326", "MMRF_2331", "MMRF_2373", "MMRF_2473",
                  "MMRF_2497", "MMRF_2505", "MMRF_2531", "MMRF_2538", "MMRF_2574",
                  "MMRF_2621", "MMRF_1338", "MMRF_1352", "MMRF_1348", "MMRF_1434",
                  "MMRF_1401", "MMRF_1484", "MMRF_1286", "MMRF_1539", "MMRF_1542",
                  "MMRF_1677", "MMRF_1722", "MMRF_1523", "MMRF_1670", "MMRF_1787",
                  "MMRF_1837", "MMRF_1565", "MMRF_1865", "MMRF_1893", "MMRF_1918",
                  "MMRF_2153", "MMRF_1755", "MMRF_2339", "MMRF_2458", "MMRF_2469",
                  "MMRF_2476", "MMRF_2457", "MMRF_2480", "MMRF_2485", "MMRF_2539",
                  "MMRF_2568", "MMRF_2579", "MMRF_2622", "MMRF_2605", "MMRF_2471"))

length(PSN_IDs[[1]])
# 185 patients but these are actually mixed patients from other places

## making a subset for these 185 patients
MMRF_PSN_labels <- filter(PSN_labels, Sample %in% PSN_IDs[[1]])


# Seeing which patients from MSSM have PSN_IDs
intersect(PSN_IDs[[1]], mssm_IDs_2[[1]])
# only 52 patients that are from mssm and have psn ids
# "MMRF_1768", "MMRF_1833", "MMRF_1839", "MMRF_1981", "MMRF_1888", 
# "MMRF_1991", "MMRF_2035", "MMRF_2057", "MMRF_2055", "MMRF_2076",
# "MMRF_2119", "MMRF_2141", "MMRF_2143", "MMRF_2180", "MMRF_2204", 
# "MMRF_2211", "MMRF_2168", "MMRF_2329", "MMRF_2291", "MMRF_2301", 
# "MMRF_2378", "MMRF_2400", "MMRF_2464", "MMRF_2436", "MMRF_2461", 
# "MMRF_2523", "MMRF_2489", "MMRF_2498", "MMRF_2516", "MMRF_2532", 
# "MMRF_2562", "MMRF_2596", "MMRF_1778", "MMRF_1823", "MMRF_1783", 
# "MMRF_1780", "MMRF_1965", "MMRF_2106", "MMRF_2122", "MMRF_1787", 
# "MMRF_1837", "MMRF_1918", "MMRF_1755", "MMRF_2458", "MMRF_2469", 
# "MMRF_2457", "MMRF_2480", "MMRF_2485", "MMRF_2539", "MMRF_2568", 
# "MMRF_2605", "MMRF_2471"

final_IDs <- list(c("MMRF_1768", "MMRF_1833", "MMRF_1839", "MMRF_1981", "MMRF_1888", 
                    "MMRF_1991", "MMRF_2035", "MMRF_2057", "MMRF_2055", "MMRF_2076",
                    "MMRF_2119", "MMRF_2141", "MMRF_2143", "MMRF_2180", "MMRF_2204", 
                    "MMRF_2211", "MMRF_2168", "MMRF_2329", "MMRF_2291", "MMRF_2301", 
                    "MMRF_2378", "MMRF_2400", "MMRF_2464", "MMRF_2436", "MMRF_2461", 
                    "MMRF_2523", "MMRF_2489", "MMRF_2498", "MMRF_2516", "MMRF_2532", 
                    "MMRF_2562", "MMRF_2596", "MMRF_1778", "MMRF_1823", "MMRF_1783", 
                    "MMRF_1780", "MMRF_1965", "MMRF_2106", "MMRF_2122", "MMRF_1787", 
                    "MMRF_1837", "MMRF_1918", "MMRF_1755", "MMRF_2458", "MMRF_2469", 
                    "MMRF_2457", "MMRF_2480", "MMRF_2485", "MMRF_2539", "MMRF_2568", 
                    "MMRF_2605", "MMRF_2471"))




## *MAIN GROUP DISTRIBUTION*
## Count unique patients in each major group (185)
pt_distribution_MG <- MMRF_PSN_labels %>%
  group_by(MajorGroup) %>%
  summarise(PatientCount = n_distinct(Sample))
print(pt_distribution_MG)

## Distribution of patients only from MSSM and have PSN labels
# Filter PSN_labels for IDs in final_IDs
filtered_PSN_labels <- PSN_labels %>%
  filter(Sample %in% final_IDs[[1]])
# Count unique patients in each major group for the filtered IDs
pt_distribution_MG_MSSM <- filtered_PSN_labels %>%
  group_by(MajorGroup) %>%
  summarise(PatientCount = n_distinct(Sample))
# Print the result
print(pt_distribution_MG_MSSM)


## *SUBGROUP DISTRIBUTION*
## Count unique patients in each subgroup
pt_distribution_sG <- MMRF_PSN_labels %>%
  group_by(SubGroup) %>%
  summarise(PatientCount = n_distinct(Sample))
print(pt_distribution_sG)

## Distribution of patients only from MSSM and have PSN labels
# Filter PSN_labels for IDs in final_IDs
filtered_PSN_labels <- PSN_labels %>%
  filter(Sample %in% final_IDs[[1]])
# Count unique patients in each major group for the filtered IDs
pt_distribution_sG_MSSM <- filtered_PSN_labels %>%
  group_by(SubGroup) %>%
  summarise(PatientCount = n_distinct(Sample))
# Print the result
print(pt_distribution_sG_MSSM)

#### Cell Population Frequency analysis ####

# Create frequency dataframe for analysis
## using lijun's df since it is the one that has the latest data
PSN_sc_lj <- filter(lijun_df, Public_ID %in% PSN_IDs[[1]])

# Column of interest
cols <- c("subcluster_V03072023", "Public_ID")

# Creating frequency columns for columns of interest
## absolute count of cells in each cell type across all samples
PSN_sc_lj_freq <- PSN_sc_lj %>%
  group_by(across(all_of(cols))) %>%
  add_count(name = "Frequency")


# Making column for percentage
## Grouping by patient and cell type
PSN_sc_lj_freq_2 <- PSN_sc_lj %>%
  group_by(Public_ID, subcluster_V03072023)

## Summarize to get counts; modify depending on data structure:
PSN_sc_lj_freq_2 <- PSN_sc_lj_freq_2 %>%
  summarise(Count_per_patient = n(), .groups = 'drop')  # Counts the number of rows for each combination of Public_ID and subcluster_V03072023

## calculating the total number of cells per patient so that I have the denominator for percentage
PSN_sc_lj_total <- PSN_sc_lj_freq_2 %>%
  group_by(Public_ID) %>%
  summarise(Total_Cells = sum(Count_per_patient))

save(PSN_sc_lj_total, file = "/data1/users/junia/MMRF/files/total_cell_count_per_patient-MMRF.RData")

## joining the three dataframes back to one
# Join the total cells data back to the frequency data
PSN_sc_lj_freq_2 <- PSN_sc_lj_freq_2 %>%
  left_join(PSN_sc_lj_total, by = "Public_ID")

# Calculating the percentage of each cell type per patient
PSN_sc_lj_freq_2 <- PSN_sc_lj_freq_2 %>%
  mutate(Percentage = (Count_per_patient / Total_Cells) * 100)


# Sanity Check
## Summing the percentages for each patient to see if they add up to 100%
percentage_sum_check <- PSN_sc_lj_freq_2 %>%
  group_by(Public_ID) %>%
  summarise(Total_Percentage = sum(Percentage))

print(percentage_sum_check)


# Joining dataframes with PSN label and percentage
PSN_percentage <- dplyr::left_join(PSN_sc_lj_freq_2, PSN_labels, by = c("Public_ID"="Sample"))
View(PSN_percentage)

save(PSN_percentage, file = "/data1/users/junia/MMRF/filesSubclusterV03072023_Count_Percentage_PSN-MMRF.RData")

# Statistical analysis 
## The Wilcoxon rank sum test is a non-parametric test used to compare two groups to determine 
## if there are statistically significant differences between their distributions. 
## This is particularly useful when the data does not meet the assumptions of parametric tests, 
## such as normal distribution.
PSN_percent_wilcox <- pairwise.wilcox.test(PSN_percentage$Percentage, 
                                           PSN_percentage$MajorGroup, 
                                           p.adjust.method = "none")
PSN_percent_wilcox


# *Making a list of cell types*
unique(PSN_percentage$subcluster_V03072023)
# 105 unique cell types part of the subluster_V03072023
subclusterV_celltypes <- list( "BEry.7", "Ery.4", "Ery.2", "Ery.7", "NkT.0", 
                               "BEry.15.0", "Myeloid.10", "Plasma.1", "BEry.4", "Plasma.0", 
                               "BEry.14", "NkT.15", "BEry.0", "NkT.1.5", "NkT.1.2", 
                               "Myeloid.12", "Myeloid.11.0", "Ery.5", "Ery.3", "Myeloid.5", 
                               "Myeloid.3", "Myeloid.0", "Myeloid.4", "Plasma.11", "Myeloid.6", 
                               "Myeloid.7", "NkT.1.1", "Ery.9", "Ery.6", "NkT.10.1", 
                               "Ery.0", "NkT.9.0", "NkT.6", "Plasma.8", "NkT.1.0", 
                               "NkT.3.2", "Myeloid.8", "Plasma.2", "BEry.11", "Myeloid.14", 
                               "Ery.1", "Plasma.10", "Full.23", "NkT.3.1", "BEry.3", 
                               "BEry.8", "NkT.8", "Plasma.4", "Ery.11", "Plasma.13", 
                               "Plasma.7", "Myeloid.9", "NkT.4", "Myeloid.1", "Plasma.18", 
                               "Plasma.6", "BEry.17", "Ery.8", "NkT.9.1", "BEry.1", 
                               "BEry.10", "NkT.1.3", "NkT.13", "BEry.9", "NkT.16", 
                               "NkT.10.0", "Plasma.5", "NkT.2.1", "NkT.5.0", "BEry.6", 
                               "NkT.12", "NkT.3.0", "NkT.2.0", "NkT.7", "NkT.5.1", 
                               "BEry.15.1", "BEry.2", "NkT.11", "Myeloid.13", "Plasma.12", 
                               "Plasma.9", "Ery.10", "BEry.16", "BEry.12", "Plasma.3", 
                               "Myeloid.2", "BEry.5", "BEry.13", "Plasma.19", "Myeloid.15", 
                               "NkT.14", "Plasma.15", "NkT.2.3", "NkT.1.4", "Myeloid.11.1", 
                               "NkT.2.2", "NkT.2.4", "Plasma.14","Plasma.16", "Myeloid.11.2", 
                               "Plasma.21", "NkT.5.2", "Plasma.22", "Plasma.23", "Plasma.17"
)



## *Function to generate plot for major group*
my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )

make_plot <- function(df, cell_type, plot_directory) {
  filtered_df <- df %>% dplyr::filter(subcluster_V03072023 == cell_type)
  
  if(nrow(filtered_df) == 0) {
    print(paste("No data for cell type: ", cell_type))
    return()
  }
  
  plot <- filtered_df %>%
    ggplot(aes(x = MajorGroup, y = Percentage, fill = MajorGroup)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(size = 1) +
    stat_compare_means(comparisons = my_comparisons, p.adjust.method = "bonferroni") +
    scale_fill_tableau("Color Blind") +
    theme_bw() +
    ggtitle(paste("Plot for ", cell_type)) +
    labs(y = "Percent cells", x = "PSN Group") +
    theme(text = element_text(size = 20))
  
  print(plot)
  
  # Construct the filename based on the cell type
  filename <- paste("Plot_for_", gsub(" ", "_", cell_type), ".png", sep="")
  filepath <- file.path(plot_directory, filename)
  
  # Save the plot
  ggsave(filepath, plot = plot, width = 10, height = 8, dpi = 1000)
}

plot_directory <- "/data1/users/junia/MMRF/files/Main_Groups"

#cell_type_id <- unique(PSN_percentage$subcluster_V03072023)

for (cell_type in subclusterV_celltypes) {
  make_plot(PSN_percentage, cell_type, plot_directory)
}


#significant percentage comparisons (p<0.05)
##Plasma.4, Plasma.6, Plasma.9, Plasma.16
##NkT.0, NkT.1.4, NkT.6, NkT.7, NkT.9.0, NkT.12, NkT.14
##BEry.1, BEry.8, BEry.10, BEry.11, BEry.12, BEry.17
##Myeloid.7, Myeloid.9


################################################################################

# Cluster analysis for significant cluster #

# Significant clusters: Plasma.4, Plasma.6, Plasma.9, Plasma.16, NkT.0, NkT.1.4, 
# NkT.6, NkT.7, NkT.9.0, NkT.12, NkT.14, BEry.1, BEry.8, BEry.10, BEry.11, BEry.12,
# BEry.17, Myeloid.7, Myeloid.9

# significant_clusters <- c( "Plasma.4", "Plasma.6", "Plasma.9", "Plasma.16", 
#                            "NkT.0", "NkT.1.4", "NkT.6", "NkT.7", "NkT.9.0", "NkT.12",
#                            "NkT.14", "BEry.1", "BEry.8", "BEry.10", "BEry.11", 
#                            "BEry.12","BEry.17", "Myeloid.7", "Myeloid.9")

significant_clusters <- c( "NkT.14", "BEry.1", "BEry.8", "BEry.10", "BEry.11", 
                           "BEry.12","BEry.17", "Myeloid.7", "Myeloid.9")

Idents(MMRF_full_obj_log) <- "subcluster_V03072023"


## RNA diff expression ##
DefaultAssay(MMRF_full_obj_log) <- 'RNA'
for (i in 1:length(significant_clusters)) {
  subset_name <- paste0("markers_adt_C", significant_clusters[i])
  markers  <- FindMarkers(MMRF_full_obj_log, ident.1 = significant_clusters[i], logfc.threshold = 0.05)
  assign(subset_name,markers)
  write.csv(markers,paste0("/data1/users/junia/MMRF/files/DE_sig_clusters/MMRF_immuneatlas_sig_clusters_RNA_markers_", significant_clusters[i], ".csv"))
}

## RNA diff with MM PSN groups 

DefaultAssay(BMMC_Tcells) <- 'RNA'

for (i in 1:length(Tcell_clusters))  {
  subset_name <- paste0("markers_RNA_C",Tcell_clusters[i],"_PFS")
  subset_group <- subset(BMMC_Tcells, ident = Tcell_clusters[i])
  Idents(subset_group) <- "PFS_36_months"
  markers <- FindMarkers(subset_group, ident.1 = "PFS > 36 Months", logfc.threshold = 0.05)
  assign(subset_name, markers)
  write.csv(markers,paste0("/data1/Immunai/JNJ_CAR_T/Cilta_Cel_Revisions/Figures_and_Plots/unsupervised_markers/RNA/BMMC_Tcell_markers_PFS_C", Tcell_clusters[i], ".csv"))
}



# ## T cells RNA volcano #
# 
# CARTcell_RNA_Diff_exp <- list(markers_RNA_C8,markers_RNA_C30)
# 
# CARTcell_RNA_Diff_exp_names <- c("markers_RNA_C8","markers_RNA_C30")
# 
# for(i in 1:length(CARTcell_RNA_Diff_exp)){
#   N_tab <- CARTcell_RNA_Diff_exp[[i]]
#   N_tab_sig <- N_tab[which(N_tab$p_val_adj < 0.05),]
#   NOlabels <- rownames(N_tab_sig)
#   volcano <- EnhancedVolcano(N_tab, rownames(N_tab), title = CARTcell_RNA_Diff_exp_names[i],
#                              x = "avg_log2FC", y = "p_val_adj", selectLab = NOlabels)
#   ggsave(filename = paste0("/data1/Immunai/JNJ_CAR_T/Cilta_Cel_Revisions/Figures_and_Plots/unsupervised_markers/RNA/CILTA_CEL_BMMC_CARTcell", 
#                            CARTcell_RNA_Diff_exp_names[i], ".pdf"), 
#          plot = volcano, width = 20, height = 10)
# }
# 
# 
# ## CARTcell RNA PFS volcano #
# 
# CARTcell_RNA_Diff_exp <- list(markers_RNA_C8_PFS,markers_RNA_C30_PFS)
# 
# CARTcell_RNA_Diff_exp_names <- c("markers_RNA_C8_PFS","markers_RNA_C30_PFS")
# 
# for(i in 1:length(CARTcell_RNA_Diff_exp_PFS)){
#   N_tab <- CARTcell_RNA_Diff_exp_PFS[[i]]
#   N_tab_sig <- N_tab[which(N_tab$p_val_adj < 0.05),]
#   NOlabels <- rownames(N_tab_sig)
#   volcano <- EnhancedVolcano(N_tab, rownames(N_tab), title = CARTcell_RNA_Diff_exp_names_PFS[i],
#                              x = "avg_log2FC", y = "p_val_adj", selectLab = NOlabels)
#   ggsave(filename = paste0("/data1/Immunai/JNJ_CAR_T/Cilta_Cel_Revisions/Figures_and_Plots/unsupervised_markers/RNA/CILTA_CEL_BMMC_PFS_CARTcell", 
#                            CARTcell_RNA_Diff_exp_names_PFS[i], ".pdf"), 
#          plot = volcano, width = 20, height = 10)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


################################################################################
# # *Function for subgroup plots*
# plot_directory_2 <- "/data1/users/junia/MMRF/files/subgroups"
# 
# comparisons_group1 <- list(c("1a", "1b"), c("1a", "1c"), c("1a", "1d"),
#                            c("1b", "1c"), c("1b", "1d"), c("1c", "1d"))
# comparisons_group2 <- list(c("2a", "2b"), c("2a", "2c"), c("2a", "2d"), c("2a", "2e"),
#                            c("2b", "2c"), c("2b", "2d"), c("2b", "2e"),
#                            c("2c", "2d"), c("2c", "2e"), c("2d", "2e"))
# comparisons_group3 <- list(c("3a", "3b"), c("3a", "3c"), c("3b", "3c"))
# 
# ### Create and save plots for each subgroup within each main group
# main_group <- (unique(PSN_percentage$MajorGroup))
# 
# 
# ### Function to create and save plots for each subgroup within each main group
# make_plot_for_subgroup <- function(df, cell_type, main_group, plot_directory) {
#   # Ensure the plot directory exists
#   if (!dir.exists(plot_directory)) {
#     dir.create(plot_directory, recursive = TRUE)
#   }
#   
#   # Select the appropriate comparisons based on the main group
#   comparisons <- switch(as.character(main_group),
#                         "1" = comparisons_group1,
#                         "2" = comparisons_group2,
#                         "3" = comparisons_group3,
#                         stop("Main group not recognized."))
#   
#   # Filter for cell type and main group
#   filtered_df <- df %>%
#     filter(subcluster_V03072023 == cell_type, MajorGroup == as.character(main_group))
#   
#   if(nrow(filtered_df) == 0) {
#     message("No data for cell type: ", cell_type, " in PSN Group ", main_group)
#     return()
#   }
#   
#   # Create the plot with the specific subgroup comparisons
#   plot <- ggplot(filtered_df, aes(x = SubGroup, y = Percentage, fill = SubGroup)) +
#     geom_boxplot(outlier.size = 0) +
#     geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
#     stat_compare_means(comparisons = comparisons, p.adjust.method = "bonferroni") +
#     scale_fill_tableau("Color Blind") +
#     theme_bw() +
#     labs(title = paste("Plot for ", cell_type, "- PSN Group ", main_group),
#          y = "Percent Cells", x = "SubGroup") +
#     theme(plot.title = element_text(size = 14),
#           axis.title = element_text(size = 12),
#           axis.text = element_text(size = 10),
#           legend.position = "bottom",
#           legend.title = element_text(size = 12),
#           legend.text = element_text(size = 10))
#   
#   # Print the plot
#   print(plot)
#   
#   # Save the plot to file
#   filename <- paste0("Plot_for_", gsub(" ", "_", cell_type), "_Group_", main_group, ".png")
#   ggsave(file.path(plot_directory, filename), plot = plot, width = 10, height = 8, dpi = 300)
# }
# 
# #### Looping through each cell type and main group to create and save plots
# for (cell_type in unique(PSN_percentage$subcluster_V03072023)) {
#   for (main_group in unique(PSN_percentage$MajorGroup)) {
#     make_plot_for_subgroup(PSN_percentage, cell_type, main_group, plot_directory_2)
#   }
# }
# 
# 
# 
# # *Getting IDs for patients that have Plasma.4 and Plasma.6*
# ## Filter the dataframe for rows where subcluster_V03072023 is 'Plasma.4'
# plasma_4_df <- PSN_percentage %>%
#   filter(subcluster_V03072023 == "Plasma.4")
# ## 139 pts with plasma.4 clusters
# plasma_4_IDs <- list(c("MMRF_1089", "MMRF_1251", "MMRF_1267", "MMRF_1286", "MMRF_1323", 
#                        "MMRF_1338", "MMRF_1352", "MMRF_1359", "MMRF_1364", "MMRF_1389", 
#                        "MMRF_1391", "MMRF_1395", "MMRF_1401", "MMRF_1408", "MMRF_1423", 
#                        "MMRF_1424", "MMRF_1433", "MMRF_1434", "MMRF_1445", "MMRF_1450", 
#                        "MMRF_1454", "MMRF_1460", "MMRF_1461", "MMRF_1484", "MMRF_1485", 
#                        "MMRF_1491", "MMRF_1497", "MMRF_1510", "MMRF_1522", "MMRF_1523", 
#                        "MMRF_1539", "MMRF_1542", "MMRF_1560", "MMRF_1565", "MMRF_1586", 
#                        "MMRF_1605", "MMRF_1613", "MMRF_1614", "MMRF_1624", "MMRF_1627", 
#                        "MMRF_1634", "MMRF_1639", "MMRF_1643", "MMRF_1648", "MMRF_1659", 
#                        "MMRF_1670", "MMRF_1677", "MMRF_1683", "MMRF_1695", "MMRF_1719", 
#                        "MMRF_1720", "MMRF_1722", "MMRF_1755", "MMRF_1768", "MMRF_1778", 
#                        "MMRF_1780", "MMRF_1787", "MMRF_1802", "MMRF_1823", "MMRF_1832", 
#                        "MMRF_1833", "MMRF_1839", "MMRF_1850", "MMRF_1893", "MMRF_1918", 
#                        "MMRF_1965", "MMRF_1967", "MMRF_1977", "MMRF_1981", "MMRF_1988", 
#                        "MMRF_2018", "MMRF_2035", "MMRF_2057", "MMRF_2062", "MMRF_2076", 
#                        "MMRF_2106", "MMRF_2122", "MMRF_2126", "MMRF_2143", "MMRF_2168", 
#                        "MMRF_2172", "MMRF_2176", "MMRF_2180", "MMRF_2187", "MMRF_2204", 
#                        "MMRF_2211", "MMRF_2231", "MMRF_2250", "MMRF_2258", "MMRF_2271", 
#                        "MMRF_2272", "MMRF_2290", "MMRF_2291", "MMRF_2297", "MMRF_2301", 
#                        "MMRF_2302", "MMRF_2326", "MMRF_2329", "MMRF_2330", "MMRF_2331", 
#                        "MMRF_2339", "MMRF_2353", "MMRF_2373", "MMRF_2378", "MMRF_2400", 
#                        "MMRF_2409", "MMRF_2422", "MMRF_2440", "MMRF_2445", "MMRF_2455", 
#                        "MMRF_2458", "MMRF_2461", "MMRF_2469", "MMRF_2471", "MMRF_2475", 
#                        "MMRF_2476", "MMRF_2480", "MMRF_2485", "MMRF_2489", "MMRF_2497", 
#                        "MMRF_2506", "MMRF_2514", "MMRF_2516", "MMRF_2523", "MMRF_2525", 
#                        "MMRF_2531", "MMRF_2532", "MMRF_2538", "MMRF_2539", "MMRF_2549", 
#                        "MMRF_2554", "MMRF_2557", "MMRF_2562", "MMRF_2574", "MMRF_2579", 
#                        "MMRF_2581", "MMRF_2595", "MMRF_2614", "MMRF_2621"))
# 
# ## Basic summary statistics
# summary_stats <- summary(plasma_4_df$Count_per_patient)
# ## Mean
# mean_count <- mean(plasma_4_df$Count_per_patient, na.rm = TRUE)
# ## Standard Deviation
# sd_count <- sd(plasma_4_df$Count_per_patient, na.rm = TRUE)
# ## Specific quantiles
# quantiles_count <- quantile(plasma_4_df$Count_per_patient, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
# 
# ## Printing the results
# print("Summary Statistics:")
# print(summary_stats)
# print(paste("Mean:", mean_count))
# print(paste("Standard Deviation:", sd_count))
# print("Quantiles (25%, 50%, 75%):")
# print(quantiles_count)
# 
# ## need to get patients that have at least 10 cells
# # Filtering plasma_4_df to include only rows where Count_per_patient >= 10
# plasma_4_df_filtered <- plasma_4_df %>%
#   filter(Count_per_patient >= 10)
# # Comes down to n=72 for patients with plasma.4 >= 10
# plasma_4_IDs_2 <- list(c("MMRF_1251", "MMRF_1286", "MMRF_1338", "MMRF_1352", "MMRF_1359", 
#                          "MMRF_1401", "MMRF_1423", "MMRF_1433", "MMRF_1445", "MMRF_1454", 
#                          "MMRF_1491", "MMRF_1497", "MMRF_1522", "MMRF_1523", "MMRF_1539", 
#                          "MMRF_1542", "MMRF_1560", "MMRF_1565", "MMRF_1586", "MMRF_1613", 
#                          "MMRF_1648", "MMRF_1670", "MMRF_1677", "MMRF_1719", "MMRF_1722", 
#                          "MMRF_1778", "MMRF_1780", "MMRF_1787", "MMRF_1832", "MMRF_1839", 
#                          "MMRF_1918", "MMRF_1965", "MMRF_1977", "MMRF_1981", "MMRF_2035", 
#                          "MMRF_2076", "MMRF_2106", "MMRF_2126", "MMRF_2143", "MMRF_2172", 
#                          "MMRF_2176", "MMRF_2187", "MMRF_2231", "MMRF_2250", "MMRF_2258", 
#                          "MMRF_2272", "MMRF_2290", "MMRF_2297", "MMRF_2301", "MMRF_2329", 
#                          "MMRF_2331", "MMRF_2339", "MMRF_2373", "MMRF_2400", "MMRF_2445", 
#                          "MMRF_2469", "MMRF_2471", "MMRF_2485", "MMRF_2523", "MMRF_2525", 
#                          "MMRF_2531", "MMRF_2532", "MMRF_2538", "MMRF_2539", "MMRF_2554", 
#                          "MMRF_2557", "MMRF_2562", "MMRF_2574", "MMRF_2579", "MMRF_2581", 
#                          "MMRF_2595", "MMRF_2621"))
# 
# ## Filter the dataframe for rows where subcluster_V03072023 is 'Plasma.6'
# plasma_6_df <- PSN_percentage %>%
#   filter(subcluster_V03072023 == "Plasma.6")
# ## 
# plasma_6_IDs <- list(c("MMRF_1089", "MMRF_1251", "MMRF_1267", "MMRF_1286", "MMRF_1323", 
#                        "MMRF_1338", "MMRF_1352", "MMRF_1359", "MMRF_1371", "MMRF_1391", 
#                        "MMRF_1395", "MMRF_1401", "MMRF_1408", "MMRF_1423", "MMRF_1424", 
#                        "MMRF_1433", "MMRF_1434", "MMRF_1436", "MMRF_1445", "MMRF_1454", 
#                        "MMRF_1461", "MMRF_1466", "MMRF_1484", "MMRF_1485", "MMRF_1497", 
#                        "MMRF_1510", "MMRF_1519", "MMRF_1522", "MMRF_1523", "MMRF_1537", 
#                        "MMRF_1538", "MMRF_1539", "MMRF_1542", "MMRF_1560", "MMRF_1565", 
#                        "MMRF_1572", "MMRF_1586", "MMRF_1613", "MMRF_1614", "MMRF_1627", 
#                        "MMRF_1634", "MMRF_1639", "MMRF_1643", "MMRF_1648", "MMRF_1659", 
#                        "MMRF_1670", "MMRF_1677", "MMRF_1719", "MMRF_1722", "MMRF_1755", 
#                        "MMRF_1768", "MMRF_1778", "MMRF_1780", "MMRF_1787", "MMRF_1823", 
#                        "MMRF_1832", "MMRF_1833", "MMRF_1837", "MMRF_1839", "MMRF_1850", 
#                        "MMRF_1865", "MMRF_1880", "MMRF_1918", "MMRF_1965", "MMRF_1977", 
#                        "MMRF_1981", "MMRF_1988", "MMRF_1991", "MMRF_2018", "MMRF_2035", 
#                        "MMRF_2055", "MMRF_2057", "MMRF_2062", "MMRF_2076", "MMRF_2106", 
#                        "MMRF_2122", "MMRF_2126", "MMRF_2143", "MMRF_2172", "MMRF_2176", 
#                        "MMRF_2183", "MMRF_2187", "MMRF_2204", "MMRF_2211", "MMRF_2232", 
#                        "MMRF_2250", "MMRF_2258", "MMRF_2272", "MMRF_2290", "MMRF_2291", 
#                        "MMRF_2297", "MMRF_2301", "MMRF_2326", "MMRF_2329", "MMRF_2330", 
#                        "MMRF_2331", "MMRF_2339", "MMRF_2353", "MMRF_2373", "MMRF_2378", 
#                        "MMRF_2400", "MMRF_2409", "MMRF_2436", "MMRF_2440", "MMRF_2445", 
#                        "MMRF_2455", "MMRF_2458", "MMRF_2461", "MMRF_2464", "MMRF_2469", 
#                        "MMRF_2471", "MMRF_2475", "MMRF_2476", "MMRF_2485", "MMRF_2505", 
#                        "MMRF_2506", "MMRF_2516", "MMRF_2523", "MMRF_2525", "MMRF_2531", 
#                        "MMRF_2532", "MMRF_2538", "MMRF_2539", "MMRF_2549", "MMRF_2554", 
#                        "MMRF_2557", "MMRF_2562", "MMRF_2568", "MMRF_2574", "MMRF_2579", 
#                        "MMRF_2581", "MMRF_2595", "MMRF_2614", "MMRF_2621"))
# 
# ## Basic summary statistics
# summary_stats <- summary(plasma_6_df$Count_per_patient)
# ## Mean
# mean_count <- mean(plasma_6_df$Count_per_patient, na.rm = TRUE)
# ## Standard Deviation
# sd_count <- sd(plasma_6_df$Count_per_patient, na.rm = TRUE)
# ## Specific quantiles
# quantiles_count <- quantile(plasma_6_df$Count_per_patient, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
# 
# ## Printing the results
# print("Summary Statistics:")
# print(summary_stats)
# print(paste("Mean:", mean_count))
# print(paste("Standard Deviation:", sd_count))
# print("Quantiles (25%, 50%, 75%):")
# print(quantiles_count)
# 
# ## need to get patients that have at least 10 cells
# # Filtering plasma_6_df to include only rows where Count_per_patient >= 10
# plasma_6_df_filtered <- plasma_6_df %>%
#   filter(Count_per_patient >= 10)
# # Comes down to n=71 for patients with plasma.6 >= 10
# plasma_6_IDs_2 <- list(c("MMRF_1251", "MMRF_1286", "MMRF_1338", "MMRF_1359", "MMRF_1395", 
#                          "MMRF_1423", "MMRF_1433", "MMRF_1434", "MMRF_1445", "MMRF_1484", 
#                          "MMRF_1510", "MMRF_1522", "MMRF_1523", "MMRF_1538", "MMRF_1539", 
#                          "MMRF_1542", "MMRF_1560", "MMRF_1565", "MMRF_1586", "MMRF_1613", 
#                          "MMRF_1648", "MMRF_1670", "MMRF_1677", "MMRF_1722", "MMRF_1778", 
#                          "MMRF_1780", "MMRF_1787", "MMRF_1839", "MMRF_1981", "MMRF_2062", 
#                          "MMRF_2126", "MMRF_2143", "MMRF_2172", "MMRF_2176", "MMRF_2187", 
#                          "MMRF_2250", "MMRF_2272", "MMRF_2290", "MMRF_2297", "MMRF_2301", 
#                          "MMRF_2329", "MMRF_2331", "MMRF_2339", "MMRF_2445", "MMRF_2455", 
#                          "MMRF_2469", "MMRF_2471", "MMRF_2485", "MMRF_2523", "MMRF_2525", 
#                          "MMRF_2531", "MMRF_2538", "MMRF_2539", "MMRF_2554", "MMRF_2562", 
#                          "MMRF_2568", "MMRF_2574", "MMRF_2579", "MMRF_2581", "MMRF_2595", 
#                          "MMRF_2621"))
# 
# # Base UMAP plot with distinct colors for each cluster
# base_umap <- DimPlot(seurat_obj, 
#                      reduction = "umap", 
#                      label = TRUE, 
#                      pt.size = 0.5, 
#                      label.size = 1.5, 
#                      raster = FALSE) + 
#   NoLegend()
# 
# View(base_umap)
# 
# 
# # Feature plot with gradient scale for GEP70_score
# feature_gradient <- FeaturePlot(seurat_obj, features = "DF_doublet_score", cols = c("lightblue", "darkblue"), raster = FALSE) +
#   scale_color_gradient(low = "yellow", high = "red") +  # Adjust colors as needed
#   theme(legend.position = "right")
# 
# # Combine the plots
# combined_plot <- base_umap + feature_gradient
# 
# # Print the combined plot
# combined_plot
