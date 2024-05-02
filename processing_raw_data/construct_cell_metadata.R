# William Pilcher
# Emory University
# Updated November 10th, 2023

# Tables are protected information - please contact the MMRF for raw tables or clinical metadata

### Build Patient and scRNA-seq Metadata Tables
# This requires the "outcomes", "visit", "per_patient", "risk", and "trt_induction" tables, along with the aliquot metadata table, and the "3/24" release of the WashU metadata.


library(ggplot2)
library(tidyverse)
library(dplyr)
library(kableExtra)
library(readxl)

## OPTIONS AND SETUP
# This will assert the counts for skerget, davies, and progression match what is previously known. If these should be different... then update the assert statements or set to false.
VALIDATE_COUNTS_FOR_RISK_AND_PROGRESSION <- TRUE

# This will load the 3/24 MD file of the object, and format the single-cell data as well as the clinical
BUILD_PER_CELL_MD_FILE <- TRUE


## PATHS TO METADATA
# Note - if changed to xlsx, change code to readxl::read_excel(file) |> as.data.frame(s)
# Change to directory with all filesf
clinicalMetadataPath <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata_1110/"

outcomes_csv <- "mmrf_commpass_IA22_active_standardized_outcomes_20230414.csv"
visit_csv <- "mmrf_commpass_IA22_active_standardized_per_patient_visit_20230414.csv"
perpat_csv <- "mmrf_commpass_IA22_active_standardized_per_patient_20230505.csv"
risk_xlsx <- "mmrf_commpass_IA21_active_standardized_risk_20231107.xlsx"
therapy_csv <- "mmrf_commpass_IA22_active_standardized_trt_induction_20230616.xlsx"

# Contains missing batch 3 MSSM samples (4 samples)
aliquot_csv <- "aliquot_metadata_20230615_2471_fix_batches_fixed.csv"
washuMD_file <- "derivative_data_WASHU_361_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_metadata_cln_V2_03242023.rds"
washuMD <- file.path(clinicalMetadataPath, washuMD_file)

## PATHS TO OUTPUT
# will save as rds and csv - dont put file ext
per_aliquot_output_name <- "1110_per_aliquot_md"
per_aliquot_output_INFO <- "1110_per_aliquot_VARIABLE_INFO.csv"

per_aliquot_output_dir <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata_1110/per_aliquot_metadata/"
dir.create(per_aliquot_output_dir)
# will save as rds and csv - dont put file ext
single_cell_output_name <- "1110_per_cell_md"
single_cell_output_INFO <- "1110_per_cell_VARIABLE_INFO.csv"

single_cell_output_dir <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata_1110/per_cell_metadata/"
dir.create(single_cell_output_dir)

## Aliquot Metadata
aliquotMD <- read.csv(file = file.path(clinicalMetadataPath, aliquot_csv), header = T)
aliquotMD <- aliquotMD %>% mutate(ID.assay = paste0(aliquot_id, ".", assay))
rownames(aliquotMD) <- aliquotMD$ID.assay
aliquotMD_scrna <- subset(aliquotMD, subset = (assay == "SCR"))
rownames(aliquotMD_scrna) <- aliquotMD_scrna$aliquot_id
visit_id_list <- unique(aliquotMD[, "d_visit_specimen_id"])
mmrfid_list <- unique(aliquotMD[, "public_id"]) # Should be around 345 (one pat may be missing)
print(paste0("Num Visits (All Atlas): ", length(visit_id_list))) # Should be arround 491 (one sample may have failed)
print(paste0("Num Patients (All Atlas): ", length(mmrfid_list)))

## Per Patient Metadata
# Demographic information, by patient
# Indexed by public_id
patMD <- read.csv(file = file.path(clinicalMetadataPath, perpat_csv), header = T)
patMD <- patMD # %>% subset(subset = (public_id %in% mmrfid_list))
rownames(patMD) <- patMD$public_id

# D_PT_mmstatus_22 is listed twice. Change to intended variable name
err_cols <- which("D_PT_mmstatus_22.1" == colnames(patMD))
colnames(patMD)[err_cols] <- "D_PT_mmstatus2_22"



## Standardized TRT (First Induction) (Curated)
firstInductionMD <- readxl::read_excel(file.path(clinicalMetadataPath, therapy_csv)) |> as.data.frame()
rownames(firstInductionMD) <- firstInductionMD$public_id
firstInductionMD <- firstInductionMD # %>% subset(subset = (public_id %in% mmrfid_list))



## Standaridized Risk
# 7/14 - now an xlsx
# Indexable by public_id (if no sample avaliable) or d_visit_specimen_id (if sample avaliable)
# atlas_subject to filter to things in our dataset Only
# Calls are now yes/no
library(readxl)
riskMD <- readxl::read_excel(file.path(clinicalMetadataPath, risk_xlsx), na = c("", "NA")) |> as.data.frame()
# riskMD <- subset(riskMD, subset = (public_id %in% mmrfid_list))

# riskMD[riskMD == "NA"] <- NA # Make 'NA' Values actually NA
# Fix string -> numeric...
callVariables <- c("nsd2_call", "ccnd3_call", "myc_call", "mafa_call", "ccnd1_call", "ccnd2_call", "maf_call", "mafb_call", "seqwgs_cp_hyperdiploid_call", "seqwgs_cp_1q21_call", "seqwgs_cp_13q14_call", "seqwgs_cp_13q34_call", "seqwgs_cp_17p13_call", "seqwgs_cp_cdkn2c_call", "seqwgs_cp_fam46c_call", "seqwgs_cp_rb1_call", "seqwgs_cp_tp53_call", "seqwgs_cp_traf3_call")
countVariables <- c("seqwgs_cp_hyperdiploid_chr_count", "seqwgs_cp_1q21", "seqwgs_cp_13q14", "seqwgs_cp_13q34", "seqwgs_cp_17p13", "seqwgs_cp_cdkn2c", "seqwgs_cp_fam46c", "seqwgs_cp_rb1", "seqwgs_cp_tp53", "seqwgs_cp_traf3")
# riskMD[countVariables] <- lapply(riskMD[countVariables], as.numeric) # Removes strings from numeric fields.

riskMD <- riskMD |> dplyr::filter(reason_for_collection == "Baseline", sample_type == "BM")


## Outcomes Metadata
# * Outcomes Metadata
#     + Contains information PER PATIENT (public_id) regarding outcomes, curated for the MMRF.
#         - Note - in 4_14 release, this is lower case. Will be upper case in document to match elsewhere
#     + Each row is indexed by public_id.
#     + Preview of table and overview of select variables below.

outcomesMD <- read.csv(file = file.path(clinicalMetadataPath, outcomes_csv), header = T)
outcomesMD$public_id <- outcomesMD$public_id %>% toupper()
outcomesMD <- outcomesMD # %>% subset(subset = (public_id %in% mmrfid_list))


library(gridExtra)
outcomesMD[, "pt_last_day_on_study_NAIS0"] <- outcomesMD[, "pt_last_day_on_study"] # For patients who are immediately lost to followup (pt_last_day_on_study is NA)
outcomesMD[is.na(outcomesMD[, "pt_last_day_on_study_NAIS0"]), "pt_last_day_on_study_NAIS0"] <- 0

# Corrections for patient with offset start day (d_lot_1_start_day = 1010)
outcomesMD <- outcomesMD %>% mutate(days_observed = (pt_last_day_on_study_NAIS0 - d_lot_1_start_day + 1))
outcomesMD <- mutate(outcomesMD, ttcpfs = (pfscdy - d_lot_1_start_day + 1))
outcomesMD <- mutate(outcomesMD, ttcos = (oscdy - d_lot_1_start_day + 1))

## Creating Progression Groups
# using !is.na(ttpfs) as this matches censpfs's definition of a progression event (includes deaths). Less strict than pdflag which p,,ots deatjs/
RP_thresh <- 548
NP_thresh <- 1460
outcomesMD <- outcomesMD %>% mutate(has_progression = !is.na(ttpfs))
outcomesMD <- outcomesMD %>% mutate(RP_threshold = (has_progression == 1 & ttpfs < RP_thresh))
outcomesMD <- outcomesMD %>% mutate(observation_threshold = (days_observed > NP_thresh))

outcomesMD <- outcomesMD %>% mutate(progression_group = dplyr::case_when(
    has_progression == 1 & ttpfs < RP_thresh ~ "RP",
    days_observed > NP_thresh & (has_progression == 0 | ttpfs >= NP_thresh) ~ "NP",
    has_progression == 1 & ttpfs >= RP_thresh & ttpfs < NP_thresh ~ "P",
    has_progression == 0 & days_observed < NP_thresh ~ "Inc",
    .default = "Inc" # setting to Inc for the one unhandled case - patient PFS data > 4 years but last observed < 4 years (likely pinned pfs to day of death after last day on study)
))


## Visit Metadata
# Contains information about each visit.
visitMD <- read.csv(file = file.path(clinicalMetadataPath, visit_csv), header = T)
visitMD <- visitMD # %>% subset(subset = (public_id %in% mmrfid_list))

# One column name is duplicated
err_lines <- which("BONE_DAYOFCOLLECT.1" == colnames(visitMD))
colnames(visitMD)[[err_lines]] <- "BONE_DAYOFCOLLECT_DUPLICATED"
visitMD_samplesOnly <- visitMD # %>% subset(subset = (d_pt_specimen_visit %in% visit_id_list))


## Selecting Variables From Metadata Objects
# Outcome Variables
pat_metadata <- outcomesMD[, c("public_id", "progression_group", "days_observed", "pdflag", "ttpfs", "censpfs", "ttcpfs", "pfscdy", "ttos", "censos", "ttcos", "oscdy", "d_pt_dc_reason", "d_pt_death_cause", "pfscdy", "d_lot_1_start_day", "pt_last_day_on_study")]
rownames(pat_metadata) <- pat_metadata$public_id

# Demographics
demog <- patMD[, c("public_id", "d_pt_race_1", "d_pt_ethnicity", "d_pt_sex", "d_dx_amm_age", "d_pt_height_cm", "d_dx_amm_weight_kg", "d_dx_amm_bmi", "d_tx_amm_regimen_1_cat")]
rownames(demog) <- demog$public_id

# ASCT, assesment
therapy_and_assesment <- patMD[, c("public_id", "d_dx_amm_iss_stage", "d_dx_amm_ecog", "d_amm_tx_asct_ever", "d_amm_tx_asct_1st", "IMWG_RISK_CLASS", "sct_elig")]
rownames(therapy_and_assesment) <- therapy_and_assesment$public_id

# First line Induction
first_induction_therapy_info <- firstInductionMD[, c("public_id", "d_tx_induction", "d_tx_induction_cat", "d_tx_induction_agents_num")]


## Creating baseline patient table
all_tables <- list(pat_metadata, demog, therapy_and_assesment, first_induction_therapy_info)
all_pat_metadata <- all_tables %>% reduce(full_join, by = "public_id")

all_pat_metadata <- all_pat_metadata %>% dplyr::mutate(has_imid = sapply(d_tx_amm_regimen_1_cat, function(x) grepl("imid", x, fixed = T)), has_pi = sapply(d_tx_amm_regimen_1_cat, function(x) grepl("pi", x, fixed = T)), has_steroid = sapply(d_tx_amm_regimen_1_cat, function(x) grepl("steroid", x, fixed = T)))

all_pat_metadata <- all_pat_metadata |> dplyr::mutate(therapy_class = dplyr::case_when(
    has_pi + has_imid + has_steroid == 3 ~ "Triplet Therapy",
    has_pi + has_imid + has_steroid == 2 ~ "Doublet Therapy",
    .default = "Singlet/Unknown"
))

risk_select <- riskMD |>
    dplyr::select(
        "public_id",
        "has_tx_data",
        "has_cn_data",
        "d_visit_specimen_id",
        "davies_based_risk",
        "skerget_based_risk", all_of(callVariables)
    ) |>
    dplyr::distinct()

all_with_risk <- dplyr::left_join(all_pat_metadata, risk_select, by = "public_id")

all_with_risk$Cohort <- "CoMMpass"

scRNAseq_MD_324 <- readRDS(washuMD)
active_aliquots <- scRNAseq_MD_324 |> dplyr::distinct(aliquot_id)
active_pat <- dplyr::left_join(active_aliquots, aliquotMD |> dplyr::filter(assay == "SCR"))

all_with_risk_IA <- all_with_risk |> dplyr::filter(public_id %in% active_pat$public_id)
all_with_risk_IA$Cohort <- "ImmuneAtlas"

all_combined <- rbind(all_with_risk_IA, all_with_risk)



library(finalfit)
library(gt)
library(gtsummary)
explanatory <- c("d_dx_amm_age", "d_dx_amm_weight_kg", "d_dx_amm_bmi", "d_pt_race_1", "d_pt_ethnicity", "d_pt_sex", "progression_group", "d_dx_amm_iss_stage", "d_dx_amm_ecog", "d_amm_tx_asct_1st", "d_tx_induction_cat", "davies_based_risk") # , "nsd2_call", "ccnd3_call", "myc_call", "mafa_call", "ccnd1_call", "ccnd2_call", "maf_call", "mafb_call", "seqwgs_cp_hyperdiploid_call", "seqwgs_cp_1q21_call", "seqwgs_cp_13q14_call", "seqwgs_cp_13q34_call", "seqwgs_cp_17p13_call", "seqwgs_cp_cdkn2c_call", "seqwgs_cp_fam46c_call", "seqwgs_cp_rb1_call", "seqwgs_cp_tp53_call", "seqwgs_cp_traf3_call")
dependent <- "Cohort"



all_combined <- all_combined |> dplyr::mutate(davies_based_risk = dplyr::case_when(
    davies_based_risk == "not_calculable" ~ NA,
    .default = davies_based_risk
))
all_combined <- all_combined |> dplyr::mutate(therapy_class = dplyr::case_when(
    therapy_class == "Singlet/Unknown" ~ NA,
    .default = therapy_class
))
all_combined <- all_combined |> dplyr::mutate(skerget_based_risk = dplyr::case_when(
    skerget_based_risk == "not_calculable" ~ NA,
    .default = skerget_based_risk
))

all_combined <- all_combined |> dplyr::mutate(d_tx_amm_regimen_1_cat = dplyr::case_when(
    d_tx_amm_regimen_1_cat == "" ~ NA,
    .default = d_tx_amm_regimen_1_cat
))
all_combined <- all_combined |> dplyr::mutate(d_pt_ethnicity = dplyr::case_when(
    d_pt_ethnicity == "" ~ NA,
    d_pt_ethnicity == "3" ~ "Unknown",
    .default = d_pt_ethnicity
))
all_combined <- all_combined |> dplyr::mutate(d_pt_race_1 = dplyr::case_when(
    d_pt_race_1 == "" ~ NA,
    .default = d_pt_race_1
))

all_combined$d_dx_amm_ecog <- factor(all_combined$d_dx_amm_ecog, levels = c(0, 1, 2, 3, 4))


# all_combined |> summary_factorlist(dependent, explanatory, p = T, add_dependent_label = T, na_to_prop = F, add_col_totals = T, na_include = T, na_to_p = F) -> t1
t1 <- all_combined |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, progression_group, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st, davies_based_risk),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", progression_group ~ "Progression Group", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received", davies_based_risk ~ "Davies Risk Score"),
        by = Cohort, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 1. Immune Atlas Cohort Information and Comparison to CoMMpass**")

t1 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table1.docx")
# knitr::kable(t1, align = c("l", "l", "r", "r", "r")) |> kableExtra::save_kable(file = "/opt/localdata/wcpilch/test_table.html", bs_theme = "paper")


t2 <- all_combined |>
    dplyr::filter(davies_based_risk %in% c("high_risk", "standard_risk"), Cohort == "ImmuneAtlas") |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, progression_group, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", progression_group ~ "Progression Group", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received"),
        by = davies_based_risk, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 2. Clinical Characteristics of High and Standard Risk Patients**")

t2 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table2.docx")



# knitr::kable(t2, align = c("l", "l", "r", "r", "r")) |> kableExtra::save_kable(file = "/opt/localdata/wcpilch/test_table_risk.html", bs_theme = "paper")
# explanatory_p <- explanatory[explanatory != "progression_group"]

t3 <- all_combined |>
    dplyr::filter(progression_group %in% c("NP", "P", "RP"), Cohort == "ImmuneAtlas") |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st, davies_based_risk),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received", davies_based_risk ~ "Davies Risk Score"),
        by = progression_group, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 3. Clinical Characteristics by Progression Group**")

t3 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table3.docx")

t3 <- all_combined |>
    dplyr::filter(progression_group %in% c("NP", "RP"), Cohort == "ImmuneAtlas") |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st, davies_based_risk),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received", davies_based_risk ~ "Davies Risk Score"),
        by = progression_group, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 3. Clinical Characteristics by Progression Group**")

t3 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table3_NPRP_ONLY.docx")

t3 <- all_combined |>
    dplyr::filter(progression_group %in% c("NP", "P", "RP"), d_tx_induction_cat %in% c("imid_pi_steroid", "chemo_imid_pi_steroid"), Cohort == "ImmuneAtlas") |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st, davies_based_risk),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received", davies_based_risk ~ "Davies Risk Score"),
        by = progression_group, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 3. Clinical Characteristics by Progression Group**")

t3 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table3_triplet_only.docx")

t3 <- all_combined |>
    dplyr::filter(progression_group %in% c("NP", "P", "RP"), davies_based_risk %in% c("standard_risk"), Cohort == "ImmuneAtlas") |>
    gtsummary::tbl_summary(
        include = c(d_dx_amm_age, d_dx_amm_weight_kg, d_dx_amm_bmi, d_pt_race_1, d_pt_ethnicity, d_pt_sex, d_dx_amm_iss_stage, d_dx_amm_ecog, d_tx_induction_cat, d_amm_tx_asct_1st),
        label = list(d_dx_amm_age ~ "Age", d_dx_amm_weight_kg ~ "Weight (kg)", d_dx_amm_bmi ~ "BMI", d_pt_race_1 ~ "Race", d_pt_ethnicity ~ "Ethnicity", d_pt_sex ~ "Sex", d_dx_amm_iss_stage ~ "ISS Stage", d_dx_amm_ecog ~ "ECOG", d_tx_induction_cat ~ "1st Line Therapy Class", d_amm_tx_asct_1st ~ "1st Line ASCT Received"),
        by = progression_group, missing = "ifany", percent = "column", missing_text = "Missing"
    ) |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    modify_footnote(
        all_stat_cols() ~ "Median (IQR) or Count (%)"
    ) |>
    bold_labels() |>
    modify_caption("**Supplementary Table 3. Clinical Characteristics by Progression Group**")

t3 |>
    as_gt() |>
    gt::gtsave(filename = "/opt/localdata/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/Supplementary_Table3_SR_ONLY.docx")



library(ggstatsplot)
library(ggpubr)
library(forestmodel)
library(gtsummary)
library(survival)
library(ggsurvfit)
library(contsurvplot)
library(ggprism)
tmp <- all_pat_metadata |> dplyr::filter(therapy_class %in% c("Triplet Therapy", "Doublet Therapy"))
tmp <- tmp |> dplyr::filter(d_lot_1_start_day < 2)
p_surv_a <- survfit2(Surv(pfscdy, censpfs) ~ therapy_class, data = tmp) %>% ggsurvfit(size = 2) +
    add_confidence_interval() +
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
        x = "Follow-up Time, Days",
        y = "Percentage Progression Free Survival"
    ) +
    add_legend_title(title = "Therapy Class") +
    coord_cartesian(xlim = c(0, 3000)) +
    add_censor_mark(size = 3) +
    add_risktable(
        risktable_height = 0.1,
        size = 3,
        theme = # increase font size of risk table title and y-axis label
            list(
                theme_risktable_default(
                    axis.text.y.size = 11,
                    plot.title.size = 11
                ),
                theme(plot.title = element_text(face = "bold"))
            )
    ) +
    theme_prism() +
    ggtitle(paste0("Progression Free Survival and Therapy Class"))
p_surv_a <- p_surv_a + ggsurvfit::add_pvalue(caption = "Log-rank {p.value}")

# Doublet Median - 656
# Triplet Median - 1325

# scale_color_manual(values = c("Low" = color_low, "High" = color_high))
p_surv_a <- ggsurvfit_build(p_surv_a)

ggsave(filename = "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/all_commpass_therapy_survival.png", plot = p_surv_a, units = "in", width = 3.2, height = 3, scale = 3)


library(ggstatsplot)
library(ggpubr)
library(forestmodel)
library(gtsummary)
library(survival)
library(ggsurvfit)
library(contsurvplot)
library(ggprism)
tmp <- all_pat_metadata |> dplyr::filter(d_lot_1_start_day < 2)
p_surv_a <- survfit2(Surv(pfscdy, censpfs) ~ 1, data = tmp) %>% ggsurvfit(size = 2) +
    add_confidence_interval() +
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    add_quantile(y_value = 0.375, color = "gray50", linewidth = 0.75) +
    add_quantile(y_value = 0.25, color = "gray50", linewidth = 0.75) +
    add_quantile(y_value = 0.625, color = "gray50", linewidth = 0.75) +
    add_quantile(y_value = 0.75, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
        x = "Follow-up Time, Days",
        y = "Percentage Progression Free Survival"
    ) +
    coord_cartesian(xlim = c(0, 3000)) +
    add_censor_mark(size = 3) +
    add_risktable(
        risktable_height = 0.1,
        size = 3,
        theme = # increase font size of risk table title and y-axis label
            list(
                theme_risktable_default(
                    axis.text.y.size = 11,
                    plot.title.size = 11
                ),
                theme(plot.title = element_text(face = "bold"))
            )
    ) +
    theme_prism() +
    ggtitle(paste0("Progression Free Survival"))
p_surv_a <- ggsurvfit_build(p_surv_a)

ggsave(filename = "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/all_commpass_survival.png", plot = p_surv_a, units = "in", width = 3.2, height = 3, scale = 3)

library(ggstatsplot)
library(ggpubr)
library(forestmodel)
library(gtsummary)
library(survival)
library(ggsurvfit)
library(contsurvplot)
library(ggprism)
tmp <- all_with_risk |> dplyr::filter(davies_based_risk %in% c("high_risk", "standard_risk"))
tmp <- tmp |> dplyr::filter(d_lot_1_start_day < 2)
p_surv_a <- survfit2(Surv(pfscdy, censpfs) ~ davies_based_risk, data = tmp) %>% ggsurvfit(size = 2) +
    add_confidence_interval() +
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
        x = "Follow-up Time, Days",
        y = "Percentage Progression Free Survival"
    ) +
    add_legend_title(title = "Cytogenetic Risk Class") +
    coord_cartesian(xlim = c(0, 3000)) +
    add_censor_mark(size = 3) +
    add_risktable(
        risktable_height = 0.1,
        size = 3,
        theme = # increase font size of risk table title and y-axis label
            list(
                theme_risktable_default(
                    axis.text.y.size = 11,
                    plot.title.size = 11
                ),
                theme(plot.title = element_text(face = "bold"))
            )
    ) +
    theme_prism() +
    ggtitle(paste0("Progression Free Survival and Cytogenetic Class"))
p_surv_a <- p_surv_a + ggsurvfit::add_pvalue(caption = "Log-rank {p.value}")

tmp <- all_with_risk |> dplyr::filter(davies_based_risk %in% c("high_risk", "standard_risk") & therapy_class == "Triplet Therapy")
tmp <- tmp |> dplyr::filter(d_lot_1_start_day < 2)
p_surv_a <- survfit2(Surv(pfscdy, censpfs) ~ davies_based_risk, data = tmp) %>% ggsurvfit(size = 2) +
    add_confidence_interval() +
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
        x = "Follow-up Time, Days",
        y = "Percentage Progression Free Survival"
    ) +
    add_legend_title(title = "Cytogenetic Risk Class") +
    coord_cartesian(xlim = c(0, 3000)) +
    add_censor_mark(size = 3) +
    add_risktable(
        risktable_height = 0.1,
        size = 3,
        theme = # increase font size of risk table title and y-axis label
            list(
                theme_risktable_default(
                    axis.text.y.size = 11,
                    plot.title.size = 11
                ),
                theme(plot.title = element_text(face = "bold"))
            )
    ) +
    theme_prism() +
    ggtitle(paste0("Progression Free Survival and Cytogenetic Class"))
p_surv_a <- p_surv_a + ggsurvfit::add_pvalue(caption = "Log-rank {p.value}")


# Doublet Median - 656
# Triplet Median - 1325

# scale_color_manual(values = c("Low" = color_low, "High" = color_high))
p_surv_a <- ggsurvfit_build(p_surv_a)

ggsave(filename = "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/all_commpass_risk_survival_triplet.png", plot = p_surv_a, units = "in", width = 3.2, height = 3, scale = 3)

p <- ggplot(tmp, aes(x = pfscdy, y = "black")) +
    geom_step(aes(y = ..y..), stat = "ecdf") +
    geom_segment(x = 0, xend = 1547, y = 0.75, yend = 0.75) +
    geom_segment(x = 1547, xend = 1547, y = 0, yend = 0.75) +
    geom_segment(x = 0, xend = 1460, y = 0.7250674, yend = 0.7250674, color = "blue") +
    geom_segment(x = 1460, xend = 1460, y = 0, yend = 0.7250674, color = "blue") +
    geom_segment(x = 0, xend = 365 * 5, y = 0.8050314, yend = 0.8050314, color = "orange") +
    geom_segment(x = 365 * 5, xend = 365 * 5, y = 0, yend = 0.8050314, color = "orange") +
    geom_segment(x = 0, xend = 365 * 6, y = 0.8050314, yend = 0.8050314, color = "red") +
    geom_segment(x = 365 * 6, xend = 365 * 6, y = 0, yend = 0.8050314, color = "red") +
    theme_prism() +
    coord_cartesian(ylim = c(0, 1), expand = F)

ggsave(filename = "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/output/all_compass_pfscdy.png", plot = p_surv_a, units = "in", width = 3.2, height = 3, scale = 3)

# Object-Specific Metadata (Simplified)
scRNAseq_MD_324 <- readRDS(washuMD)
md <- scRNAseq_MD_324

md$DO_NOT_USE_OLD_OBJECT_PUBLIC_ID <- md$public_id
md$DO_NOT_USE_OLD_IMMUNE_ANALYSIS_ID <- md$immune_analysis_id
md$DO_NOT_USE_OLD_CYTOGENETIC_RISK_1Q <- md$MMRF_Cytogenetic_Risk_1q
md <- md |>
    dplyr::select(aliquot_id, sample_id, sample_aliquot_id, DO_NOT_USE_OLD_OBJECT_PUBLIC_ID, DO_NOT_USE_OLD_IMMUNE_ANALYSIS_ID, DO_NOT_USE_OLD_CYTOGENETIC_RISK_1Q, visit_type, Study_Site, Batch, siteXbatch) |>
    dplyr::distinct()

# Add missing visit type for MMRF_2471_1
md$visit_type[is.na(md$visit_type)] <- "baseline_diagnosis"


aliquotMD_scrna <- aliquotMD_scrna |> dplyr::select("aliquot_id", "d_visit_specimen_id", "public_id", "collection_event", "batch")
md <- dplyr::left_join(md, aliquotMD_scrna, by = "aliquot_id")

all_pat_metadata_select <- all_pat_metadata |> dplyr::select(
    "public_id",
    "d_pt_race_1",
    "d_pt_ethnicity",
    "d_pt_sex",
    "d_dx_amm_age",
    "d_dx_amm_ecog",
    "d_dx_amm_iss_stage",
    "d_dx_amm_bmi",
    "d_pt_height_cm",
    "d_dx_amm_weight_kg",
    "d_tx_induction",
    "d_tx_induction_cat",
    "d_amm_tx_asct_1st",
    "d_amm_tx_asct_ever",
    "progression_group",
    "pt_last_day_on_study",
    "days_observed",
    "d_lot_1_start_day",
    "censpfs",
    "pfscdy",
    "ttcpfs",
    "censos",
    "oscdy",
    "ttcos",
    "pdflag",
    "ttpfs"
)

md <- dplyr::left_join(md, all_pat_metadata_select, by = "public_id")

visit_select <- visitMD_samplesOnly |>
    dplyr::select("VISIT", "VJ_INTERVAL", "VISITDY", "d_pt_specimen_visit") |>
    distinct()

fixcol <- which(colnames(visit_select) == "d_pt_specimen_visit")
colnames(visit_select)[[fixcol]] <- "d_visit_specimen_id"
md <- dplyr::left_join(md, visit_select, by = "d_visit_specimen_id")

riskMD$visit_in_risk_MD <- "yes"

# risk_select$davies_based_risk[is.na(risk_select$davies_based_risk)] <- "not_calculable"
# risk_select$skerget_based_risk[is.na(risk_select$skerget_based_risk)] <- "not_calculable"


md <- dplyr::left_join(md, risk_select, by = "d_visit_specimen_id")
md$has_tx_data[is.na(md$has_tx_data)] <- "no"
md$has_cn_data[is.na(md$has_cn_data)] <- "no"
md$davies_based_risk[is.na(md$davies_based_risk)] <- "no_risk_data"
md$skerget_based_risk[is.na(md$skerget_based_risk)] <- "no_risk_data"

if (VALIDATE_COUNTS_FOR_RISK_AND_PROGRESSION) {
    ## Verify Numbers
    baseline <- md |> dplyr::filter(VJ_INTERVAL == "Baseline")
    davies_table <- table(baseline$davies_based_risk)
    # stopifnot(davies_table["high_risk"] == 123)
    # stopifnot(davies_table["standard_risk"] == 108)

    # 123 HR, 108 SR, 3 NC

    skerget_table <- table(baseline$skerget_based_risk)
    # stopifnot(skerget_table["high_risk"] == 61)
    # stopifnot(skerget_table["standard_risk"] == 168)
    # 61 HR, 168 SR, 5 NC

    progression_table <- table(baseline$progression_group)
    stopifnot(progression_table["NP"] == 83)
    stopifnot(progression_table["RP"] == 67)
    # 83 NP, 67 RP, 71 P, 42 Inc

    numAliquots <- length(unique(md$aliquot_id))
    stopifnot(numAliquots == 361)
    # 361

    numPat <- length(unique(md$public_id))
    stopifnot(numPat == 263)
    # 263

    numVisits <- length(unique(md$d_visit_specimen_id))
    stopifnot(numVisits == 358)
    # 358
}

# These differ by one due to one visit/public_id being NA in the original metadata.
md$DO_NOT_USE_OLD_OBJECT_PUBLIC_ID[is.na(md$DO_NOT_USE_OLD_OBJECT_PUBLIC_ID)] <- "PUBLIC_ID_MISSING"
md$DO_NOT_USE_OLD_IMMUNE_ANALYSIS_ID[is.na(md$DO_NOT_USE_OLD_IMMUNE_ANALYSIS_ID)] <- "VISIT_ID_MISSING"
print("Should be 1 False, 360 True (One patient's pat id missing previously)")
table(md$public_id == md$DO_NOT_USE_OLD_OBJECT_PUBLIC_ID)
print("Should be 1 False, 360 True (One patient's pat id missing previously)")
table(md$d_visit_specimen_id == md$DO_NOT_USE_OLD_IMMUNE_ANALYSIS_ID)

# Reorder Columns
md <- md |> dplyr::select(
    aliquot_id, sample_id, d_visit_specimen_id, public_id, Study_Site, Batch, siteXbatch, collection_event, visit_type, VISIT, VJ_INTERVAL, VISITDY, # Sample Info
    d_pt_sex, d_pt_race_1, d_pt_ethnicity, d_dx_amm_age, d_dx_amm_ecog, d_dx_amm_iss_stage, d_dx_amm_bmi, d_pt_height_cm, d_dx_amm_weight_kg, # Demographics
    d_tx_induction, d_tx_induction_cat, d_amm_tx_asct_1st, progression_group, ttpfs, pdflag, censpfs, pfscdy, ttcpfs, censos, oscdy, ttcos, d_lot_1_start_day, days_observed, pt_last_day_on_study, # Progression
    has_tx_data, has_cn_data, davies_based_risk, skerget_based_risk, nsd2_call, ccnd3_call, myc_call, mafa_call, ccnd1_call, ccnd2_call, maf_call, mafb_call, seqwgs_cp_hyperdiploid_call, seqwgs_cp_1q21_call, seqwgs_cp_13q14_call, seqwgs_cp_13q34_call, seqwgs_cp_17p13_call, seqwgs_cp_cdkn2c_call, seqwgs_cp_fam46c_call, seqwgs_cp_rb1_call, seqwgs_cp_tp53_call, seqwgs_cp_traf3_call
)

# write.csv(md, file.path(per_aliquot_output_dir, paste0(per_aliquot_output_name, ".csv")))
# saveRDS(md, file.path(per_aliquot_output_dir, paste0(per_aliquot_output_name, ".rds")))

## Create Single Cell RNA Sequencing MD
scRNA_MD <- scRNAseq_MD_324
scRNA_MD$cellname <- rownames(scRNA_MD)
## Reorganize Subcluster Factors
# Note: These are from the Jan 20th 2023 release of subcluster metadata
sorting <- scRNA_MD$subcluster_V03072023
valuesToSort <- unique(sorting)
sortDF <- data.frame(valuesToSort)
sortDF <- sortDF %>% separate(valuesToSort, c("comp", "parent", "sub"), fill = "right") # This will give a harmless warning for the instances without a subcluster
sortDF$string <- valuesToSort
sortDF$sub[is.na(sortDF$sub)] <- 0
sortDF <- arrange(sortDF, comp, as.numeric(parent), as.numeric(sub))
fullNameOrder <- sortDF$string
scRNA_MD$subcluster_V03072023 <- factor(scRNA_MD$subcluster_V03072023, levels = fullNameOrder)
stopifnot(all(!is.na(scRNA_MD$subcluster_V03072023))) # if this is NA, your factors are messed up

# Split by '.' to create 'compartment' and 'parent' variables for subclusters
scRNA_MD$subcluster_V03072023_compartment <- unlist(lapply(strsplit(as.character(scRNA_MD$subcluster_V03072023), "\\."), "[[", 1))
parent_cluster <- unlist(lapply(strsplit(as.character(scRNA_MD$subcluster_V03072023), "\\."), "[[", 2))
scRNA_MD$subcluster_V03072023_parent <- paste0(scRNA_MD$subcluster_V03072023_compartment, ".", parent_cluster)

# convert 'compartment' and 'parent' variables into sorted factors.
sorting <- scRNA_MD$subcluster_V03072023_parent
valuesToSort <- unique(sorting)
sortDF <- data.frame(valuesToSort)
sortDF <- sortDF %>% separate(valuesToSort, c("comp", "parent")) # This will give a harmless warning for the instances without a subcluster
sortDF$string <- valuesToSort
sortDF <- arrange(sortDF, comp, as.numeric(parent))
fullNameOrder <- sortDF$string
scRNA_MD$subcluster_V03072023_parent <- factor(scRNA_MD$subcluster_V03072023_parent, levels = fullNameOrder)

sorting <- scRNA_MD$subcluster_V03072023_compartment
valuesToSort <- unique(sorting)
fullNameOrder <- sort(valuesToSort)
scRNA_MD$subcluster_V03072023_compartment <- factor(scRNA_MD$subcluster_V03072023_compartment, levels = fullNameOrder)

## Set Doublet clusters
# Base off of 'subcluster_cellcount_doublet' as of 7/17
# Recommended to remove 'doublet' and 'possible doublet'

isDoublet <- c(
    "BEry.13", "BEry.14", "BEry.5", "BEry.8",
    "Ery.10", "Ery.11", "Ery.7",
    "Myeloid.14", "Myeloid.7", "Myeloid.9",
    "NkT.15", "NkT.16",
    "Plasma.10", "Plasma.11"
)
isPlasma_B <- c("Plasma.13", "Plasma.22", "Plasma.23")
isPossibleDoulbet <- c("Ery.9", "Plasma.3", "Plasma.9")
isPossibleSinglet <- c("Plasma.7", "Plasma.8")

# check cluster names
all(isDoublet %in% scRNA_MD$subcluster_V03072023)
all(isPlasma_B %in% scRNA_MD$subcluster_V03072023)
all(isPossibleDoulbet %in% scRNA_MD$subcluster_V03072023)
all(isPossibleSinglet %in% scRNA_MD$subcluster_V03072023)

scRNA_MD$OLD_DOUBLET_PRED <- scRNA_MD$doublet_pred

scRNA_MD <- scRNA_MD |> dplyr::mutate(doublet_pred = dplyr::case_when(
    subcluster_V03072023 %in% isDoublet ~ "dblet_cluster",
    subcluster_V03072023 %in% isPlasma_B ~ "Plasma_B",
    subcluster_V03072023 %in% isPossibleDoulbet ~ "poss_dblet_cluster",
    subcluster_V03072023 %in% isPossibleSinglet ~ "poss_singlet_cluster",
    .default = "singlet"
))

print("doublet_pred previously incorrectly ommitted Plasma.10, leading to a difference of 3221 cells")
table(scRNA_MD$doublet_pred, scRNA_MD$OLD_DOUBLET_PRED)

scRNA_MD <- scRNA_MD |> dplyr::select(
    aliquot_id,
    sample_id,
    nCount_RNA,
    nFeature_RNA,
    RNA_snn_res.0.5,
    percent.mt,
    ambient_rho,
    predicted_cell_type,
    predicted_cell_type_probability,
    scrublet_doublet_score,
    scrublet_pred_dbl,
    DF_doublet_score,
    DF_pred_dbl,
    pegasus_doublet_score,
    pegasus_pred_dbl,
    doublet_pred,
    general_ct,
    compartment,
    seurat_clusters,
    subcluster_V03072023,
    subcluster_V03072023_compartment,
    subcluster_V03072023_parent,
    cellname
)

scRNA_MD <- dplyr::left_join(scRNA_MD, md)

rownames(scRNA_MD) <- scRNA_MD$cellname

## Checking values with known counts
TEST_VALUES <- scRNA_MD |>
    dplyr::filter(VJ_INTERVAL == "Baseline") |>
    dplyr::distinct(aliquot_id, d_visit_specimen_id, public_id, progression_group, davies_based_risk, skerget_based_risk)
print("Should be 263...")
nrow(TEST_VALUES)

print("Should be 42 Inc, 83 NP, 71 P, 67 RP")
table(TEST_VALUES$progression_group)

print("Should be 123 high_risk, 108 standard_risk")
table(TEST_VALUES$davies_based_risk)

print("Should be 61 high_risk, 167 standard_risk, 5 not calculable, 29 no data")
table(TEST_VALUES$skerget_based_risk)

print("Should be 263...")
length(unique(TEST_VALUES$public_id))

print("Should have no true values...")
table(is.na(TEST_VALUES$public_id))

print("Should be 358...")
length(unique(scRNA_MD$d_visit_specimen_id))

print("Should have no true values...")
table(is.na(scRNA_MD$d_visit_specimen_id))

print("Should be 361...")
length(unique(scRNA_MD$aliquot_id))

print("Should have no true values...")
table(is.na(scRNA_MD$aliquot_id))

# write.csv(scRNA_MD, file.path(single_cell_output_dir, paste0(single_cell_output_name, ".csv")))
# saveRDS(scRNA_MD, file.path(single_cell_output_dir, paste0(single_cell_output_name, ".rds")))

### Documenting all variables in object

quickTemplate <- function(dataSource, linkedby, dataFile, lastUpdated, notes = NULL) {
    if (is.null(notes)) {
        notes <- "See Dictionary/Wiki"
    }
    return(list(
        "Origin" = dataSource,
        "Linked by" = linkedby,
        "File" = dataFile,
        "lastUpdated" = lastUpdated,
        "notes" = notes
    ))
}


aliquotMD_date <- "6/15/23"
riskMD_date <- "11/10/23"
perPatMD_date <- "5/05/23"
therapyMD_date <- "6/16/23"
outcomesMD_date <- "4/14/23"
visitMD_date <- "4/14/23"
info <- list(
    "aliquot_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Aliquot ID in standard MMRF format derived from seurat object. Used to link other variables"),
    "sample_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Internal seurat object name containing the aliquotID and the public id SM indicates spike-in moues samples"),
    "sample_aliquot_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Internal seurat object name containing the aliquotID and the public id, with no indicator for spike-in"),
    "d_visit_specimen_id" = quickTemplate("Aliquot Metadata", "aliquot_id", aliquot_csv, aliquotMD_date, "Formerly visit_id or immmune_analysis_id. Linked from aliquot metadata. Sample identifier (from which multiple aliquots can be derived). Used to link risk and visit metadata"),
    "public_id" = quickTemplate("Aliquot Metadata", "aliquot_id", aliquot_csv, aliquotMD_date, "Public_id in standard MMRF Format identifying patients, linked from aliquot metadata table. Aliquot metadata table manually corrected to incoroporate four missing MSSM samples. Used to link outcome, treatment, and demographic metadata"),
    "Study_Site" = quickTemplate("Seurat Metadata", "aliquot_id", washuMD_file, "N/A", "Processing site from seurat metadata"),
    "Batch" = quickTemplate("Seurat Metadata", "aliquot_id", washuMD_file, "N/A", "Batch number from seurat metadata. This has the correct batch number information."),
    "siteXbatch" = quickTemplate("Seurat Metadata", "aliquot_id", washuMD_file, "N/A", "Concatenation of Study_Site and Batch for harmonization"),
    "collection_event" = quickTemplate("Aliquot Metadata", "d_visit_specimen_id", aliquot_csv, aliquotMD_date, "Reason for sample collection from aliquot metadata. Note one patient has a 'Screening' field instead of 'Baseline"),
    "visit_type" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Original visit_type from Seurat Metadata. Provenance unknown. Does differ for certain values"),
    "VISIT" = quickTemplate("Visit Metadata", "d_visit_specimen_id", visit_csv, visitMD_date, "Interval as an integer (0=baseline, -1 = EOS/discontiue)"),
    "VJ_INTERVAL" = quickTemplate("Visit Metadata", "d_visit_specimen_id", visit_csv, visitMD_date, "Visit Interval as a string"),
    "VISITDY" = quickTemplate("Visit Metadata", "d_visit_specimen_id", visit_csv, visitMD_date, "Day of Visit"),
    "d_pt_sex" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_pt_race_1" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_pt_ethnicity" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_dx_amm_age" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_dx_amm_ecog" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_dx_amm_iss_stage" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_dx_amm_bmi" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_pt_height_cm" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_dx_amm_weight_kg" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "d_tx_induction" = quickTemplate("1st Induction", "public_id", therapy_csv, therapyMD_date),
    "d_tx_induction_cat" = quickTemplate("1st Induction", "public_id", therapy_csv, therapyMD_date),
    "d_amm_tx_asct_1st" = quickTemplate("Per Patient", "public_id", perpat_csv, perPatMD_date),
    "progression_group" = quickTemplate("(Derived) Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Progression group. RP: Progression within 18 months. NP: No progression while under observation for 4 years. P: Progression observed between 18mo and 4yrs. Inc: No progression observed, but not observed for four years. Definition for 'progression' matches that used for censpfs."),
    "ttpfs" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Time to Progression. Includes progression events and deaths. NA means no progression event and patient is not (known to be) dead"),
    "pdflag" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "pdflag variable. Unlike censpfs, appears to omit deaths not associated with disease progression."),
    "censpfs" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date),
    "pfscdy" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Day of censoring, not time to"),
    "ttcpfs" = quickTemplate("(Derived) Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Use for survival. ttcpfs = pfscdy - d_lot_1_start_day + 1. Corrects for patient with offset start day."),
    "censos" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date),
    "oscdy" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Day of censoring, not time to"),
    "ttcos" = quickTemplate("(Derived) Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Use for survival. ttcos = oscdy - d_lot_1_start_day + 1. Corrects for patient with offset start day."),
    "d_lot_1_start_day" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Start day variable. 1 for all patients, except for one, who's value is 1010. Used to correct for 'days' for this patient"),
    "days_observed" = quickTemplate("(Derived) Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Correction for patient with offset start day: pt_last_day_on_study - d_lot_1_start_day + 1"),
    "pt_last_day_on_study" = quickTemplate("Outcomes", "public_id", outcomes_csv, outcomesMD_date, "Uncorrected value for days_observed - last day on study."),
    "has_tx_data" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Set to 'no' if visit not in risk table as well"),
    "has_cn_data" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Set to 'no' if visit not in risk table as well"),
    "davies_based_risk" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "See dictionairy/wiki for definition. This variable is TIME SENSITIVE and may differ between baseline and follow-ups. Follow-ups may not always have risk information. Values with no metadata in risk table will have value 'no_risk_data'"),
    "skerget_based_risk" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "See dictionairy/wiki for definition. This variable is TIME SENSITIVE and may differ between baseline and follow-ups. Follow-ups may not always have risk information. Values with no metadata in risk table will have value 'no_risk_data'"),
    "nsd2_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Manually corrected naming error in column"),
    "ccnd3_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Manually corrected naming error in column"),
    "myc_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "mafa_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "ccnd1_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Manually corrected naming error in column"),
    "ccnd2_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "Manually corrected naming error in column"),
    "maf_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "mafb_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_hyperdiploid_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date, "TRUE/FALSE swapped to yes/no to match others"),
    "seqwgs_cp_1q21_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_13q14_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_13q34_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_17p13_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_cdkn2c_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_fam46c_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_tp53_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date),
    "seqwgs_cp_traf3_call" = quickTemplate("Risk Metadata", "d_visit_specimen_id", risk_xlsx, riskMD_date)
)

df <- lapply(info, data.frame) |> bind_rows()
rownames(df) <- names(info)
df$variable_name <- rownames(df)
# write.csv(df, file.path(per_aliquot_output_dir, per_aliquot_output_INFO))

info_singleCell <- list(
    "aliquot_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Aliquot ID in standard MMRF format derived from seurat object. Used to link other variables"),
    "sample_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Internal seurat object name containing the aliquotID and the public id SM indicates spike-in moues samples"),
    "sample_aliquot_id" = quickTemplate("Seurat Metadata", "N/A", washuMD_file, "N/A", "Internal seurat object name containing the aliquotID and the public id, with no indicator for spike-in"),
    "nCount_RNA" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Number unique UMIs per Cell"),
    "nFeature_RNA" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Number unique Genes per Cell"),
    "RNA_snn_res.0.5" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Initial Clustering Results on Full Object"),
    "percent.mt" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Percent of UMIs mapped to mitochondrial transcripts per cell"),
    "ambient_rho" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "(?Unsure) Related to Ambient RNA"),
    "predicted_cell_type" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Automatic Cell Type Annotation"),
    "predicted_cell_type_probability" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Probability cell is of cell type"),
    "scrublet_doublet_score" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet score computed via Scrublet"),
    "scurublet_pred_dbl" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet call computed via Scrublet (False, True, or NA)"),
    "DF_doublet_score" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet score computed via Doublet Finder"),
    "DF_pred_dbl" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet call computed via Doublet Finder (Singlet, Doublet)"),
    "pegasus_doublet_score" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet score computed via Pegasus"),
    "pegasus_pred_dbl" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Doublet call via Pegasus (False, NA. No Trues present in column)"),
    "doublet_pred" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Group Consensus on clusters. Filter 'dblet_cluster' and 'poss_dblet_cluster'. Other values: 'singlet', 'poss_singlet_cluster', 'Plasma_B'. See: https://docs.google.com/spreadsheets/d/1VaosVSKbHF1A4I-8LHXsPYF1BleYYzrP1BSxNaIap8w/edit#gid=627540304"),
    "general_ct" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "(Undocumented) General Cell Type description of cluster"),
    "compartment" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Major Compartment cell type belongs to. Values are: B_ery (B + Erythroblast), Erythrocyte, Fibroblast, Myeloid, NK_T, Plasma"),
    "seurat_clusters" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Initial Seurat Clusters as Computed on Full Object"),
    "subcluster_V03072023" = quickTemplate("Seurat Metadata", "Cell Name", paste0(washuMD_file, " and https://mmrfinc.sharepoint.com/:u:/r/sites/IAbioinformatics/Shared%20Documents/General/v1_subclustering/WP_Cluster_Metadata/cluster_metadata_WP_1_20_23.rds?csf=1&web=1&e=8usgKc"), "1/20/23", "Factor of the Subcluster Name for each cell, in format of 'Compartment'.'Parent'.'Subcluster'. Use for most analyses."),
    "subcluster_V03072023_compartment" = quickTemplate("Seurat Metadata", "Cell Name", paste0(washuMD_file, " and https://mmrfinc.sharepoint.com/:u:/r/sites/IAbioinformatics/Shared%20Documents/General/v1_subclustering/WP_Cluster_Metadata/cluster_metadata_WP_1_20_23.rds?csf=1&web=1&e=8usgKc"), "1/20/23", "Factor of the Compartment for each cell. Uses the same naming scheme as 'subcluster_V03072023' for consistency"),
    "subcluster_V03072023_parent" = quickTemplate("Seurat Metadata", "Cell Name", paste0(washuMD_file, " and https://mmrfinc.sharepoint.com/:u:/r/sites/IAbioinformatics/Shared%20Documents/General/v1_subclustering/WP_Cluster_Metadata/cluster_metadata_WP_1_20_23.rds?csf=1&web=1&e=8usgKc"), "1/20/23", "Factor of the original per-compartment cluster for each cell, in the format of 'Compartment'.'Parent'"),
    "cellname" = quickTemplate("Seurat Metadata", "Cell Name", washuMD_file, "N/A", "Cell Identifier. Rownames of seurat metadata object.")
)

info_filt <- info
info_filt[c("aliquot_id", "sample_aliquot_id", "sample_id")] <- NULL
allMD_info <- c(info_singleCell, info_filt)

df <- lapply(allMD_info, data.frame) |> bind_rows()
rownames(df) <- names(allMD_info)
df$variable_name <- rownames(df)
# write.csv(df, file.path(single_cell_output_dir, single_cell_output_INFO))

OLD_MD <- "/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/data/metadata_84_outputs/"
OLD_PER_ALIQUOT <- readRDS(file.path(OLD_MD, "717_per_aliquot_md_v8423.rds"))
OLD_PER_CELL <- readRDS(file.path(OLD_MD, "717_per_cell_md_v8423.rds"))

aliquot_MISMATCH <- list()
aliquot_EXTRA <- list()
idx_mismatch <- 1
idx_extra <- 1

for (feature in colnames(scRNA_MD)) {
    if (feature %in% colnames(OLD_PER_CELL)) {
        test_feat_new <- na.omit(scRNA_MD[, feature])
        test_feat_old <- na.omit(OLD_PER_CELL[, feature])
        if (!all(test_feat_new == test_feat_old)) {
            aliquot_MISMATCH[[idx_mismatch]] <- feature
            idx_mismatch <- idx_mismatch + 1
        }
        # stopifnot(all(test_feat_new == test_feat_old))
    } else {
        aliquot_EXTRA[[idx_extra]] <- feature
        idx_extra <- idx_extra + 1
    }
}

percell_MISMATCH <- list()
percell_EXTRA <- list()
idx_mismatch <- 1
idx_extra <- 1

feature <- "maf_call"

for (feature in colnames(md)) {
    if (feature %in% colnames(OLD_PER_ALIQUOT)) {
        test_feat_new <- na.omit(md[, feature])
        test_feat_old <- na.omit(OLD_PER_ALIQUOT[, feature])
        if (!all(test_feat_new == test_feat_old)) {
            percell_MISMATCH[[idx_mismatch]] <- feature
            idx_mismatch <- idx_mismatch + 1
        }
    } else {
        percell_EXTRA[[idx_extra]] <- feature
        idx_extra <- idx_extra + 1
    }
}
