source('/opt/megaseq-data/sarthak/projects/MMRF/scripts/20240214_Survival_load.R')
CUTOFF = 90

## MYRED: 4 compartments - Myeloid
setwd('/opt/megaseq-data/sarthak/projects/MMRF/')
merged_MY = subset(merged, subset = subcluster_V03072023_compartment %in% c('Myeloid'))
merged_MY = subset(merged_MY, subset = doublet_pred %in% c('singlet','poss_singlet'))
MY_compartments_interest = c('Myeloid.5', 'Myeloid.2',  'Myeloid.8', 'Myeloid.12') #  (M2 Macro)
merged_MY = subset(merged_MY, subset = subcluster_V03072023 %in% c(MY_compartments_interest))
merged_MY

# Important files
sig_name = file.path('outs',CUTOFF,'MYRED_auc.csv')

## Module score approach
color_low <- "#2c90c3"
color_high <- "#5a221b"

## AUC value approach
sig <- read.csv(sig_name, check.names = FALSE)
rownames(sig) <- sig$Cell
colnames(sig) <- make.names(colnames(sig))
merged_MY <- AddMetaData(merged_MY, sig)

# > colnames(sig)
#  [1] "Cell"     "ATF3..."  "BRCA1..." "E2F1..."  "E2F7..."  "E2F8..." 
#  [7] "ETS1..."  "FOS..."   "FOSL2..." "IKZF2..." "IRF1..."  "IRF2..." 
# [13] "IRF7..."  "IRF9..."  "KLF9..."  "STAT1..."

subset_patMD <- merged_MY@meta.data |>
    dplyr::filter(d_lot_1_start_day < 2) |>
    dplyr::distinct(public_id, progression_group, censpfs, pfscdy, d_pt_sex, d_dx_amm_age, d_dx_amm_iss_stage, siteXbatch, d_tx_induction_cat, d_amm_tx_asct_1st)

merged_MYaverageMod <- merged_MY@meta.data %>%
    dplyr::group_by(public_id
    ) %>%
    dplyr::summarise(
        BRCA1 = mean(!!sym(paste0(colnames(sig)[3]))),
        E2F1 = mean(!!sym(paste0(colnames(sig)[4]))),
        E2F8 =  mean(!!sym(paste0(colnames(sig)[6]))),
        IRF2 =  mean(!!sym(paste0(colnames(sig)[12]))),
        IRF7 = mean(!!sym(paste0(colnames(sig)[13]))),
        IRF9 =  mean(!!sym(paste0(colnames(sig)[14]))),
        STAT1 =  mean(!!sym(paste0(colnames(sig)[16]))),
        IRF1 =  mean(!!sym(paste0(colnames(sig)[11]))),
        KLF9 =  mean(!!sym(paste0(colnames(sig)[15])))

    )

patMD_averageMod <- dplyr::left_join(subset_patMD, merged_MYaverageMod, by = "public_id")

featList <- list(
    "BRCA1",
    "E2F1",
    "E2F8",
    "IRF2",
    "IRF7",
    "IRF9",
    "STAT1",
    "IRF1",
    "KLF9"
)

#library(ggstatsplot)
library(ggpubr)
library(forestmodel)
library(gtsummary)
library(survival)
library(ggsurvfit)
library(contsurvplot)
library(ggprism)
library('pammtools')
output_surv <- file.path(output_F3, "Myeloids_AUC_4_Fig5E",CUTOFF)
dir.create(output_surv, recursive = T, showWarnings = F)

surv_stats <- data.frame()
source('/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/source/Helper_Functions/pdf_and_png.R')
for (feat in c(featList)) {
    tmp <- patMD_averageMod

    tmp$survVar <- tmp[, feat]
    c1 <- survfit2(Surv(pfscdy, censpfs) ~ survVar, data = tmp)
    c1 <- survMisc::cutp(c1)
    med_val <- c1[[1]][[1]][1]

    tmp <- tmp %>% mutate(median_cutoff = case_when(
        .data[[feat]] < med_val ~ "Low",
        .data[[feat]] >= med_val ~ "High"
    ))

    tmp$time <- tmp$pfscdy
    p_naive <- survfit2(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp) %>% ggsurvfit(size = 2) +
        labs(
            x = "Days",
            y = "Progression Free Survival"
        ) +
        add_censor_mark(size = 3) +
        # add_risktable() +
        theme_prism() +
        ggtitle(paste0(feat)) +
        scale_color_manual(values = c("Low" = color_low, "High" = color_high))
    
    print(paste0(feat, ':LogP:', survdiff(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)$p))
    LRP = survdiff(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)$p
    tmp$median_cutoff <- factor(tmp$median_cutoff, levels = c('Low', 'High'))
    cox_naive <- coxph(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)
    HR = summary(cox_naive)$coef[2]
    HRP = summary(cox_naive)$coef[5]
    surv_stats <- rbind(surv_stats, cbind(feat, LRP, HR, HRP))
    cox_cont_sxb <- coxph(Surv(pfscdy, censpfs) ~ survVar + siteXbatch, data = tmp, x = TRUE)

    p_area <- plot_surv_area(
        time = "pfscdy",
        status = "censpfs",
        variable = "survVar",
        data = tmp,
        model = cox_cont_sxb,
        start_color = color_low,
        end_color = color_high
    ) + ggtitle(paste0(feat))



    pval <- summary(cox_cont_sxb)$coefficients["survVar", "Pr(>|z|)"]
    HR <- summary(cox_cont_sxb)$coefficients["survVar", "exp(coef)"]

    pval <- round(pval, digits = 4)

    annotations <- data.frame(
        xpos = c(-Inf),
        ypos = c(-Inf),
        annotateText = c(
            paste0("HR: ", round(HR, digits = 3), "\n", "CoxPH P.Value: ", scales::pvalue(pval))
        ),
        hjustvar = c(0),
        vjustvar = c(0)
    )
    p_naive_cox <- p_naive #+ geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText, fontface = "bold"), color = c("black"), size = 6)
    p_naive_cox <- ggsurvfit_build(p_naive_cox)

    pdf_and_png(p_naive_cox, output_surv, paste0(feat, "_surv"), pdfWidth = W_DBL / (3*1.1), pdfHeight = H_FULL /1.5, scale = 3)

    annotations_b <- data.frame(
        xpos = c(Inf),
        ypos = c(Inf),
        annotateText = c(
            paste0("HR: ", round(HR, digits = 3), "\nCoxPH P.Value: ", scales::pvalue(pval))
        ),
        hjustvar = c(1),
        vjustvar = c(1)
    )

    p_area <- p_area + geom_text(data = annotations_b, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText, fontface = "bold"), color = c("black"), size = 6)

    combined <- p_naive_cox + p_area

    #pdf_and_png(p_area, output_surv, paste0(feat, "_area_surv"), pdfWidth = W_DBL / 3, pdfHeight = H_FULL, scale = 3)

    #pdf_and_png(combined, output_surv, paste0(feat, "_median_and_area"), pdfWidth = 2 * W_DBL / 3, pdfHeight = H_FULL, scale = 3)
}

# sobj$NK_Focus_List1 <- (sobj$NK_Focus_List1  - mean(sobj$NK_Focus_List1 )) / sd(sobj$NK_Focus_List1 )
write.csv(surv_stats, file.path(output_surv,'surv_stats_cutp.csv'))

library(RColorBrewer)
library(ggprism)
merged_MY@meta.data = merged_MY@meta.data |> dplyr::mutate(cellID_short = dplyr::case_when(cellID_short == "CD14+Mono_INF" ~ "CD14+Mono_IFN", .default = cellID_short))

Idents(merged_MY) <- 'cellID_short'


plot_feature <- function(obj, reg, filePath= output_surv){
    TF = gsub("...","",reg, fixed = T)
    title = paste0(TF, "_Regulon")
    filename = file.path(output_surv, paste0("FeaturePlot_SCENIC_",TF,"_Regulon.pdf"))
    ggsave(filename,
        FeaturePlot(obj, 
                    features = reg,
                    min.cutoff = 'q1',
                    max.cutoff = 'q99',
                    pt.size = 2, order = T, label = T, label.size = 6, reduction = "umap.sub", ) +
            labs(x = "UMAP1", y = "UMAP2", title = title) +
            scale_color_gradientn(colors = rev(brewer.pal(n = 10, "Spectral"))) +
            theme_prism() +
            scale_y_continuous(limits = c(-6,5)) + 
            scale_x_continuous(limits = c(-8,12)) + 
            theme(legend.title = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.position = c(0.9, 0.2),
                legend.background = element_rect(fill = "white", linewidth = 0)),
        width = 6,
        height = 6)
}

for(i in 2:length(colnames(sig))){
    reg = colnames(sig)[i]
    plot_feature(obj = merged_MY,reg = reg)
}


# output_surv <- file.path(output_F3, "Myeloids_AUC_4_Fig5E",CUTOFF,"MEDIAN")
# dir.create(output_surv, recursive = T, showWarnings = F)
# surv_stats <- data.frame()
# source('/opt/pmedseq-data/wcpilch/Projects/MMRF-ImmuneAtlas-Annotation/source/Helper_Functions/pdf_and_png.R')
# for (feat in c(featList)) {
#     tmp <- patMD_averageMod

#     tmp$survVar <- tmp[, feat]
#     med_val <- median(tmp$survVar)

#     tmp <- tmp %>% mutate(median_cutoff = case_when(
#         .data[[feat]] < med_val ~ "Low",
#         .data[[feat]] >= med_val ~ "High"
#     ))

#     tmp$time <- tmp$pfscdy
#     p_naive <- survfit2(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp) %>% ggsurvfit(size = 2) +
#         labs(
#             x = "Days",
#             y = "Progression Free Survival"
#         ) +
#         add_censor_mark(size = 3) +
#         # add_risktable() +
#         theme_prism() +
#         ggtitle(paste0(feat)) +
#         scale_color_manual(values = c("Low" = color_low, "High" = color_high))
    
#     print(paste0(feat, ':LogP:', survdiff(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)$p))
#     LRP = survdiff(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)$p
#     tmp$median_cutoff <- factor(tmp$median_cutoff, levels = c('Low', 'High'))
#     cox_naive <- coxph(Surv(pfscdy, censpfs) ~ median_cutoff, data = tmp)
#     HR = summary(cox_naive)$coef[2]
#     HRP = summary(cox_naive)$coef[5]
#     surv_stats <- rbind(surv_stats, cbind(feat, LRP, HR, HRP))
#     cox_cont_sxb <- coxph(Surv(pfscdy, censpfs) ~ survVar + siteXbatch, data = tmp, x = TRUE)

#     p_area <- plot_surv_area(
#         time = "pfscdy",
#         status = "censpfs",
#         variable = "survVar",
#         data = tmp,
#         model = cox_cont_sxb,
#         start_color = color_low,
#         end_color = color_high
#     ) + ggtitle(paste0(feat))



#     pval <- summary(cox_cont_sxb)$coefficients["survVar", "Pr(>|z|)"]
#     HR <- summary(cox_cont_sxb)$coefficients["survVar", "exp(coef)"]

#     pval <- round(pval, digits = 4)

#     annotations <- data.frame(
#         xpos = c(-Inf),
#         ypos = c(-Inf),
#         annotateText = c(
#             paste0("HR: ", round(HR, digits = 3), "\n", "CoxPH P.Value: ", scales::pvalue(pval))
#         ),
#         hjustvar = c(0),
#         vjustvar = c(0)
#     )
#     p_naive_cox <- p_naive #+ geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText, fontface = "bold"), color = c("black"), size = 6)
#     p_naive_cox <- ggsurvfit_build(p_naive_cox)

#     pdf_and_png(p_naive_cox, output_surv, paste0(feat, "_surv"), pdfWidth = W_DBL / (3*1.1), pdfHeight = H_FULL /1.5, scale = 3)

#     annotations_b <- data.frame(
#         xpos = c(Inf),
#         ypos = c(Inf),
#         annotateText = c(
#             paste0("HR: ", round(HR, digits = 3), "\nCoxPH P.Value: ", scales::pvalue(pval))
#         ),
#         hjustvar = c(1),
#         vjustvar = c(1)
#     )

#     p_area <- p_area + geom_text(data = annotations_b, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText, fontface = "bold"), color = c("black"), size = 6)

#     combined <- p_naive_cox + p_area

#     #pdf_and_png(p_area, output_surv, paste0(feat, "_area_surv"), pdfWidth = W_DBL / 3, pdfHeight = H_FULL, scale = 3)

#     #pdf_and_png(combined, output_surv, paste0(feat, "_median_and_area"), pdfWidth = 2 * W_DBL / 3, pdfHeight = H_FULL, scale = 3)
# }
# write.csv(surv_stats, file.path(output_surv,'surv_stats.csv'))
library(ggplot2)
cutp_stats <- read.csv('plots/Survival/FIG3/Myeloids_AUC_4_Fig5E/90/surv_stats_cutp.csv')
auc_Z <- read.csv('outs/90/MYRED_aucZMean_cluster.csv')
auc_Z <- reshape::melt(auc_Z, id=c("regulon"))

auc_Z <- auc_Z |> filter(regulon %in% c("BRCA1(+)","E2F1(+)","E2F8(+)","KLF9(+)","IRF1(+)","IRF2(+)","IRF7(+)","IRF9(+)","STAT1(+)"))
auc_Z$regulon <- factor(auc_Z$regulon, levels = c("BRCA1(+)","E2F1(+)","E2F8(+)","KLF9(+)","IRF1(+)","IRF2(+)","IRF7(+)","IRF9(+)","STAT1(+)"))
auc_Z$variable <- factor(auc_Z$variable, levels = c("Myeloid.2","Myeloid.8", "Myeloid.12", "Myeloid.5" ))
pdf(file.path(output_surv,"geom_tile.pdf"), height = 4, width = 4)
ggplot(auc_Z, aes(x = variable, y = regulon)) + 
  geom_tile(colour="black", size = 0.25, aes(fill = value)) +
  scale_fill_gradientn(colours = c("#603F8B", "#F8F8FF", "#F3E5AB","#FF2E2E"),
                       na.value = "white"#,
                       #values = c(0.75, 1, 1, 1.75)#,
                       # = c(0, 1) +
                        ) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

auc_Z$cluster <- plyr::revalue(auc_Z$variable, c(
    "Myeloid.2" = "CD14+Mono_INF",
    "Myeloid.8" = "CD14+Mono_pro-inflam", 
    "Myeloid.12" = "Macro/Mono", 
    "Myeloid.5" = "GMP"
))

pdf(file.path(output_surv,"geom_tile_cluster.pdf"), height = 4, width = 4)
ggplot(auc_Z, aes(x = cluster, y = regulon)) + 
  geom_tile(colour="black", size = 0.25, aes(fill = value)) +
  scale_fill_gradientn(colours = c("#603F8B", "#F8F8FF", "#F3E5AB","#FF2E2E"),
                       na.value = "white"#,
                       #values = c(0.75, 1, 1, 1.75)#,
                       # = c(0, 1) +
                        ) +
  theme(axis.text.y = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, size = 8))
dev.off()



cutp_stats$X = rep('HR',nrow(cutp_stats))
cutp_stats$Y = rep('p.value',nrow(cutp_stats))

cutp_stats$feat <- factor(cutp_stats$feat, levels = c("BRCA1","E2F1","E2F8","KLF9","IRF1","IRF2","IRF7","IRF9","STAT1"))


pdf(file.path(output_surv,"geom_tile_stats.pdf"), height = 4, width = 2.5)
ggplot(cutp_stats, aes(x = X, y = feat)) + 
  geom_tile(colour="black", size = 0.25, aes( fill =  HR)) +
  scale_fill_gradientn(colours = c("#603F8B", "#F8F8FF", "#F3E5AB","#FF2E2E"),
                       na.value = "white"#, 
                       #values = c(0.75, 0.90, 1.10, 1.75)#,
                       # = c(0, 1) +
                        )
dev.off()

pdf(file.path(output_surv,"geom_tile_stats_pval.pdf"), height = 3, width = 2)
ggplot(cutp_stats, aes(x = Y, y = feat)) + 
  geom_tile(colour="black", size = 0.25, aes( fill =  HRP_cat)) +
  scale_fill_manual(colours = c("#cb416b", "#FFB6C1","#D3D3D3","#FFB6C1"))
dev.off()

pdf(file.path(output_surv,"geom_tile_stats_pval.pdf"), height = 4, width = 2.5)
ggplot(cutp_stats, aes(x = Y, y = feat)) + 
  geom_tile(colour="black", size = 0.25, aes( fill =  HRP_cat)) +
  scale_fill_manual(values = c("#cb416b", "#FFB6C1", "#D3D3D3","#A9A9A9")
                        )
dev.off()



cutp_stats$HRP_cat <- cut(cutp_stats$HRP,
                       breaks=c(0,0.01,0.05,0.1,1),
                       labels=c('< 0.01', '< 0.05', '< 0.10', '> 0.10'))



targets <- read.csv('outs/90/MYRED_regulon_target.csv')
targets <- targets |> filter(Regulon %in% c("BRCA1(+)","E2F1(+)","E2F8(+)","KLF9(+)","IRF1(+)","IRF2(+)","IRF7(+)","IRF9(+)","STAT1(+)"))

targets$Target <- gsub("dict_keys([","",targets$Target, fixed = TRUE)
targets$Target <- gsub("])","",targets$Target, fixed = TRUE)
targets$list = strsplit(targets$Target, ',')
data = data.frame()
for(i in 1:nrow(targets)){
    data = rbind(data, cbind(targets$Regulon[i],paste0(targets$list[[i]], collapse = ' '),length(targets$list[[i]])))
}
data$V2 <- gsub("'  '"," ", data$V2)
data$V2 <- gsub("'","", data$V2)

regulons = targets$list
names(regulons) = targets$Regulon

table_R <- crossprod(table(stack(regulons)))
dat2 <-
 table_R %>%
  as_tibble()
row.names(dat2) <- colnames(dat2)
dat2 <-  dat2 |> 
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")
pdf(file.path(output,"co_occurance.pdf"), height = 5, width = 6)

  ggplot(dat2, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "#f5f5dc", high = "#000080")
dev.off()
