require(Seurat)
require(ggplot2)
require(SeuratWrappers)
require(RColorBrewer)
require(ggprism)
require(dplyr)
require(limma)
require(edgeR)

# DEG and Visualization Functions
source("pdf_and_png.R")
source("__DotPlot_custom.R")

# Code for generaing dot plots, violin plots, and feature plots
# Used primarily for Figure 2.
plotMarkers <- function(
    seurat_object, markerList, plotLists = list(
        "Vln", "Feat",
        "Dot", "Tile"
    ), group.by = "seurat_clusters", vlnWidth = 15, vlnHeight = 10, featWidth = 20,
    featHeight = 20, dotBaseWidth = 4, dotBaseHeight = 2, dotScaleWidth = 0.35, dotScaleHeight = 0.35, facetScaleWidth = NA, facetScaleHeight = NA, idents.grouping = NA, dotScale = T, nameMod = "", output = "output",
    umapSlot = NULL, DOT_SIZE_X = 1, DOT_SIZE_Y = 1, SAVE_DOT_RDS = F, CLUSTER_COLORS = NULL, DOT_OVERRIDE_H = NULL, DOT_OVERRIDE_W = NULL) {
    if (is.na(facetScaleHeight)) {
        facetScaleHeight <- dotScaleHeight
    }
    if (is.na(facetScaleWidth)) {
        facetScaleWidth <- dotScaleWidth
    }

    outputDir <- file.path(output, nameMod)
    dir.create(outputDir, recursive = T)

    plotOuts <- list()

    if (is.null(group.by)) {
        # Best way to deal with this??
        group.by <- "TMP_IDENT_FIELD"
        seurat_object$TMP_IDENT_FIELD <- Idents(seurat_object)
    }

    if ("Vln" %in% plotLists) {
        pa <- VlnPlot(seurat_object, pt.size = 0, features = markerList, group.by = group.by)

        numIdents <- length(unique(seurat_object@meta.data[[group.by]]))
        vlnWidth <- numIdents * 10 / 18 * min(length(markerList), 3)
        vlnHeight <- ceiling(length(markerList) / 3) * 3
        plotOuts[["Vln"]] <- pdf_and_png(pa, outputDir, paste0(nameMod, "_VIOLIN"),
            pdfWidth = vlnWidth,
            pdfHeight = vlnHeight
        )
        # pdf(file.path(outputDir, paste0(nameMod, '_VIOLIN.pdf')), width =
        # vlnWidth, height = vlnHeight) print(pa) dev.off()
    }

    if ("Feat" %in% plotLists) {
        if (is.null(umapSlot)) {
            pb <- FeaturePlot(seurat_object, features = markerList, reduction = umapSlot)
            plotOuts[["Feat"]] <- pdf_and_png(pb, outputDir, paste0(nameMod, "_FEATUREPLOT"),
                pdfWidth = featWidth,
                pdfHeight = featHeight
            )
            # pdf(file.path(outputDir, paste0(nameMod, '_FEATUREPLOT.pdf')),
            # width = featWidth, height = featHeight) print(pb) dev.off()
        } else {
            plotOuts[["Feat"]] <- list()
            for (red in umapSlot) {
                pb <- FeaturePlot(seurat_object, features = markerList, reduction = red)
                plotOuts[["Feat"]][[red]] <- pdf_and_png(pb, outputDir, paste0(nameMod, "_FEATUREPLOT_", red),
                    pdfWidth = featWidth, pdfHeight = featHeight
                )
                # pdf(file.path(outputDir, paste0(nameMod, '_FEATUREPLOT_',
                # red, '.pdf')), width = featWidth, height = featHeight)
                # print(pb) dev.off()
            }
        }
    }
    # scale_colour_gradientn(colors = c(
    #             "cyan",
    #             "gray", "yellow", "red"
    #         ))

    # Color scheme from pan-cancer tcell atlas
    if ("Dot" %in% plotLists || "Tile" %in% plotLists) {
        facets_vert <- length(names(idents.grouping))
        facets_horiz <- length(names(markerList))
        if (facets_vert > 0) {
            facets_vert <- facets_vert - 1
        }
        if (facets_horiz > 0) {
            facets_horiz <- facets_horiz - 1
        }
        Idents(seurat_object) <- seurat_object@meta.data[[group.by]] # bad color behavior with group.by set
        print(levels(seurat_object))
        if (is.list(idents.grouping) || any(!is.na(names(idents.grouping))) || !is.null(CLUSTER_COLORS)) {
            pc <- customDotPlot(seurat_object, features = markerList, col.min = -2, col.max = 2, idents.grouping = idents.grouping, scale = dotScale, USE_PATCHWORK_INSTEAD_OF_FACETS = T, CLUSTER_COLORS = CLUSTER_COLORS, DO_TILE = F)
        } else {
            pc <- DotPlot(seurat_object, features = markerList, col.min = -2, col.max = 2, scale = dotScale)
        }
        if (is.null(CLUSTER_COLORS)) {
            if ("Tile" %in% plotLists) {
                pc_t <- pc + geom_tile(aes(fill = avg.exp.scaled), colour = "grey") +
                    # geom_point(aes(size = pct.exp), shape = 21, colour = "grey", stroke = 0.5) +
                    scale_fill_gradient2(
                        low =
                            "#5CACDB",
                        mid = "white",
                        high = "#EA7FA3"
                    ) + xlab("") + ylab("") + theme(axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5, hjust = 1, size = rel(DOT_SIZE_X),
                    ), axis.text.y = element_text(size = rel(DOT_SIZE_Y))) + guides(size = guide_legend(override.aes = list(shape = 21, colour = "grey", fill = "white")))
            }
            pc <- pc +
                geom_point(aes(size = pct.exp), shape = 21, colour = "grey", stroke = 0.5) +
                theme(panel.grid = element_line(
                    color = "grey",
                    size = 0.25,
                    linetype = "solid"
                )) +
                scale_colour_gradient2(
                    low =
                        "#5CACDB",
                    mid = "white",
                    high = "#EA7FA3"
                ) + xlab("") + ylab("") +
                theme(
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5, hjust = 1, size = rel(DOT_SIZE_X)
                    ),
                    axis.text.y = element_text(size = rel(DOT_SIZE_Y)),
                    panel.border = element_rect(color = "grey", fill = NA, size = 0.25),
                    strip.clip = "off",
                    panel.spacing.x = unit(facetScaleWidth, "in"),
                    panel.spacing.y = unit(facetScaleHeight, "in")
                ) +
                guides(size = guide_legend(override.aes = list(shape = 21, colour = "grey", fill = "white")))

            pc_order <- DotPlot(seurat_object, features = markerList, cluster.idents = T) +
                scale_colour_gradientn(colors = c(
                    "#60BDE6",
                    "white", "#EA7FA3"
                )) +
                xlab("") + ylab("") + theme(axis.text.x = element_text(
                    angle = 90,
                    vjust = 0.5, hjust = 1
                ))
        } else {
            if ("Tile" %in% plotLists) {
                pc_t <- customDotPlot(seurat_object, features = markerList, col.min = -2, col.max = 2, idents.grouping = idents.grouping, scale = dotScale, USE_PATCHWORK_INSTEAD_OF_FACETS = T, CLUSTER_COLORS = CLUSTER_COLORS, DO_TILE = T)
            }
        }
        numIdents <- length(unique(seurat_object@meta.data[[group.by]]))
        dotWidthAuto <- dotBaseWidth + dotScaleWidth * length(unlist(markerList)) + facetScaleWidth * facets_horiz # was 0.75
        dotHeightAuto <- dotBaseHeight + dotScaleHeight * numIdents + facetScaleHeight * facets_vert
        print(dotWidthAuto)
        print(dotHeightAuto)
        if (!is.null(DOT_OVERRIDE_H)) {
            dotHeightAuto <- DOT_OVERRIDE_H
        }
        if (!is.null(DOT_OVERRIDE_W)) {
            dotWidthAuto <- DOT_OVERRIDE_W
        }
        plotOuts[["Dot"]] <- pdf_and_png(pc, outputDir, paste0(nameMod, "_DOTPLOT_LEGEND"),
            pdfWidth = dotWidthAuto,
            pdfHeight = dotHeightAuto,
            SAVE_RDS = SAVE_DOT_RDS
        )

        plotOuts[["Dot"]][["NoLegend"]] <- pdf_and_png(pc & theme(legend.position = "none"), outputDir, paste0(nameMod, "_DOTPLOT"),
            pdfWidth = dotWidthAuto,
            pdfHeight = dotHeightAuto,
            SAVE_RDS = SAVE_DOT_RDS
        )

        plotOuts[["Dot"]][["botlegend"]] <- pdf_and_png(pc & theme(legend.position = "bottom"), outputDir, paste0(nameMod, "_DOTPLOT_BOTLEGEND"),
            pdfWidth = dotWidthAuto,
            pdfHeight = dotHeightAuto,
            SAVE_RDS = SAVE_DOT_RDS
        )
        plotOuts[["Dot"]][["Plot"]] <- pc

        if ("Tile" %in% plotLists) {
            plotOuts[["Tile"]] <- pdf_and_png(pc_t, outputDir, paste0(nameMod, "_TILE"),
                pdfWidth = dotWidthAuto,
                pdfHeight = dotHeightAuto,
                SAVE_RDS = SAVE_DOT_RDS
            )
        }
        # pdf(file.path(outputDir, paste0(nameMod, '_DOTPLOT.pdf')), width =
        # dotWidthAuto, height = dotHeightAuto) print(pc) dev.off()
        # plotOuts[["Dot_Ordered"]] <- pdf_and_png(pc_order, outputDir, paste0(nameMod, "_DOTPLOT_ORDERED"),
        #     pdfWidth = dotWidthAuto,
        #     pdfHeight = dotHeightAuto,
        #     SAVE_RDS = SAVE_DOT_RDS
        # )
        # pdf(file.path(outputDir, paste0(nameMod, '_DOTPLOT_ORDERED.pdf')),
        # width = 2 + 0.75 * length(markerList), height = 2 + 0.4 *
        # numIdents) print(pc_order) dev.off()
    }
    return(plotOuts)
}

# Modifies the output of presto so that it matches Seurat's Find Markers.
convert_presto_formatting <- function(prestoOutput) {
    df <- data.frame(p_val = prestoOutput[, "pval"], avg_log2FC = prestoOutput[
        ,
        "logFC"
    ], pct.1 = prestoOutput[, "pct_in"] / 100, pct.2 = prestoOutput[
        ,
        "pct_out"
    ] / 100, p_val_adj = prestoOutput[, "padj"], cluster = prestoOutput[
        ,
        "group"
    ], gene = prestoOutput[, "feature"], prestoStatistic = prestoOutput[
        ,
        "statistic"
    ], prestoAUC = prestoOutput[, "auc"])
    rownames(df) <- make.names(df$gene, unique = T)
    df <- dplyr::arrange(df, cluster, p_val_adj, desc(avg_log2FC))
    return(df)
}

# Filter the output of presto
filterPresto <- function(mkrtest, only.pos = T, padj_max = 0.01, logfc.thresh = 0.25) {
    if (only.pos == T) {
        mkrtest <- mkrtest %>%
            dplyr::filter(.data$padj <= padj_max & .data$logFC >= logfc.thresh)
    } else {
        mkrtest <- mkrtest %>%
            dplyr::filter(.data$padj <= padj_max & .data$logFC >= logfc.thresh ||
                .data$logFC <= logfc.thresh)
    }
}

# Find Cluster Markers using either presto directly (presto::wilcoxauc) or the Seurat-Wrappers implementation of presto
FindAllMarkers_presto <- function(
    seurat_object, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25,
    padj_threshold = 0.01, USE_ORIGINAL_PRESTO = F) {
    if (USE_ORIGINAL_PRESTO) {
        print("CALCULATING MARKERS WITH PRESTO")
        # Seurat report the log of the arithmetic mean of the counts, presto
        # reports the geometric mean of the log-normed data. Discrepency in
        # log-fc
        mkrs <- presto::wilcoxauc(seurat_object, pct_in_min = min.pct, padj_max = padj_threshold)
        mkrs <- filterPresto(mkrs,
            only.pos = only.pos, logfc.thresh = logfc.threshold,
            padj_max = 1
        )
        mkrs <- convert_presto_formatting(mkrs)
    } else {
        print("CALCULATING MARKERS WITH SEURAT-WRAPPERS PRESTO")
        mkrs <- RunPrestoAll(seurat_object,
            only.pos = only.pos, min.pct = min.pct,
            logfc.threshold = logfc.threshold
        )
    }
    return(mkrs)
}

# Find the top markers with respect to a provided covariate using FindAllMarkers, and generate a heatmap and dot plot. If markers alerady exist in output directy, re-use those markers.
markers_and_heatmap <- function(
    seurat_object, logfc.threshold = 0.25, min.pct = 0.1,
    heatmapWidth = 25, heatmapHeight = 20, n_markers = 5, nameMod = "", output = "output",
    findMarkersVariable = NULL, displayVariable = NULL, USEPRESTO = T, USE_ORIGINAL_PRESTO = F,
    FORCE_OVERWRITE = F, SCALE_FEATURES_IN_HEATMAP = T, geneText = 18, labelText = 24) {
    outputDir <- file.path(output, nameMod)
    dir.create(outputDir, recursive = T)
    markersName <- file.path(outputDir, paste0(nameMod, "_markers.csv"))
    print(markersName)
    print(file.exists(markersName))
    mapal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)
    if (!is.null(findMarkersVariable)) {
        Idents(seurat_object) <- seurat_object@meta.data[[findMarkersVariable]]
    }
    if (!FORCE_OVERWRITE && file.exists(markersName)) {
        print("RELOADING MARKERS")
        mkrs <- read.csv(markersName, row.names = 1)
        print(head(mkrs))
    } else {
        if (USEPRESTO) {
            mkrs <- tryCatch(
                FindAllMarkers_presto(seurat_object,
                    only.pos = T,
                    min.pct = min.pct, logfc.threshold = logfc.threshold, USE_ORIGINAL_PRESTO = USE_ORIGINAL_PRESTO
                ),
                error = function(cond) {
                    message("FIND MARKERS PRESTO ERRROR:")
                    message(cond)
                    message("")
                    return(NULL)
                }
            )
        } else {
            mkrs <- tryCatch(FindAllMarkers(seurat_object,
                only.pos = T, min.pct = min.pct,
                logfc.threshold = logfc.threshold
            ), error = function(cond) {
                message("FIND MARKERS ERRROR:")
                message(cond)
                message("")
                return(NULL)
            })
        }
        write.csv(mkrs, file.path(outputDir, paste0(nameMod, "_markers.csv")))
    }

    # FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = min.pct,
    # logfc.threshold = logfc.threshold)

    if (!is.null(mkrs)) {
        tryCatch(
            {
                topn <- mkrs %>%
                    group_by(cluster) %>%
                    top_n(n = n_markers, wt = avg_log2FC)
                numUniqueGenes <- length(unique(topn$gene))
                heatmapHeight <- 5 + numUniqueGenes / 4
                if (!is.null(displayVariable)) {
                    Idents(seurat_object) <- seurat_object@meta.data[[displayVariable]]
                }
                if (SCALE_FEATURES_IN_HEATMAP & DefaultAssay(seurat_object) == "RNA") {
                    seurat_object <- ScaleData(seurat_object, features = topn$gene)
                }
                p <- DoHeatmap(seurat_object, features = topn$gene, angle = 0, hjust = 0.5) +
                    theme(axis.text.y = element_text(size = geneText), text = element_text(size = labelText)) +
                    scale_fill_gradientn(colours = rev(mapal))

                #            browser()
                pdf_and_png(p, outputDir, paste0(nameMod, "_MKR_HEATMAP"),
                    pdfWidth = heatmapWidth,
                    pdfHeight = heatmapHeight
                )
                # pdf(file.path(outputDir, paste0(nameMod, '_MKR_HEATMAP.pdf')),
                # width = heatmapWidth, height = heatmapHeight) print(p)
                # dev.off() png(file.path(outputDir, paste0(nameMod,
                # '_MKR_HEATMAP.png')), width = 72*heatmapWidth, height =
                # 72*heatmapHeight) print(p) dev.off()
                plotMarkers(seurat_object,
                    markerList = unique(topn$gene), nameMod = paste0(
                        nameMod,
                        "_HEATMAP_DERIVED"
                    ), output = outputDir, group.by = displayVariable,
                    plotList = c("Dot")
                )
            },
            error = function(cond) {
                message("FIND MARKERS ERRROR - DISPLAY:")
                message(cond)
                message("")
                return(NULL)
            }
        )
    }

    return(list("mkrs" = mkrs, "heatmap" = p))
}

# Render a heatmap from a set of markers.
heatmap_from_markers <- function(
    seurat_object, markerFile, displayVariable, output,
    name, n_markers = 5, heatmapWidth = NULL, heatmapHeight = 20) {
    mkrs <- read.csv(markerFile, row.names = 1)
    mapal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)
    topn <- mkrs %>%
        group_by(cluster) %>%
        top_n(n = n_markers, wt = avg_log2FC)
    numUniqueGenes <- length(unique(topn$gene))
    if (is.null(heatmapWidth)) {
        heatmapHeight <- 5 + numUniqueGenes / 2
    }
    if (!is.null(displayVariable)) {
        Idents(seurat_object) <- seurat_object@meta.data[[displayVariable]]
    }
    p <- DoHeatmap(seurat_object, features = topn$gene, angle = 0, hjust = 0.5) +
        theme(axis.text.y = element_text(size = 12), text = element_text(size = 12)) +
        scale_fill_gradientn(colours = rev(mapal))
    pdf_and_png(p, output, paste0(name, "_MKR_HEATMAP"),
        pdfWidth = heatmapWidth,
        pdfHeight = heatmapHeight
    )

    plotMarkers(seurat_object,
        markerList = unique(topn$gene), nameMod = paste0(
            name,
            "_HEATMAP_DERIVED_SUBOBJECT"
        ), output = output, group.by = displayVariable,
        plotList = c("Dot")
    )
}

library(EnhancedVolcano)

# Find the top markers, and generate a volcano plot
# (Not Used)
markers_and_volcano <- function(
    seurat_object, logfc.threshold = 0.5, pCutoff = 10e-50,
    volcanoWidth = 9, volcanoHeight = 8, nameMod = "", output = "output",
    findMarkersVariable = NULL, positive_ident, negative_ident = NULL, USEPRESTO = T, USE_ORIGINAL_PRESTO = F,
    FORCE_OVERWRITE = F) {
    outputDir <- file.path(output, nameMod)
    dir.create(outputDir, recursive = T)
    markersName <- file.path(outputDir, paste0(nameMod, "_markers_unfiltered.csv"))
    print(markersName)
    print(file.exists(markersName))
    mapal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)
    if (!is.null(findMarkersVariable)) {
        Idents(seurat_object) <- seurat_object@meta.data[[findMarkersVariable]]
    }
    if (!FORCE_OVERWRITE && file.exists(markersName)) {
        print("RELOADING MARKERS")
        mkrs <- read.csv(markersName, row.names = 1)
        print(head(mkrs))
    } else {
        # if (USEPRESTO) {
        #     mkrs <- tryCatch(
        #         FindAllMarkers_presto(seurat_object,
        #             only.pos = T,
        #             min.pct = min.pct, logfc.threshold = logfc.threshold, USE_ORIGINAL_PRESTO = USE_ORIGINAL_PRESTO
        #         ),
        #         error = function(cond) {
        #             message("FIND MARKERS PRESTO ERRROR:")
        #             message(cond)
        #             message("")
        #             return(NULL)
        #         }
        #     )
        # } else {
        mkrs <- tryCatch(RunPresto(seurat_object,
            only.pos = F, ident.1 = positive_ident, ident.2 = negative_ident, min.pct = 0.05,
            logfc.threshold = 0.0
        ), error = function(cond) {
            message("FIND MARKERS ERRROR:")
            message(cond)
            message("")
            return(NULL)
        })
        # }
        write.csv(mkrs, file.path(outputDir, paste0(nameMod, "_markers_unfiltered.csv")))
    }

    # FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = min.pct,
    # logfc.threshold = logfc.threshold)

    if (!is.null(mkrs)) {
        tryCatch(
            {
                topn <- mkrs %>%
                    top_n(n = 10, wt = avg_log2FC)
                bottomn <- mkrs %>% top_n(n = 10, wt = -1 * avg_log2FC)

                topn_label <- unique(c(rownames(topn), rownames(bottomn)))
                positive_label <- positive_ident
                negative_label <- negative_ident
                if (is.null(negative_ident)) {
                    negative_label <- "all others"
                }
                p <- mkrs |> EnhancedVolcano(
                    lab = rownames(mkrs),
                    x = "avg_log2FC",
                    y = "p_val_adj",
                    title = paste0(positive_label, "v", negative_label),
                    selectLab = topn_label,
                    # xlim = c(-1.5, 1.5),
                    pCutoff = pCutoff,
                    FCcutoff = logfc.threshold,
                    pointSize = 3.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.75,
                    max.overlaps = 20,
                    lengthConnectors = unit(0.007, "npc"),
                    labSize = 5.0
                )
                #            browser()
                pdf_and_png(p, outputDir, paste0(nameMod, "_VOLCANO"),
                    pdfWidth = volcanoWidth,
                    pdfHeight = volcanoHeight
                )
                # pdf(file.path(outputDir, paste0(nameMod, '_MKR_HEATMAP.pdf')),
                # width = heatmapWidth, height = heatmapHeight) print(p)
                # dev.off() png(file.path(outputDir, paste0(nameMod,
                # '_MKR_HEATMAP.png')), width = 72*heatmapWidth, height =
                # 72*heatmapHeight) print(p) dev.off()
            },
            error = function(cond) {
                message("FIND MARKERS ERRROR - DISPLAY:")
                message(cond)
                message("")
                return(NULL)
            }
        )
    }

    return(list("mkrs" = mkrs, "heatmap" = p))
}

# Not Used
make_all_contrasts <- function(group, delim = "_vs_") {
    suppressMessages(require(limma))

    # / ensure that group levels are unique
    group <- sort(unique(as.character(group)))

    # / make all combinations
    cb <- combn(group, 2, FUN = function(x) {
        paste0(x[1], "-", x[2])
    })

    # / make contrasts
    contrasts <- limma::makeContrasts(contrasts = cb, levels = group)
    colnames(contrasts) <- gsub("-", delim, colnames(contrasts))

    return(contrasts)
}

# (Not used - see relevant differential expression scripts) - Perform differential expression using limma.
doLimma <- function(seurat_object = NULL, variable, groupA, groupB, output = ".", nameMod = "limma", covariates = "siteXbatch", reload_existing_markers = T, expr = NULL, md = NULL) {
    lmFit_str <- sprintf("lmFit_%s_%s", variable, covariates)
    lmContr_str <- sprintf("Contr_%s_%s-%s", lmFit_str, groupA, groupB)

    subDir <- file.path(output, nameMod)
    dir.create(subDir, showWarnings = F, recursive = T)

    fit_all_path <- file.path(output, nameMod, paste0(lmFit_str, ".rds"))
    fit_final_path <- file.path(output, nameMod, paste0(lmContr_str, ".rds"))
    fit_csv_path <- file.path(output, nameMod, paste0(lmContr_str, ".csv"))

    if (!is.null(seurat_object)) {
        md <- seurat_object@meta.data
    }

    if (file.exists(fit_final_path) && reload_existing_markers) {
        fit_ebayes <- readRDS(fit_final_path)
    } else {
        if (file.exists(fit_all_path) && reload_existing_markers) {
            fit_all <- readRDS(fit_all_path)
            rm(seurat_object)
            gc()
        } else {
            print("Expr Matrix")
            if (!is.null(seurat_object)) {
                expr <- as.matrix(GetAssayData(seurat_object, slot = "data", assay = "RNA"))
                rm(seurat_object)
                gc()
            }
            print("Model Matrix")
            mm <- model.matrix(as.formula(sprintf("~ 0 + %s + %s", variable, covariates)), data = md)

            print("fit")
            fit_all <- lmFit(expr, mm)
            rm(expr)
            gc()
            saveRDS(fit_all, fit_all_path)
        }

        contrastA <- paste0(variable, groupA)
        contrastB <- paste0(variable, groupB)
        contrast_mm <- makeContrasts(contrasts = sprintf("%s - %s", contrastA, contrastB), levels = colnames(mm))

        fit_contr <- contrasts.fit(fit_all, contrast_mm)

        fit_ebayes <- eBayes(fit_contr, trend = T, robust = T)

        saveRDS(fit_ebayes, fit_final_path)
    }

    fit_tt <- topTable(fit_ebayes, coef = 1, adjust = "BH", number = Inf, sort.by = "t", resort.by = "logFC")

    write.csv(fit_tt, fit_csv_path)

    return(list("model" = fit_ebayes, "gene" = fit_tt))
}



# Not Used - Render a volcano plot from limma
renderVolcano <- function(limma_model, groupA, groupB, xlim) {
    forVolcano <- topTable(limma_model, coef = 1, adjust = "BH", n = Inf, sort.by = "t", resort.by = "logFC")
    forVolcano$gene <- rownames(forVolcano)
    forVolcano <- dplyr::filter(forVolcano, !grepl("MT-", gene))

    nGenes <- length(rownames(forVolcano))
    genelist <- rownames(forVolcano)[c(1:15, (nGenes - 14):nGenes)]

    volcano <- forVolcano |> EnhancedVolcano(
        lab = rownames(forVolcano),
        x = "logFC",
        y = "adj.P.Val",
        title = sprintf("%s - %s", groupA, groupB),
        xlim = c(-xlim, xlim),
        selectLab = genelist,
        pCutoff = 10e-50,
        FCcutoff = 0.12,
        pointSize = 3.0,
        drawConnectors = TRUE,
        widthConnectors = 0.75,
        max.overlaps = 20,
        lengthConnectors = unit(0.007, "npc"),
        labSize = 4.0
    )
    return(volcano)
}
