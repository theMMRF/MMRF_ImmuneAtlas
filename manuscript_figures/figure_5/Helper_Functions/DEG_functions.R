require(Seurat)
require(ggplot2)
require(SeuratWrappers)
require(RColorBrewer)
require(ggprism)
require(dplyr)
require(limma)
require(edgeR)
# library(meta)
# DEG and Visualization Functions
source("pdf_and_png.R")
source("__DotPlot_custom.R")

# Takes two lmFit outputs with non-overlapping gene sets and combines them.
merge_limma_lmFit_outputs <- function(lmFitA, lmFitB) {
    lmFitA$coefficients <- rbind(lmFitA$coefficients, lmFitB$coefficients)
    lmFitA$df.residual <- c(lmFitA$df.residual, lmFitB$df.residual)
    lmFitA$sigma <- c(lmFitA$sigma, lmFitB$sigma)
    lmFitA$stdev.unscaled <- rbind(lmFitA$stdev.unscaled, lmFitB$stdev.unscaled)
    lmFitA$Amean <- c(lmFitA$Amean, lmFitB$Amean)
    return(lmFitA)
}


# Computes DEGs using limma-trend on the scRNA-seq matrix. Corrects for a single batch covariate if included (batchVariable)
# If groupA and groupB are specified - only performs the requested contrast. Otherwise, returns all requested pairs of contrasts
# IF the variable is intended to be treated as a continuous value, specify 'treatVariableAsContinuous'
# Output is a named list, with fields
#   $lmFit - all contrast models
#   $fullModel - fit prior to contrasts
#   $topTable - list of topTable outputs, one per contrast. List is named by the contrast (%var%groupA-%var%groupB)
#   $mergedTable - all topTables merged into one large table. contrast column, along with groupA and groupB columns, indicate what comparison the row corresponds to

# Inputs
#   seurat_object - (seurat object)
#   variable - what variable are you interested in comparing. If a factor, contrasts will be computed in order of the levels
#   groupA, groupB - If specified, the specific values in variable to perform a contrast in-between.
#   treatVariableAsContinuous - if variable is a continuous value (e.g. age), fits the model directly with that as a continuous covariate instead of a factor
#   batchVariable - (Currently one value), if Specified, what to adjust for in the model
#   subsetVariableGroupList - a list of values in variable to subset 'variable' to prior to fitting the model. Set to 'AUTO' to automatically subset to cells with values of 'groupA' and 'groupB'.
#   batch_features - Memory management, how many genes should a model be fit to at a time. -1, or a value above the number of genes you have will fit things in one batch. Increase to improve speed, Decrease if memory is a concern
#   RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES - In the output tables, convert logFC, P, and the adjusted P value to the seurat standard naming pattern.
#   DISCARD_CONTRAST_MODELS - If true, deletes the fitted contrast models. These contain residuals, and may take up extra memory. Output will contain NULL
#   DISCARD_FULL_MODELS - If true, deletes the full fitted model when done. Output will contain NULL
#   EXPERIMENTAL_META_ANALYSES - NYI
doLimma <- function(seurat_object, variable = NULL, groupA = NULL, groupB = NULL, treatVariableAsContinuous = F, batchVariable = NULL, CORRECT_WITH_INTERACTION_TERM = FALSE, INTERACTION_TERM_BATCH_LEVELS = NULL, subsetVariableGroupList = NULL, batch_features = 6000, RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES = FALSE, DISCARD_CONTRAST_MODELS = FALSE, DISCARD_FULL_MODELS = FALSE, EXPERIMENTAL_META_ANALYSES = FALSE) {
    if (is.null(variable)) {
        seurat_object$IDENTS <- Idents(seurat_object)
        variable <- "IDENTS"
    }

    if (!treatVariableAsContinuous) {
        Idents(seurat_object) <- seurat_object@meta.data[, variable]
    }

    cellList <- !is.na(seurat_object@meta.data[, variable])
    if (!is.null(batchVariable)) {
        cellList <- cellList & !is.na(seurat_object@meta.data[, batchVariable])
    }

    if (!is.null(subsetVariableGroupList)) {
        if (subsetVariableGroupList == "AUTO") {
            subsetVariableGroupList <- c(groupA, groupB)
        }
        cellList <- cellList & (seurat_object@meta.data[, variable] %in% subsetVariableGroupList)
    }

    seurat_object <- subset(seurat_object, cells = which(cellList))
    if (batch_features == -1) {
        batch_features <- nrow(seurat_object)
    }

    expr <- Seurat::GetAssayData(seurat_object, slot = "data", assay = "RNA")
    md <- seurat_object@meta.data

    model_string <- sprintf("~ 0 + %s", variable)

    if (!is.null(batchVariable)) {
        if (CORRECT_WITH_INTERACTION_TERM) {
            if (is.null(INTERACTION_TERM_BATCH_LEVELS)) {
                if (is.factor(md[[batchVariable]])) {
                    INTERACTION_TERM_BATCH_LEVELS <- levels(md[[batchVariable]])
                } else {
                    INTERACTION_TERM_BATCH_LEVELS <- unique(md[[batchVariable]])
                }
            }
            md[[batchVariable]] <- factor(md[[batchVariable]], levels = INTERACTION_TERM_BATCH_LEVELS)

            model_string <- sprintf("%s * %s", model_string, batchVariable)
        } else {
            model_string <- sprintf("%s + %s", model_string, batchVariable)
        }
    }
    design <- model.matrix(as.formula(model_string), data = md)
    colnames(design) <- make.names(colnames(design))

    lmFit_tmp <- NULL
    for (idx_start in seq(from = 1, to = nrow(expr), by = batch_features)) {
        start <- idx_start
        end <- min(idx_start + batch_features - 1, nrow(expr))
        print(sprintf("Gene IDx %d to %d of %d", start, end, nrow(expr)))

        expr_subset <- as.matrix(expr[start:end, ])
        fit_subset <- limma::lmFit(expr_subset, design)

        if (is.null(lmFit_tmp)) {
            lmFit_tmp <- fit_subset
        } else {
            lmFit_tmp <- merge_limma_lmFit_outputs(lmFit_tmp, fit_subset)
            # lmFit_tmp$coefficients <- rbind(lmFit_tmp$coefficients, fit_subset$coefficients)
            # lmFit_tmp$df.residual <- c(lmFit_tmp$df.residual, fit_subset$df.residual)
            # lmFit_tmp$sigma <- c(lmFit_tmp$sigma, fit_subset$sigma)
            # lmFit_tmp$stdev.unscaled <- rbind(lmFit_tmp$stdev.unscaled, fit_subset$stdev.unscaled)
            # lmFit_tmp$Amean <- c(lmFit_tmp$Amean, fit_subset$Amean)
        }
        rm(fit_subset)
        gc()
    }

    lmfit_contrasts_out <- list()
    contrasts_lists <- c()
    variable_pairs <- list()
    toptable_list <- list()



    if (!treatVariableAsContinuous) {
        if (!is.null(groupA) && !is.null(groupB)) {
            if (CORRECT_WITH_INTERACTION_TERM) {
                # Probably can refactor to reuse this block
                nBatch <- length(INTERACTION_TERM_BATCH_LEVELS)

                termA <- sprintf("%s%s", variable, groupA)
                termB <- sprintf("%s%s", variable, groupB)

                # interactTermA <- ""
                interactTermB <- ""

                for (batchIDX in 2:nBatch) {
                    batch_factor <- INTERACTION_TERM_BATCH_LEVELS[[batchIDX]]
                    if (batchIDX != 2) {
                        # interactTermA <- sprintf("%s+", interactTermA)
                        interactTermB <- sprintf("%s + ", interactTermB)
                    }
                    # interactTermA <- sprintf("%s%s%s.%s%s", interactTermA, variable, groupA, batchVariable, batch_factor)
                    interactTermB <- sprintf("%s%s%s.%s%s", interactTermB, variable, groupB, batchVariable, batch_factor)
                }

                termA_full <- termA # sprintf("%s + (%s)/%d", termA, interactTermA, nBatch)
                termB_full <- sprintf("%s + (%s)/%d", termB, interactTermB, nBatch)

                contrast_str <- c(sprintf("(%s) - (%s)", termA_full, termB_full))
            } else {
                contrast_str <- c(sprintf("%s%s - %s%s", variable, groupA, variable, groupB))
            }
            contrasts_lists <- c(contrast_str)
            variable_pairs[[contrast_str]] <- c(groupA, groupB)
        } else {
            if (is.factor(seurat_object@meta.data[, variable])) {
                all_levels <- levels(droplevels(seurat_object@meta.data[, variable]))
            } else {
                all_levels <- unique(seurat_object@meta.data[, variable])
            }
            if (length(all_levels) < 2) {
                stop(sprintf("Variable passed to perform differential expression with has %d levels or unique values. Need at least 2.", length(all_levels)))
            }
            for (i in 1:(length(all_levels) - 1)) {
                for (j in (i + 1):length(all_levels)) {
                    if (CORRECT_WITH_INTERACTION_TERM) {
                        groupA <- all_levels[[i]]
                        groupB <- all_levels[[j]]

                        nBatch <- length(INTERACTION_TERM_BATCH_LEVELS)


                        termA <- sprintf("%s%s", variable, groupA)
                        termB <- sprintf("%s%s", variable, groupB)

                        interactTermA <- ""
                        interactTermB <- ""

                        for (batchIDX in 2:nBatch) {
                            batch_factor <- INTERACTION_TERM_BATCH_LEVELS[[batchIDX]]
                            if (batchIDX != 2) {
                                interactTermA <- sprintf("%s+", interactTermA)
                                interactTermB <- sprintf("%s+", interactTermB)
                            }
                            interactTermA <- sprintf("%s%s%s.%s%s", interactTermA, variable, groupA, batchVariable, batch_factor)
                            interactTermB <- sprintf("%s%s%s.%s%s", interactTermB, variable, groupB, batchVariable, batch_factor)
                        }


                        if (i == 1) {
                            termA_full <- termA # sprintf("%s + (%s)/%d", termA, interactTermA, nBatch)
                        } else {
                            termA_full <- sprintf("%s + (%s)/%d", termA, interactTermA, nBatch)
                        }
                        termB_full <- sprintf("%s + (%s)/%d", termB, interactTermB, nBatch)

                        contrast_str <- sprintf("%s - %s", termA_full, termB_full)
                    } else {
                        contrast_str <- sprintf("%s%s-%s%s", variable, all_levels[[i]], variable, all_levels[[j]])
                    }
                    contrasts_lists <- c(contrasts_lists, contrast_str)
                    variable_pairs[[contrast_str]] <- c(all_levels[[i]], all_levels[[j]])
                }
            }
        }
        for (x in contrasts_lists) {
            print(sprintf("Computing Contrast %s", x))
            contrast <- limma::makeContrasts(contrasts = x, levels = design) # NP is positive, RP is negative

            lmFit_cont <- limma::contrasts.fit(lmFit_tmp, contrast)
            TMPmd <- limma::eBayes(lmFit_cont, trend = TRUE, robust = TRUE)
            toptable_list[[x]] <- topTable(TMPmd, coef = 1, adjust = "BH", n = Inf, sort.by = "t", resort.by = "t")
            toptable_list[[x]]$log2FC <- log2(exp(toptable_list[[x]]$logFC))
            if (!DISCARD_CONTRAST_MODELS) {
                lmfit_contrasts_out[[x]] <- TMPmd
            } else {
                rm(TMPmd)
            }

            seurat_style_foldchange <- Seurat::FoldChange(seurat_object, ident.1 = variable_pairs[[x]][[1]], ident.2 = variable_pairs[[x]][[2]], group.by = variable)
            colnames(seurat_style_foldchange)[[1]] <- "SEURAT_ESTIMATED_avg_log2FC"

            seurat_style_foldchange$gene <- rownames(seurat_style_foldchange)
            toptable_list[[x]]$gene <- rownames(toptable_list[[x]])
            toptable_list[[x]] <- dplyr::left_join(toptable_list[[x]], seurat_style_foldchange, by = "gene") |> dplyr::relocate(gene, logFC, pct.1, pct.2, t, P.Value, adj.P.Val, SEURAT_ESTIMATED_avg_log2FC)
            toptable_list[[x]]$group1 <- variable_pairs[[x]][[1]]
            toptable_list[[x]]$group2 <- variable_pairs[[x]][[2]]
            toptable_list[[x]]$contrast_string <- x
            toptable_list[[x]]$variable <- variable

            if (RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES) {
                toptable_list[[x]] <- toptable_list[[x]] |> dplyr::rename(avg_log2FC = "log2FC", p_val = "P.Value", p_val_adj = "p_val_adj")
            }
        }
        if (EXPERIMENTAL_META_ANALYSES) {
            # todo
            print("not yet implemented")
        }
    } else {
        x <- variable
        TMPmd <- limma::eBayes(lmFit_cont, trend = TRUE, robust = TRUE)
        toptable_list[[x]] <- topTable(TMPmd, coef = 1, adjust = "BH", n = Inf, sort.by = "t", resort.by = "t")
        toptable_list[[x]]$log2FC <- log2(exp(toptable_list[[x]]$logFC))

        if (!DISCARD_CONTRAST_MODELS) {
            lmfit_contrasts_out[[x]] <- TMPmd
        } else {
            rm(TMPmd)
        }

        toptable_list[[x]]$variable <- variable
        if (RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES) {
            toptable_list[[x]] <- toptable_list[[x]] |> dplyr::rename(avg_log2FC = "log2FC", p_val = "P.Value", p_val_adj = "p_val_adj")
        }
    }
    if (DISCARD_FULL_MODELS) {
        lmFit_tmp <- NULL
    }
    if (DISCARD_CONTRAST_MODELS) {
        lmfit_contrasts_out <- NULL
    }

    SINGLE_TABLE <- NULL
    for (tmp_table in toptable_list) {
        if (is.null(SINGLE_TABLE)) {
            SINGLE_TABLE <- tmp_table
        } else {
            SINGLE_TABLE <- rbind(SINGLE_TABLE, tmp_table)
        }
    }
    return(list("lmFit" = lmfit_contrasts_out, "topTable" = toptable_list, "mergedTable" = SINGLE_TABLE, "fullModel" = lmFit_tmp))
}


# Computes DEGs using limma-trend on the scRNA-seq matrix. Corrects for a single batch covariate if included (batchVariable)
# If groupA and groupB are specified - only performs the requested contrast. Otherwise, returns all requested pairs of contrasts
# IF the variable is intended to be treated as a continuous value, specify 'treatVariableAsContinuous'
# Output is a named list, with fields
#   $lmFit - all contrast models
#   $topTable - list of topTable outputs, one per level List is named by the levels of the varaible.
#   $mergedTable - all topTables merged into one large table. contrast column, along with groupA and groupB columns, indicate what comparison the row corresponds to

# Inputs
#   seurat_object - (seurat object)
#   variable - what variable are you interested in comparing. Will perform all 1 v all comparisons for the variable of interest.
#   VARIABLE_OMIT_NA - Should NA values be treated as 'Other' values in 1 v All comparisons, or should they be removed? Defaults to TRUE (removing NA values in variable)
#   treatVariableAsContinuous - if variable is a continuous value (e.g. age), fits the model directly with that as a continuous covariate instead of a factor
#   batchVariable - (Currently one value), if Specified, what to adjust for in the model
#   subsetVariableGroupList - a list of values in variable to subset 'variable' to prior to fitting the model. Set to 'AUTO' to automatically subset to cells with values of 'groupA' and 'groupB'.
#   batch_features - Memory management, how many genes should a model be fit to at a time. -1, or a value above the number of genes you have will fit things in one batch. Increase to improve speed, Decrease if memory is a concern
#   RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES - In the output tables, convert logFC, P, and the adjusted P value to the seurat standard naming pattern.
#   DISCARD_MODELS - Tells doLimma to discard any models to save memory.

doLimmaAll <- function(seurat_object, variable = NULL, VARIABLE_OMIT_NA = TRUE, batchVariable = "siteXbatch", subsetVariableGroupList = NULL, batch_features = 6000, RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES = FALSE, DISCARD_MODELS = TRUE) {
    if (is.null(variable)) {
        seurat_object@meta.data$IDENTS <- Idents(seurat_object)
        variable <- "IDENTS"
    }

    if (VARIABLE_OMIT_NA) {
        cellList <- !is.na(seurat_object@meta.data[, variable])
    } else {
        cellList <- rep(TRUE, n = ncol(seurat_object))
    }
    if (!is.null(batchVariable)) {
        cellList <- cellList & !is.na(seurat_object@meta.data[, batchVariable])
    }


    if (!is.null(subsetVariableGroupList)) {
        if (subsetVariableGroupList == "AUTO") {
            subsetVariableGroupList <- c(groupA, groupB)
        }
        cellList <- cellList & (seurat_object@meta.data[, variable] %in% subsetVariableGroupList)
    }

    seurat_object <- subset(seurat_object, cells = which(cellList))


    if (is.factor(seurat_object@meta.data[, variable])) {
        all_levels <- levels(droplevels(seurat_object@meta.data[, variable]))
    } else {
        all_levels <- unique(seurat_object@meta.data[, variable])
    }
    if (length(all_levels) < 2) {
        stop("Need 2 or more unique values in %s (have %d)", variable, length(all_levels))
    }

    topTable_outs <- list()
    model_outs <- list()
    SINGLE_TABLE <- NULL
    for (i in all_levels) {
        print(sprintf("Computing DEGs for Ident %s versus all others...", i))

        seurat_object@meta.data$cellname <- rownames(seurat_object@meta.data)
        seurat_object@meta.data <- seurat_object@meta.data |> dplyr::mutate(tmp_covariate = dplyr::case_when(
            !!sym(variable) == i ~ i,
            .default = sprintf("ALL_BUT_%s", i)
        ))
        rownames(seurat_object@meta.data) <- seurat_object$cellname

        tmp_outs <- doLimma(seurat_object, variable = "tmp_covariate", groupA = i, groupB = sprintf("ALL_BUT_%s", i), treatVariableAsContinuous = F, batchVariable = batchVariable, subsetVariableGroupList = NULL, batch_features = batch_features, RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES = RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES, DISCARD_CONTRAST_MODELS = DISCARD_MODELS, DISCARD_FULL_MODELS = DISCARD_MODELS)

        if (!DISCARD_MODELS) {
            model_outs[[i]] <- tmp_outs[["lmFit"]][[1]]
        }
        topTable_outs[[i]] <- tmp_outs[["topTable"]][[1]]
        topTable_outs[[i]]$cluster <- i

        if (is.null(SINGLE_TABLE)) {
            SINGLE_TABLE <- topTable_outs[[i]]
        } else {
            SINGLE_TABLE <- rbind(SINGLE_TABLE, topTable_outs[[i]])
        }
    }
    if (DISCARD_MODELS) {
        model_outs <- NULL
    }
    return(list("lmFit" = model_outs, "topTable" = topTable_outs, "mergedTable" = SINGLE_TABLE))
}


volcano_from_limma_table <- function(limmaTopTable, groupA, groupB, gene_lists = NULL, base_color = "grey", base_alpha = 0.5, highlight_alpha = 1, color_lists = NULL, shape_lists = NULL, pval_thresh = 10e-50, logFC_thresh = 0.1, plotTitle = "DEG", overlay_comparisons = T, groupA_up_color = "red", groupB_up_color = "blue", xlim = c(-0.5, 0.5), ylim = c(-5, 325), pointSize = 3.0, widthConnectors = 0.75, max.overlaps = 20, lengthConnectors = unit(0.007, "npc"), labSize = 4) {
    if (!(groupA %in% unique(limmaTopTable$group1) || groupA %in% unique(limmaTopTable$group2))) {
        stop(sprintf("%s not found in table...", groupA))
    }
    if (!(groupB %in% unique(limmaTopTable$group1) || groupB %in% unique(limmaTopTable$group2))) {
        stop(sprintf("%s not found in table...", groupB))
    }

    if (!is.null(xlim)) {
        X_EXTREME_VALUE <- max(abs(limmaTopTable$logFC))

        xlim <- c(-X_EXTREME_VALUE * 1.2, X_EXTREME_VALUE * 1.2)
    }

    if (!is.null(ylim)) {
        YMAX_OVERRIDE <- max(-log10(limmaTopTable$adj.P.Val + 10e-300))

        YMAX_OFFSET <- (YMAX_OVERRIDE) + (YMAX_OVERRIDE) * (25 / 305)

        ylim <- c(-5, YMAX_OFFSET)
    }

    limmaTopTable <- limmaTopTable |> dplyr::filter(group1 %in% c(groupA, groupB) & group2 %in% c(groupA, groupB))
    group1table <- unique(limmaTopTable$group1)
    group2table <- unique(limmaTopTable$group2)

    if (length(group1table) > 1 || length(group2table) > 1) {
        stop("bad assumption on table filtering. check")
    }

    if (group1table == groupB && group2table == groupA) {
        limmaTopTable$logFC <- -limmaTopTable$logFC
        limmaTopTable$t <- -limmaTopTable$t
        tmp1 <- limmaTopTable$pct.1
        tmp2 <- limmaTopTable$pct.2
        limmaTopTable$pct.2 <- tmp1
        limmaTopTable$pct.1 <- tmp2
    } else if (group1table == groupA && group2table == groupB) {

    } else {
        stop("Table groups do not match what was passed?? Tbl1: %s, Tbl2: %s, gA: %s, gB: %s", group1table, group2table, groupA, groupB)
    }


    if (!is.null(groupA_up_color)) {
        groupA_up_color <- "red"
    }
    if (!is.null(groupB_up_color)) {
        groupB_up_color <- "blue"
    }

    limmaTopTable$avg_log2FC <- limmaTopTable$log2FC
    limmaTopTable$p_val_adj <- limmaTopTable$adj.P.Val
    limmaTopTable$p_val <- limmaTopTable$P



    if (is.null(gene_lists)) {
        UP_GENES <- limmaTopTable |> dplyr::filter(p_val_adj < pval_thresh & avg_log2FC > logFC_thresh)
        DOWN_GENES <- limmaTopTable |> dplyr::filter(p_val_adj < pval_thresh & avg_log2FC < -logFC_thresh)

        up_str <- paste0(groupA, " UP")
        down_str <- paste0(groupA, " DOWN")

        gene_lists <- list("UP" = UP_GENES$gene, "DOWN" = DOWN_GENES$gene)
        names(gene_lists) <- c(up_str, down_str)

        color_lists <- list("UP" = groupA_up_color, "DOWN" = groupB_up_color)
        names(color_lists) <- c(up_str, down_str)
    }

    keyvals.color <- rep(base_color, nrow(limmaTopTable))
    names(keyvals.color) <- rep("Gene", nrow(limmaTopTable))

    keyvals.shape <- rep(1, nrow(limmaTopTable))
    names(keyvals.shape) <- rep("Gene", nrow(limmaTopTable))

    selectLab_df <- limmaTopTable |> dplyr::filter(p_val_adj < pval_thresh & abs(avg_log2FC) > logFC_thresh)
    selectLab <- selectLab_df$gene

    labels.color <- rep(base_color, length(selectLab))

    if (!is.null(color_lists)) {
        names(keyvals.color) <- rep("Other", length(names(keyvals.shape)))
        for (gene_groupings in names(gene_lists)) {
            keyvals.color[limmaTopTable$gene %in% gene_lists[[gene_groupings]]] <- color_lists[[gene_groupings]]

            names(keyvals.color)[limmaTopTable$gene %in% gene_lists[[gene_groupings]]] <- gene_groupings

            labels.color[selectLab %in% gene_lists[[gene_groupings]]] <- color_lists[[gene_groupings]]
        }
    }

    if (!is.null(shape_lists)) {
        names(keyvals.shape) <- rep("Other", length(names(keyvals.shape)))
        for (gene_groupings in names(gene_lists)) {
            keyvals.shape[limmaTopTable$gene %in% gene_lists[[gene_groupings]]] <- shape_lists[[gene_groupings]]
            names(keyvals.shape)[limmaTopTable$gene %in% gene_lists[[gene_groupings]]] <- gene_groupings
        }
    }

    highlight_genes <- c()
    for (genes in gene_lists) {
        highlight_genes <- c(highlight_genes, intersect)
    }


    keyvals.alpha <- rep(0.5, nrow(limmaTopTable))
    for (gene_groupings in names(gene_lists)) {
        keyvals.alpha[limmaTopTable$gene %in% gene_lists[gene_groupings]] <- 1
    }

    p <- limmaTopTable |> EnhancedVolcano(
        lab = limmaTopTable$gene,
        x = "avg_log2FC",
        y = "p_val_adj",
        title = plotTitle,
        selectLab = selectLab,
        labCol = labels.color,
        colCustom = keyvals.color,
        colAlpha = keyvals.alpha,
        caption = "",
        xlim = xlim,
        ylim = ylim,
        pCutoff = pval_thresh,
        FCcutoff = logFC_thresh,
        pointSize = pointSize,
        drawConnectors = T,
        widthConnectors = widthConnectors,
        max.overlaps = max.overlaps,
        lengthConnectors = lengthConnectors,
        labSize = labSize
    ) +
        theme_prism() +
        theme(plot.caption = element_blank(), plot.subtitle = element_blank())

    if (overlay_comparisons) {
        annotations <- data.frame(
            xpos = c(-Inf, Inf),
            ypos = c(Inf, Inf),
            annotateText = c(
                sprintf("%s Enriched", groupB),
                sprintf("%s Enriched", groupA)
            ),
            hjustvar = c(-.05, 1),
            vjustvar = c(0.9, 0.9)
        )
        p <- p + geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText, fontface = "bold"), color = c(groupB_up_color, groupA_up_color), size = 6)
    }

    return(p)
}

library("org.Hs.eg.db")
library(clusterProfiler)
library(ReactomePA)

pathwayAnalysis_From_Limma <- function(limmaTopTable, groupA, groupB, sortListBy = "logFC", PATHWAY = "reactome", pvalueCutoff = 1.0, pvalueCutoff_FOR_DISPLAY = 0.2, color_hi = "red", color_lo = "blue", ENSEMBLE_ID_COLUMN = NULL) {
    if (!(groupA %in% unique(limmaTopTable$group1) || groupA %in% unique(limmaTopTable$group2))) {
        stop(sprintf("%s not found in table...", groupA))
    }
    if (!(groupB %in% unique(limmaTopTable$group1) || groupB %in% unique(limmaTopTable$group2))) {
        stop(sprintf("%s not found in table...", groupB))
    }

    limmaTopTable <- limmaTopTable |> dplyr::filter(group1 %in% c(groupA, groupB) & group2 %in% c(groupA, groupB))
    group1table <- unique(limmaTopTable$group1)
    group2table <- unique(limmaTopTable$group2)

    if (group1table == groupB && group2table == groupA) {
        limmaTopTable$logFC <- -limmaTopTable$logFC
        limmaTopTable$t <- -limmaTopTable$t
        tmp1 <- limmaTopTable$pct.1
        tmp2 <- limmaTopTable$pct.2
        limmaTopTable$pct.2 <- tmp1
        limmaTopTable$pct.1 <- tmp2
    } else if (group1table == groupA && group2table == groupB) {

    } else {
        stop("Table groups do not match what was passed?? Tbl1: %s, Tbl2: %s, gA: %s, gB: %s", group1table, group2table, groupA, groupB)
    }

    limmaTopTable <- limmaTopTable |> dplyr::filter(group1 %in% c(groupA, groupB) & group2 %in% c(groupA, groupB))
    m_markers <- limmaTopTable |> arrange(pick(sortListBy))
    gene_list <- unlist(m_markers$gene)
    symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    colnames(symbols) <- c("gene", "ENTREZID")
    m_markers <- dplyr::left_join(symbols, m_markers, by = "gene") # Drops unmapped
    rownames(m_markers) <- m_markers$ENTREZID
    logfc_m <- m_markers[, sortListBy]
    names(logfc_m) <- rownames(m_markers)
    logfc_m <- logfc_m |> sort(decreasing = T)

    if (PATHWAY == "reactome") {
        GSEA_OUT <- gsePathway(logfc_m,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = "BH",
            verbose = FALSE
        )
    } else if (PATHWAY == "H") {
        m_df <- msigdbr(species = "Homo sapiens")
        m_H <- msigdbr(species = "Homo sapiens", category = "H") %>%
            dplyr::select(gs_name, gene_symbol)
        names(logfc_m) <- m_markers$gene

        GSEA_OUT <- GSEA(logfc_m, TERM2GENE = m_H, pvalueCutoff = pvalueCutoff)
    }

    filt <- GSEA_OUT |> dplyr::filter(p.adjust < pvalueCutoff_FOR_DISPLAY)
    results <- NULL
    if (nrow(filt) > 0) {
        results <- filt@result |> arrange(NES)
        # results <- cd8_gsea@result |> arrange(NES)
        results$Description_wrap <- stringr::str_wrap(results$Description, width = 45)
        results <- results |> dplyr::mutate(signed.p.adjust = -log10(p.adjust) * sign(NES))
        results$Description_wrap <- factor(results$Description_wrap, levels = results$Description_wrap)

        p2 <- results |>
            arrange(NES) |>
            ggplot() +
            geom_bar(aes(y = Description_wrap, x = NES, fill = signed.p.adjust), stat = "identity") +
            scale_fill_gradient2(
                low =
                    color_lo,
                mid = "white",
                high = color_hi
            ) +
            geom_vline(xintercept = 0, color = "black") +
            theme_prism() +
            theme(axis.text.y = element_text(size = 10)) +
            xlab("Normalized Enrichment Score") +
            ylab("Significantly Enriched Pathways") +
            labs(fill = "Sign(NES)*\nLog10(FDR)") +
            theme(legend.title = element_text())
    } else {
        p2 <- NULL
    }
    return(list("GSEA_out" = GSEA_OUT, "Filtered_GSEA_out" = results, "Figure" = p2, "Entries_in_figure" = nrow(filt)))
}


# Code for generaing dot plots, violin plots, and feature plots
# Used primarily for Figure 2.
plotMarkers <- function(seurat_object, markerList, plotLists = list(
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
    findMarkersVariable = NULL, displayVariable = NULL, method = "Presto", USEPRESTO = NULL, USE_ORIGINAL_PRESTO = F, limmaBatchVar = NULL, limmaBatchFeatures = 6000,
    FORCE_OVERWRITE = F, SCALE_FEATURES_IN_HEATMAP = T, geneText = 18, labelText = 24) {
    if (USEPRESTO) {
        method <- "Presto"
    }

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
        if (method == "Presto") {
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
        } else if (method == "Seurat") {
            mkrs <- tryCatch(FindAllMarkers(seurat_object,
                only.pos = T, min.pct = min.pct,
                logfc.threshold = logfc.threshold
            ), error = function(cond) {
                message("FIND MARKERS ERRROR:")
                message(cond)
                message("")
                return(NULL)
            })
        } else if (method == "Limma") {
            mkrs <- tryCatch(doLimmaAll(seurat_object, DISCARD_MODELS = T, RETURN_TOPTABLE_WITH_SEURAT_STYLE_NAMES = T, batch_features = limmaBatchFeatures, batchVariable = limmaBatchVar)$mergedTable)
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
