library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(presto)
library(SeuratWrappers)
library(ggprism)
library(slingshot)
library(SingleCellExperiment)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(Matrix)
library(grid)
library(gridGraphics)
library(tradeSeq)
library(clustree)
library(harmony)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(patchwork)
library(scater)
library(ggbeeswarm)
# Load with chdir = T !!!!

"%ni%" <- Negate("%in%")

source("pdf_and_png.R")
source("DEG_functions.R")

# Given a Seurat Object...
## Find Variable Feature
## Scale Data
## Compute PCA
## Correct with Harmony
## Cluster at Varying Resolution
## Find top cluster markers
# Save .RDS file with object, along with heatmaps of top markers at each resolution
# Used for subclustering clusters (e.g., NkT.2 -> NkT.2.0, NkT.2.1, ... etc.)
subset_and_markers <- function(
    object, name, res = 0.1, output = output, rebuildObject = F,
    rebuildMarkers = F, SKIP_MARKERS = F, FORCE_RECLUSTER = F, DIMPLOT_ADDITIONAL_DISPLAY = NULL,
    harmonyVariable = "siteXbatch", dims = 30, OLD_PC_BEHAVIOR = F, nfeatures = 2000) {
    # NOW DEFAULTING TO siteXbatch
    output <- file.path(output, name)
    dir.create(output, recursive = T)
    rdsOutput <- file.path(output, paste0(name, ".rds"))

    print(rebuildObject)
    print(file.exists(rdsOutput))
    pcDims <- dims
    if (OLD_PC_BEHAVIOR) {
        print("USING OLD PC BEHAVIOR - HARMONY COMPUTES ON 50PCs THEN IS SUBSET")
        pcDims <- NULL
    }
    changesOccurred <- FALSE
    # Batch corretion and UMAP generation. Batch correction only uses harmony. Skip if object already exists, unless user has specified for object to be recreated. Unaware of object changes
    # NOTE - could change to check parameters of existing object?
    if (rebuildObject == T || !file.exists(rdsOutput)) {
        if (is.null(harmonyVariable)) {
            validateHarmony <- 0
        } else {
            validateHarmony <- length(unique(object@meta.data[, harmonyVariable]))
        }
        if (validateHarmony > 1) {
            print("REBUILDING OBJECT")
            changesOccurred <- TRUE
            object <- object %>%
                FindVariableFeatures(nfeatures = nfeatures) %>%
                ScaleData() %>%
                RunPCA(npcs = pcDims) %>%
                RunHarmony(dims.use = 1:dims, group.by.vars = harmonyVariable) %>%
                FindNeighbors(dims = 1:dims, reduction = "harmony") %>%
                RunUMAP(
                    reduction = "harmony", reduction.name = "umapsub", redution.key = "UMAPSUB_",
                    dims = 1:dims
                )
        } else {
            changesOccurred <- TRUE
            print("REBUILDING OBJECT - HARMONY SKIPPED")
            object <- object %>%
                FindVariableFeatures() %>%
                ScaleData() %>%
                RunPCA() %>%
                FindNeighbors(dims = 1:dims) %>%
                RunUMAP(
                    reduction.name = "umapsub", redution.key = "UMAPSUB_",
                    dims = 1:dims
                )
        }
    } else {
        print("RELOADING OBJECT")
        object <- readRDS(rdsOutput)
    }

    # Cluster generation and cluster tree
    # change harmony validation check
    # If new clusters are being requested, generate only those clusters. Resave object if new clusters are generated. Generate clustertree plots.
    # NOTE: If dims changes, user must explicilty pass the 'FORCE_RECLUSTER' flag to recompute clusters.
    # NOTE: Could be made more efficient. Filter out existing clusters from the 'res' vector, then pass all into find clusters. Will run resolutions in parallel threads.
    if (!is.null(res)) {
        for (i in res) {
            if (!(paste0("RNA_snn_res.", i) %in% colnames(object@meta.data)) || FORCE_RECLUSTER || i == 0.5) {
                print(paste0("NEW cluster res ", i))
                changesOccurred <- TRUE
                object <- object %>%
                    FindClusters(dims = 1:dims, res = i, reduction = "harmony")
                object[[paste0("seurat_clusters_", i)]] <- object$seurat_clusters
            }
        }
        if (changesOccurred) {
            paste0("Cluster Object updated, overwritting previous")
            saveRDS(object, rdsOutput)
        }


        ct <- clustree(object)
        pdf_and_png(ct, output, paste0(name, "CLUSTTREE"), pdfWidth = 8, pdfHeight = 11)

        cto_label <- clustree_overlay(object,
            red_dim = "umap", x_value = "umap1",
            y_value = "umap2", label_nodes = TRUE
        )
        cto <- clustree_overlay(object,
            red_dim = "umap", x_value = "umap1", y_value = "umap2",
            label_nodes = FALSE
        )
        pdf_and_png(cto_label, output, paste0(name, "CLUSTTREE_OVERLAY_LABEL"),
            pdfWidth = 8, pdfHeight = 11
        )
        # pdf(file.path(output, paste0(name, 'CLUSTTREE_OVERLAY_LABEL.pdf')),
        # width = 11) print(cto_label) dev.off()
        pdf_and_png(cto, output, paste0(name, "CLUSTTREE_OVERLAY"),
            pdfWidth = 8,
            pdfHeight = 11
        )
        # pdf(file.path(output, paste0(name, 'CLUSTTREE_OVERLAY.pdf')), width
        # = 11) print(cto) dev.off()

        ctosub_label <- clustree_overlay(object,
            red_dim = "umapsub", x_value = "umapsub1",
            y_value = "umapsub2", label_nodes = TRUE
        )
        ctosub <- clustree_overlay(object,
            red_dim = "umapsub", x_value = "umapsub1",
            y_value = "umapsub2", label_nodes = FALSE
        )
        pdf_and_png(ctosub_label, output, paste0(name, "CLUSTTREE_OVERLAY_SUB_LABEL"),
            pdfWidth = 8, pdfHeight = 11
        )
        pdf_and_png(ctosub, output, paste0(name, "CLUSTTREE_OVERLAY_SUB"),
            pdfWidth = 8,
            pdfHeight = 11
        )

        for (i in res) {
            p <- DimPlot(object,
                group.by = paste0("seurat_clusters_", i), label = T,
                repel = T
            )


            pdf_and_png(p, output, paste0(name, "_", i), pdfWidth = 11)

            # pdf(file.path(output, paste0(name, '_', i, '.pdf')), width =
            # 11) print(p) dev.off()
            p <- DimPlot(object,
                group.by = paste0("seurat_clusters_", i), reduction = "umapsub",
                label = T, repel = T
            )
            pdf_and_png(p, output, paste0(name, "_", i, "_NEWUMAP"), pdfWidth = 11)
            # pdf(file.path(output, paste0(name, '_', i, '_NEWUMAP.pdf')),
            # width = 11) print(p) dev.off()
        }
    } else {
        if (changesOccurred) {
            print("Object built with no clustering... saving")
            saveRDS(object, rdsOutput)
        }
    }

    # Generates a DimPlot by an additional group.by variable. Skipped by default
    # NOTE: Allow group.by to be a vector
    if (!is.null(DIMPLOT_ADDITIONAL_DISPLAY)) {
        p <- DimPlot(object,
            group.by = DIMPLOT_ADDITIONAL_DISPLAY, label = T,
            repel = T
        )
        pdf_and_png(p, output, paste0(name), pdfWidth = 11)
        p <- DimPlot(object,
            group.by = DIMPLOT_ADDITIONAL_DISPLAY, label = T,
            repel = T, reduction = "umapsub"
        )
        pdf_and_png(p, output, paste0(name, "_NEWUMAP"), pdfWidth = 11)
    }

    # Avoid futures - parallel in parallel is awkward.
    # Generates markers for each cluster resolution. Displays a heatmap. Can be skipped.
    futures.list <- list()
    counter <- 1
    if (!SKIP_MARKERS) {
        for (i in res) {
            markerOutput <- file.path(output, paste0(name, "_", res), paste0(
                name,
                "_", i, "_MKR_HEATMAP.pdf"
            ))

            Idents(object) <- object@meta.data[[paste0("seurat_clusters_", i)]]
            if ((rebuildMarkers || !file.exists(markerOutput))) {
                if (length(unique(Idents(object))) > 1) {
                    print("REBUILDING MARKERS")
                    futures.list[[counter]] <- markers_and_heatmap(object,
                        n_markers = 10,
                        nameMod = paste0(name, "_", i), output = output
                    )
                    counter <- counter + 1
                }
            }
        }
    }
    # for (futureMkrs in futures.list) { tmpVal <- value(futureMkrs, signal =
    # F) print(tmpVal) }
    return(object)
}

# Perform trajectory analysis on a given seurat object.
runSlingshot <- function(
    seurat_object, groupingVar = "seurat_clusters", output = ".",
    nameMod = "", slingshotDim = "HARMONY", dims = NULL, visualizationDim = "UMAP", startCluster = NULL, endCluster = NULL,
    CREATE_SUBFOLDER_WITH_NAMEMOD = T, REBUILD_TRAJECTORY = F, DO_RF_GENES = F) {
    if (CREATE_SUBFOLDER_WITH_NAMEMOD) {
        output <- file.path(output, nameMod)
        dir.create(output, showWarnings = F)
    }
    print("CREATING SCE")
    Idents(seurat_object) <- seurat_object@meta.data[[groupingVar]]

    seurat_object.sce <- as.SingleCellExperiment(seurat_object) # change to seurat_object...seeing if its a seurat_clusters issue
    # Select number of dims
    if (!is.null(dims)) {
        reducedDim(
            seurat_object.sce,
            slingshotDim
        ) <- reducedDim(
            seurat_object.sce,
            slingshotDim
        )[, dims]
    }
    cl <- seurat_object.sce$ident
    if (file.exists(file.path(output, paste0("sce_", nameMod, ".rds"))) & !REBUILD_TRAJECTORY) {
        print("RELOADING TRAJECTORY")
        seurat_object.sce <- readRDS(file.path(output, paste0(
            "sce_", nameMod,
            ".rds"
        )))
        # seurat_object.lnes <-
        # readRDS(file.path(output,paste0('lnes_',nameMod,'.rds')))
    } else {
        print("SLINGSHOT START")
        seurat_object.sce <- slingshot(seurat_object.sce,
            clusterLabels = "ident",
            reducedDim = slingshotDim, allow.breaks = T, start.clus = startCluster, end.clus = endCluster
        )
        print("LINEAGES START")
        print("DONE")
        saveRDS(seurat_object.sce, file.path(output, paste0("sce_", nameMod, ".rds")))
        # saveRDS(seurat_object.lnes,file.path(output,paste0('lnes_',nameMod,'.rds')))
    }
    sds <- SlingshotDataSet(seurat_object.sce)
    print("generating .lnes")


    for (reduc in visualizationDim) {
        # seurat_object.crvs[[reduc]] <- getCurves(seurat_object.lnes[[reduc]])
        print("colors")
        my_color <- createPalette(length(unique(seurat_object.sce$ident)), c(
            "#010101",
            "#FF0000"
        ), M = 1000)
        names(my_color) <- unique(as.character(seurat_object.sce$ident))
        # p <- ggplot(slingshot_df,
        # aes(x=slingPseudotime_1,y=ident,colour=ident)) +
        # geom_quasirandom(groupOnX=F) + theme_classic() + xlab('First
        # Slingshot Pseudotime') + ylab('cell type') + ggtitle('Cells ordered
        # by slingshot pseudotime') + scale_colour_manual(values=my_color)
        # pdf_and_png(p,output,'slingshotB_test',pdfWidth=11)
        print("plot")
        # browser()

        X <- reducedDim(seurat_object.sce, reduc) |> as.data.frame()
        colnames(X)[1:2] <- c("Dim1", "Dim2")

        clusterLabels <- slingClusterLabels(sds)
        connectivity <- slingMST(sds, as.df = T)
        clusters <- unique(connectivity$Cluster)
        nclus <- length(clusters)
        centers <- t(vapply(clusters, function(clID) {
            w <- clusterLabels[, clID]
            return(apply(X, 2, weighted.mean, w = w))
        }, rep(0, ncol(X)))) # gets center of UMAP points
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
            drop = FALSE
        ]

        centers <- centers |> as.data.frame()
        centers$Cluster <- rownames(centers)

        p <- ggplot(X, aes(x = Dim1, y = Dim2)) +
            geom_point(aes(fill = seurat_object.sce$ident), col = "grey70", shape = 21) +
            theme_classic()

        mst_with_centers <- dplyr::left_join(connectivity, centers, by = "Cluster")

        p <- p + geom_point(data = mst_with_centers, size = 4) +
            geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
            geom_label_repel(data = centers, aes(label = Cluster))

        pdf_and_png(p, output, paste0("TRAJECTORY_", nameMod, "_LINES_OVERLAY_", reduc, ".pdf"), pdfHeight = 7, pdfWidth = 11, splitDirectories = T)

        p <- ggplot(data = mst_with_centers, aes(x = Dim1, y = Dim2)) +
            geom_point(data = mst_with_centers, size = 4) +
            geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
            geom_label_repel(data = centers, aes(label = Cluster)) +
            theme_classic()

        pdf_and_png(p, output, paste0("TRAJECTORY_", nameMod, "_LINES_NOPOINTS_", reduc, ".pdf"), pdfHeight = 7, pdfWidth = 11, splitDirectories = T)


        lineages <- unique(mst_with_centers$Lineage)
        for (i in lineages) {
            mst_subset <- mst_with_centers |> dplyr::filter(Lineage == i)

            p <- ggplot(X, aes(x = Dim1, y = Dim2)) +
                geom_point(aes(fill = seurat_object.sce$ident), col = "grey70", shape = 21) +
                theme_classic()

            mst_subset$clustLabel <- paste0(mst_subset$Order, ": ", mst_subset$Cluster)

            p <- p + geom_point(data = mst_subset, size = 4) +
                geom_path(data = mst_subset %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = mst_subset, aes(label = clustLabel))

            pdf_and_png(p, output, paste0("TRAJECTORY_", nameMod, "_LINES_OVERLAY_", reduc, "_lineage", i, ".pdf"), pdfHeight = 7, pdfWidth = 11, splitDirectories = T)

            p <- ggplot(data = mst_with_centers, aes(x = Dim1, y = Dim2)) +
                geom_point(data = mst_subset, size = 4) +
                geom_path(data = mst_subset %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = mst_subset, aes(label = clustLabel)) +
                theme_classic()

            pdf_and_png(p, output, paste0("TRAJECTORY_", nameMod, "_LINES_NOPOINTS_", reduc, "_lineage", i, ".pdf"), pdfHeight = 7, pdfWidth = 11, splitDirectories = T)
        }

        # Can't do principal curves without some harmony -> umap mapping
        print("INIT PLOTS DONE")

        nc <- 3
        pt <- slingPseudotime(seurat_object.sce)
        nms <- 1:ncol(pt) # colnames(pt)
        nr <- ceiling(length(nms) / nc)
        pal <- viridis(100, end = 0.95)


        par(mfrow = c(nr, nc))

        # getContour_ident <- function(X_pt, targetIdent) {
        #     x <- X_pt |> dplyr::filter(ident == targetIdent)

        #     kd <- ks::kde(x[, c("Dim1", "Dim2")], gridsize = rep(401, 2), compute.cont = TRUE)
        #     contour_95 <- with(kd, contourLines(
        #         x = eval.points[[1]], y = eval.points[[2]],
        #         z = estimate, levels = cont["5%"]
        #     )[[1]])
        #     contour_95 <- data.frame(contour_95)
        # }
        # getContour_allident <- function(X_pt) {
        #     ident_list <- unique(X_pt$ident)
        #     out <- lapply(ident_list, function(x) {
        #         tmp <- getContour_ident(X_pt, x)
        #         tmp$ident <- x
        #         return(tmp)
        #     })
        #     out <- do.call(rbind, out)
        # }
        X_id <- X
        X_id$ident <- seurat_object.sce$ident
        # allContours <- getContour_allident(X_id)
        for (i in nms) {
            pt_tmp <- pt
            colnames(pt_tmp)[i] <- "Pseudotime"
            X_pt <- cbind(X_id, pt_tmp)

            X_pt <- X_pt |> arrange(!is.na(Pseudotime), Pseudotime)

            ptest <- ggplot(X_pt, aes(x = Pseudotime, y = ident)) +
                geom_jitter(aes(fill = ident, color = ident), alpha = 0.5) +
                geom_boxplot(aes(fill = ident, group = ident), color = "black", alpha = 0.5) +
                theme_classic()

            pdf_and_png(ptest, output, paste0("PSEUDOTIME_", nameMod, "_lin_", i, "_per_cluster_pseudotime"), pdfWidth = 8, pdfHeight = 2 + 0.5 * length(unique(X_pt$ident)), splitDirectories = T)

            p <- ggplot(X_pt, aes(x = Dim1, y = Dim2)) +
                geom_point(aes(fill = Pseudotime), color = "grey70", stroke = 0, shape = 21) +
                theme_classic() +
                scale_fill_gradientn(colours = pal)
            # geom_bag(data = X_id, aes(group = ident, color = ident), alpha=0, prop=0.80)

            p_base <- p + geom_point(data = mst_with_centers, size = 4) +
                geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = centers, aes(label = Cluster))

            pdf_and_png(p_base, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)
            mst_subset <- mst_with_centers |> dplyr::filter(Lineage == i)
            mst_subset$clustLabel <- paste0(mst_subset$Order, ": ", mst_subset$Cluster)

            p_alt <- p + geom_point(data = mst_subset, size = 4) +
                geom_path(data = mst_subset %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = mst_subset, aes(label = clustLabel))

            pdf_and_png(p_alt, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i, "_REMOVEOTHERLINEAGES"), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)

            X_pt <- X_pt[!is.na(X_pt$Pseudotime), ]
            p <- ggplot(X_pt, aes(x = Dim1, y = Dim2)) +
                geom_point(aes(fill = Pseudotime), color = "grey70", stroke = 0, shape = 21) +
                theme_classic() +
                scale_fill_gradientn(colours = pal)
            # geom_bag(data = X_id, aes(group = ident, color = ident), alpha=0, prop=0.80)

            p_base <- p + geom_point(data = mst_with_centers, size = 4) +
                geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = centers, aes(label = Cluster))

            pdf_and_png(p_base, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i, "_NAOMIT"), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)

            p_alt <- p + geom_point(data = mst_subset, size = 4) +
                geom_path(data = mst_subset %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = mst_subset, aes(label = clustLabel))

            pdf_and_png(p_alt, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i, "_NAOMIT", "_REMOVEOTHERLINEAGES"), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)
        }

        pt <- slingPseudotime(seurat_object.sce, na = F)

        for (i in nms) {
            pt_tmp <- pt
            colnames(pt_tmp)[i] <- "Pseudotime"
            X_pt <- cbind(X_id, pt_tmp)

            X_pt <- X_pt |> arrange(!is.na(Pseudotime), Pseudotime)

            ptest <- ggplot(X_pt, aes(x = Pseudotime, y = ident)) +
                geom_jitter(aes(fill = ident, color = ident), alpha = 0.5) +
                geom_boxplot(aes(fill = ident, group = ident), color = "black", alpha = 0.5) +
                theme_classic()

            pdf_and_png(ptest, output, paste0("PSEUDOTIME_", nameMod, "_", i, "per_cluster_pseudotime_NAINTERP"), pdfWidth = 8, pdfHeight = 2 + 0.5 * length(unique(X_pt$ident)), splitDirectories = T)

            p <- ggplot(X_pt, aes(x = Dim1, y = Dim2)) +
                geom_point(aes(fill = Pseudotime), color = "grey70", stroke = 0, shape = 21) +
                theme_classic() +
                scale_fill_gradientn(colours = pal)
            # geom_bag(data = X_id, aes(group = ident, color = ident), alpha=0, prop=0.80)

            p_base <- p + geom_point(data = mst_with_centers, size = 4) +
                geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = centers, aes(label = Cluster))

            pdf_and_png(p_base, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i, "_NAINTERP"), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)
            mst_subset <- mst_with_centers |> dplyr::filter(Lineage == i)
            mst_subset$clustLabel <- paste0(mst_subset$Order, ": ", mst_subset$Cluster)

            p_alt <- p + geom_point(data = mst_subset, size = 4) +
                geom_path(data = mst_subset %>% arrange(Order), aes(group = Lineage), size = 2) +
                geom_label_repel(data = mst_subset, aes(label = clustLabel))

            pdf_and_png(p_alt, output, paste0("PSEUDOTIME_", nameMod, "_", reduc, "lin", i, "_NAINTERP", "_REMOVEOTHERLINEAGES"), pdfWidth = 11, pdfHeight = 7, splitDirectories = T)
        }
    }
    # print("HVF EXTRACT")
    # dynverse.. https://bustools.github.io/BUS_notebooks_R/slingshot.html
    # Get top highly variable genes
    if (DO_RF_GENES) {
        top_hvg <- HVFInfo(seurat_object) %>%
            mutate(., bc = rownames(.)) %>%
            arrange(desc(variance.standardized)) %>%
            top_n(300, variance.standardized) %>%
            pull(bc)
        # Prepare data for random forest
        dat_use <- t(GetAssayData(seurat_object, slot = "data")[top_hvg, ])
        top_genes_list <- list()
        val_results_list <- list()
        truth_list <- list()
        estimate_list <- list()
        var_imp_list <- list()
        for (k in 1:ncol(slingPseudotime(seurat_object.sce))) {
            dat_use_df <- cbind(slingPseudotime(seurat_object.sce)[, k], dat_use) # FROM SITE: - but only one curve - will vary, may split out to other block or loop through curves? Do curve 2, so 2nd columnn
            colnames(dat_use_df)[1] <- "pseudotime"
            dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[, 1]), ])

            dat_split <- initial_split(dat_use_df)
            dat_train <- training(dat_split)
            dat_val <- testing(dat_split)
            # random forest, based on tidymodels
            model_cache_string <- paste0("tidymodel_", nameMod, "_", k)
            if (!file.exists(file.path(output, model_cache_string)) && !REBUILD_TRAJECTORY) {
                print("TRAINING")

                model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
                    set_engine("ranger", importance = "impurity", num.threads = 16) %>%
                    fit(pseudotime ~ ., data = dat_train)
                saveRDS(model, file.path(output, model_cache_string))
                print("TRIANING DONE")
            } else {
                print("RELOADING")

                model <- readRDS(file.path(output, model_cache_string))
            }
            val_results <- dat_val %>%
                dplyr::mutate(estimate = predict(model, .[, -1]) %>%
                    pull()) %>%
                dplyr::select(truth = pseudotime, estimate)
            print(metrics(data = val_results, truth, estimate))

            var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
            top_genes <- names(var_imp)[1:12]
            # Convert to gene symbol top_gene_name <-
            # gns$gene_name[match(top_genes, gns$gene)]
            top_gene_name <- top_genes
            for (reduc in visualizationDim) {
                print("PLOTTING")

                for (i in seq_along(top_genes)) {
                    colors <- pal[cut(dat_use[, top_genes[i]], breaks = 100)]
                    plot(reducedDim(seurat_object.sce, reduc),
                        col = colors, pch = 16,
                        cex = 0.5, main = top_gene_name[i]
                    )
                    lines(SlingshotDataSet(seurat_object.lnes[[reduc]]),
                        lwd = 2, col = "black",
                        type = "lineages"
                    )
                }
                grid.echo()
                p <- grid.grab()
                pdf_and_png(p, output, paste0(
                    "SLINGSHOT_", nameMod, "_differentially_expressed_genes_",
                    reduc, "_CURVE-", k
                ), pdfHeight = 8, pdfWidth = 12)
            }
            var_imp_list[[i]] <- var_imp
            top_genes_list[[i]] <- top_genes
            val_results_list[[i]] <- val_results
        }
    } else {
        var_imp_list <- NULL
        top_genes_list <- NULL
        val_results_list <- NULL
    }
    thingsToReturn <- list(
        TrajectoryObject = seurat_object.sce,
        numLineages = nms,
        AllGenes = var_imp_list,
        TopGenes = top_genes_list,
        val_results = val_results_list
    )

    return(thingsToReturn)
}



library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)

# (Not Used) Perform FGSEA between two idents
doFGSEA <- function(object, pathways = c("H", "C2", "C5", "C7"), species = "Homo sapiens", group.by = NULL, ident.1, ident.2 = NULL, nameMod = "FGSEA", sortingMethod = "pval", output, RERUN_MARKERS = F) {
    if (!is.null(group.by)) {
        Idents(object) <- object@meta.data[, group.by]
    }
    subDir <- file.path(output, nameMod)
    dir.create(subDir)
    if (!RERUN_MARKERS && file.exists(file.path(subDir, paste0("sorted_markers_", nameMod, ".rds")))) {
        ranks <- readRDS(file.path(subDir, paste0("sorted_markers_", nameMod, ".rds")))
    } else {
        Deg <- RunPresto(object, ident.1 = ident.1, ident.2 = ident.2, min.pct = 0.1, logfc.threshold = 0)
        Deg <- Deg |> dplyr::filter(p_val_adj < 1.0) # to eliminate ties?
        if (sortingMethod == "pval") {
            ranks <- -1 * sign(Deg$avg_log2FC) * log(Deg$p_val_adj + .Machine$double.eps)
        } else {
            ranks <- Deg$avg_log2FC
        }
        names(ranks) <- rownames(Deg)

        write.csv(ranks, file.path(subDir, paste0("sorted_markers_", nameMod, ".csv")))
        write.csv(Deg, file.path(subDir, paste0("full_markers_", nameMod, ".csv")))
        saveRDS(ranks, file.path(subDir, paste0("sorted_markers_", nameMod, ".rds")))
    }
    for (pathway in pathways) {
        print(pathway)
        if (pathway != "reactome") {
            if (pathway == "C2CP") {
                print("C2CP Override")

                m_df <- msigdbr(species = species, category = "C2")
                m_df <- dplyr::filter(m_df, stringr::str_detect(gs_subcat, "^CP"))
            } else if (pathway == "C2CGP") {
                print("C2CGP Override")

                m_df <- msigdbr(species = species, category = "C2", subcategory = "CGP")
            } else if (pathway == "C7") {
                print("C7 Override")
                m_df <- msigdbr(species = species, category = "C7", subcategory = "IMMUNESIGDB")
            } else if (pathway == "C3TFT") {
                print("C3 Override")
                m_df <- msigdbr(species = species, category = "C3", subcategory = "TFT:GTRD")
            } else {
                print("downloading")
                m_df <- msigdbr(species = species, category = pathway)
            }

            fgsea_sets <- split(x = m_df$gene_symbol, f = m_df$gs_name)
            print("sets obtained")
            fgseaRes <- fgsea(
                pathways = fgsea_sets,
                stats = ranks,
                minSize = 15,
                maxSize = 500,
                nperm = 100000
            )
        } else if (pathway == "reactome") {
            tmp <- mapIds(org.Hs.eg.db,
                keys = names(ranks),
                column = "ENTREZID", keytype = "SYMBOL"
            )
            react_ranks <- ranks
            names(react_ranks) <- tmp
            fgseaRes <- fgsea(
                pathways = reactomePathways(names(react_ranks)),
                stats = react_ranks,
                minSize = 15,
                maxSize = 500,
                nperm = 100000
            )
        }
        print("gsea done")
        saveRDS(fgseaRes, file.path(subDir, paste0("fgseaRes_", pathway, "_.rds")))
        fgseaResTidy <- fgseaRes %>%
            as_tibble() %>%
            arrange(desc(NES))
        write.csv(unnest_wider(as.data.frame(fgseaResTidy), "leadingEdge", names_sep = "."), file.path(subDir, paste0("fgseaResTidy_", pathway, "_.csv")))


        p <- ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n = 20), aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = NES < 7.5)) +
            coord_flip() +
            labs(
                x = "Pathway", y = "Normalized Enrichment Score",
                title = paste0(pathway, " pathways NES from GSEA")
            ) +
            theme_minimal()

        pdf_and_png(p, subDir, paste0(nameMod, "_GSEA_", pathway), pdfWidth = 8, pdfHeight = 8)


        testCatch <- tryCatch({
            test <- fgseaResTidy |>
                filter(size >= 6) |>
                filter(pval < 0.05)

            if (nrow(test) > 19) {
                fgseaRes_display <- test
            } else {
                fgseaRes_display <- fgseaResTidy
            }
            p <- fgseaRes_display |>
                arrange(abs(NES)) |>
                dplyr::slice(1:10, (n() - 9):n()) |>
                ggplot(aes(reorder(pathway, NES), NES)) + # %>% head(n= 20)
                geom_col(aes(fill = -log10(pval))) +
                coord_flip() +
                labs(
                    x = pathway, y = "Normalized Enrichment Score",
                    title = paste0(pathway, " fgsea findmarkerss for ", nameMod)
                ) +
                theme_minimal()
            pdf_and_png(p, subDir, paste0("ALT_", nameMod, "_GSEA_", pathway), pdfWidth = 8, pdfHeight = 8)
        })
    }
}
