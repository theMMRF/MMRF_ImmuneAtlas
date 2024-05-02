source("general.R", chdir = TRUE)
library(slingshot)
library(ggplot2)
library(dplyr)
library(viridis)

getClusterPositions <- function(sce, reduc = "UMAP") {
    sds <- SlingshotDataSet(sce)
    X <- reducedDim(sce, reduc) |> as.data.frame()
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

    return(list("centers" = centers, "connectivity" = connectivity))
}


plotLineages <- function(sce, reduc = "UMAP", selectLineage = NULL, pseudotimeLineage = NULL, label_map = NULL) {
    X <- reducedDim(sce, reduc) |> as.data.frame()
    colnames(X)[1:2] <- c("Dim1", "Dim2")


    out <- getClusterPositions(sce, reduc = reduc)

    centers <- out[["centers"]]
    connectivity <- out[["connectivity"]]

    if (is.null(pseudotimeLineage)) {
        p <- ggplot(X, aes(x = Dim1, y = Dim2)) +
            geom_point(aes(fill = sce$ident), col = "grey70", shape = 21) +
            theme_classic()
    } else {
        X_id <- X
        X_id$ident <- sce$ident
        if (pseudotimeLineage == "AVG") {
            pt <- slingAvgPseudotime(sce) |> as.data.frame()
            colnames(pt) <- "Pseudotime"
        } else {
            pt <- slingPseudotime(sce)
            colnames(pt_tmp)[pseudotimeLineage] <- "Pseudotime"
        }
        nms <- 1:ncol(pt) # colnames(pt)
        pal <- viridis(100, end = 0.95)

        pt_tmp <- pt
        X_pt <- cbind(X_id, pt_tmp)

        X_pt <- X_pt |> arrange(!is.na(Pseudotime), Pseudotime)

        p <- ggplot(X_pt, aes(x = Dim1, y = Dim2)) +
            geom_point(aes(fill = Pseudotime), col = "grey70", shape = 21) +
            theme_classic() +
            scale_fill_gradientn(colours = pal)
    }
    mst_with_centers <- dplyr::left_join(connectivity, centers, by = "Cluster")

    if (!is.null(selectLineage)) {
        mst_with_centers <- mst_with_centers |> dplyr::filter(Lineage == selectLineage)
    }


    if (!is.null(label_map)) {
        mst_with_centers$ident <- mst_with_centers$Cluster
        mst_with_centers <- dplyr::left_join(mst_with_centers, label_map, by = "ident")
        mst_with_centers$Cluster <- mst_with_centers$label
    }

    if (is.null(selectLineage)) {
        centers <- mst_with_centers |> dplyr::distinct(Cluster, Dim1, Dim2)
        p <- p + geom_point(data = mst_with_centers, size = 4) +
            geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
            geom_label_repel(data = centers, aes(label = paste0(Cluster)))

        p <- p + ggtitle("All Lineages")
    } else {
        p <- p + geom_point(data = mst_with_centers, size = 4) +
            geom_path(data = mst_with_centers %>% arrange(Order), aes(group = Lineage), size = 2) +
            geom_label_repel(data = mst_with_centers, aes(label = paste0(Order, ": ", Cluster)))

        p <- p + ggtitle("Pseudotime: ", selectLineage)
    }
    return(p)
}


plotGenePseudotime <- function(sce, i, geneName, scaledPseudo = F, i2 = NULL, i3 = NULL, color = NULL) {
    if (i == "ALL") {
        pt <- sce$avgpseudo
    } else {
        pt <- slingPseudotime(sce)
    }
    if (is.null(color)) {
        color <- "blue"
    }
    if (scaledPseudo) {
        for (i in 1:ncol(pt)) {
            pt[, i] <- pt[, i] / max(pt[, i], na.rm = T)
        }
    }
    w <- slingCurveWeights(sce)
    p <- ggplot()
    colors <- list("blue", "red")
    idx <- 0
    if (i == "ALL") {
        df <- data.frame(pseudotime = pt, weight = rep(1, length(pt)))
    } else {
        df <- data.frame(pseudotime = pt[, i], weight = w[, i])
    }
    for (geneid in geneName) {
        y <- logcounts(sce)[geneid, , drop = FALSE][1, ]
        df[[geneid]] <- y
    }
    for (geneid in geneName) {
        idx <- idx + 1

        p <- p +
            # geom_point(col = "grey70") +
            geom_smooth(data = df, aes(x = pseudotime, y = !!sym(geneid), weight = weight, color = geneName[[idx]]), method = "gam", formula = y ~ s(x))
        if (!is.null(i2)) {
            df2 <- data.frame(pseudotime = pt[, i2], weight = w[, i2], gene = y)
            p <- p + geom_smooth(data = df2, method = "gam", formula = y ~ s(x), aes(weight = weight, color = paste0("Pseudotime ", i2)))
        }
        if (!is.null(i3)) {
            df3 <- data.frame(pseudotime = pt[, i3], weight = w[, i3], gene = y)
            p <- p + geom_smooth(data = df3, method = "gam", formula = y ~ s(x), aes(weight = weight, color = paste0("Pseudotime ", i3)))
        }
    }
    p <- p + theme_classic() +
        ylab("Gene Expression") +
        xlab(paste0("Pseudotime"))
    if (is.null(i2)) {
        i2 <- ""
    }
    if (is.null(i3)) {
        i3 <- ""
    }
    pdf_and_png(p, output = output, filename = paste0(geneName, "_", i, i2, i3), pdfWidth = 6, pdfHeight = 4)
    return(p)
}

library(SingleCellExperiment)
plotOrderedDensity <- function(sce, reduc = "UMAP", lineage = "AVG", group.by = "ident", split.by = NULL, allow_na_psuedotime = TRUE, output = ".", nameMod = "", pdfWidth = 10, pdfHeight = 12) {
    X <- reducedDim(sce, reduc) |> as.data.frame()
    colnames(X)[1:2] <- c("Dim1", "Dim2")

    X$ident <- colData(sce)[[group.by]]

    if (lineage == "AVG") {
        pt <- slingAvgPseudotime(sce) |> as.data.frame()
        colnames(pt) <- "Pseudotime"
    } else {
        pt <- slingPseudotime(sce, na = allow_na_psuedotime)
        colnames(pt)[lineage] <- "Pseudotime"
    }

    if (!is.null(split.by)) {
        X_id$split.by <- colData(sce)[[split.by]]
    } else {
        X_id <- X
    }
    X_pt <- cbind(X_id, pt)

    X_pt <- X_pt |> arrange(!is.na(Pseudotime), Pseudotime)

    test <- X_pt |>
        group_by(ident) |>
        summarise(median = median(Pseudotime, na.rm = T)) |>
        arrange(median)

    X_pt$ident <- factor(X_pt$ident, levels = test$ident)
    ptest <- ggplot(X_pt, aes(x = Pseudotime)) +
        geom_density(aes(fill = ident, group = ident), color = "black", alpha = 0.5) +
        theme_prism() +
        facet_grid(ident ~ .)

    pdf_and_png(ptest, output, paste0("ordered_lin_", lineage, "_plot"), pdfWidth = pdfWidth, pdfHeight = pdfHeight)

    return(list("ordering" = test$ident, "plot" = ptest))
}
