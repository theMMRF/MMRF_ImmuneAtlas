library(cowplot)



# https://stackoverflow.com/questions/48531257/ggplot-geom-point-how-to-set-font-of-custom-plotting-symbols
point_with_family <- function(layer, family) {
    old_geom <- layer$geom
    new_geom <- ggproto(
        NULL, old_geom,
        draw_panel = function(self, data, panel_params, coord, na.rm = FALSE) {
            pts <- ggproto_parent(GeomPoint, self)$draw_panel(
                data, panel_params, coord,
                na.rm = na.rm
            )
            pts$gp$fontfamily <- family
            pts
        },
        draw_key = function(self, data, params, size) {
            pts <- ggproto_parent(GeomPoint, self)$draw_key(
                data, params, size
            )
            pts$gp$fontfamily <- family
            pts
        }
    )
    layer$geom <- new_geom
    layer
}

# Modification of Seurat's dotplot object to 1) Allow for grouping of cell types and 2) Replace facets with patchworks so that coord_fixed() could be used.
# This does not touch any of the code for computation of z-scores, % expression, and otherwise. The modifications only lie in the actual rendering of the plot
# Called by plotMarkers()

customDotPlot <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA,
    idents.grouping = NA,
    USE_PATCHWORK_INSTEAD_OF_FACETS = F,
    CLUSTER_COLORS = NULL,
    DO_TILE = F) {
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
    scale.func <- switch(
        EXPR = scale.by,
        "size" = scale_size,
        "radius" = scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
        feature.groups <- unlist(x = sapply(
            X = 1:length(features),
            FUN = function(x) {
                return(rep(x = names(x = features)[x], each = length(features[[x]])))
            }
        ))
        if (any(is.na(x = feature.groups))) {
            warning(
                "Some feature groups are unnamed.",
                call. = FALSE,
                immediate. = TRUE
            )
        }
        features <- unlist(x = features)
        names(x = feature.groups) <- features
    }

    idents.groups <- NULL
    if (is.list(idents.grouping) | any(!is.na(names(idents.grouping)))) {
        idents.groups <- unlist(x = sapply(
            X = 1:length(idents.grouping),
            FUN = function(x) {
                return(rep(x = names(x = idents.grouping)[x], each = length(idents.grouping[[x]])))
            }
        ))
        if (any(is.na(x = idents.grouping))) {
            warning(
                "Some feature groups are unnamed.",
                call. = FALSE,
                immediate. = TRUE
            )
        }
        idents.grouping <- unlist(x = idents.grouping)
        names(x = idents.groups) <- idents.grouping
    }

    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

    data.features <- FetchData(object = object, vars = features, cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)[cells, drop = TRUE]
    } else {
        object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }

    id.levels <- levels(x = data.features$id)
    if (is.list(idents.grouping) | any(!is.na(names(idents.grouping)))) {
        id.levels <- unname(idents.grouping)
    }
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
        if (split.colors) {
            if (length(x = unique(x = splits)) > length(x = cols)) {
                stop("Not enough colors for the number of groups")
            }
            cols <- cols[1:length(x = unique(x = splits))]
            names(x = cols) <- unique(x = splits)
        }
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(
        X = unique(x = data.features$id),
        FUN = function(ident) {
            data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
            avg.exp <- apply(
                X = data.use,
                MARGIN = 2,
                FUN = function(x) {
                    return(mean(x = expm1(x = x)))
                }
            )
            pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
            return(list(avg.exp = avg.exp, pct.exp = pct.exp))
        }
    )
    names(x = data.plot) <- unique(x = data.features$id)
    if (cluster.idents) {
        mat <- do.call(
            what = rbind,
            args = lapply(X = data.plot, FUN = unlist)
        )
        mat <- scale(x = mat)
        id.levels <- id.levels[hclust(d = dist(x = mat))$order]
    }
    data.plot <- lapply(
        X = names(x = data.plot),
        FUN = function(x) {
            data.use <- as.data.frame(x = data.plot[[x]])
            data.use$features.plot <- rownames(x = data.use)
            data.use$id <- x
            return(data.use)
        }
    )
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    ngroup <- length(x = levels(x = data.plot$id))
    if (ngroup == 1) {
        scale <- FALSE
        warning(
            "Only one identity present, the expression values will be not scaled",
            call. = FALSE,
            immediate. = TRUE
        )
    } else if (ngroup < 5 & scale) {
        warning(
            "Scaling data with a low number of groups may produce misleading results",
            call. = FALSE,
            immediate. = TRUE
        )
    }
    avg.exp.scaled <- sapply(
        X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min, max = col.max)
            } else {
                data.use <- log1p(x = data.use)
            }
            return(data.use)
        }
    )
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (split.colors) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(
        x = data.plot$features.plot,
        levels = features
    )
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (split.colors) {
        splits.use <- vapply(
            X = as.character(x = data.plot$id),
            FUN = gsub,
            FUN.VALUE = character(length = 1L),
            pattern = paste0(
                "^((",
                paste(sort(x = levels(x = object), decreasing = TRUE), collapse = "|"),
                ")_)"
            ),
            replacement = "",
            USE.NAMES = FALSE
        )
        data.plot$colors <- mapply(
            FUN = function(color, value) {
                return(colorRampPalette(colors = c("grey", color))(20)[value])
            },
            color = cols[splits.use],
            value = avg.exp.scaled
        )
    }
    color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
        data.plot$feature.groups <- factor(
            x = feature.groups[data.plot$features.plot],
            levels = unique(x = feature.groups)
        )
    }
    if (!is.null(x = idents.groups)) {
        data.plot$idents.groups <- factor(
            x = idents.groups[data.plot$id],
            levels = unique(x = idents.groups)
        )
    }
    if (!USE_PATCHWORK_INSTEAD_OF_FACETS) {
        plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) +
            geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
            guides(size = guide_legend(title = "Percent Expressed")) +
            labs(
                x = "Features",
                y = ifelse(test = is.null(x = split.by), yes = "Identity", no = "Split Identity")
            ) +
            theme_cowplot()
        if (!is.null(x = feature.groups) && !is.null(x = idents.groups)) {
            plot <- plot + facet_grid(
                facets = idents.groups ~ feature.groups,
                scales = "free",
                space = "free",
                switch = "y"
            ) + theme(
                panel.spacing = unit(x = 1, units = "lines"),
                strip.background = element_blank()
            )
        } else if (!is.null(x = feature.groups)) {
            plot <- plot + facet_grid(
                facets = ~feature.groups,
                scales = "free_x",
                space = "free_x",
                switch = "y"
            ) + theme(
                panel.spacing = unit(x = 1, units = "lines"),
                strip.background = element_blank()
            )
        } else if (!is.null(x = idents.groups)) {
            plot <- plot + facet_grid(
                facets = idents.groups ~ .,
                scales = "free_y",
                space = "free_y",
            ) + theme(
                panel.spacing = unit(x = 1, units = "lines"),
                strip.background = element_blank()
            )
        }
        if (split.colors) {
            plot <- plot + scale_color_identity()
        } else if (length(x = cols) == 1) {
            plot <- plot + scale_color_distiller(palette = cols)
        } else {
            plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
        }
        if (!split.colors) {
            plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
        }
    } else {
        plotList <- list()
        plotNum <- 1
        scale.min <- NA
        scale.max <- NA
        last_ident <- idents.groups[[length(idents.groups)]]
        first_ident <- idents.groups[[1]] # length(feature.groups)]]

        first_feature <- feature.groups[[1]] # length(feature.groups)]]
        last_feature <- feature.groups[[length(feature.groups)]]

        if (!is.null(CLUSTER_COLORS)) {
            feature.groups_customfacet <- c("CUSTOMCOLOR", feature.groups)
            data.plot <- data.plot |> dplyr::mutate(TILE_COLOR = CLUSTER_COLORS[id], color_x_plot = 0)
        } else {
            feature.groups_customfacet <- feature.groups
        }

        if (!DO_TILE) {
            for (y in unique(idents.groups)) {
                for (x in unique(feature.groups_customfacet)) {
                    if (x == "CUSTOMCOLOR") {
                        data.plot.tmp <- data.plot |> dplyr::filter(idents.groups == y)

                        plotList[[plotNum]] <- ggplot(data = data.plot.tmp, mapping = aes_string(x = "color_x_plot", y = "id")) +
                            geom_raster(fill = "white", mapping = aes_string(x = "color_x_plot", y = "id")) +
                            point_with_family(geom_point(shape = "▶", mapping = aes(color_x_plot, y = id, color = id), size = 10, show.legend = F), "DejaVu Sans") +
                            scale_color_manual(values = CLUSTER_COLORS) +
                            # geom_text(label = "▶", size = 10, color = "black", mapping = aes_string(x = "color_x_plot", y = "id")) +
                            # geom_text(label = "▶", size = 9, color = data.plot.tmp$TILE_COLOR, mapping = aes_string(x = "color_x_plot", y = "id"), vjust=0, hjust=0) +
                            theme_cowplot() +
                            # scale_x_discrete(expand = c(0, 0)) +
                            # scale_y_discrete(expand = c(0, 0)) +
                            theme(panel.grid = element_blank(), panel.border = element_blank()) +
                            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                            theme(plot.margin = margin(t = 1, r = 0, b = 1, l = 1)) +
                            theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
                            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
                            coord_fixed()

                        plotNum <- plotNum + 1
                    } else {
                        data.plot.tmp <- data.plot |> dplyr::filter(feature.groups == x & idents.groups == y)
                        plotList[[plotNum]] <- ggplot(data = data.plot.tmp, mapping = aes_string(x = "features.plot", y = "id")) +
                            geom_point(mapping = aes_string(size = "pct.exp", fill = color.by), shape = 21, colour = "grey", stroke = 0.5) +
                            theme_cowplot() +
                            theme(panel.grid = element_line(
                                color = "grey",
                                size = 0.25,
                                linetype = "solid"
                            ), panel.border = element_rect(color = "grey", fill = NA, size = 0.25)) +
                            theme(
                                axis.text.x = element_text(
                                    angle = 90,
                                    vjust = 0.5, hjust = 1
                                )
                            ) +
                            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                            theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)) +
                            scale.func(name = "Percent Expression", labels = c("0%", "25%", "50%", "75%", "100%"), breaks = c(0, 25, 50, 75, 100), range = c(0, dot.scale), limits = c(0, 100)) +
                            scale_fill_gradient2(
                                low =
                                    "#5CACDB",
                                mid = "white",
                                high = "#EA7FA3",
                                limits = range(c(-2, 2)),
                                guide = guide_colorbar(title = "Average Expression")
                                # limits = range(c(min(data.plot[, color.by]), max(data.plot[, color.by])))
                            ) +
                            guides(size = guide_legend(title = "Percent Expressed")) +
                            coord_fixed()
                        if (y != last_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
                        }

                        if (x != first_feature) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
                        }
                        if (!is.null(CLUSTER_COLORS)) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
                            if (x == first_feature) {
                                plotList[[plotNum]] <- plotList[[plotNum]] + theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 0))
                            }
                        }

                        if (x == last_feature && y == first_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(idents.groups ~ feature.groups)
                        } else if (x == last_feature) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(idents.groups ~ .)
                        } else if (y == first_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(. ~ feature.groups)
                        }
                        plotList[[plotNum]] <- plotList[[plotNum]] + theme(strip.clip = "off")
                        plotNum <- plotNum + 1
                    }
                }
            }
            plot <- patchwork::wrap_plots(plotList, ncol = length(unique(feature.groups_customfacet)), byrow = T, guides = "collect", axes = "collect", axes_titles = "collect")
        }


        if (DO_TILE) {
            for (y in unique(idents.groups)) {
                for (x in unique(feature.groups_customfacet)) {
                    if (x == "CUSTOMCOLOR") {
                        data.plot.tmp <- data.plot |> dplyr::filter(idents.groups == y)

                        plotList[[plotNum]] <- ggplot(data = data.plot.tmp, mapping = aes_string(x = "color_x_plot", y = "id")) +
                            geom_raster(fill = "white", mapping = aes_string(x = "color_x_plot", y = "id")) +
                            point_with_family(geom_point(shape = "▶", mapping = aes(color_x_plot, y = id, color = id), size = 10, show.legend = F), "DejaVu Sans") +
                            scale_color_manual(values = CLUSTER_COLORS) +
                            # geom_text(label = "▶", size = 10, color = "black", mapping = aes_string(x = "color_x_plot", y = "id")) +
                            # geom_text(label = "▶", size = 9, color = data.plot.tmp$TILE_COLOR, mapping = aes_string(x = "color_x_plot", y = "id"), vjust=0, hjust=0) +
                            theme_cowplot() +
                            # scale_x_discrete(expand = c(0, 0)) +
                            # scale_y_discrete(expand = c(0, 0)) +
                            theme(panel.grid = element_blank(), panel.border = element_blank()) +
                            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                            theme(plot.margin = margin(t = 1, r = 0, b = 1, l = 1)) +
                            theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
                            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
                            coord_fixed()

                        plotNum <- plotNum + 1
                    } else {
                        data.plot.tmp <- data.plot |> dplyr::filter(feature.groups == x & idents.groups == y)
                        plotList[[plotNum]] <- ggplot(data = data.plot.tmp, mapping = aes_string(x = "features.plot", y = "id")) +
                            geom_tile(mapping = aes_string(fill = color.by), colour = "grey") +
                            theme_cowplot() +
                            theme(panel.grid = element_line(
                                color = "grey",
                                size = 0.25,
                                linetype = "solid"
                            ), panel.border = element_rect(color = "grey", fill = NA, size = 0.25)) +
                            theme(
                                axis.text.x = element_text(
                                    angle = 90,
                                    vjust = 0.5, hjust = 1
                                )
                            ) +
                            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                            theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)) +
                            scale_fill_gradient2(
                                low =
                                    "#5CACDB",
                                mid = "white",
                                high = "#EA7FA3",
                                limits = range(c(-2, 2)),
                                guide = guide_colorbar(title = "Average Expression")
                                # limits = range(c(min(data.plot[, color.by]), max(data.plot[, color.by])))
                            ) +
                            coord_fixed()
                        if (y != last_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
                        }

                        if (x != first_feature) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
                        }
                        if (!is.null(CLUSTER_COLORS)) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
                            if (x == first_feature) {
                                plotList[[plotNum]] <- plotList[[plotNum]] + theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 0))
                            }
                        }

                        if (x == last_feature && y == first_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(idents.groups ~ feature.groups)
                        } else if (x == last_feature) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(idents.groups ~ .)
                        } else if (y == first_ident) {
                            plotList[[plotNum]] <- plotList[[plotNum]] + facet_grid(. ~ feature.groups)
                        }
                        plotList[[plotNum]] <- plotList[[plotNum]] + theme(strip.clip = "off")
                        plotNum <- plotNum + 1
                    }
                }
            }
            plot <- patchwork::wrap_plots(plotList, ncol = length(unique(feature.groups_customfacet)), byrow = T, guides = "collect", axes = "collect", axes_titles = "collect")
        }
        # pdf_and_png(test, "output", "test.png", pdfWidth = 24, pdfHeight = 12)
        # pdf_and_png(test, "output", "testb.png", pdfWidth = 18, pdfHeight = 8)

        # test2 <- test / test
        # pdf_and_png(test2, "output", "vtest", pdfWidth = 18, pdfHeight = 8)
    }
    return(plot)
}
