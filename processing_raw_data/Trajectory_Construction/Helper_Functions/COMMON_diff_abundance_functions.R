library(dplyr)
library(ggplot2)
library(effectsize)

# Implementation of various methods to perform differential abundance testing

# Computes the wilcoxon rank sum test of all clusters between two groups. Returns a log2FC and a ranked-biserial correlation coefficient
clusterWilcoxonTest <- function(counts, clusters_filter = NULL, fitted_variable = "progression_group", groupA = "NP", groupB = "RP") {
    # Variable mode - tbd

    p_sig_thresh <- 0.05
    padj_sig_thresh <- 0.1

    if (!is.null(clusters_filter)) {
        counts <- counts |> dplyr::filter(cluster %in% clusters_filter)
    }

    counts <- counts |> dplyr::filter(!!sym(fitted_variable) %in% c(groupA, groupB))
    cluster_size <- counts |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(cluster_size = sum(Counts))

    empty_clusters <- cluster_size |> dplyr::filter(cluster_size == 0)

    counts <- dplyr::left_join(counts, cluster_size)

    # remove empty clusters
    counts <- counts |>
        dplyr::filter(cluster_size != 0) |>
        droplevels()

    counts <- counts |>
        dplyr::group_by(public_id) |>
        dplyr::summarise(pat_total = sum(Counts)) |>
        dplyr::right_join(counts) |>
        dplyr::filter(pat_total != 0) |> # remove patients with zero cells in group
        dplyr::mutate(prop = Counts / pat_total)

    counts[[fitted_variable]] <- factor(counts[[fitted_variable]], levels = c(groupA, groupB))
    df_p_val <- counts |>
        rstatix::group_by(cluster) %>%
        rstatix::wilcox_test(as.formula(sprintf("prop ~ %s", fitted_variable)), detailed = T)

    counts_GA <- counts |> dplyr::filter(!!sym(fitted_variable) == groupA)
    counts_GB <- counts |> dplyr::filter(!!sym(fitted_variable) == groupB)

    rank_bis_list <- list()
    for (i in unique(df_p_val$cluster)) {
        tmp_GA <- counts_GA |> dplyr::filter(cluster == i)
        tmp_GB <- counts_GB |> dplyr::filter(cluster == i)
        rank_bis_list[[i]] <- effectsize::rank_biserial(tmp_GA$prop, tmp_GB$prop)$r_rank_biserial
    }

    rank_bis_df <- data.frame(cluster = names(rank_bis_list), rank_biserial = unlist(rank_bis_list))

    df_p_val <- dplyr::left_join(df_p_val, rank_bis_df, by = "cluster")

    means <- counts |>
        dplyr::group_by(cluster, !!sym(fitted_variable)) |>
        dplyr::summarise(mean = mean(prop))

    means_A <- means |> dplyr::filter(!!sym(fitted_variable) == groupA)
    means_B <- means |> dplyr::filter(!!sym(fitted_variable) == groupB)

    means_A <- means_A |> dplyr::select(cluster, mean)
    means_B <- means_B |> dplyr::select(cluster, mean)

    colnames(means_A)[2] <- "mean_A"
    colnames(means_B)[2] <- "mean_B"

    df_p_val <- dplyr::left_join(df_p_val, means_A) |> dplyr::left_join(means_B)

    df_p_val$fold_change <- (df_p_val$mean_A) / (df_p_val$mean_B)
    df_p_val$fold_change_for_display <- (df_p_val$mean_A + 1e-10) / (df_p_val$mean_B + 1e-10)

    df_p_val$log2_effect_size <- log2(df_p_val$fold_change_for_display)

    df_p_val$padj <- p.adjust(df_p_val$p, method = "BH")

    df_p_val$xlabel <- "Log2 Mean Fold Change"
    df_p_val$ylabel <- "P Value, Wilcoxon"
    df_p_val$fitted_variable <- fitted_variable
    df_p_val$groupA <- groupA
    df_p_val$groupB <- groupB
    df_p_val$covariates <- "None"
    df_p_val$adjustment <- "BH"
    df_p_val$p_sig_thresh <- p_sig_thresh
    df_p_val$padj_sig_thresh <- padj_sig_thresh

    df_p_val$sign <- sign(df_p_val$rank_biserial)
    df_p_val <- df_p_val |> dplyr::mutate(direction = dplyr::case_when(
        sign < 0 & p < p_sig_thresh ~ sprintf("%s (SIG)", groupB),
        sign >= 0 & p < p_sig_thresh ~ sprintf("%s (SIG)", groupA),
        sign < 0 ~ sprintf("%s", groupB),
        sign > 0 ~ sprintf("%s", groupA),
        sign == 0 ~ "Equal Mean",
        .default = "UNHANDELED CASE"
    ))

    df_p_val <- df_p_val |> dplyr::mutate(direction_adj = dplyr::case_when(
        sign < 0 & padj < padj_sig_thresh ~ sprintf("%s (SIG)", groupB),
        sign >= 0 & padj < padj_sig_thresh ~ sprintf("%s (SIG)", groupA),
        sign < 0 ~ sprintf("%s", groupB),
        sign > 0 ~ sprintf("%s", groupA),
        sign == 0 ~ "Equal Mean",
        .default = "UNHANDELED CASE"
    ))

    df_p_val <- dplyr::mutate(df_p_val, signed_pval = sign * p)
    # if (nrow(empty_clusters) > 0) {
    #     # empty_cluster_df <- data.frame(cluster = empty_clusters$cluster, .y. = "prop", group1 = "EMPTY", group2 = "EMPTY", n1 = 0, n2 = 0, statistic = 0, p = 1, mean_A = 0, mean_B = 0, fold_change = 0, fold_change_for_display = 0, sign = 0, direction = "No Cells", signed_pval = 1)

    #     df_p_val <- rbind(df_p_val, empty_cluster_df)
    # }
    return(list("df_p_val" = df_p_val, "counts_with_props" = counts))
}

# (Not used) Computes differential abundance using a beta regression model with fixed precision, correcting for siteXbatch.

betaReg_AllClusters <- function(counts, md, clusterlist, groupVar) {
    prop <- counts / rowSums(counts)
    empty_clusters <- colSums(counts) == 0
    prop <- prop |> as.data.frame.matrix()
    prop <- prop + 0.000001
    prop$public_id <- rownames(prop)

    patRatios <- dplyr::left_join(md, prop)
    breg_results <- c()
    breg_mdl <- list()

    for (i in clusterlist) {
        if (empty_clusters[[i]]) {
            breg_mdl[[i]] <- NULL
            breg_results[[i]] <- 0
        } else {
            breg_mdl[[i]] <- betareg(as.formula(sprintf("%s ~ %s + siteXbatch | 1", i, groupVar)), data = patRatios)
            breg_results[i] <- summary(breg_mdl[[i]])$coef$mean[2, 1]
        }
    }
    return(list("models" = breg_mdl, "effect_sizes" = breg_results))
}

library(betareg)
# (Not used) Computes differential abundance using a beta regression model. Precision parameter can be customized

clusterBetaregTest <- function(counts, clusters_filter = NULL, fitted_variable = "progression_group", groupA = "NP", groupB = "RP", precision_model = "1", contrast_codes = NULL, factor_levels = NULL) {
    # Variable mode - tbd

    p_sig_thresh <- 0.05
    padj_sig_thresh <- 0.1

    if (!is.null(clusters_filter)) {
        counts <- counts |> dplyr::filter(cluster %in% clusters_filter)
    }

    if (is.null(contrast_codes)) {
        counts <- counts |> dplyr::filter(!!sym(fitted_variable) %in% c(groupA, groupB))
    }

    if (!is.null(factor_levels)) {
        counts <- counts |> dplyr::filter(!!sym(fitted_variable) %in% factor_levels)
    }
    cluster_size <- counts |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(cluster_size = sum(Counts))


    empty_clusters <- cluster_size |> dplyr::filter(cluster_size == 0)

    counts <- dplyr::left_join(counts, cluster_size)

    # remove empty clusters
    counts <- counts |>
        dplyr::filter(cluster_size != 0) |>
        droplevels()

    counts <- counts |>
        dplyr::group_by(public_id) |>
        dplyr::summarise(pat_total = sum(Counts)) |>
        dplyr::right_join(counts) |>
        dplyr::filter(pat_total != 0) |> # remove patients with zero cells in group
        dplyr::mutate(prop = Counts / pat_total)
    # browser()
    breg_mdl <- list()
    breg_effect_size <- c()
    breg_pval <- c()

    for (clust in unique(counts$cluster)) {
        tmp <- counts |> dplyr::filter(cluster == clust)
        tmp <- tmp |> dplyr::mutate(trnsf = ((prop * (pat_total - 1) + 0.5) / (pat_total)))

        if (!is.null(contrast_codes)) {
            tmp[, fitted_variable] <- factor(tmp[, fitted_variable], levels = factor_levels)
            colnames(codes) <- c("NP-RP", "P-RP", "NP-P")
            contrasts(tmp[, fitted_variable]) <- codes
        }
        breg_mdl[[clust]] <- betareg(as.formula(sprintf("%s ~ %s + siteXbatch | %s", "trnsf", fitted_variable, precision_model)), data = tmp)
        breg_effect_size[clust] <- summary(breg_mdl[[clust]])$coef$mean[2, 1]
        breg_pval[clust] <- summary(breg_mdl[[clust]])$coef$mean[2, 4]
    }

    n_group_a <- sum(counts[, fitted_variable] == groupA)
    n_group_b <- sum(counts[, fitted_variable] == groupB)


    df_p_val <- data.frame("cluster" = names(breg_pval), ".y." = "trnsf", group1 = groupA, group2 = groupB, n1 = n_group_a, n2 = n_group_b, p = breg_pval, log2_effect_size = breg_effect_size)

    # df_p_val <- counts |>
    #     rstatix::group_by(cluster) %>%
    #     rstatix::wilcox_test(as.formula(sprintf("prop ~ %s", fitted_variable)))

    means <- counts |>
        dplyr::group_by(cluster, !!sym(fitted_variable)) |>
        dplyr::summarise(mean = mean(prop))

    means_A <- means |> dplyr::filter(!!sym(fitted_variable) == groupA)
    means_B <- means |> dplyr::filter(!!sym(fitted_variable) == groupB)

    means_A <- means_A |> dplyr::select(cluster, mean)
    means_B <- means_B |> dplyr::select(cluster, mean)

    colnames(means_A)[2] <- "mean_A"
    colnames(means_B)[2] <- "mean_B"

    df_p_val <- dplyr::left_join(df_p_val, means_A) |> dplyr::left_join(means_B)

    # df_p_val$fold_change <- (df_p_val$mean_A) / (df_p_val$mean_B)
    # df_p_val$fold_change_for_display <- (df_p_val$mean_A + 1e-10) / (df_p_val$mean_B + 1e-10)
    df_p_val$fold_change <- exp(df_p_val$log2_effect_size)
    df_p_val$fold_change_for_display <- exp(df_p_val$log2_effect_size)
    # df_p_val$log2_effect_size <- log2(df_p_val$fold_change_for_display)

    df_p_val$padj <- p.adjust(df_p_val$p, method = "BH")

    df_p_val$xlabel <- "Log Effect Size"
    df_p_val$ylabel <- "P Value, BetaReg"
    df_p_val$fitted_variable <- fitted_variable
    df_p_val$groupA <- groupA
    df_p_val$groupB <- groupB
    df_p_val$covariates <- "None"
    df_p_val$adjustment <- "BH"
    df_p_val$p_sig_thresh <- p_sig_thresh
    df_p_val$padj_sig_thresh <- padj_sig_thresh

    df_p_val$sign <- sign(log10(df_p_val$fold_change_for_display))
    df_p_val <- df_p_val |> dplyr::mutate(direction = dplyr::case_when(
        sign < 0 & p < p_sig_thresh ~ sprintf("%s (SIG)", groupB),
        sign >= 0 & p < p_sig_thresh ~ sprintf("%s (SIG)", groupA),
        sign < 0 ~ sprintf("%s", groupB),
        sign > 0 ~ sprintf("%s", groupA),
        sign == 0 ~ "Equal Mean",
        .default = "UNHANDELED CASE"
    ))

    df_p_val <- df_p_val |> dplyr::mutate(direction_adj = dplyr::case_when(
        sign < 0 & padj < padj_sig_thresh ~ sprintf("%s (SIG)", groupB),
        sign >= 0 & padj < padj_sig_thresh ~ sprintf("%s (SIG)", groupA),
        sign < 0 ~ sprintf("%s", groupB),
        sign > 0 ~ sprintf("%s", groupA),
        sign == 0 ~ "Equal Mean",
        .default = "UNHANDELED CASE"
    ))

    df_p_val <- dplyr::mutate(df_p_val, signed_pval = sign * p)
    # if (nrow(empty_clusters) > 0) {
    #     # empty_cluster_df <- data.frame(cluster = empty_clusters$cluster, .y. = "prop", group1 = "EMPTY", group2 = "EMPTY", n1 = 0, n2 = 0, statistic = 0, p = 1, mean_A = 0, mean_B = 0, fold_change = 0, fold_change_for_display = 0, sign = 0, direction = "No Cells", signed_pval = 1)

    #     df_p_val <- rbind(df_p_val, empty_cluster_df)
    # }
    return(list("df_p_val" = df_p_val, "counts_with_props" = counts))
}

# Simplifies the output of the wilcoxon model code for output to a .csv
simplify_output_wilcoxon <- function(results) {
    return(results |> dplyr::select(cluster, fold_change, p, direction))
}


# Plotting a scatterplot of wilcoxon p value and effect size
plotWilcoxonResults <- function(df_p_val) {
    # df_p_val <- df_p_val |> dplyr::mutate(fold_change = dplyr::case_when(
    #     fold_change > 2^2.5 ~ 2^2.5,
    #     fold_change < 2^-2.5 ~ 2^-2.5,
    #     .default = fold_change
    # ))

    p <- ggplot(df_p_val) +
        geom_point(aes(x = rank_biserial, y = -log10(p + 1e-200), fill = direction, color = direction, alpha = 0.5), show.legend = F) +
        geom_label_repel(aes(x = rank_biserial, y = -log10(p + 1e-200), label = cluster, fill = direction), max.overlap = 15) +
        geom_vline(xintercept = 0, color = "black", size = 1) +
        geom_hline(yintercept = -log10(0.05 + 1e-200), color = "black", size = 1) +
        xlim(-2.5, 2.5) +
        ggtitle("Diff Abundance - Wilcoxon") +
        theme_prism()

    return(p)
}

# Plots a box plot for a cluster

plotClusters <- function(counts, grouping_var = NULL, groupA = NULL, groupB = NULL, groups = NULL, scale_cluster_plots = TRUE, ylabel = "Average Cluster Proportion", xlabel = "") {
    if (!is.null(grouping_var)) {
        if (xlabel == "") {
            xlabel <- grouping_var
        }
        if (!is.null(groupA) && !is.null(groupB)) {
            groups <- c(groupA, groupB)
        }
        if (!is.null(groups)) {
            counts[[grouping_var]] <- factor(counts[[grouping_var]], levels = c(groupA, groupB))
        } else {
            groups <- unique(counts[[grouping_var]])
        }
    }

    cluster_size <- counts |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(cluster_size = sum(Counts))

    counts <- dplyr::left_join(counts, cluster_size)


    counts <- counts |>
        dplyr::group_by(public_id) |>
        dplyr::summarise(pat_total = sum(Counts)) |>
        dplyr::right_join(counts) |>
        dplyr::filter(pat_total != 0) |> # remove patients with zero cells in group
        dplyr::mutate(prop = Counts / pat_total)

    output_lists <- list()
    p_lists <- list()
    for (clustID in unique(counts$cluster)) {
        tmp <- counts |> dplyr::filter(cluster == clustID)

        prop_quantiles <- quantile(tmp$prop, c(0, 0.05, 0.25, 0.75, 0.95, 1))

        iqr <- prop_quantiles[["75%"]] - prop_quantiles[["25%"]]
        full_range <- prop_quantiles[["100%"]] - prop_quantiles[["0%"]]

        iqr_ratio <- iqr / full_range

        if (iqr_ratio < 0.25) {
            y_limits <- c(max(0, prop_quantiles[["25%"]] - iqr * 1.5), min(1, prop_quantiles[["75%"]] + 1.5 * iqr))
        } else {
            y_limits <- c(prop_quantiles[["100%"]], prop_quantiles[["0%"]])
        }

        if (!is.null(grouping_var)) {
            p <- ggplot(tmp, aes(x = !!sym(grouping_var), fill = !!sym(grouping_var), color = !!sym(grouping_var), group = !!sym(grouping_var), y = prop))
        } else {
            p <- ggplot(tmp, aes(x = 0, y = prop))
        }

        p <- p +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter() +
            ylab(ylabel) +
            xlab(xlabel) +
            theme_prism() +
            scale_fill_prism("floral") +
            scale_color_prism("floral") +
            ylim(y_limits) +
            ggtitle(clustID)

        p_lists[[clustID]] <- p
    }

    return(list("Plots" = p_lists, "Outputs" = output_lists))
}

# output_lists[[clustID]] <- pdf_and_png(p, output_subdir, paste0(nameMod, "_", clustID), pdfWidth = 6, pdfHeight = 3.5)
# Perform all wilcoxon tests and generate plots

quickWilcoxAndPlots <- function(counts, grouping_var, groupA, groupB) {
    wilcox_out_list <- clusterWilcoxonTest(counts, fitted_variable = grouping_var, groupA = groupA, groupB = groupB)
    wilcox_out <- wilcox_out_list[["df_p_val"]]
    counts_with_prop <- wilcox_out_list[["counts_with_props"]]

    plotWilcoxOut <- plotWilcoxonResults(wilcox_out)
    plotClusterOut <- plotClusters(counts, grouping_var = grouping_var, groupA = groupA, groupB = groupB)
    return(list("wilcox_out" = wilcox_out, "plotWilcoxOut" = plotWilcoxOut, "plotClusterOut" = plotClusterOut, "counts_with_props" = counts_with_prop))
}

# Perform all beta regression tests and generate plots

quickBetaAndPlots <- function(counts, grouping_var, groupA, groupB, precision_model = "1", contrast_codes = NULL, factor_levels = NULL) {
    wilcox_out_list <- clusterBetaregTest(counts, fitted_variable = grouping_var, groupA = groupA, groupB = groupB, precision_model = precision_model, contrast_codes = contrast_codes, factor_levels = factor_levels)
    wilcox_out <- wilcox_out_list[["df_p_val"]]
    counts_with_prop <- wilcox_out_list[["counts_with_props"]]

    plotWilcoxOut <- plotWilcoxonResults(wilcox_out)
    plotClusterOut <- plotClusters(counts, grouping_var = grouping_var, groupA = groupA, groupB = groupB)
    return(list("wilcox_out" = wilcox_out, "plotWilcoxOut" = plotWilcoxOut, "plotClusterOut" = plotClusterOut, "counts_with_props" = counts_with_prop))
}


library(DirichletReg)
# Helper functions to grab the relevant coefficients for certain models
HelperFun <- function(fit) {
    u <- summary(fit)
    pvals <- u$coef.mat[grep("Condition", rownames(u$coef.mat), invert = FALSE), 4]
    v <- names(pvals)

    pvals <- matrix(pvals, ncol = length(u$varnames))
    rownames(pvals) <- gsub("Condition", "", v[1:nrow(pvals)])
    colnames(pvals) <- u$varnames

    return(pvals)
}

# From Yered
HelperFun2 <- function(fit) {
    u <- summary(fit)
    # p-values
    pvals <- u$coef.mat[grep("progression_groupRP", rownames(u$coef.mat), invert = FALSE), 4]
    v <- names(pvals)
    # test statistic
    zstat <- u$coef.mat[grep("progression_groupRP", rownames(u$coef.mat), invert = FALSE), 3]

    res <- data.frame(
        cluster = u$varnames,
        z.ratio = zstat,
        p.value = pvals
    )

    rownames(res) <- u$varnames
    return(
        res
    )
}

HelperFun2_risk <- function(fit) {
    u <- summary(fit)
    # p-values
    pvals <- u$coef.mat[grep("davies_based_risk", rownames(u$coef.mat), invert = FALSE), 4]
    v <- names(pvals)
    # test statistic
    zstat <- u$coef.mat[grep("davies_based_risk", rownames(u$coef.mat), invert = FALSE), 3]

    res <- data.frame(
        cluster = u$varnames,
        z.ratio = zstat,
        p.value = pvals
    )

    rownames(res) <- u$varnames
    return(
        res
    )
}

library(tidyr)

# Reformatting function to prepare for dirichlet regression
reformat_df_to_table <- function(counts) {
    cluster_rownames <- unique(counts$cluster)
    counts <- counts |>
        # dplyr::select(-isDoublet, -isLQ, -isBcell, -isImmune, -isNonMalignantImmune, -nkt_subtype, -subcluster_V03072023_compartment) |> # Variables that depend on cluster which I added
        dplyr::group_by(public_id) |>
        tidyr::spread(cluster, Counts)

    return(list("table_df" = counts, "cluster_rownames" = cluster_rownames))
}

library(mvtnorm)

# return the mean and quantiles of a vector
mean.quant <- function(x, probs) {
    a <- mean(x, na.rm = T)
    b <- t(quantile(x, probs = probs, na.rm = T))
    d <- c(a, b)
    names(d) <- c("mean", "lower", "upper")
    return(d)
}

# Giving a fit Dirichlet model corrected for a set of batch covariates, construct an estimated 95% CI of the fitted Dirichlet mean by sampling from the mean and covariance of the various parameters.
# Can be used even in cases where a given covariate shows large differences in precision.
# (p values - not use. Fold Changes used for Figure 5A)
dirichlet_pvalue_simulation <- function(
    dr_fit_common,
    reformated_counts,
    cluster_colnames,
    fitted_variable = "progression_group",
    groupA = "NP",
    groupB = "RP",
    sim_per_siteXbatch = 50000,
    OVERRIDE_RETURN_MEAN_ONLY = F) {
    set.seed(42)

    # Group B is the reference variable

    # Extract model mean and covariance
    means <- unlist(coef(dr_fit_common))
    vc <- vcov(dr_fit_common)

    # Get the number of siteXbatch variables in the model
    siteXbatch_ref <- sort(unique(reformated_counts$siteXbatch))
    Study_Site_ref <- sort(unique(reformated_counts$Study_Site))
    Batch_ref <- sort(unique(reformated_counts$Batch))

    siteXbatch_all <- siteXbatch_ref
    siteXbatch_ref <- siteXbatch_ref[2:length(siteXbatch_ref)]

    siteXbatch_replicated <- c()
    for (i in siteXbatch_all) {
        siteXbatch_replicated <- c(siteXbatch_replicated, rep(i, sim_per_siteXbatch))
    }

    # add one for the base case
    numSiteXBatch <- length(siteXbatch_ref) + 1

    # Number of simulated model covariates to generate.
    # Draw parameters for the dirichlet model from a multivariate normal distribution based on the mean and variance of the coefficients.

    # FALSE: Just generate N independent models. TRUE: Reuse fitted models across siteXbatches
    if (!OVERRIDE_RETURN_MEAN_ONLY) {
        N <- sim_per_siteXbatch #* (length(siteXbatch_ref) + 1)
        rnd <- rmvnorm(N, mean = means, sigma = vc)
    } else {
        # browser()
        N <- 10
        rnd <- rmvnorm(N, mean = means, sigma = matrix(0, nrow = nrow(vc), ncol = ncol(vc))) # Not clean but allows for code reuse. All rows should be the same.
    }
    modobj <- dr_fit_common # calls best model fit object `modobj`
    # generate linear predictors for each randomization at each of 39 original observations
    dm <- do.call(cbind, modobj$X) # is the design matrix, don't need the last column


    # Grab the fitted variable names
    intercept_col <- grep("Intercept", colnames(dm))
    siteXbatch_all_col <- grep("siteXbatch", colnames(dm))
    Batch_all_col <- grep("Batch", colnames(dm))
    Study_Site_all_col <- grep("Study_Site", colnames(dm))

    progression_group_call_a <- grep(paste0(fitted_variable, groupA), colnames(dm))
    progression_group_call_b <- grep(paste0(fitted_variable, groupB), colnames(dm))

    # Get the variable numbers for the siteXbatch variables.
    siteXbatch_calls <- list()
    for (i in siteXbatch_ref) {
        siteXbatch_calls[[i]] <- grep(i, colnames(dm))
    }

    Batch_calls <- list()
    for (i in Batch_ref) {
        Batch_calls[[i]] <- grep(i, colnames(dm))
    }

    Study_Site_calls <- list()
    for (i in Study_Site_ref) {
        Study_Site_calls[[i]] <- grep(i, colnames(dm))
    }

    covariate_variable <- ""
    if (length(siteXbatch_all_col) != 0) {
        covariate_variable <- "siteXbatch"
        covariate_columns <- siteXbatch_all
        num_covariate_columns <- numSiteXBatch

        # generate a 12 x nVariables design matrix.
        dm_custom <- array(0, dim = c(1 + length(siteXbatch_ref), ncol(dm)))

        # Set intercept columns to 1
        dm_custom[, intercept_col] <- 1

        # Each of the 12 rows corresponds to a different siteXbatch. Base case: All 0. Other cases: set corresponding siteXbatch variable to one.
        for (i in 1:(length(siteXbatch_ref) + 1)) {
            dm_custom[i, siteXbatch_all_col] <- 0.0
            if (i > 1) {
                select_cols <- siteXbatch_calls[[i - 1]]
                dm_custom[i, select_cols] <- 1.0
            }
        }

        # Copy the design matrix - one with progression_group_call == true, one with progression_group_call == false
        dm_custom_np <- dm_custom
        dm_custom_rp <- dm_custom

        # if one is the reference case, it is already all zeros.
        if (length(progression_group_call_a) > 0) {
            dm_custom_np[, progression_group_call_a] <- 1.0 # set groupA to 1.0 in groupA
            dm_custom_rp[, progression_group_call_a] <- 0.0 # set groupA to 0.0 in groupB
        }

        if (length(progression_group_call_b) > 0) {
            dm_custom_np[, progression_group_call_b] <- 0.0 # set groupB to 0.0 in groupA
            dm_custom_rp[, progression_group_call_b] <- 1.0 # set groupB to 1.0 in groupB
        }
    } else if (length(Study_Site_all_col) != 0 || length(Batch_all_col) != 0) {
        if (length(Study_Site_all_col) != 0 && length(Batch_all_col) != 0) {
            covariate_variable <- "siteXbatch"

            covariate_columns <- siteXbatch_all
            num_covariate_columns <- numSiteXBatch


            # Create design matrix that gives sxb-styled results

            # generate a 12 x nVariables design matrix.
            dm_custom <- array(0, dim = c(num_covariate_columns, ncol(dm)))

            # Set intercept columns to 1
            dm_custom[, intercept_col] <- 1

            # Each of the 12 rows corresponds to a different siteXbatch. Base case: All 0. Other cases: set corresponding siteXbatch variable to one.
            # splite siteX batch into site and batch number
            # get the cols for each
            for (i in 1:(length(siteXbatch_ref) + 1)) {
                dm_custom[i, Study_Site_all_col] <- 0.0
                dm_custom[i, Batch_all_col] <- 0.0

                siteXbatch_pseudo <- siteXbatch_all[[i]]
                SITE <- strsplit(siteXbatch_pseudo, "_")[[1]][[1]]
                BATCH <- strsplit(siteXbatch_pseudo, "_")[[1]][[2]] |> substr(2, 2)
                BATCH <- paste0("Batch_", BATCH)

                select_cols_1 <- Batch_calls[[BATCH]]
                select_cols_2 <- Study_Site_calls[[SITE]]
                if (length(select_cols_1) > 0) {
                    dm_custom[i, select_cols_1] <- 1.0
                }
                if (length(select_cols_2) > 0) {
                    dm_custom[i, select_cols_2] <- 1.0
                }
            }
        } else {
            if (length(Study_Site_all_col) != 0) {
                covariate_variable <- "Study_Site"

                covariate_columns <- Study_Site_ref
                num_covariate_columns <- length(Study_Site_Ref)


                # Create design matrix that gives sxb-styled results

                # generate a 12 x nVariables design matrix.
                dm_custom <- array(0, dim = c(num_covariate_columns, ncol(dm)))

                # Set intercept columns to 1
                dm_custom[, intercept_col] <- 1

                # Each of the 12 rows corresponds to a different siteXbatch. Base case: All 0. Other cases: set corresponding siteXbatch variable to one.
                # splite siteX batch into site and batch number
                # get the cols for each
                for (i in 1:(length(covariate_columns))) {
                    dm_custom[i, siteXbatch_all_col] <- 0.0
                    select_cols_1 <- Study_Site_calls[[i]]
                    if (length(select_cols_1) > 0) {
                        dm_custom[i, select_cols_1] <- 1.0
                    }
                }
            } else {
                covariate_variable <- "Batch"

                covariate_columns <- Batch_ref
                num_covariate_columns <- length(Batch_ref)


                # Create design matrix that gives sxb-styled results

                # generate a 12 x nVariables design matrix.
                dm_custom <- array(0, dim = c(num_covariate_columns, ncol(dm)))

                # Set intercept columns to 1
                dm_custom[, intercept_col] <- 1

                # Each of the 12 rows corresponds to a different siteXbatch. Base case: All 0. Other cases: set corresponding siteXbatch variable to one.
                # splite siteX batch into site and batch number
                # get the cols for each
                for (i in 1:(length(covariate_columns))) {
                    dm_custom[i, Batch_all_col] <- 0.0
                    select_cols_1 <- Batch_calls[[i]]
                    if (length(select_cols_1) > 0) {
                        dm_custom[i, select_cols_1] <- 1.0
                    }
                }
            }
        }

        # Copy the design matrix - one with progression_group_call == true, one with progression_group_call == false
        dm_custom_np <- dm_custom
        dm_custom_rp <- dm_custom

        # if one is the reference case, it is already all zeros.
        if (length(progression_group_call_a) > 0) {
            dm_custom_np[, progression_group_call_a] <- 1.0 # set groupA to 1.0 in groupA
            dm_custom_rp[, progression_group_call_a] <- 0.0 # set groupA to 0.0 in groupB
        }

        if (length(progression_group_call_b) > 0) {
            dm_custom_np[, progression_group_call_b] <- 0.0 # set groupB to 0.0 in groupA
            dm_custom_rp[, progression_group_call_b] <- 1.0 # set groupB to 1.0 in groupB
        }
    } else {
        covariate_variable <- "NONE"
        covariate_columns <- "NONE"
        num_covariate_columns <- 1

        dm_custom <- array(0, dim = c(num_covariate_columns, ncol(dm)))

        # Set intercept columns to 1
        dm_custom[, intercept_col] <- 1

        dm_custom_np <- dm_custom
        dm_custom_rp <- dm_custom

        # if one is the reference case, it is already all zeros.
        if (length(progression_group_call_a) > 0) {
            dm_custom_np[, progression_group_call_a] <- 1.0 # set groupA to 1.0 in groupA
            dm_custom_rp[, progression_group_call_a] <- 0.0 # set groupA to 0.0 in groupB
        }

        if (length(progression_group_call_b) > 0) {
            dm_custom_np[, progression_group_call_b] <- 0.0 # set groupB to 0.0 in groupA
            dm_custom_rp[, progression_group_call_b] <- 1.0 # set groupB to 1.0 in groupB
        }
    }

    n <- nrow(reformated_counts)

    per_cluster_linear_predictors_NP <- list()
    per_cluster_linear_predictors_RP <- list()

    colnames(dm_custom_np) <- colnames(dm)
    colnames(dm_custom_rp) <- colnames(dm)

    # Compute the per-cluster alpha values for all N simulations, for NP and RP, for all design matrix permutations.
    # Log alpha:
    # Per cluster: N model samples x num_sxb covariates.
    # Rowsum to get 'averaged' linear predictor per sxb.

    # This will work for "subcluster_V03072023" as long as we don't have any cases where "Comp.A" and "Comp.A.B" are valid clusters (which shouldn't be the case)
    for (i in cluster_colnames) {
        ncols_grabbed <- length(grep(paste0("^", i, "\\."), colnames(rnd)))
        # if (ncols_grabbed > 15) {
        #     browser()
        #     stopifnot(ncols_grabbed < 16)
        # }
        per_cluster_linear_predictors_NP[[i]] <- rnd[, grep(paste0("^", i, "\\."), colnames(rnd))] %*% t(dm_custom_np)[grep(paste0("^", i, "\\."), colnames(rnd)), ]
        per_cluster_linear_predictors_RP[[i]] <- rnd[, grep(paste0("^", i, "\\."), colnames(rnd))] %*% t(dm_custom_rp)[grep(paste0("^", i, "\\."), colnames(rnd)), ]
    }

    alphas_averaged <- array(NA, dim = c(length(cluster_colnames), N), dimnames = list(cluster_colnames, c(1:N)))

    alphas_output_NP_averaged <- alphas_averaged
    alphas_output_RP_averaged <- alphas_averaged

    # browser()
    # Gives same output as above
    if (dim(per_cluster_linear_predictors_NP[[1]])[[2]] > 1) {
        for (j in cluster_colnames) {
            alphas_output_NP_averaged[j, ] <- exp(rowMeans(per_cluster_linear_predictors_NP[[j]][, ]))
            alphas_output_RP_averaged[j, ] <- exp(rowMeans(per_cluster_linear_predictors_RP[[j]][, ]))
        }
    } else {
        for (j in cluster_colnames) {
            alphas_output_NP_averaged[j, ] <- exp(per_cluster_linear_predictors_NP[[j]])
            alphas_output_RP_averaged[j, ] <- exp(per_cluster_linear_predictors_RP[[j]])
        }
    }
    # GroupA_v_B, Cluster, Model permutation
    mu_output_averaged <- array(NA, dim = c(2, length(cluster_colnames), N), dimnames = list(c("RP", "NP"), cluster_colnames, c(1:N)))

    test_denom_NP <- colSums(alphas_output_NP_averaged)
    test_denom_RP <- colSums(alphas_output_RP_averaged)

    # next time check all required operations before setting matrix dimensions
    mu_output_averaged[1, , ] <- t(t(alphas_output_NP_averaged[, ]) / test_denom_NP)
    mu_output_averaged[2, , ] <- t(t(alphas_output_RP_averaged[, ]) / test_denom_RP)

    # GroupA_v_B, Cluster, Model permutation, sxb covariate

    diff_output <- mu_output_averaged[1, , ] - mu_output_averaged[2, , ]
    diff_ratio <- mu_output_averaged[1, , ] / mu_output_averaged[2, , ]

    mean_output <- rowMeans(diff_output)
    mean_ratio <- rowMeans(diff_ratio)
    # browser()

    p_NP <- rowSums(diff_output < 0) / N
    p_RP <- rowSums(diff_output > 0) / N
    p_all_2side <- pmin(2 * p_NP, 2 * p_RP)



    IQR_avg <- apply(diff_output[, ], 1, mean.quant, probs = c(0.025, 0.975))

    if (dim(per_cluster_linear_predictors_NP[[1]])[[2]] > 1) {
        alphas_per_sxb <- array(NA, dim = c(length(cluster_colnames), N, num_covariate_columns), dimnames = list(cluster_colnames, c(1:N), covariate_columns))

        alphas_output_NP_sxb <- alphas_per_sxb
        alphas_output_RP_sxb <- alphas_per_sxb

        for (j in cluster_colnames) {
            alphas_output_NP_sxb[j, , ] <- exp(per_cluster_linear_predictors_NP[[j]][, ]) # exp(RMF_lp[i, ]) / denom
            alphas_output_RP_sxb[j, , ] <- exp(per_cluster_linear_predictors_RP[[j]][, ]) # exp(RMF_lp[i, ]) / denom
        }
        # browser()


        mu_output_perSiteXBatch <- array(NA, dim = c(2, length(cluster_colnames), N, num_covariate_columns), dimnames = list(c("RP", "NP"), cluster_colnames, c(1:N), covariate_columns))
        denom_test_NP <- apply(alphas_output_NP_sxb, MARGIN = c(2, 3), sum)
        denom_test_RP <- apply(alphas_output_RP_sxb, MARGIN = c(2, 3), sum)

        mu_output_perSiteXBatch[1, , , ] <- sweep(alphas_output_NP_sxb, MARGIN = c(2, 3), denom_test_NP, FUN = "/")
        mu_output_perSiteXBatch[2, , , ] <- sweep(alphas_output_RP_sxb, MARGIN = c(2, 3), denom_test_RP, FUN = "/")



        diff_output_perSiteXbatch <- mu_output_perSiteXBatch[1, , , ] - mu_output_perSiteXBatch[2, , , ]
        mean_output_perSitexbatch <- colMeans(diff_output_perSiteXbatch, dim = 1)
        diff_ratio_perSiteXbatch <- mu_output_perSiteXBatch[1, , , ] / mu_output_perSiteXBatch[2, , , ]

        p_NP_persxb <- colSums(aperm(diff_output_perSiteXbatch > 0, c(2, 1, 3))) / N
        p_RP_persxb <- colSums(aperm(diff_output_perSiteXbatch < 0, c(2, 1, 3))) / N

        sig_NP_count <- rowSums(p_NP_persxb < 0.025)
        sig_RP_count <- rowSums(p_RP_persxb > 0.025)

        p_all_2side_persxb <- pmin(2 * p_NP_persxb, 2 * p_RP_persxb)
        p_all_2side_persxb_count <- rowSums(p_all_2side_persxb < 0.05)

        diff_output_indp_sxb <- diff_output_perSiteXbatch
        curDim <- dim(diff_output_indp_sxb)
        N_indpsxb <- curDim[[2]] * curDim[[3]]
        dim(diff_output_indp_sxb) <- c(curDim[[1]], curDim[[2]] * curDim[[3]])
        rownames(diff_output_indp_sxb) <- rownames(diff_output)
        p_NP_indp_sxb <- rowSums(diff_output_indp_sxb > 0) / N_indpsxb
        p_RP_indp_sxb <- rowSums(diff_output_indp_sxb < 0) / N_indpsxb
        p_all_2side_indp_sxb <- pmin(2 * p_NP_indp_sxb, 2 * p_RP_indp_sxb)

        names(p_NP_indp_sxb) <- names(p_NP)
        names(p_RP_indp_sxb) <- names(p_NP)

        diff_ratio_indp_sxb <- diff_ratio_perSiteXbatch
        curDim <- dim(diff_ratio_indp_sxb)
        dim(diff_ratio_indp_sxb) <- c(curDim[[1]], curDim[[2]] * curDim[[3]])
        dimnames(diff_ratio_indp_sxb) <- list(cluster_colnames, 1:(curDim[[2]] * curDim[[3]]))

        mu_output_indp_sxb <- mu_output_perSiteXBatch
        curDim <- dim(mu_output_indp_sxb)
        N_indpsxb <- curDim[[3]] * curDim[[4]]
        dim(mu_output_indp_sxb) <- c(2, curDim[[2]], curDim[[3]] * curDim[[4]])
        dimnames(mu_output_indp_sxb) <- list(c("NP", "RP"), cluster_colnames, 1:(curDim[[3]] * curDim[[4]]))




        IQR_indp_sxb <- apply(diff_output_indp_sxb[, ], 1, mean.quant, probs = c(0.025, 0.975))
        IQR_per_sxb <- array(NA, dim = c(3, length(cluster_colnames), dim(diff_output_perSiteXbatch)[[3]]), dimnames = list(c("mean", "lower", "upper"), cluster_colnames, covariate_columns))


        for (i in 1:dim(diff_output_perSiteXbatch)[[3]]) {
            IQR_per_sxb[, , i] <- apply(diff_output_perSiteXbatch[, , i], 1, mean.quant, probs = c(0.025, 0.975))
        }

        colnames(IQR_indp_sxb) <- colnames(IQR_avg)

        outputs_persxb <- list(
            "groupA" = groupA,
            "groupB" = groupB,
            "fitted_variable" = fitted_variable,
            "covariate_variable" = covariate_variable,
            "estimation_mode" = "per_sxb",
            "diff_output" = diff_output_perSiteXbatch,
            "diff_ratio" = diff_ratio_perSiteXbatch,
            "mu_groupA" = mu_output_perSiteXBatch[1, , , ],
            "mu_groupB" = mu_output_perSiteXBatch[2, , , ],
            "p_2sided" = p_all_2side_persxb,
            "nsig_2sided" = p_all_2side_persxb_count,
            "p_groupA" = p_NP_persxb,
            "p_groupB" = p_RP_persxb,
            "nsig_groupA" = sig_NP_count,
            "nsig_groupB" = sig_RP_count,
            "IQR" = IQR_per_sxb
        )

        outputs_indp_sxb <- list(
            "groupA" = groupA,
            "groupB" = groupB,
            "fitted_variable" = fitted_variable,
            "covariate_variable" = covariate_variable,
            "estimation_mode" = "indp_sxb",
            "mu_groupA" = mu_output_indp_sxb[1, , ],
            "mu_groupB" = mu_output_indp_sxb[2, , ],
            "diff_output" = diff_output_indp_sxb,
            "diff_ratio" = diff_ratio_indp_sxb,
            "p_2sided" = p_all_2side_indp_sxb,
            "p_groupA" = p_NP_indp_sxb,
            "p_groupB" = p_RP_indp_sxb,
            "IQR" = IQR_indp_sxb
        )

        # outputs_all <- list(
        #     "groupA" = groupA,
        #     "groupB" = groupB,
        #     "fitted_variable" = fitted_variable,
        #     "covariate_variable" = covariate_variable,
        #     "estimation_mode" = "average_covariates",
        #     "diff_output" = diff_output,
        #     "diff_ratio" = diff_ratio,
        #     "mu_groupA" = mu_output_averaged[1, , ],
        #     "mu_groupB" = mu_output_averaged[2, , ],
        #     "p_2sided" = p_all_2side,
        #     "p_groupA" = p_NP,
        #     "p_groupB" = p_RP,
        #     "IQR_avg" = IQR_avg
        # )

        outputs_all <- outputs_indp_sxb
    } else {
        outputs_all <- list(
            "groupA" = groupA,
            "groupB" = groupB,
            "fitted_variable" = fitted_variable,
            "covariate_variable" = covariate_variable,
            "diff_output" = diff_output,
            "diff_ratio" = diff_ratio,
            "mu_groupA" = mu_output_averaged[1, , ],
            "mu_groupB" = mu_output_averaged[2, , ],
            "p_groupA" = p_NP,
            "p_groupB" = p_RP,
            "p_2sided" = p_all_2side,
            "IQR_avg" = IQR_avg
        )
    }


    return(outputs_all)
    # output <- array(NA, dim = c(n, length(cluster_colnames), N), dimnames = list(c(1:n), cluster_colnames, c(1:N)))
    # output.diff <- array(NA, dim = c(66, 15, N))


    # spp <- reformated_counts |> dplyr::select(public_id, progression_group)
    # # Need alpha link function



    # quant_nkt20 <- data.frame(reformated_counts, t(nkt20))

    # nkt6 <- apply(output[, "NkT.6", ], 1, mean.quant, probs = c(0.025, 0.975))

    # quant_nkt6 <- data.frame(reformated_counts, t(nkt6))

    # nkt3 <- apply(output[, "NkT.3.0", ], 1, mean.quant, probs = c(0.025, 0.975))

    # quant_nkt3 <- data.frame(reformated_counts, t(nkt3))
}

# Fits a dirichlet regression model for various conditions.
# Optionally - derive additional metrics by estimating the distribution of the fitted cluster proportions conditioned on uniform sampling of the batch covariate.
# Only used for Figure 5a
cluster_dirichlet_test <- function(counts, clusters_filter = NULL, variable_mode = "PFS", BC_mode = "None", DO_SIMULATION_STEP = F, SIMULATION_OVERRIDE_RETURN_MEAN_ONLY = T) {
    # Variable mode - tbd
    if (!is.null(clusters_filter)) {
        counts <- counts |> dplyr::filter(cluster %in% clusters_filter)
    }

    # DR_data requires a different formatting
    if (variable_mode == "PFS_NO_P" && BC_mode == "siteXbatch") {
        counts <- counts |> dplyr::filter(siteXbatch != "MAYO_B3") # Only one NP or RP sample
    }
    outputs <- reformat_df_to_table(counts)
    reformated_counts <- outputs[["table_df"]]
    cluster_colnames <- outputs[["cluster_rownames"]]
    rownames(reformated_counts) <- reformated_counts$public_id
    offset <- 0.0
    transformation_magnitude <- 10^-100
    # browser()
    if (variable_mode == "PFS" || variable_mode == "PFS_NO_P") {
        groupA <- "NP"
        groupB <- "RP"
        fitted_variable <- "progression_group"
        if (variable_mode == "PFS") {
            reformated_counts <- reformated_counts |> dplyr::filter(progression_group %in% c("NP", "RP", "P"))
            reformated_counts$progression_group <- factor(reformated_counts$progression_group, levels = c("NP", "RP", "P"))
        } else {
            reformated_counts <- reformated_counts |> dplyr::filter(progression_group %in% c("NP", "RP"))
            reformated_counts$progression_group <- factor(reformated_counts$progression_group, levels = c("NP", "RP"))
        }
        # Adding 1 for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)
        # browser()
        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ progression_group, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ progression_group + siteXbatch, # siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ progression_group + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }
        # out <- HelperFun2(dr_fit_common)
    } else if (variable_mode == "Risk") {
        groupA <- "standard_risk"
        groupB <- "high_risk"
        fitted_variable <- "davies_based_risk"
        reformated_counts <- reformated_counts |> dplyr::filter(davies_based_risk %in% c("high_risk", "standard_risk"))
        reformated_counts$davies_based_risk <- factor(reformated_counts$davies_based_risk, levels = c("standard_risk", "high_risk"))

        # Adding 1 for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)

        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ davies_based_risk, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ davies_based_risk + siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ davies_based_risk + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }

        # out <- HelperFun2_risk(dr_fit_common)
    } else if (variable_mode == "pfs_risk_SR") {
        groupA <- "NP.SR"
        groupB <- "RP.SR"
        fitted_variable <- "pfs_risk"
        if (variable_mode == "pfs_risk_SR") {
            reformated_counts <- reformated_counts |> dplyr::filter(progression_group %in% c("NP", "RP", "P"), davies_based_risk == "standard_risk")
            reformated_counts$pfs_risk <- factor(reformated_counts$pfs_risk, levels = c("NP.SR", "RP.SR", "Prog.SR"))
        }
        # Adding a small offset for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)
        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk + siteXbatch, # siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }
        # out <- HelperFun2(dr_fit_common)
    } else if (variable_mode == "pfs_risk_HR") {
        groupA <- "NP.HR"
        groupB <- "RP.HR"
        fitted_variable <- "pfs_risk"
        if (variable_mode == "pfs_risk_HR") {
            reformated_counts <- reformated_counts |> dplyr::filter(pfs_risk %in% c("NP.HR", "RP.HR", "Prog.HR", "NP.SR", "RP.SR", "Prog.SR"))
            reformated_counts$pfs_risk <- factor(reformated_counts$pfs_risk, levels = c("NP.HR", "RP.HR", "Prog.HR", "NP.SR", "RP.SR", "Prog.SR"))
        }
        # Adding a small offset for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)
        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk + siteXbatch, # siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_risk + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }
        # out <- HelperFun2(dr_fit_common)
    } else if (variable_mode == "pfs_Therapy_T") {
        groupA <- "NP.T"
        groupB <- "RP.T"
        fitted_variable <- "pfs_therapy"
        if (variable_mode == "pfs_Therapy_T") {
            reformated_counts <- reformated_counts |> dplyr::filter(pfs_therapy %in% c("NP.T", "RP.T", "Prog.T", "NP.D", "RP.D", "Prog.D"))
            reformated_counts$pfs_risk <- factor(reformated_counts$pfs_therapy, levels = c("NP.T", "RP.T", "Prog.T", "NP.D", "RP.D", "Prog.D"))
        }
        # Adding a small offset for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)
        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy + siteXbatch, # siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }
        # out <- HelperFun2(dr_fit_common)
    } else if (variable_mode == "pfs_Therapy_D") {
        groupA <- "NP.D"
        groupB <- "RP.D"
        fitted_variable <- "pfs_therapy"
        if (variable_mode == "pfs_Therapy_D") {
            reformated_counts <- reformated_counts |> dplyr::filter(pfs_therapy %in% c("NP.T", "RP.T", "Prog.T", "NP.D", "RP.D", "Prog.D"))
            reformated_counts$pfs_risk <- factor(reformated_counts$pfs_therapy, levels = c("NP.D", "RP.D", "Prog.D", "NP.T", "RP.T", "Prog.T"))
        }
        # Adding a small offset for zero entries
        cell_proportions <- DR_data(reformated_counts[, cluster_colnames] + offset, trafo = transformation_magnitude)
        if (BC_mode == "None") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "siteXbatch") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy + siteXbatch, # siteXbatch, # Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        } else if (BC_mode == "site_batch_split") {
            dr_fit_common <- DirichReg(
                cell_proportions ~ pfs_therapy + Study_Site + Batch,
                reformated_counts,
                model = "common"
            )
        }
        # out <- HelperFun2(dr_fit_common)
    }


    print("Starting Simulation Stage")
    outputs_all <- dirichlet_pvalue_simulation(dr_fit_common, reformated_counts, cluster_colnames, fitted_variable, groupA, groupB, sim_per_siteXbatch = 50000, OVERRIDE_RETURN_MEAN_ONLY = SIMULATION_OVERRIDE_RETURN_MEAN_ONLY)
    prop_with_var <- cbind(reformated_counts |> dplyr::select(public_id, progression_group, davies_based_risk, siteXbatch), cell_proportions)
    prop_with_var <- prop_with_var |> gather(variable, value, -public_id, -progression_group, -davies_based_risk, -siteXbatch)
    colnames(prop_with_var) <- c("public_id", "progression_group", "davies_based_risk", "siteXbatch", "cluster", "proportion")
    # do_all_dirichlet_plots_per(outputs_all, prop_with_var = prop_with_var, output = output_file, nameMod = "test_cd8")
    print("Simulation Stage Complete")

    muGroupA <- rowMeans(outputs_all[["mu_groupA"]])
    muGroupB <- rowMeans(outputs_all[["mu_groupB"]])

    outputs_all$fold_change <- muGroupA / muGroupB
    outputs_all$log2FC <- log2(outputs_all$fold_change)


    muGroupA <- outputs_all[["mu_groupA"]]
    muGroupB <- outputs_all[["mu_groupB"]]

    outputs_all$all_fold_change <- muGroupA / muGroupB
    outputs_all$all_log2FC <- log2(outputs_all$all_fold_change)
    outputs_all$average_log2FC <- rowMeans(outputs_all$all_log2FC)

    out_list <- list(
        # "out" = out,
        "model" = dr_fit_common,
        "simulation_out" = outputs_all,
        "prop_with_var" = prop_with_var
    )

    return(out_list)
}
