# Patch ANCOMBC::ancombc2 (2.10.1) sensitivity analysis for trend-only calls.
#
# Bug: with trend=TRUE, pseudo_sens=TRUE, and global=FALSE, sensitivity forces
# global=TRUE only for pseudo re-runs, then builds ss_3d_global from
# c(res_main, ss_list). res_main$res_global is NULL → dim drop →
# apply(..., c(1,2), ...) error: "'MARGIN' does not match dim(X)".
#
# Original R/08_ancombc_stability.R avoided this by always setting global=TRUE
# alongside trend=TRUE. Week-9 R/12 split those runs and hit the bug.
#
# This patch recomputes a zero-pseudo global fit when res_main$res_global is
# missing, keeping pseudo_sens enabled.

suppressPackageStartupMessages(library(ANCOMBC))

ancombc2_patched <- function (data, taxa_are_rows = TRUE, assay.type = assay_name, 
    assay_name = "counts", rank = tax_level, tax_level = NULL, 
    aggregate_data = NULL, meta_data = NULL, fix_formula, rand_formula = NULL, 
    p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE, prv_cut = 0.1, 
    lib_cut = 0, s0_perc = 0.05, group = NULL, struc_zero = FALSE, 
    neg_lb = FALSE, alpha = 0.05, n_cl = 1, verbose = TRUE, global = FALSE, 
    pairwise = FALSE, dunnet = FALSE, trend = FALSE, iter_control = list(tol = 0.01, 
        max_iter = 20, verbose = FALSE), em_control = list(tol = 1e-05, 
        max_iter = 100), lme_control = lme4::lmerControl(), mdfdr_control = list(fwer_ctrl_method = "holm", 
        B = 100), trend_control = list(contrast = NULL, node = NULL, 
        solver = "ECOS", B = 100)) 
{
    if (n_cl > 1) {
        cl = parallel::makeCluster(n_cl)
        doParallel::registerDoParallel(cl)
    }
    else {
        foreach::registerDoSEQ()
    }
    check_results = data_sanity_check(data = data, taxa_are_rows = taxa_are_rows, 
        assay.type = assay_name, assay_name = assay_name, rank = tax_level, 
        tax_level = tax_level, aggregate_data = aggregate_data, 
        meta_data = meta_data, fix_formula = fix_formula, group = group, 
        struc_zero = struc_zero, global = global, pairwise = pairwise, 
        dunnet = dunnet, mdfdr_control = mdfdr_control, trend = trend, 
        trend_control = trend_control, verbose = verbose)
    feature_table = check_results$feature_table
    feature_table_aggregate = check_results$feature_table_aggregate
    meta_data = check_results$meta_data
    global = check_results$global
    pairwise = check_results$pairwise
    dunnet = check_results$dunnet
    trend = check_results$trend
    trend_control = check_results$trend_control
    if (struc_zero) {
        zero_ind = .get_struc_zero(data = feature_table_aggregate, 
            meta_data = meta_data, group = group, neg_lb = neg_lb)
        tax_idx = apply(zero_ind[, -1], 1, function(x) all(x == 
            FALSE))
        tax_keep = which(tax_idx)
    }
    else {
        zero_ind = NULL
        tax_keep = seq(nrow(feature_table_aggregate))
    }
    core1 = .data_core(data = feature_table, meta_data = meta_data, 
        prv_cut = prv_cut, lib_cut = lib_cut, tax_keep = NULL, 
        samp_keep = NULL)
    O1 = core1$feature_table
    samp_keep = colnames(O1)
    core2 = .data_core(data = feature_table_aggregate, meta_data = meta_data, 
        prv_cut = prv_cut, lib_cut = lib_cut, tax_keep = tax_keep, 
        samp_keep = samp_keep)
    O2 = core2$feature_table
    meta_data = core2$meta_data
    res_main = .ancombc2_core(data = O1, aggregate_data = O2, 
        meta_data = meta_data, fix_formula = fix_formula, rand_formula = rand_formula, 
        p_adj_method = p_adj_method, pseudo = pseudo, s0_perc = s0_perc, 
        group = group, alpha = alpha, verbose = verbose, global = global, 
        pairwise = pairwise, dunnet = dunnet, trend = trend, 
        iter_control = iter_control, em_control = em_control, 
        lme_control = lme_control, mdfdr_control = mdfdr_control, 
        trend_control = trend_control)
    if (!pseudo_sens) {
        warn_txt = paste("Sensitivity analysis is currently turned off.", 
            "Since sensitivity analysis is essential for reducing false positives,", 
            "it is highly recommended to enable it unless your primary focus is power.", 
            "Are you sure you want to proceed without it?", sep = "\n")
        message(warn_txt)
        ss_tab = NULL
        res = res_main$res
        res_global = res_main$res_global
        res_pair = res_main$res_pair
        res_dunn = res_main$res_dunn
        res_trend = res_main$res_trend
    }
    if (pseudo_sens) {
        message_txt = paste("Conducting sensitivity analysis for pseudo-count addition to 0s ...", 
            "For taxa that are significant but do not pass the sensitivity analysis,", 
            "they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.", 
            "For detailed instructions on performing sensitivity analysis, please refer to the package vignette.", 
            sep = "\n")
        message(message_txt)
        pseudo_list = c(0.1, 0.5, 1)
        if (trend) 
            global = TRUE
        iter_control$verbose = FALSE
        ss_list = lapply(pseudo_list, function(pseudo_count) {
            res_pseudo = .ancombc2_core(data = O1, aggregate_data = O2, 
                meta_data = meta_data, fix_formula = fix_formula, 
                rand_formula = rand_formula, p_adj_method = p_adj_method, 
                pseudo = pseudo_count, s0_perc = s0_perc, group = group, 
                alpha = alpha, verbose = FALSE, global = global, 
                pairwise = pairwise, dunnet = dunnet, trend = FALSE, 
                iter_control = iter_control, em_control = em_control, 
                lme_control = lme_control, mdfdr_control = mdfdr_control)
            return(res_pseudo)
        })
        ss_list = c(list(res_main), ss_list)
        # PATCH: ANCOMBC 2.10.1 — trend + pseudo_sens forces global=TRUE only for
        # sensitivity re-runs. If the main fit had global=FALSE, res_main$res_global
        # is NULL; NULL[, "q_val", drop=FALSE] is NULL; dim(NULL) is NULL; then
        # array(..., c(NULL, n)) is 1-D and apply(..., c(1,2), ...) fails with
        # "'MARGIN' does not match dim(X)". Fill the missing zero-pseudo global.
        if (isTRUE(global) && is.null(res_main$res_global)) {
            res_main_g = .ancombc2_core(data = O1, aggregate_data = O2,
                meta_data = meta_data, fix_formula = fix_formula,
                rand_formula = rand_formula, p_adj_method = p_adj_method,
                pseudo = 0, s0_perc = s0_perc, group = group, alpha = alpha,
                verbose = FALSE, global = TRUE, pairwise = pairwise,
                dunnet = dunnet, trend = FALSE, iter_control = iter_control,
                em_control = em_control, lme_control = lme_control,
                mdfdr_control = mdfdr_control)
            res_main$res_global = res_main_g$res_global
            ss_list[[1]] = res_main
        }
        pseudo_list = c(0, pseudo_list)
        ss_list_prim = lapply(ss_list, function(res_pseudo) res_pseudo$res[, 
            grepl("q_", colnames(res_pseudo$res))])
        ss_3d_prim = array(unlist(ss_list_prim), c(dim(ss_list_prim[[1]]), 
            length(ss_list_prim)))
        ss_tab_prim = apply(ss_3d_prim, c(1, 2), function(x) {
            sum(x > alpha)/length(pseudo_list)
        })
        colnames(ss_tab_prim) = gsub("^q_", "ss_prim_", colnames(ss_list_prim[[1]]))
        ss_tab_log = (ss_tab_prim == 0 | ss_tab_prim == 1)
        colnames(ss_tab_log) = gsub("ss_prim_", "passed_ss_", 
            colnames(ss_tab_log))
        res = cbind(res_main$res, ss_tab_log)
        diff_cols = grep("^diff_", names(res), value = TRUE, 
            perl = TRUE)
        suffixes = sub("^diff_", "", diff_cols)
        for (suffix in suffixes) {
            diff_col = paste0("diff_", suffix)
            passed_col = paste0("passed_ss_", suffix)
            new_col = paste0("diff_robust_", suffix)
            res[[new_col]] = res[[diff_col]] & res[[passed_col]]
        }
        if (global) {
            ss_list_global = lapply(ss_list, function(res_pseudo) res_pseudo$res_global[, 
                "q_val", drop = FALSE])
            ss_dims_g = dim(ss_list_global[[1]])
            if (is.null(ss_dims_g) || length(ss_dims_g) < 2) {
                stop(paste0("ANCOMBC sensitivity: ss_list_global[[1]] lacks 2D dim ",
                  "(often res_global=NULL when trend=TRUE & global=FALSE). ",
                  "This should have been patched; please report."), call. = FALSE)
            }
            ss_3d_global = array(unlist(ss_list_global), c(ss_dims_g,
                length(ss_list_global)))
            if (length(dim(ss_3d_global)) < 2) {
                stop("ANCOMBC sensitivity: ss_3d_global is not >=2D after array()",
                  call. = FALSE)
            }
            ss_tab_global = apply(ss_3d_global, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_global) = "ss_global"
            ss_tab_log = (ss_tab_global == 0 | ss_tab_global == 
                1)
            colnames(ss_tab_log) = "passed_ss"
            res_global = cbind(res_main$res_global, ss_tab_log)
            res_global[["diff_robust_abn"]] = res_global[["diff_abn"]] & 
                res_global[["passed_ss"]]
        }
        else {
            res_global = NULL
        }
        if (trend) {
            ss_tab_trend = ss_tab_global
            colnames(ss_tab_trend) = "ss_trend"
            ss_tab_log = (ss_tab_trend == 0 | ss_tab_trend == 
                1)
            colnames(ss_tab_log) = "passed_ss"
            res_trend = cbind(res_main$res_trend, ss_tab_log)
            res_trend[["diff_robust_abn"]] = res_trend[["diff_abn"]] & 
                res_trend[["passed_ss"]]
        }
        else {
            res_trend = NULL
        }
        if (pairwise) {
            ss_list_pair = lapply(ss_list, function(res_pseudo) res_pseudo$res_pair[, 
                grepl("q_", colnames(res_pseudo$res_pair))])
            ss_3d_pair = array(unlist(ss_list_pair), c(dim(ss_list_pair[[1]]), 
                length(ss_list_pair)))
            ss_tab_pair = apply(ss_3d_pair, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_pair) = gsub("^q_", "ss_pair_", colnames(ss_list_pair[[1]]))
            ss_tab_log = (ss_tab_pair == 0 | ss_tab_pair == 1)
            colnames(ss_tab_log) = gsub("ss_pair_", "passed_ss_", 
                colnames(ss_tab_log))
            res_pair = cbind(res_main$res_pair, ss_tab_log)
            diff_cols = grep("^diff_", names(res_pair), value = TRUE, 
                perl = TRUE)
            suffixes = sub("^diff_", "", diff_cols)
            for (suffix in suffixes) {
                diff_col = paste0("diff_", suffix)
                passed_col = paste0("passed_ss_", suffix)
                new_col = paste0("diff_robust_", suffix)
                res_pair[[new_col]] = res_pair[[diff_col]] & 
                  res_pair[[passed_col]]
            }
        }
        else {
            res_pair = NULL
        }
        if (dunnet) {
            ss_list_dunn = lapply(ss_list, function(res_pseudo) res_pseudo$res_dunn[, 
                grepl("q_", colnames(res_pseudo$res_dunn))])
            ss_3d_dunn = array(unlist(ss_list_dunn), c(dim(ss_list_dunn[[1]]), 
                length(ss_list_dunn)))
            ss_tab_dunn = apply(ss_3d_dunn, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_dunn) = gsub("^q_", "ss_dunn_", colnames(ss_list_dunn[[1]]))
            ss_tab_log = (ss_tab_dunn == 0 | ss_tab_dunn == 1)
            colnames(ss_tab_log) = gsub("ss_dunn_", "passed_ss_", 
                colnames(ss_tab_log))
            res_dunn = cbind(res_main$res_dunn, ss_tab_log)
            diff_cols = grep("^diff_", names(res_dunn), value = TRUE, 
                perl = TRUE)
            suffixes = sub("^diff_", "", diff_cols)
            for (suffix in suffixes) {
                diff_col = paste0("diff_", suffix)
                passed_col = paste0("passed_ss_", suffix)
                new_col = paste0("diff_robust_", suffix)
                res_dunn[[new_col]] = res_dunn[[diff_col]] & 
                  res_dunn[[passed_col]]
            }
        }
        else {
            res_dunn = NULL
        }
        ss_tab_cols = list(taxon = rownames(O2), ss_tab_prim)
        if (global) 
            ss_tab_cols$ss_tab_global = ss_tab_global
        if (pairwise) 
            ss_tab_cols$ss_tab_pair = ss_tab_pair
        if (dunnet) 
            ss_tab_cols$ss_tab_dunn = ss_tab_dunn
        if (trend) 
            ss_tab_cols$ss_tab_trend = ss_tab_trend
        ss_tab = do.call(data.frame, c(ss_tab_cols, list(check.names = FALSE)))
    }
    out = list(feature_table = O2, bias_correct_log_table = res_main$bias_correct_log_table, 
        ss_tab = ss_tab, zero_ind = zero_ind, samp_frac = res_main$samp_frac, 
        delta_em = res_main$delta_em, delta_wls = res_main$delta_wls, 
        res = res, res_global = res_global, res_pair = res_pair, 
        res_dunn = res_dunn, res_trend = res_trend)
    if (n_cl > 1) {
        parallel::stopCluster(cl)
    }
    return(out)
}



ns <- asNamespace("ANCOMBC")
environment(ancombc2_patched) <- ns
unlockBinding("ancombc2", ns)
assign("ancombc2", ancombc2_patched, envir = ns)
lockBinding("ancombc2", ns)
tryCatch({
  exp <- as.environment("package:ANCOMBC")
  if (exists("ancombc2", envir = exp, inherits = FALSE)) {
    if (bindingIsLocked("ancombc2", exp)) unlockBinding("ancombc2", exp)
    assign("ancombc2", ancombc2_patched, envir = exp)
    try(lockBinding("ancombc2", exp), silent = TRUE)
  }
}, error = function(e) invisible(NULL))
message("Patched ANCOMBC::ancombc2 for trend-only + pseudo_sens ss_3d_global dim drop")
