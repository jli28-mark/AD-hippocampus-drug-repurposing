# =============================================================================
# Four-study AD hippocampal gene-level random-effects meta-analysis.
# Input: four_studies_common_meta_input_long.csv from the SE harmonization step.
# Only the main analysis is fitted. Source files are never modified.
#
# Statistical model: yi ~ N(mu, vi + tau^2); vi = SE^2.
# REML estimates between-study variance tau^2.
# Standard Wald Z-test (test = "z") is used for hypothesis testing across
# multi-cohort transcriptomic families, avoiding the df=3 Knapp-Hartung penalty.
# BH correction uses the entire planned common-gene family before any
# effect-size or direction filter. A failed fit retains NA statistics.
#
# Candidate signature criteria (prespecified):
# Meta_FDR < 0.05, abs(Meta_log2FC) >= 0.58, >= 3/4 signs matching the pooled sign.
# Study significance is NOT required. Opposite study-FDR evidence is flagged,
# not silently removed. All FDR-significant genes are saved separately.
#
# Documentation:
# https://wviechtb.github.io/metafor/reference/misc-recs.html
# https://wviechtb.github.io/metafor/reference/rma.uni.html
# https://wviechtb.github.io/metafor/reference/predict.rma.html
# =============================================================================

run_four_study_meta_analysis <- function(
    base_dir = "D:/My Desktop/drug/3.Signature_integration",
    input_file = file.path(base_dir, "four_studies_common_meta_input_long.csv"),
    test_method = "z",
    n_workers = 2L,
    fdr_cutoff = 0.05,
    fc_cutoff = 0.58,
    min_direction_support = 3L,
    chunk_size = 200L,
    heatmap_genes_per_direction = 20L,
    forest_n = 6L,
    png_dpi = 600L,
    install_missing = TRUE) {

  studies <- c("GSE278723", "GSE173955", "GSE48350", "GSE5281")
  stopifnot(
    length(test_method) == 1L, test_method %in% c("z", "knha"),
    length(n_workers) == 1L, is.finite(n_workers), n_workers >= 1,
    n_workers == as.integer(n_workers), length(fdr_cutoff) == 1L,
    is.finite(fdr_cutoff), fdr_cutoff > 0, fdr_cutoff < 1,
    length(fc_cutoff) == 1L, is.finite(fc_cutoff), fc_cutoff > 0,
    length(min_direction_support) == 1L, min_direction_support %in% 1:4,
    length(chunk_size) == 1L, is.finite(chunk_size), chunk_size >= 1,
    length(heatmap_genes_per_direction) == 1L,
    is.finite(heatmap_genes_per_direction), heatmap_genes_per_direction >= 1,
    length(forest_n) == 1L, is.finite(forest_n), forest_n >= 0,
    length(png_dpi) == 1L, is.finite(png_dpi), png_dpi >= 72
  )
  chunk_size <- as.integer(chunk_size)
  forest_n <- as.integer(forest_n)
  heatmap_genes_per_direction <- as.integer(heatmap_genes_per_direction)
  required <- c("data.table", "metafor", "ggplot2", "ggrepel")
  
  message("[1/7] Checking packages and creating a new output folder ...")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && install_missing) {
    utils::install.packages(missing, repos = "https://cloud.r-project.org")
  }
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these packages first: ",
    paste(missing, collapse = ", "), call. = FALSE)
  if (!dir.exists(base_dir)) stop("Folder not found: ", base_dir, call. = FALSE)
  if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)

  output_root <- file.path(base_dir, "random_effects_meta")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = output_root)
  dir.create(out, recursive = TRUE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  plots_dir <- file.path(out, "plots")
  dir.create(plots_dir)
  
  csv <- function(x, name) data.table::fwrite(x, file.path(out, name),
    na = "NA", quote = "auto", bom = TRUE)
  write_log <- function(...) {
    s <- paste0(...)
    message(s)
    cat(s, "\n", file = file.path(out, "run.log"), append = TRUE)
  }
  
  checks <- list()
  check <- function(name, ok, details = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      Check = name, Passed = isTRUE(ok), Details = details)
    if (!isTRUE(ok)) {
      csv(do.call(rbind, checks), "input_validation_checks.csv")
      stop(name, ": ", details, " [see ", out, "]", call. = FALSE)
    }
    invisible(TRUE)
  }
  
  settings <- list(
    analysis = "Four-study gene-level random-effects meta-analysis",
    studies = studies, contrast = "AD minus Control/NC",
    method = "REML", test = test_method, confidence_level = 95,
    multiple_testing = "BH over the entire input common-gene family",
    FDR_cutoff = fdr_cutoff, pooled_abs_log2FC_cutoff = fc_cutoff,
    min_direction_support = as.integer(min_direction_support),
    opposite_study_FDR = "reported, not a primary candidate exclusion",
    study_significance_requirement = "none",
    heterogeneity_filter = "none", study_weight_filter = "none",
    requested_workers = as.integer(n_workers), chunk_size = chunk_size,
    heatmap_genes_per_direction = heatmap_genes_per_direction,
    forest_n = forest_n, png_dpi = as.integer(png_dpi),
    input_file = normalizePath(input_file, winslash = "/"),
    input_md5 = unname(tools::md5sum(input_file)),
    started = as.character(Sys.time()))
  saveRDS(settings, file.path(out, "analysis_settings.rds"))
  dput(settings, file = file.path(out, "analysis_settings.R"))
  csv(data.frame(File = settings$input_file, MD5 = settings$input_md5),
    "input_checksums.csv")
  check("input_snapshot", file.copy(input_file, file.path(out, basename(input_file))),
    "A copy of the exact input is retained inside the new run folder.")
  csv(data.frame(Package = required,
    Version = vapply(required, function(p) as.character(utils::packageVersion(p)), character(1))),
    "package_versions.csv")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  
  writeLines(c("run_four_study_meta_analysis <-",
    deparse(run_four_study_meta_analysis, width.cutoff = 100L)),
    file.path(out, "analysis_function_used.R"))

  write_log("[2/7] Validating study identities, gene alignment and paired effects/SE ...")
  dat <- data.table::fread(input_file, data.table = FALSE, check.names = FALSE,
    colClasses = c(Dataset = "character", HGNC_ID = "character", gene_symbol = "character"),
    na.strings = c("", "NA", "NaN"), showProgress = FALSE)
  names(dat) <- trimws(sub("^\ufeff", "", names(dat)))
  required_cols <- c("Dataset", "HGNC_ID", "gene_symbol", "log2FC", "SE",
    "SE_method", "Model_coefficient")
  check("required_columns", all(required_cols %in% names(dat)),
    paste("Required:", paste(required_cols, collapse = ", ")))
  for (nm in c("Dataset", "HGNC_ID", "gene_symbol", "SE_method", "Model_coefficient")) {
    dat[[nm]] <- trimws(as.character(dat[[nm]]))
  }
  for (nm in c("log2FC", "SE")) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  for (nm in c("P_value", "FDR")) {
    if (!nm %in% names(dat)) dat[[nm]] <- NA_real_
    original <- dat[[nm]]
    dat[[nm]] <- suppressWarnings(as.numeric(original))
    check(paste0("numeric_", nm),
      !any(!is.na(original) & is.na(dat[[nm]])),
      paste(nm, "must contain numeric probabilities or NA."))
    check(paste0("range_", nm),
      all(is.na(dat[[nm]]) | (is.finite(dat[[nm]]) & dat[[nm]] >= 0 & dat[[nm]] <= 1)),
      paste(nm, "must be within [0, 1] when supplied."))
  }
  check("nonempty_identifiers",
    all(!is.na(dat$Dataset) & nzchar(dat$Dataset) &
        !is.na(dat$HGNC_ID) & nzchar(dat$HGNC_ID) &
        !is.na(dat$gene_symbol) & nzchar(dat$gene_symbol)),
    "Study, HGNC ID and gene symbol must be present on every row.")
  check("exactly_four_expected_studies", setequal(unique(dat$Dataset), studies),
    paste("Expected:", paste(studies, collapse = ", ")))
  check("one_row_per_gene_study",
    !anyDuplicated(paste(dat$Dataset, dat$HGNC_ID, sep = "|")),
    "Duplicate gene-study records must be resolved upstream.")
  good <- is.finite(dat$log2FC) & is.finite(dat$SE) & dat$SE > 0 &
    is.finite(dat$SE^2) & dat$SE^2 > 0
  if (any(!good)) csv(dat[!good, , drop = FALSE], "invalid_effect_SE_input_rows.csv")
  check("finite_effect_and_positive_finite_SE", all(good),
    "No effect/SE cap, floor or imputation is applied.")
  if ("variance" %in% names(dat)) {
    supplied_v <- suppressWarnings(as.numeric(dat$variance))
    check("variance_matches_SE_squared",
      all(is.finite(supplied_v) & supplied_v > 0 &
        abs(supplied_v - dat$SE^2) <=
          1e-8 * pmax(supplied_v, dat$SE^2, .Machine$double.eps)),
      "The model uses SE squared as sampling variance.")
  }
  if ("locus_group" %in% names(dat)) {
    check("protein_coding_annotation",
      all(!is.na(dat$locus_group) & dat$locus_group == "protein-coding gene"),
      "The harmonized input must contain protein-coding genes only.")
  }
  expected_se <- c(
    GSE278723 = "limma_ordinary_SE_sigma_times_stdev_unscaled",
    GSE173955 = "DESeq2_Wald_model_SE_unshrunken",
    GSE48350 = "limma_ordinary_SE_sigma_times_stdev_unscaled",
    GSE5281 = "limma_ordinary_SE_sigma_times_stdev_unscaled")
  expected_coef <- c(GSE278723 = "diagnosisAD",
    GSE173955 = "Group: AD - Control",
    GSE48350 = "AD_vs_Control", GSE5281 = "GroupAD")
  check("native_model_SE_types",
    all(!is.na(dat$SE_method) & dat$SE_method == unname(expected_se[dat$Dataset])),
    "The supplied SE must come from the audited main models.")
  check("AD_minus_control_coefficients",
    all(!is.na(dat$Model_coefficient) &
      dat$Model_coefficient == unname(expected_coef[dat$Dataset])),
    "All four study effects must retain the audited AD-minus-control direction.")
  symbol_per_id <- vapply(split(dat$gene_symbol, dat$HGNC_ID),
    function(z) length(unique(z)), integer(1))
  id_per_symbol <- vapply(split(dat$HGNC_ID, dat$gene_symbol),
    function(z) length(unique(z)), integer(1))
  check("consistent_HGNC_symbol_pairs",
    all(symbol_per_id == 1L) && all(id_per_symbol == 1L),
    "Each HGNC ID and approved symbol must identify one common gene.")
  genes <- sort(unique(dat$HGNC_ID))
  n_genes <- length(genes)
  check("nonempty_common_universe", n_genes > 0L, "At least one common gene is required.")
  aligned <- setNames(lapply(studies, function(s) {
    z <- dat[dat$Dataset == s, , drop = FALSE]
    z[match(genes, z$HGNC_ID), , drop = FALSE]
  }), studies)
  check("four_studies_per_gene",
    nrow(dat) == 4L * n_genes &&
      all(vapply(aligned, function(z) !anyNA(z$HGNC_ID) &&
        identical(z$HGNC_ID, genes), logical(1))),
    "Every planned gene must have exactly one paired effect/SE in each study.")
  gene_info <- aligned[[1L]][, c("HGNC_ID", "gene_symbol"), drop = FALSE]
  Y <- do.call(cbind, lapply(aligned, function(z) z$log2FC))
  V <- do.call(cbind, lapply(aligned, function(z) z$SE^2))
  S <- sqrt(V)
  F <- do.call(cbind, lapply(aligned, function(z) z$FDR))
  colnames(Y) <- colnames(V) <- colnames(S) <- colnames(F) <- studies
  rownames(Y) <- rownames(V) <- rownames(S) <- rownames(F) <- genes
  csv(do.call(rbind, checks), "input_validation_checks.csv")
  csv(data.frame(Dataset = studies,
    N_genes = vapply(aligned, nrow, integer(1)),
    SE_min = apply(S, 2L, min), SE_median = apply(S, 2L, stats::median),
    SE_max = apply(S, 2L, max),
    Missing_study_FDR = colSums(is.na(F))), "input_summary.csv")
  settings$n_common_genes <- n_genes
  write_log("Input: ", n_genes, " common genes x 4 studies = ", nrow(dat), " rows.")

  # 拟合单基因 REML 随机效应模型
  meta_fit_one <- function(yi, vi, study_names, test_type = "z") {
    notes <- character()
    attempts <- 1L
    do_fit <- function(control) {
      withCallingHandlers(
        tryCatch(metafor::rma.uni(yi = yi, vi = vi, method = "REML",
          test = test_type, level = 95, control = control),
          error = function(e) e),
        warning = function(w) {
          notes <<- c(notes, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    }
    m <- do_fit(list(maxiter = 1000L))
    first_error <- ""
    if (inherits(m, "error")) {
      first_error <- conditionMessage(m)
      attempts <- 2L
      m <- do_fit(list(maxiter = 10000L, stepadj = 0.5))
    }
    answer <- data.frame(
      Meta_log2FC = NA_real_, Meta_SE = NA_real_,
      Meta_SE_inverse_variance = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      Meta_stat = NA_real_, DF = if (test_type == "knha") 3L else Inf,
      Meta_P = NA_real_, Tau2 = NA_real_,
      I2_percent = NA_real_, Q = NA_real_, Q_P = NA_real_,
      PI_low = NA_real_, PI_high = NA_real_, Fit_status = "failed",
      Fit_attempts = attempts, Fit_warning = paste(unique(notes), collapse = " | "),
      First_fit_error = first_error, Fit_error = "", stringsAsFactors = FALSE)
    for (s in study_names) answer[[paste0(s, "_Weight_percent")]] <- NA_real_
    if (inherits(m, "error")) {
      answer$Fit_error <- conditionMessage(m)
      return(answer)
    }
    values <- c(m$b[1L], m$se[1L], m$ci.lb[1L], m$ci.ub[1L],
      m$zval[1L], m$pval[1L], m$tau2, m$I2)
    if (length(values) != 8L || any(!is.finite(values)) ||
        m$se[1L] <= 0 || m$pval[1L] < 0 || m$pval[1L] > 1 ||
        m$tau2 < 0 || m$ci.lb[1L] > m$ci.ub[1L]) {
      answer$Fit_status <- "non_estimable_inference"
      answer$Fit_error <- "Non-finite or invalid pooled estimates/inference."
      return(answer)
    }
    w <- 1 / (vi + m$tau2)
    w <- w / sum(w)
    if (!all(is.finite(w)) || any(w <= 0) ||
        abs(sum(w) - 1) > 1e-8 ||
        abs(sum(w * yi) - as.numeric(m$b[1L])) >
          1e-7 * max(1, abs(as.numeric(m$b[1L])))) {
      answer$Fit_status <- "invalid_weight_check"
      answer$Fit_error <- "Random-effects weights failed numerical verification."
      return(answer)
    }
    prediction <- withCallingHandlers(
      tryCatch(stats::predict(m, level = 95), error = function(e) e),
      warning = function(w) {
        notes <<- c(notes, paste0("Prediction: ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    answer$Meta_log2FC <- as.numeric(m$b[1L])
    answer$Meta_SE <- as.numeric(m$se[1L])
    answer$Meta_SE_inverse_variance <- sqrt(1 / sum(1 / (vi + m$tau2)))
    answer$CI_low <- as.numeric(m$ci.lb[1L])
    answer$CI_high <- as.numeric(m$ci.ub[1L])
    answer$Meta_stat <- as.numeric(m$zval[1L])
    answer$Meta_P <- as.numeric(m$pval[1L])
    answer$Tau2 <- as.numeric(m$tau2)
    answer$I2_percent <- as.numeric(m$I2)
    answer$Q <- as.numeric(m$QE)
    answer$Q_P <- as.numeric(m$QEp)
    if (!inherits(prediction, "error")) {
      answer$PI_low <- as.numeric(prediction$pi.lb[1L])
      answer$PI_high <- as.numeric(prediction$pi.ub[1L])
    } else {
      notes <- c(notes, paste0("Prediction unavailable: ", conditionMessage(prediction)))
    }
    answer$Fit_status <- "ok"
    answer$Fit_warning <- paste(unique(notes), collapse = " | ")
    for (j in seq_along(study_names)) {
      answer[[paste0(study_names[j], "_Weight_percent")]] <- 100 * w[j]
    }
    answer
  }
  environment(meta_fit_one) <- baseenv()

  test_desc <- if (test_method == "z") "standard Wald Z inference" else "Knapp-Hartung inference"
  write_log("[3/7] Fitting per-gene REML models with ", test_desc, " ...")
  chunks <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / chunk_size))
  available <- parallel::detectCores(logical = TRUE)
  if (is.na(available)) available <- 2L
  workers <- min(as.integer(n_workers), max(1L, available - 1L), length(chunks))
  cl <- NULL
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (workers > 1L) {
    cl <- tryCatch({
      candidate <- parallel::makePSOCKcluster(workers)
      candidate
    }, error = function(e) {
      write_log("Parallel startup failed; using one process: ", conditionMessage(e))
      NULL
    })
    if (!is.null(cl)) {
      worker_setup <- function(paths) {
        .libPaths(paths)
        Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
          MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
        suppressPackageStartupMessages(library(metafor))
        NULL
      }
      environment(worker_setup) <- baseenv()
      setup <- tryCatch({
        parallel::clusterCall(cl, worker_setup, .libPaths())
        exports <- new.env(parent = baseenv())
        exports$meta_Y <- Y
        exports$meta_V <- V
        exports$meta_studies <- studies
        exports$meta_test_method <- test_method
        exports$meta_fit_one <- meta_fit_one
        parallel::clusterExport(cl,
          c("meta_Y", "meta_V", "meta_studies", "meta_test_method", "meta_fit_one"), envir = exports)
        TRUE
      }, error = function(e) {
        write_log("Worker setup failed; using one process: ", conditionMessage(e))
        FALSE
      })
      if (!setup) {
        parallel::stopCluster(cl)
        cl <- NULL
      }
    }
  }
  if (is.null(cl)) workers <- 1L
  settings$actual_workers <- workers
  write_log("Processes: ", workers, "; common genes: ", n_genes, ".")
  
  meta_fit_batch <- function(ii) {
    lapply(ii, function(i) meta_fit_one(meta_Y[i, ], meta_V[i, ], meta_studies, test_type = meta_test_method))
  }
  environment(meta_fit_batch) <- .GlobalEnv
  records <- vector("list", n_genes)
  start <- proc.time()[["elapsed"]]
  last_checkpoint <- 0L
  for (first in seq.int(1L, length(chunks), by = workers)) {
    jobs <- seq.int(first, min(length(chunks), first + workers - 1L))
    if (!is.null(cl)) {
      payload <- parallel::parLapply(cl, chunks[jobs], meta_fit_batch)
    } else {
      payload <- lapply(chunks[jobs], function(ii) {
        lapply(ii, function(i) meta_fit_one(Y[i, ], V[i, ], studies, test_type = test_method))
      })
    }
    for (j in seq_along(jobs)) records[chunks[[jobs[j]]]] <- payload[[j]]
    completed <- max(unlist(chunks[jobs], use.names = FALSE))
    write_log("  Fitted ", completed, "/", n_genes, " genes; elapsed ",
      round(proc.time()[["elapsed"]] - start, 1), " s.")
    if (completed - last_checkpoint >= 2000L || completed == n_genes) {
      saveRDS(list(completed = completed, gene_order = genes,
        fits = as.data.frame(data.table::rbindlist(records[seq_len(completed)],
          use.names = TRUE, fill = TRUE)), settings = settings),
        file.path(out, "fit_checkpoint.rds"))
      last_checkpoint <- completed
    }
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  fit_table <- as.data.frame(data.table::rbindlist(records, use.names = TRUE, fill = TRUE))
  check("one_output_per_planned_gene", nrow(fit_table) == n_genes,
    "Failed fits are retained rather than silently dropped.")
  r <- cbind(gene_info, fit_table)
  ok <- r$Fit_status == "ok"
  r$Meta_FDR <- NA_real_
  if (any(ok)) r$Meta_FDR[ok] <- stats::p.adjust(r$Meta_P[ok],
    method = "BH", n = n_genes)

  write_log("[4/7] Applying full-family BH correction, then assessing direction support ...")
  r$N_studies <- 4L
  r$N_UP <- rowSums(Y > 0)
  r$N_DOWN <- rowSums(Y < 0)
  r$N_zero_effect <- rowSums(Y == 0)
  pooled_sign <- sign(r$Meta_log2FC)
  same <- sign(Y) == matrix(pooled_sign, nrow = n_genes, ncol = 4L)
  opposite <- sign(Y) == -matrix(pooled_sign, nrow = n_genes, ncol = 4L) &
    matrix(pooled_sign != 0, nrow = n_genes, ncol = 4L)
  same <- same & matrix(pooled_sign != 0, nrow = n_genes, ncol = 4L)
  study_sig <- is.finite(F) & F < fdr_cutoff
  r$N_same_direction <- rowSums(same)
  r$N_opposite_direction <- rowSums(opposite)
  r$N_study_FDR_significant <- rowSums(study_sig)
  r$N_same_direction_FDR_significant <- rowSums(same & study_sig)
  r$N_opposite_direction_FDR_significant <- rowSums(opposite & study_sig)
  r$N_study_FDR_missing <- rowSums(is.na(F))
  r$Opposite_FDR_evidence <- ifelse(!ok, "Not_estimable",
    ifelse(r$N_opposite_direction_FDR_significant > 0, "Present",
      ifelse(rowSums(opposite & is.na(F)) > 0, "Unknown_missing_FDR", "Not_observed")))
  r$Opposite_FDR_datasets <- vapply(seq_len(n_genes), function(i) {
    if (!ok[i]) return(NA_character_)
    paste(studies[which(opposite[i, ] & study_sig[i, ])], collapse = ";")
  }, character(1))
  r$Pooled_direction <- ifelse(!ok, "Not_estimable",
    ifelse(pooled_sign > 0, "UP", ifelse(pooled_sign < 0, "DOWN", "ZERO")))
  r$Direction_support <- ifelse(!ok, "Not_estimable",
    paste0(r$N_same_direction, "/4"))
  r$Meta_FDR_significant <- ok & !is.na(r$Meta_FDR) & r$Meta_FDR < fdr_cutoff
  r$Meta_FDR_and_FC <- r$Meta_FDR_significant & abs(r$Meta_log2FC) >= fc_cutoff
  r$Signature_candidate <- r$Meta_FDR_and_FC &
    !is.na(r$N_same_direction) & r$N_same_direction >= min_direction_support
  r$Signature_4of4 <- r$Signature_candidate &
    !is.na(r$N_same_direction) & r$N_same_direction == 4L
  r$Volcano_group <- ifelse(!ok, "Not_estimable",
    ifelse(r$Meta_FDR_and_FC, r$Pooled_direction, "Not_DEG"))
  weight_cols <- paste0(studies, "_Weight_percent")
  W <- as.matrix(r[, weight_cols, drop = FALSE])
  r$Max_study_weight_percent <- NA_real_
  r$Max_weight_dataset <- NA_character_
  if (any(ok)) {
    r$Max_study_weight_percent[ok] <- apply(W[ok, , drop = FALSE], 1L, max)
    r$Max_weight_dataset[ok] <- studies[max.col(W[ok, , drop = FALSE], ties.method = "first")]
  }
  for (j in seq_along(studies)) {
    for (v in c("log2FC", "SE", "P_value", "FDR")) {
      r[[paste0(studies[j], "_", v)]] <- aligned[[j]][[v]]
    }
  }
  rank_order <- order(!ok, r$Meta_FDR, r$Meta_P, -abs(r$Meta_log2FC),
    r$HGNC_ID, na.last = TRUE)
  r <- r[rank_order, , drop = FALSE]
  rownames(r) <- NULL
  r$Meta_rank <- NA_integer_
  fitted_rows <- which(r$Fit_status == "ok")
  r$Meta_rank[fitted_rows] <- seq_along(fitted_rows)
  candidates <- r[r$Signature_candidate, , drop = FALSE]
  significant <- r[r$Meta_FDR_significant, , drop = FALSE]
  diagnostics <- r[, c("HGNC_ID", "gene_symbol", "Fit_status", "Fit_attempts",
    "Fit_warning", "First_fit_error", "Fit_error"), drop = FALSE]
  summary <- data.frame(Stage = c("Input_common_genes", "Successful_models",
    "Failed_or_nonestimable_models", "Models_with_warnings", "Optimizer_retries",
    "Meta_FDR_pass", "Meta_FDR_and_FC_pass", "Signature_candidates",
    "Signature_UP", "Signature_DOWN", "Signature_4of4_subset",
    "Signature_with_opposite_study_FDR_evidence"),
    N_genes = c(n_genes, sum(r$Fit_status == "ok"), sum(r$Fit_status != "ok"),
      sum(nzchar(r$Fit_warning)), sum(r$Fit_attempts > 1L),
      nrow(significant), sum(r$Meta_FDR_and_FC), nrow(candidates),
      sum(candidates$Pooled_direction == "UP"),
      sum(candidates$Pooled_direction == "DOWN"), sum(r$Signature_4of4),
      sum(candidates$Opposite_FDR_evidence == "Present")))
  direction_summary <- do.call(rbind, lapply(
    c("All_fitted", "Meta_FDR_pass", "Meta_FDR_and_FC_pass", "Signature_candidates"),
    function(stage) {
      use <- switch(stage,
        All_fitted = r$Fit_status == "ok",
        Meta_FDR_pass = r$Meta_FDR_significant,
        Meta_FDR_and_FC_pass = r$Meta_FDR_and_FC,
        Signature_candidates = r$Signature_candidate)
      data.frame(Stage = stage, Direction_support = paste0(0:4, "/4"),
        N_genes = as.integer(table(factor(r$N_same_direction[use], levels = 0:4))))
    }))

  write_log("[5/7] Saving full results, candidate tables and study weights ...")
  csv(r, "meta_all_genes.csv")
  csv(significant, "meta_FDR_significant.csv")
  csv(candidates, "signature_candidates.csv")
  csv(candidates[candidates$Pooled_direction == "UP", , drop = FALSE], "signature_UP.csv")
  csv(candidates[candidates$Pooled_direction == "DOWN", , drop = FALSE], "signature_DOWN.csv")
  csv(r[r$Signature_4of4, , drop = FALSE], "signature_4of4_subset.csv")
  csv(candidates[candidates$Opposite_FDR_evidence == "Not_observed", , drop = FALSE],
    "signature_no_observed_opposite_study_FDR.csv")
  csv(summary, "analysis_summary.csv")
  csv(direction_summary, "direction_support_summary.csv")
  csv(diagnostics, "fit_diagnostics.csv")
  csv(diagnostics[diagnostics$Fit_status != "ok", , drop = FALSE], "failed_models.csv")
  long_weight <- do.call(rbind, lapply(studies, function(s) {
    z <- aligned[[s]]
    k <- match(z$HGNC_ID, r$HGNC_ID)
    data.frame(Dataset = s, HGNC_ID = z$HGNC_ID, gene_symbol = z$gene_symbol,
      log2FC = z$log2FC, SE = z$SE, Study_FDR = z$FDR,
      Weight_percent = r[[paste0(s, "_Weight_percent")]][k],
      Meta_log2FC = r$Meta_log2FC[k], Meta_FDR = r$Meta_FDR[k],
      Signature_candidate = r$Signature_candidate[k],
      Fit_status = r$Fit_status[k], stringsAsFactors = FALSE)
  }))
  csv(long_weight, "study_effects_weights_long.csv")
  weight_summary <- do.call(rbind, lapply(studies, function(s) {
    w <- long_weight$Weight_percent[long_weight$Dataset == s &
      long_weight$Fit_status == "ok"]
    data.frame(Dataset = s, N_genes = length(w),
      Min_percent = if (length(w)) min(w) else NA_real_,
      Median_percent = if (length(w)) stats::median(w) else NA_real_,
      Mean_percent = if (length(w)) mean(w) else NA_real_,
      Max_percent = if (length(w)) max(w) else NA_real_)
  }))
  csv(weight_summary, "study_weight_summary.csv")
  R <- stats::cor(Y, method = "spearman")
  csv(data.frame(Dataset = rownames(R), R, check.names = FALSE),
    "study_log2FC_spearman.csv")
  check("retained_input_identity",
    identical(settings$input_md5, unname(tools::md5sum(input_file))),
    "The input must not change while the analysis is running.")
  csv(do.call(rbind, checks), "input_validation_checks.csv")
  settings$completed_models <- sum(r$Fit_status == "ok")
  settings$BH_family_size <- n_genes
  settings$completed_model_fitting <- as.character(Sys.time())
  dput(settings, file = file.path(out, "analysis_settings.R"))
  saveRDS(settings, file.path(out, "analysis_settings.rds"))
  saveRDS(list(results = r, candidates = candidates, summary = summary,
    settings = settings, study_effects = Y, study_variances = V,
    study_FDR = F, gene_information = gene_info,
    input_checks = do.call(rbind, checks)),
    file.path(out, "meta_analysis_objects.rds"))
  print(summary, row.names = FALSE)
  if (!nrow(candidates)) {
    write_log("No genes meet the candidate criteria; empty tables are saved. ",
      "The statistical model and thresholds are unchanged.")
  }

  write_log("[6/7] Saving publication-style PDF and PNG figures ...")
  colour <- c(DOWN = "#3B7EA1", Not_DEG = "#B8BDC3", UP = "#D55E00")
  paper <- function() {
    ggplot2::theme_classic(base_size = 11, base_family = "sans") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.4),
        axis.ticks = ggplot2::element_line(linewidth = 0.35),
        axis.text = ggplot2::element_text(colour = "black"),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 9),
        plot.caption = ggplot2::element_text(size = 8, hjust = 0),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        plot.margin = ggplot2::margin(9, 12, 9, 9))
  }
  plot_log <- list()
  save_plot <- function(p, stem, width = 7, height = 6) {
    for (extension in c("pdf", "png")) {
      destination <- file.path(plots_dir, paste0(stem, ".", extension))
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      status <- tryCatch({
        if (extension == "pdf") {
          device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
          ggplot2::ggsave(destination, plot = p, device = device,
            width = width, height = height, units = "in", bg = "white",
            limitsize = FALSE)
        } else {
          ggplot2::ggsave(destination, plot = p, device = "png", dpi = png_dpi,
            width = width, height = height, units = "in", bg = "white",
            limitsize = FALSE)
        }
        list(Status = "saved", Error = "")
      }, error = function(e) list(Status = "failed", Error = conditionMessage(e)))
      plot_log[[length(plot_log) + 1L]] <<- data.frame(
        File = destination, Status = status$Status, Error = status$Error)
      if (status$Status == "failed") write_log("Plot export failed: ",
        destination, ": ", status$Error)
    }
    invisible(p)
  }
  blank_plot <- function(title, note) {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = note, size = 3.5) + ggplot2::labs(title = title) +
      ggplot2::theme_void(base_family = "sans") +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
  plot_data <- r[r$Fit_status == "ok", , drop = FALSE]
  if (nrow(plot_data)) {
    plot_data$neglog10_FDR <- -log10(pmax(plot_data$Meta_FDR, .Machine$double.xmin))
    plot_data$Volcano_group <- factor(plot_data$Volcano_group,
      levels = c("DOWN", "Not_DEG", "UP"))
    plot_data <- plot_data[order(plot_data$Volcano_group != "Not_DEG"), , drop = FALSE]
    label_data <- do.call(rbind, lapply(c("DOWN", "UP"), function(s) {
      head(r[r$Volcano_group == s, , drop = FALSE], 5L)
    }))
    label_data$neglog10_FDR <- -log10(pmax(label_data$Meta_FDR, .Machine$double.xmin))
    counts <- table(plot_data$Volcano_group)
    legend_labels <- c(
      DOWN = paste0("Down (", counts["DOWN"], ")"),
      Not_DEG = paste0("Other tested genes (", counts["Not_DEG"], ")"),
      UP = paste0("Up (", counts["UP"], ")"))
    model_label <- if (test_method == "z") "REML / Wald Z-test" else "REML / Knapp-Hartung"
    volcano <- ggplot2::ggplot(plot_data,
      ggplot2::aes(x = Meta_log2FC, y = neglog10_FDR)) +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
        colour = "grey45", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed",
        colour = "grey45", linewidth = 0.35) +
      ggplot2::geom_point(ggplot2::aes(colour = Volcano_group),
        size = 1.05, alpha = 0.65, shape = 16) +
      ggplot2::scale_colour_manual(values = colour, labels = legend_labels, drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 2, alpha = 1))) +
      ggplot2::labs(title = "AD hippocampal transcriptome meta-analysis",
        subtitle = paste0(model_label, "; 4 studies; ", nrow(plot_data), " fitted genes"),
        x = expression(paste("Pooled ", log[2], " fold change (AD / control)")),
        y = expression(-log[10]("Meta FDR")),
        caption = paste0("Colours: BH FDR < ", fdr_cutoff,
          " and |pooled log2FC| >= ", fc_cutoff, ".\n",
          "Direction support is assessed after multiple-testing correction.")) + paper()
    if (nrow(label_data)) {
      volcano <- volcano + ggrepel::geom_text_repel(data = label_data,
        ggplot2::aes(label = gene_symbol), size = 3, seed = 42,
        min.segment.length = 0, segment.size = 0.25, segment.colour = "grey50",
        box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE)
    }
  } else {
    volcano <- blank_plot("AD hippocampal transcriptome meta-analysis",
      "No models have estimable inference.\nSee failed_models.csv.")
  }
  save_plot(volcano, "meta_volcano", 7, 6.2)

  direction_summary$Stage <- factor(direction_summary$Stage,
    levels = c("All_fitted", "Meta_FDR_pass", "Meta_FDR_and_FC_pass", "Signature_candidates"),
    labels = c("All fitted genes", "Meta FDR significant",
      "Meta FDR and effect cutoff", "Candidate signature"))
  direction_plot <- ggplot2::ggplot(direction_summary,
    ggplot2::aes(x = Direction_support, y = N_genes)) +
    ggplot2::geom_col(fill = "#3B7EA1", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = N_genes), vjust = -0.3, size = 3) +
    ggplot2::facet_wrap(~Stage, scales = "free_y", ncol = 2) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = "Study direction support",
      x = "Studies matching the pooled direction", y = "Number of genes",
      caption = "Signs use all study estimates, regardless of single-study significance.") +
    paper() + ggplot2::theme(strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
  save_plot(direction_plot, "direction_support", 7.4, 6.1)

  if (nrow(plot_data)) {
    heterogeneity <- ggplot2::ggplot(plot_data, ggplot2::aes(x = I2_percent)) +
      ggplot2::geom_histogram(binwidth = 5, boundary = 0,
        fill = "#B8BDC3", colour = "white", linewidth = 0.2) +
      ggplot2::labs(title = "Between-study heterogeneity",
        subtitle = "All estimable common genes",
        x = expression(paste(I^2, " (%)")), y = "Number of genes",
        caption = "I2 is reported descriptively; no heterogeneity cutoff is used.") + paper()
    weight_data <- long_weight[long_weight$Fit_status == "ok", , drop = FALSE]
    weight_data$Dataset <- factor(weight_data$Dataset, levels = studies)
    weight_plot <- ggplot2::ggplot(weight_data,
      ggplot2::aes(x = Dataset, y = Weight_percent)) +
      ggplot2::geom_boxplot(fill = "#DBE7EF", width = 0.6, outlier.shape = NA,
        linewidth = 0.4) +
      ggplot2::coord_cartesian(ylim = c(0, 100)) +
      ggplot2::labs(title = "Study weights across fitted genes",
        x = NULL, y = "Random-effects weight (%)",
        caption = "Weights are proportional to 1 / (SE^2 + tau^2).\n",
        subtitle = "Weights describe precision; they do not establish independent replication.") +
      paper()
  } else {
    heterogeneity <- blank_plot("Between-study heterogeneity", "No estimable models.")
    weight_plot <- blank_plot("Study weights", "No estimable models.")
  }
  save_plot(heterogeneity, "heterogeneity_I2", 6.4, 4.7)
  save_plot(weight_plot, "study_weights", 7, 4.7)
  cor_data <- expand.grid(Dataset_x = studies, Dataset_y = studies, stringsAsFactors = FALSE)
  cor_data$Spearman <- as.vector(R)
  cor_data$Dataset_x <- factor(cor_data$Dataset_x, levels = studies)
  cor_data$Dataset_y <- factor(cor_data$Dataset_y, levels = rev(studies))
  cor_plot <- ggplot2::ggplot(cor_data,
    ggplot2::aes(x = Dataset_x, y = Dataset_y, fill = Spearman)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Spearman)), size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#3B7EA1", mid = "white", high = "#D55E00",
      midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Study effect-size correlations",
      subtitle = paste(n_genes, "common genes; unweighted Spearman correlation"),
      x = NULL, y = NULL) + paper() +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank())
  save_plot(cor_plot, "study_effect_correlations", 6.7, 5.5)

  display_pool <- if (nrow(candidates)) candidates else r[r$Fit_status == "ok", , drop = FALSE]
  display_source <- if (nrow(candidates)) "Candidate signature genes" else
    "Exploratory top-ranked genes: no candidate signature"
  chosen <- do.call(rbind, lapply(c("UP", "DOWN"), function(s) {
    head(display_pool[display_pool$Pooled_direction == s, , drop = FALSE],
      heatmap_genes_per_direction)
  }))
  if (nrow(chosen)) {
    hm <- long_weight[long_weight$HGNC_ID %in% chosen$HGNC_ID, , drop = FALSE]
    hm$Dataset <- factor(hm$Dataset, levels = studies)
    hm$gene_symbol <- factor(hm$gene_symbol, levels = rev(chosen$gene_symbol))
    csv(hm, "displayed_gene_effects.csv")
    limit <- max(abs(hm$log2FC))
    if (!is.finite(limit) || limit == 0) limit <- 1
    hm_plot <- ggplot2::ggplot(hm,
      ggplot2::aes(x = Dataset, y = gene_symbol, fill = log2FC)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
      ggplot2::scale_fill_gradient2(low = "#3B7EA1", mid = "white", high = "#D55E00",
        midpoint = 0, limits = c(-limit, limit), name = "log2FC") +
      ggplot2::labs(title = "Study-level gene effects", subtitle = display_source,
        x = NULL, y = NULL,
        caption = "Genes are ranked by meta FDR within each pooled direction.\n",
        tag = NULL) + paper() +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 8))
    save_plot(hm_plot, "selected_gene_effects_heatmap", 6.7,
      max(4, 2.2 + 0.19 * nrow(chosen)))
  } else {
    save_plot(blank_plot("Study-level gene effects", "No estimable nonzero pooled effects."),
      "selected_gene_effects_heatmap", 6.7, 4)
  }

  selected <- head(display_pool, forest_n)
  forest_data_list <- list()
  pooled_label <- if (test_method == "z") "Pooled (Wald Z)" else "Pooled (KH)"
  if (nrow(selected)) for (i in seq_len(nrow(selected))) {
    z <- selected[i, , drop = FALSE]
    f <- long_weight[long_weight$HGNC_ID == z$HGNC_ID, , drop = FALSE]
    f <- f[match(studies, f$Dataset), , drop = FALSE]
    fd <- data.frame(Label = c(studies, pooled_label),
      Estimate = c(f$log2FC, z$Meta_log2FC),
      Low = c(f$log2FC - stats::qnorm(0.975) * f$SE, z$CI_low),
      High = c(f$log2FC + stats::qnorm(0.975) * f$SE, z$CI_high),
      Type = c(rep("Study", 4), "Pooled"),
      Weight_percent = c(f$Weight_percent, NA_real_),
      HGNC_ID = z$HGNC_ID, gene_symbol = z$gene_symbol)
    forest_data_list[[i]] <- fd
    fd$Label <- factor(fd$Label, levels = rev(c(studies, pooled_label)))
    forest <- ggplot2::ggplot(fd, ggplot2::aes(x = Estimate, y = Label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
        colour = "grey55", linewidth = 0.35) +
      ggplot2::geom_segment(ggplot2::aes(x = Low, xend = High, yend = Label),
        linewidth = 0.5, colour = "grey30") +
      ggplot2::geom_point(ggplot2::aes(shape = Type, colour = Type), size = 2.8) +
      ggplot2::scale_shape_manual(values = c(Study = 16, Pooled = 18)) +
      ggplot2::scale_colour_manual(values = c(Study = "#3B7EA1", Pooled = "#D55E00")) +
      ggplot2::labs(title = paste0(z$gene_symbol, " (", z$HGNC_ID, ")"),
        subtitle = paste0("Meta FDR = ", formatC(z$Meta_FDR, format = "g", digits = 3),
          "; I2 = ", round(z$I2_percent, 1), "%; direction support = ", z$Direction_support),
        x = expression(paste(log[2], " fold change (AD / control)")), y = NULL,
        caption = paste0("Study bars: estimate +/- 1.96 SE; pooled bar: 95% confidence interval.\n",
          if (z$Signature_candidate) "Meets the candidate signature criteria." else
            "Exploratory display; does not meet the candidate signature criteria.")) +
      paper() + ggplot2::theme(legend.position = "none")
    safe_name <- gsub("[^A-Za-z0-9_.-]", "_",
      paste(z$gene_symbol, z$HGNC_ID, sep = "_"))
    save_plot(forest, file.path("forest", safe_name), 7, 4.2)
  }
  if (length(forest_data_list)) csv(do.call(rbind, forest_data_list), "forest_plot_data.csv")
  plot_exports <- do.call(rbind, plot_log)
  csv(plot_exports, "plot_export_log.csv")

  write_log("[7/7] Saving interpretation rules and final reproducibility records ...")
  rules <- c(
    "Primary analysis: four-study, intercept-only random-effects meta-analysis of log2FC.",
    "Input effects are AD minus Control/NC, paired with the native audited model SE.",
    "Each common gene contributes exactly four estimates. No study DEG filter is applied.",
    "Sampling variance = SE squared; random-effects weight = 1/(SE squared + tau squared).",
    paste0("Tau squared: REML. Hypothesis testing and pooled 95% CIs: ",
      ifelse(test_method == "z", "Standard Wald Z-test (standard normal distribution)",
             "Knapp-Hartung adjustment (t-distribution, df=3)"), "."),
    "BH adjustment uses the entire planned common-gene family before effect-size/direction filters.",
    "Failed fits remain in meta_all_genes.csv with NA P/FDR and count toward BH family size.",
    "Meta_FDR_significant identifies statistical evidence for a nonzero average study effect.",
    paste0("Candidate signature: Meta_FDR < ", fdr_cutoff, ", abs(pooled log2FC) >= ",
      fc_cutoff, ", at least ", min_direction_support, "/4 matching signs."),
    "A zero study effect does not support either direction.",
    "Single-study statistical significance is not required for directional support.",
    "Study FDR values are preserved from each original differential analysis.",
    "Opposite-direction study FDR evidence is flagged, not a primary candidate exclusion.",
    "signature_no_observed_opposite_study_FDR.csv is an explicitly labelled additional subset.",
    "The 4/4 subset is reported separately and is not required for the primary candidate list.",
    "I2, tau2 and maximum study weight are reported; no cutoff on these is used.",
    "Direction support and study weights do not prove that a result is free of single-study influence.",
    "Forest study intervals use normal approximations from the input SE.",
    "Forest pooled intervals use the model confidence interval.",
    "Study estimates are treated as independent; this script cannot verify donor overlap.",
    "Different tissue sampling and covariate models remain sources of between-study differences.",
    "Documentation: https://wviechtb.github.io/metafor/reference/rma.uni.html")
  writeLines(rules, file.path(out, "README_analysis_rules.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  n_bad_plots <- sum(plot_exports$Status == "failed")
  write_log("Finished. Candidate signature: ", nrow(candidates), " genes (",
    sum(candidates$Pooled_direction == "UP"), " UP; ",
    sum(candidates$Pooled_direction == "DOWN"), " DOWN).")
  write_log("Models without estimable inference: ", sum(r$Fit_status != "ok"),
    "; failed figure exports: ", n_bad_plots, ".")
  write_log("Results: ", out)
  if (n_bad_plots) warning("Numerical results are saved, but some figures failed. ",
    "See plot_export_log.csv in ", out, call. = FALSE)
  invisible(list(output_dir = out, all_genes = r, meta_significant = significant,
    signature = candidates, summary = summary, direction_summary = direction_summary,
    study_weights = weight_summary, plot_exports = plot_exports,
    settings = settings, volcano = volcano))
}

# 1. 检查输入文件在不在
print(file.exists("D:/My Desktop/drug/3.Signature_integration/four_studies_common_meta_input_long.csv"))

# 2. 检查函数是否已成功加载进内存
print(exists("run_four_study_meta_analysis"))

# 3. 用单线程安全模式启动
if (exists("run_four_study_meta_analysis")) {
  meta_results <- run_four_study_meta_analysis(n_workers = 1L)
}