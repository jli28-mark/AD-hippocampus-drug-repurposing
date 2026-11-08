# GSE5281 hippocampal main analysis only: 13 Control and 10 AD.
# Input: positive, unlogged GEO probe-set expression intensities (not CEL files).
# log2 once -> joint quantile normalization -> unique annotation -> mean probes/gene.
# Model: Group + centered continuous age + sex. Contrast: AD minus Control.
# No sample-exclusion sensitivity analysis; no protein-coding filter.
# GSM119628: use 85 YEARS as confirmed by the user, retaining the GEO unit discrepancy.
# source("D:/My Desktop/drug/1.DEG/GSE5281/GSE5281_DEG_age_sex.R", encoding="UTF-8")
# GSE5281_results <- run_GSE5281()

run_GSE5281 <- function(
    base_dir = "D:/My Desktop/drug/1.DEG/GSE5281",
    clinical_file = NULL,
    matrix_file = file.path(base_dir, "GSE5281_series_matrix.txt"),
    fdr_cutoff = 0.05, fc_cutoff = 0.58, n_threads = 4L, png_dpi = 600L) {
  stopifnot(fdr_cutoff > 0, fdr_cutoff < 1, fc_cutoff > 0)
  message("[1/8] Checking packages ...")
  cran <- c("data.table", "readxl", "ggplot2", "ggrepel")
  bioc <- c("limma", "AnnotationDbi", "hgu133plus2.db")
  for (p in cran) if (!requireNamespace(p, quietly = TRUE))
    utils::install.packages(p, repos = "https://cloud.r-project.org")
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
  for (p in c(cran, bioc)) if (!requireNamespace(p, quietly = TRUE)) stop("Package unavailable: ", p)
  if (is.null(clinical_file)) {
    candidates <- file.path(base_dir, c("sample.xlsx", "sample(1).xlsx"))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates)) stop("Set clinical_file to the clinical XLSX path.")
    clinical_file <- candidates[1L]
  }
  if (!file.exists(clinical_file) || !file.exists(matrix_file)) stop("Input file not found.")
  parent <- file.path(base_dir, "adjusted_results")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  csv <- function(d, name) utils::write.csv(d, file.path(out, name), row.names = FALSE, na = "")
  check <- function(ok, text) if (!isTRUE(ok)) stop(text, " [", out, "]", call. = FALSE)
  clean <- function(z) trimws(as.character(z))
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(n_threads)
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)

  message("[2/8] Reading the 23-sample clinical cohort ...")
  meta <- as.data.frame(readxl::read_excel(clinical_file, col_types = "text"))
  names(meta) <- clean(names(meta))
  check(!anyDuplicated(names(meta)), "Duplicated clinical column names.")
  required <- c("Sample_ID", "Group", "Sex", "Age", "Region")
  check(all(required %in% names(meta)), paste("Required columns:", paste(required, collapse = ", ")))
  meta <- meta[!apply(meta, 1L, function(z) all(is.na(z) | clean(z) == "")), , drop = FALSE]
  csv(meta, "clinical_input_snapshot.csv")
  meta$Sample_ID <- clean(meta$Sample_ID)
  gm <- c(control = "Control", ctrl = "Control", nc = "Control", ad = "AD")
  sm <- c(male = "M", m = "M", female = "F", f = "F")
  meta$Group <- unname(gm[tolower(clean(meta$Group))])
  meta$Sex <- unname(sm[tolower(clean(meta$Sex))])
  meta$Age <- suppressWarnings(as.numeric(clean(meta$Age)))
  check(!anyNA(meta[, required]) && !anyDuplicated(meta$Sample_ID) &&
    all(is.finite(meta$Age) & meta$Age > 0 & meta$Age <= 120) &&
    all(tolower(clean(meta$Region)) == "hippocampus"), "Invalid clinical fields.")
  controls <- paste0("GSM", 119628:119640)
  cases <- paste0("GSM", 238799:238808)
  check(setequal(meta$Sample_ID[meta$Group == "Control"], controls) &&
    setequal(meta$Sample_ID[meta$Group == "AD"], cases), "Unexpected hippocampal cohort or diagnosis.")
  meta <- meta[match(c(controls, cases), meta$Sample_ID), , drop = FALSE]
  check(meta$Age[meta$Sample_ID == "GSM119628"] == 85,
        "GSM119628 must have Age=85 (years), as confirmed for this analysis.")
  meta$Group <- factor(meta$Group, levels = c("Control", "AD"))
  meta$Sex <- factor(meta$Sex, levels = c("M", "F"))
  meta$Age_c <- meta$Age - mean(meta$Age)
  rownames(meta) <- meta$Sample_ID
  csv(meta, "clinical_used.csv")
  csv(data.frame(Sample_ID = "GSM119628", GEO_original = "Age: 85 days",
    Analysis_age_years = 85, Reason = "User-confirmed interpretation: years; presumed GEO unit error; no independent documentary confirmation claimed."),
    "age_annotation_decision.csv")
  cohort <- do.call(rbind, lapply(levels(meta$Group), function(g) {
    z <- meta[meta$Group == g, , drop = FALSE]
    data.frame(Group = g, N = nrow(z), Female = sum(z$Sex == "F"), Male = sum(z$Sex == "M"),
      Age_mean = mean(z$Age), Age_SD = stats::sd(z$Age), Age_min = min(z$Age), Age_max = max(z$Age))
  }))
  print(cohort); csv(cohort, "cohort_summary.csv")

  message("[3/8] Reading selected probe intensities; log2 once and joint quantile normalization ...")
  # Stream through the file to support extracted tables and full GEO preambles/footers.
  locate_table <- function(path) {
    con <- file(path, "rt"); on.exit(close(con))
    n <- 0L; header <- NA_integer_; footer <- NA_integer_
    repeat {
      s <- readLines(con, n = 1000L, warn = FALSE)
      if (!length(s)) break
      if (is.na(header)) {
        j <- which(grepl('^\\s*"?ID_REF"?\\s*\\t', sub("^\ufeff", "", s), perl = TRUE))
        if (length(j)) header <- n + j[1L]
      }
      j <- which(grepl('^\\s*"?!series_matrix_table_end"?\\s*$', s, perl = TRUE))
      if (length(j)) footer <- n + j[1L]
      n <- n + length(s)
      if (!is.na(footer)) break
    }
    check(!is.na(header), "No ID_REF table header found.")
    list(skip = header - 1L, nrows = if (is.na(footer)) n - header else footer - header - 1L)
  }
  loc <- locate_table(matrix_file)
  tab <- data.table::fread(matrix_file, skip = loc$skip, nrows = loc$nrows,
    select = c("ID_REF", meta$Sample_ID), colClasses = list(character = "ID_REF"),
    check.names = FALSE, data.table = FALSE)
  check(all(c("ID_REF", meta$Sample_ID) %in% names(tab)), "Missing expression columns.")
  probes <- clean(tab$ID_REF)
  check(!anyNA(probes) && all(nzchar(probes)) && !anyDuplicated(probes), "Invalid probe IDs.")
  check(all(vapply(tab[meta$Sample_ID], is.numeric, logical(1))), "Expression columns must be numeric.")
  raw <- as.matrix(tab[meta$Sample_ID]); rownames(raw) <- probes; rm(tab)
  check(all(is.finite(raw)) && all(raw > 0), "Expected finite, strictly positive linear intensities.")
  check(identical(colnames(raw), rownames(meta)), "Sample order mismatch.")
  before <- log2(raw)
  norm <- limma::normalizeBetweenArrays(before, method = "quantile")
  csv(data.frame(Sample_ID = meta$Sample_ID, Raw_median = apply(raw, 2L, stats::median),
    Log2_median = apply(before, 2L, stats::median), Normalized_median = apply(norm, 2L, stats::median)),
    "normalization_summary.csv")

  message("[4/8] Unique probe-to-gene mapping and averaging normalized log2 probes ...")
  db <- getExportedValue("hgu133plus2.db", "hgu133plus2.db")
  valid_keys <- intersect(probes, AnnotationDbi::keys(db, keytype = "PROBEID"))
  mapping <- unique(AnnotationDbi::select(db, keys = valid_keys, keytype = "PROBEID",
                                        columns = c("SYMBOL", "ENTREZID")))
  csv(mapping, "probe_mapping_raw.csv")
  good <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL) &
                    !is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE]
  nc <- table(good$PROBEID)
  single <- names(nc)[nc == 1L]
  good <- good[good$PROBEID %in% single & !grepl("^AFFX", good$PROBEID, ignore.case = TRUE), , drop = FALSE]
  good <- good[order(match(good$PROBEID, probes)), , drop = FALSE]
  check(nrow(good) > 100L, "Too few uniquely annotated probes.")
  csv(good, "probe_mapping_retained.csv")
  csv(data.frame(Probe_ID = probes, Retained = probes %in% good$PROBEID), "probe_filter_audit.csv")
  gene <- limma::avereps(norm[match(good$PROBEID, probes), , drop = FALSE], ID = good$SYMBOL)
  variance <- apply(gene, 1L, stats::var)
  retained <- is.finite(variance) & variance > 0
  csv(data.frame(Gene = rownames(gene), Variance = variance, Retained = retained), "gene_filter_audit.csv")
  gene <- gene[retained, , drop = FALSE]
  message("Retained ", nrow(good), " probes; testing ", nrow(gene), " genes. No protein-coding filter.")

  message("[5/8] Main limma model: AD versus Control, adjusted for age and sex ...")
  design <- stats::model.matrix(~ Group + Age_c + Sex, data = meta)
  check(qr(design)$rank == ncol(design) && nrow(design) > ncol(design), "Design is rank deficient or has no residual degrees of freedom.")
  csv(data.frame(Sample_ID = rownames(design), design, check.names = FALSE), "design_matrix.csv")
  fit <- limma::lmFit(gene, design)
  fit <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
  res <- limma::topTable(fit, coef = "GroupAD", number = Inf, adjust.method = "BH", sort.by = "P")
  res <- data.frame(Gene = rownames(res), res, row.names = NULL, check.names = FALSE)
  res$log2FC <- res$logFC; res$FDR <- res$adj.P.Val
  res$Regulation <- ifelse(!is.na(res$FDR) & res$FDR < fdr_cutoff & res$logFC >= fc_cutoff, "UP",
    ifelse(!is.na(res$FDR) & res$FDR < fdr_cutoff & res$logFC <= -fc_cutoff, "DOWN", "Not_DEG"))
  res$N_probes <- as.integer(table(good$SYMBOL)[res$Gene])
  csv(res, "GSE5281_all_genes.csv")
  csv(res[res$Regulation != "Not_DEG", ], "GSE5281_DEG.csv")
  for (g in c("UP", "DOWN")) {
    csv(res[res$Regulation == g, ], paste0("GSE5281_", g, ".csv"))
    writeLines(res$Gene[res$Regulation == g], file.path(out, paste0("GSE5281_", g, "_genes.txt")))
  }
  summary <- data.frame(Regulation = c("DOWN", "Not_DEG", "UP"),
    N = as.integer(table(factor(res$Regulation, levels = c("DOWN", "Not_DEG", "UP")))))
  print(summary); csv(summary, "DEG_summary.csv")

  save_plot <- function(p, stem, width = 6.5, height = 5.3) {
    if (capabilities("cairo")) ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), p,
      width = width, height = height, device = grDevices::cairo_pdf, bg = "white") else
      ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height,
        device = grDevices::pdf, useDingbats = FALSE, bg = "white")
    ggplot2::ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
      dpi = png_dpi, bg = "white")
  }
  theme_paper <- function() ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8), plot.margin = ggplot2::margin(10, 12, 10, 10))
  message("[6/8] Saving volcano PDF and 600-dpi PNG ...")
  vp <- res[is.finite(res$logFC) & is.finite(res$FDR) & res$FDR >= 0 & res$FDR <= 1, ]
  vp$y <- -log10(pmax(vp$FDR, .Machine$double.xmin))
  vp$Regulation <- factor(vp$Regulation, levels = c("DOWN", "Not_DEG", "UP"))
  vp <- vp[order(vp$Regulation != "Not_DEG"), ] # grey points underneath significant points
  labels <- do.call(rbind, lapply(c("UP", "DOWN"), function(g) {
    z <- vp[vp$Regulation == g, ]; head(z[order(z$FDR, -abs(z$logFC)), ], 5L)
  }))
  legend_labels <- setNames(paste0(c("Down", "Not DEG", "Up"), " (n = ", summary$N, ")"), summary$Regulation)
  p <- ggplot2::ggplot(vp, ggplot2::aes(logFC, y)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", linewidth = 0.35, colour = "grey50") +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", linewidth = 0.35, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = Regulation), size = 1.1, alpha = 0.65, shape = 16) +
    ggplot2::scale_colour_manual(values = c(DOWN = "#3B7EA1", Not_DEG = "#B8BDC3", UP = "#D55E00"),
      breaks = names(legend_labels), labels = legend_labels, drop = FALSE) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = Gene), seed = 5281,
      size = 3, box.padding = 0.4, max.overlaps = Inf, min.segment.length = 0,
      segment.colour = "grey55", segment.linewidth = 0.3, show.legend = FALSE) +
    ggplot2::labs(title = "GSE5281 | Hippocampus", subtitle = "AD versus control; adjusted for age and sex",
      x = expression(log[2]~"fold change (AD / Control)"), y = expression(-log[10]~FDR)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.2, alpha = 1))) + theme_paper()
  save_plot(p, "GSE5281_volcano")

  message("[7/8] Saving normalization, PCA and sample-correlation QC ...")
  cols <- ifelse(meta$Group == "AD", "#D55E00", "#0072B2")
  grDevices::pdf(file.path(out, "QC_distributions.pdf"), width = 11, height = 8, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(2, 1), mar = c(7, 4, 3, 1))
    graphics::boxplot(as.data.frame(before), outline = FALSE, las = 2, cex.axis = 0.6, col = cols,
      main = "Before joint quantile normalization", ylab = "Log2 expression")
    graphics::boxplot(as.data.frame(norm), outline = FALSE, las = 2, cex.axis = 0.6, col = cols,
      main = "After joint quantile normalization", ylab = "Log2 expression")
  }, finally = grDevices::dev.off())
  idx <- head(order(apply(gene, 1L, stats::var), decreasing = TRUE), 2000L)
  pc <- stats::prcomp(t(gene[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2) # all components in the denominator
  scores <- data.frame(Sample_ID = meta$Sample_ID, Group = meta$Group, Sex = meta$Sex,
    Age = meta$Age, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  csv(scores, "PCA_scores.csv")
  csv(data.frame(PC = seq_along(pct), Percent_variance = pct), "PCA_variance.csv")
  pp <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, colour = Group, shape = Sex)) +
    ggplot2::geom_point(size = 2.7) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), size = 2.4, seed = 5281, max.overlaps = Inf,
      show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(Control = "#0072B2", AD = "#D55E00")) +
    ggplot2::labs(title = "GSE5281 | PCA", subtitle = "Up to 2,000 most variable genes; no covariates regressed out",
      x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) + theme_paper()
  save_plot(pp, "QC_PCA", width = 7, height = 5.5)
  cor_mat <- stats::cor(gene, method = "pearson")
  utils::write.csv(cor_mat, file.path(out, "sample_correlations.csv"))
  grDevices::pdf(file.path(out, "QC_sample_correlation.pdf"), width = 9, height = 9, useDingbats = FALSE)
  tryCatch(stats::heatmap(cor_mat, symm = TRUE, scale = "none", margins = c(9, 9),
    ColSideColors = cols, RowSideColors = cols, cexRow = 0.65, cexCol = 0.65,
    main = "Pearson correlation of normalized gene expression"), finally = grDevices::dev.off())

  message("[8/8] Saving objects and analysis record ...")
  settings <- list(contrast = "AD-Control", model = "~Group + Age_c + Sex", FDR = fdr_cutoff,
    abs_log2FC = fc_cutoff, normalization = "log2 once; joint probe-level quantile normalization of all 23 samples",
    annotation = "hgu133plus2.db; one complete SYMBOL/ENTREZID mapping per probe; mean log2 probes per symbol",
    eBayes_trend = FALSE, eBayes_robust = FALSE, protein_coding_filter = FALSE,
    sample_exclusions = "None; all 13 Control and 10 AD retained", age_GSM119628 = "85 years; user-confirmed unit interpretation",
    n_threads = n_threads, png_dpi = png_dpi)
  dput(settings, file = file.path(out, "analysis_settings.R"))
  versions <- data.frame(Package = c(cran, bioc), Version = vapply(c(cran, bioc), function(z) as.character(utils::packageVersion(z)), character(1)))
  csv(versions, "package_versions.csv")
  csv(data.frame(File = c(clinical_file, matrix_file), MD5 = unname(tools::md5sum(c(clinical_file, matrix_file)))), "input_checksums.csv")
  writeLines(c("GSE5281 main analysis only: 13 Control / 10 AD hippocampal neuronal samples.",
    "Starting material: submitted positive linear probe-set expression values, not CEL files.",
    "GEO AD records specify MAS5.0 target scaling to 150; control processing fields are missing.",
    "Additional quantile normalization is a step of this reanalysis, not asserted to be the author's procedure.",
    "No background correction, RMA, TMM, voom, ComBat, or sample-exclusion sensitivity analysis was applied.",
    "GSM119628: GEO says 85 days; analysis uses 85 years at user confirmation. Original file is unchanged.",
    "GSM119631 (102 years) retained. No age threshold or age winsorization.",
    "Cell-type wording is not used as a covariate; these are selected neuronal samples, not unsorted bulk tissue.",
    "No batch variable supplied; GSM prefixes are not treated as experimental batches.",
    "Normalization assumes most probes do not undergo a global unidirectional shift.",
    "All tested genes enter one BH correction; DEG labels also require the specified estimated effect size.",
    "Full gene-level results include noncoding genes; downstream protein-coding harmonization is separate.",
    paste("Tested genes:", nrow(gene)), paste("Residual degrees of freedom:", nrow(design) - ncol(design))),
    file.path(out, "analysis_record.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  saveRDS(list(meta = meta, design = design, gene_expression = gene, fit = fit, results = res,
    settings = settings, probe_mapping = good, PCA = pc), file.path(out, "GSE5281_analysis_objects.rds"))
  message("Finished. Results: ", normalizePath(out, winslash = "/"))
  invisible(list(output_dir = normalizePath(out, winslash = "/"), all_genes = res,
    DEG_summary = summary, clinical = meta, settings = settings, volcano = p))
}
