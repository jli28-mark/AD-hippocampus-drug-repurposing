# GSE173955: NCBI counts -> age/sex-adjusted DESeq2 -> all results, SE and plots.
# Uses NCBI GRCh38.p13 counts, not the author's hg19/StringTie result table.
# Positive log2FC = AD minus Control. Wald MLE effects; no lfcShrink applied.
# source("D:/My Desktop/drug/1.DEG/GSE173955/GSE173955_DESeq2_age_sex.R", encoding="UTF-8")
# GSE173955_results <- run_GSE173955()

run_GSE173955 <- function(
    base_dir = "D:/My Desktop/drug/1.DEG/GSE173955",
    n_workers = 2L, fdr_cutoff = 0.05, fc_cutoff = 0.58, png_dpi = 600L) {
  stopifnot(fdr_cutoff > 0, fdr_cutoff < 1, fc_cutoff > 0, n_workers >= 1)
  message("[1/8] Checking packages ...")
  cran <- c("ggplot2", "ggrepel")
  bioc <- c("DESeq2", "BiocParallel")
  for (p in cran) if (!requireNamespace(p, quietly = TRUE))
    utils::install.packages(p, repos = "https://cloud.r-project.org")
  missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  for (p in c(cran, bioc)) if (!requireNamespace(p, quietly = TRUE)) stop("Package unavailable: ", p)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(base_dir, "NCBI_counts_download")
  dir.create(cache, showWarnings = FALSE)
  parent <- file.path(base_dir, "adjusted_results")
  dir.create(parent, showWarnings = FALSE)
  out <- tempfile(paste0("DESeq2_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, " [", out, "]", call. = FALSE)
  old_timeout <- getOption("timeout"); options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  count_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
    "&acc=GSE173955&file=GSE173955_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  annot_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
    "&file=Human.GRCh38.p13.annot.tsv.gz")
  download_cached <- function(url, name) {
    dest <- file.path(cache, name)
    if (!file.exists(dest)) {
      tmp <- tempfile("download_", tmpdir = cache)
      on.exit(unlink(tmp))
      status <- utils::download.file(url, tmp, mode = "wb")
      check(status == 0 && file.info(tmp)$size > 100, paste("Download failed:", name))
      con <- gzfile(tmp, "rt")
      first <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
      check(length(first) == 1L && grepl("GeneID", first, fixed = TRUE),
        paste("Downloaded file is not the expected GeneID table:", name))
      check(file.rename(tmp, dest), paste("Cannot save", name))
    }
    dest
  }
  read_gz <- function(path, ...) {
    con <- gzfile(path, "rt"); on.exit(close(con))
    utils::read.delim(con, check.names = FALSE, quote = "", comment.char = "",
      stringsAsFactors = FALSE, ...)
  }
  message("[2/8] Downloading/reusing NCBI counts and gene annotations ...")
  count_file <- download_cached(count_url, "GSE173955_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  # NCBI can return an HTML verification page for the shared annotation URL.
  # Try a series-scoped URL, then the shared URL. Never accept HTML as annotation.
  annotation_file <- NULL
  annotation_source <- NULL
  annotation_urls <- c(paste0(annot_url, "&acc=GSE173955"), annot_url)
  for (u in annotation_urls) {
    candidate <- tryCatch(download_cached(u, "Human.GRCh38.p13.annot.tsv.gz"),
      error = function(e) {message("NCBI annotation attempt failed: ", conditionMessage(e)); NULL})
    if (!is.null(candidate)) {
      annotation_file <- candidate
      annotation_source <- u
      break
    }
  }
  if (is.null(annotation_file)) {
    message("Using org.Hs.eg.db as the GeneID-to-symbol annotation fallback. Statistical models are unchanged.")
    fallback_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")
    missing_fallback <- fallback_pkgs[!vapply(fallback_pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_fallback)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(missing_fallback, ask = FALSE, update = FALSE)
    }
    for (p in fallback_pkgs) check(requireNamespace(p, quietly = TRUE), paste("Missing fallback package", p))
    db <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene_ids <- AnnotationDbi::keys(db, keytype = "ENTREZID")
    gene_symbols <- AnnotationDbi::mapIds(db, keys = gene_ids, keytype = "ENTREZID",
      column = "SYMBOL", multiVals = "asNA")
    fallback <- data.frame(GeneID = gene_ids, Symbol = unname(gene_symbols[gene_ids]))
    annotation_file <- file.path(out, "orgHs_GeneID_annotation_fallback.tsv.gz")
    con <- gzfile(annotation_file, "wt")
    tryCatch(utils::write.table(fallback, con, sep = "\t", quote = FALSE,
      row.names = FALSE, na = ""), finally = close(con))
    annotation_source <- paste0("Bioconductor org.Hs.eg.db v", utils::packageVersion("org.Hs.eg.db"),
      "; current ENTREZID-to-SYMBOL mapping; ambiguous symbols set to NA; not the original NCBI GRCh38.p13 annotation")
  }
  csv(data.frame(File = c(count_file, annotation_file), URL = c(count_url, annot_url),
    Actual_source = c(count_url, annotation_source),
    MD5 = unname(tools::md5sum(c(count_file, annotation_file)))), "input_sources_checksums.csv")
  writeLines(annotation_source, file.path(out, "gene_annotation_source.txt"))

  message("[3/8] Matching the 18 clinical records to GSM columns ...")
  meta <- data.frame(
    Sample_ID = paste0("GSM", 5283449:5283466),
    Group = c(rep("AD", 8), rep("Control", 10)),
    Sex = c("F", "F", "F", "F", "M", "M", "M", "F",
            "F", "F", "F", "M", "M", "F", "F", "M", "M", "M"),
    Age = c(88, 95, 95, 100, 99, 83, 90, 84, 87, 80, 84, 77, 55, 72, 78, 83, 80, 74),
    stringsAsFactors = FALSE)
  # Explicit GSM mapping from the supplied GEO sample records; no title substring matching.
  meta$Group <- factor(meta$Group, levels = c("Control", "AD"))
  meta$Sex <- factor(meta$Sex, levels = c("F", "M"))
  meta$Age_c <- meta$Age - mean(meta$Age)
  rownames(meta) <- meta$Sample_ID
  csv(meta, "clinical_used.csv")
  cohort <- do.call(rbind, lapply(levels(meta$Group), function(g) {
    z <- meta[meta$Group == g, ]
    data.frame(Group = g, N = nrow(z), Female = sum(z$Sex == "F"), Male = sum(z$Sex == "M"),
      Age_mean = mean(z$Age), Age_SD = stats::sd(z$Age), Age_min = min(z$Age), Age_max = max(z$Age))
  }))
  csv(cohort, "cohort_summary.csv"); print(cohort)
  design <- stats::model.matrix(~ Age_c + Sex + Group, meta)
  check(qr(design)$rank == ncol(design) && nrow(design) > ncol(design), "Invalid model design.")
  csv(data.frame(Sample_ID = meta$Sample_ID, design, check.names = FALSE), "design_matrix.csv")
  tab <- read_gz(count_file)
  check("GeneID" %in% names(tab) && !anyDuplicated(names(tab)), "Invalid count-table header.")
  available <- setdiff(names(tab), "GeneID")
  check(all(meta$Sample_ID %in% available), "Some clinical GSM IDs are absent from counts.")
  csv(data.frame(Count_column = available, Selected = available %in% meta$Sample_ID), "sample_selection.csv")
  ids <- as.character(tab$GeneID)
  check(!anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids), "Missing or duplicated GeneID.")
  check(all(vapply(tab[meta$Sample_ID], is.numeric, logical(1))), "Nonnumeric count columns.")
  counts <- as.matrix(tab[meta$Sample_ID]); rownames(counts) <- ids; rm(tab)
  check(all(is.finite(counts)) && all(counts >= 0) &&
    all(abs(counts - round(counts)) < 1e-8) && max(counts) <= .Machine$integer.max,
    "Expected finite nonnegative integer counts, not normalized expression.")
  storage.mode(counts) <- "integer"
  check(identical(colnames(counts), rownames(meta)), "Sample order mismatch.")
  check(all(colSums(counts) > 0), "A sample has zero total counts.")
  csv(data.frame(Sample_ID = meta$Sample_ID, Group = meta$Group,
    Total_counts = colSums(counts), Detected_genes = colSums(counts > 0)), "sample_count_QC.csv")

  message("[4/8] Filtering low counts and fitting age/sex-adjusted DESeq2 Wald models ...")
  min_samples <- min(table(meta$Group)) # 8 samples, mirroring the supplied GEO2R filter
  keep <- rowSums(counts >= 10L) >= min_samples
  csv(data.frame(GeneID = ids, Samples_count_ge10 = rowSums(counts >= 10L), Retained = keep), "gene_filter_audit.csv")
  check(sum(keep) > 100L, "Too few genes after count filtering.")
  message("Counts: ", nrow(counts), " genes x ", ncol(counts), " samples; retained ", sum(keep), " genes.")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
    colData = meta, design = ~ Age_c + Sex + Group)
  if (n_workers > 1L) {
    # Windows-compatible PSOCK workers; use n_workers=1 if memory is limited.
    bp <- BiocParallel::SnowParam(workers = as.integer(n_workers), type = "SOCK",
      progressbar = TRUE, RNGseed = 173955)
    on.exit(BiocParallel::bpstop(bp), add = TRUE)
  } else bp <- BiocParallel::SerialParam()
  dds <- DESeq2::DESeq(dds, test = "Wald", sfType = "poscount", betaPrior = FALSE,
    parallel = n_workers > 1L, BPPARAM = bp)
  r <- DESeq2::results(dds, contrast = c("Group", "AD", "Control"),
    alpha = fdr_cutoff, pAdjustMethod = "BH", independentFiltering = TRUE,
    parallel = n_workers > 1L, BPPARAM = bp)
  print(DESeq2::resultsNames(dds))
  # No top-N truncation; no effect shrinkage. Missing P/padj values remain explicit.
  res <- data.frame(GeneID = rownames(r), as.data.frame(r), row.names = NULL)

  message("[5/8] Annotating and saving all genes, model SE and DEG lists ...")
  annot <- read_gz(annotation_file, colClasses = "character")
  check(all(c("GeneID", "Symbol") %in% names(annot)), "Annotation lacks GeneID or Symbol.")
  annot$GeneID <- as.character(annot$GeneID)
  check(!anyDuplicated(annot$GeneID), "Repeated GeneID in NCBI annotation; inspect before joining.")
  j <- match(res$GeneID, annot$GeneID)
  res$gene_symbol <- as.character(annot$Symbol[j])
  if ("Description" %in% names(annot)) res$Description <- annot$Description[j]
  res$logFC <- res$log2FoldChange
  res$SE <- res$lfcSE
  res$P.Value <- res$pvalue
  res$adj.P.Val <- res$padj
  res$variance <- res$SE^2
  res$Meta_eligible <- is.finite(res$logFC) & is.finite(res$SE) & res$SE > 0 & is.finite(res$pvalue)
  res$status <- "Not_DEG"
  res$status[is.na(res$padj)] <- "Not_evaluable_FDR"
  res$status[!is.na(res$padj) & res$padj < fdr_cutoff & res$logFC >= fc_cutoff] <- "UP"
  res$status[!is.na(res$padj) & res$padj < fdr_cutoff & res$logFC <= -fc_cutoff] <- "DOWN"
  res <- res[order(res$padj, res$pvalue, na.last = TRUE), ]
  csv(res, "GSE173955_all_genes.csv")
  csv(res[res$status %in% c("UP", "DOWN"), ], "GSE173955_DEG.csv")
  # GeneID is the source identity. Do not silently average rows sharing a symbol.
  csv(res[res$Meta_eligible, c("GeneID", "gene_symbol", "logFC", "SE", "variance", "pvalue", "padj")],
    "GSE173955_meta_input_unharmonized.csv")
  for (g in c("UP", "DOWN")) csv(res[res$status == g, ], paste0("GSE173955_", g, ".csv"))
  summary <- as.data.frame(table(factor(res$status, levels = c("DOWN", "Not_DEG", "Not_evaluable_FDR", "UP"))))
  names(summary) <- c("Status", "N_genes")
  csv(summary, "DEG_summary.csv"); print(summary)
  csv(data.frame(GeneID = res$GeneID[!res$Meta_eligible],
    logFC = res$logFC[!res$Meta_eligible], SE = res$SE[!res$Meta_eligible],
    P_value = res$pvalue[!res$Meta_eligible]), "meta_ineligible_records.csv")

  save_plot <- function(p, stem, width = 6.5, height = 5.3) {
    if (capabilities("cairo")) ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), p,
      width = width, height = height, device = grDevices::cairo_pdf, bg = "white") else
      ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height,
        device = grDevices::pdf, useDingbats = FALSE, bg = "white")
    ggplot2::ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
      dpi = png_dpi, bg = "white")
  }
  paper <- function() ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8))
  message("[6/8] Saving volcano PDF and PNG ...")
  v <- res[is.finite(res$logFC) & is.finite(res$padj), ]
  v$y <- -log10(pmax(v$padj, .Machine$double.xmin))
  v$status <- factor(v$status, levels = c("DOWN", "Not_DEG", "UP"))
  v$Label <- ifelse(is.na(v$gene_symbol) | !nzchar(v$gene_symbol), v$GeneID, v$gene_symbol)
  v <- v[order(v$status != "Not_DEG"), ]
  labels <- do.call(rbind, lapply(c("UP", "DOWN"), function(g) {
    z <- v[v$status == g, ]; z <- z[order(z$padj, -abs(z$logFC)), ]
    utils::head(z[!duplicated(z$Label), ], 5L)
  }))
  sizes <- table(factor(v$status, levels = c("DOWN", "Not_DEG", "UP")))
  p <- ggplot2::ggplot(v, ggplot2::aes(logFC, y)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", linewidth = 0.35, colour = "grey50") +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", linewidth = 0.35, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = status), size = 1.1, alpha = 0.65) +
    ggplot2::scale_colour_manual(values = c(DOWN = "#3B7EA1", Not_DEG = "#B8BDC3", UP = "#D55E00"),
      breaks = names(sizes), labels = paste0(c("Down", "Not DEG", "Up"), " (n = ", as.integer(sizes), ")"), drop = FALSE) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = Label), seed = 173955,
      size = 3, box.padding = 0.4, max.overlaps = Inf, min.segment.length = 0,
      segment.color = "grey55", segment.size = 0.3, show.legend = FALSE) +
    ggplot2::labs(title = "GSE173955 | Hippocampus", subtitle = "AD versus control; adjusted for age and sex",
      x = expression(log[2]~"fold change (AD / Control)"), y = expression(-log[10]~FDR),
      caption = paste(sum(!is.finite(res$padj) | !is.finite(res$logFC)), "rows without finite FDR/log2FC omitted from plot; retained in all-gene table.")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2.2))) + paper()
  save_plot(p, "GSE173955_volcano")

  message("[7/8] Saving sample QC, dispersion and PCA ...")
  # VST is used only for QC; count-model fitting above always uses untransformed counts.
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  e <- SummarizedExperiment::assay(vsd)
  qc <- data.frame(Sample_ID = meta$Sample_ID, Group = meta$Group, Age = meta$Age, Sex = meta$Sex,
    Size_factor = DESeq2::sizeFactors(dds))
  csv(qc, "sample_normalization_QC.csv")
  grDevices::pdf(file.path(out, "QC_dispersions.pdf"), width = 6, height = 5, useDingbats = FALSE)
  tryCatch(DESeq2::plotDispEsts(dds, main = "GSE173955"), finally = grDevices::dev.off())
  grDevices::pdf(file.path(out, "QC_VST_distributions.pdf"), width = 10, height = 5, useDingbats = FALSE)
  tryCatch({
    graphics::par(mar = c(7, 4, 3, 1))
    graphics::boxplot(as.data.frame(e), outline = FALSE, las = 2, cex.axis = 0.65,
      col = ifelse(meta$Group == "AD", "#D55E00", "#0072B2"), ylab = "VST expression",
      main = "GSE173955 | Sample distributions")
  }, finally = grDevices::dev.off())
  idx <- utils::head(order(apply(e, 1, stats::var), decreasing = TRUE), 2000L)
  pc <- stats::prcomp(t(e[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  scores <- data.frame(meta, PC1 = pc$x[, 1], PC2 = pc$x[, 2], row.names = NULL)
  csv(scores, "PCA_scores.csv")
  csv(data.frame(PC = seq_along(pct), Percent_variance = pct), "PCA_variance.csv")
  pp <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, colour = Group, shape = Sex)) +
    ggplot2::geom_point(size = 2.7) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), size = 2.4, seed = 173955,
      max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(Control = "#0072B2", AD = "#D55E00")) +
    ggplot2::labs(title = "GSE173955 | PCA", subtitle = "VST; up to 2,000 variable genes; covariates not regressed out",
      x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) + paper()
  save_plot(pp, "QC_PCA", width = 7, height = 5.5)
  utils::write.csv(stats::cor(e), file.path(out, "sample_correlations.csv"))

  message("[8/8] Saving model and analysis provenance ...")
  settings <- list(source = "NCBI GRCh38.p13 RNA-seq counts", model = "~ Age_c + Sex + Group",
    contrast = c("Group", "AD", "Control"), test = "Wald", sfType = "poscount", betaPrior = FALSE,
    low_count_filter = paste("counts >=10 in at least", min_samples, "samples"),
    independentFiltering = TRUE, pAdjustMethod = "BH", FDR = fdr_cutoff, abs_log2FC = fc_cutoff,
    lfcShrink = FALSE, n_workers = n_workers, protein_coding_filter = FALSE)
  dput(settings, file = file.path(out, "analysis_settings.R"))
  writeLines(c("GSE173955 reanalysis: all 8 AD and 10 Control samples; continuous age and sex adjusted.",
    "Sample metadata supplied from GEO records; no age cutoff; the 55-year Control and oldest donors retained.",
    "NCBI GRCh38.p13 counts differ from the author's hg19/StringTie pipeline; this is a new analysis.",
    "DESeq2 Wald maximum-likelihood log2FC and corresponding lfcSE; no LRT-based SE approximation.",
    "Default DESeq2 outlier handling and independent filtering retained; missing P/padj remain in full output.",
    "Rows with valid effect, positive finite SE and finite P can enter later meta input even if padj is NA.",
    "Do not combine these SE with the author's old logFC; use the paired estimates in this output.",
    "GeneID is the row identity. Symbols are annotations; duplicate symbols are not silently collapsed.",
    "No protein-coding filter, HGNC harmonization, meta-analysis or RRA is performed here.",
    "For downstream work, harmonize this new all-gene table and preserve its SE column; do not mix old/new study estimates.",
    "No RIN/PMI/batch covariates supplied. No donor-block term: one sample per clinical record.",
    "PCA and sample correlations use VST values only; they are not covariate-adjusted residuals.",
    paste("Actual gene annotation source:", annotation_source)),
    file.path(out, "analysis_record.txt"))
  all_pkgs <- c(cran, bioc)
  if (grepl("org.Hs.eg.db", annotation_source, fixed = TRUE))
    all_pkgs <- unique(c(all_pkgs, "AnnotationDbi", "org.Hs.eg.db"))
  csv(data.frame(Package = all_pkgs, Version = vapply(all_pkgs,
    function(z) as.character(utils::packageVersion(z)), character(1))), "package_versions.csv")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  saveRDS(list(dds = dds, metadata = meta, design = design, results = res, settings = settings),
    file.path(out, "GSE173955_analysis_objects.rds"))
  message("Finished. Results: ", normalizePath(out, winslash = "/"))
  invisible(list(output_dir = normalizePath(out, winslash = "/"), all_genes = res,
    DEG_summary = summary, clinical = meta, volcano = p))
}
