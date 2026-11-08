# =============================================================================
# GSE48350 hippocampus: AD versus Control
# Prepared for sexage.xlsx uploaded on 2026-10-21.
# Run in R/RStudio:
# source("D:/My Desktop/drug/1.DEG/GSE48350/GSE48350_DEG_age_sex.R",
#        encoding = "UTF-8")
# GSE48350_results <- run_GSE48350()
#
# Main: age >= 60 in BOTH groups, ~ 0 + Group + Age_c + Sex.
# Sensitivity: all 62 hippocampal samples, same model and tested genes.
# Input: previously normalised, positive LINEAR expression ratios, NOT counts.
# Log2 once. No extra RMA, quantile normalisation, TMM, voom, or pseudocount.
# Probe handling: exclude ambiguous mappings; mean of log2 probes per symbol.
# DEGs: BH FDR < 0.05 AND absolute estimated log2FC >= 0.58.
# Positive logFC means higher expression in AD.
# The fold-change rule is a post-test filter, not a TREAT hypothesis test.
#
# Sources:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11882
# https://bioconductor.org/packages/limma/
# https://bioconductor.org/packages/AnnotationDbi/
# =============================================================================

run_GSE48350 <- function(
    base_dir = "D:/My Desktop/drug/1.DEG/GSE48350",
    clinical_file = file.path(base_dir, "sexage.xlsx"),
    matrix_file = file.path(base_dir, "GSE48350_series_matrix.txt"),
    output_dir = NULL,
    min_age = 60,
    fdr_cutoff = 0.05,
    log2fc_cutoff = 0.58,
    run_all_age_sensitivity = TRUE,
    install_missing = TRUE,
    n_threads = 4L,
    n_labels_per_direction = 5L,
    png_dpi = 600L,
    save_rds = TRUE) {

  stopifnot(length(min_age) == 1L, is.finite(min_age),
            length(fdr_cutoff) == 1L, fdr_cutoff > 0, fdr_cutoff < 1,
            length(log2fc_cutoff) == 1L, log2fc_cutoff > 0,
            is.finite(n_threads), n_threads >= 1,
            is.finite(n_labels_per_direction), n_labels_per_direction >= 0)
  if (!file.exists(clinical_file)) stop("Clinical file not found: ", clinical_file)
  if (!file.exists(matrix_file)) stop("Expression file not found: ", matrix_file)
  if (is.null(output_dir)) {
    output_dir <- file.path(base_dir, "adjusted_results",
                            paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  if (dir.exists(output_dir)) {
    stop("Output directory already exists; choose a new output_dir: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  log_file <- file.path(output_dir, "run_log.txt")
  say <- function(...) {
    x <- paste0(..., collapse = "")
    message(x)
    cat(x, "\n", file = log_file, append = TRUE)
  }
  check <- function(ok, message) {
    if (!isTRUE(ok)) {
      cat("ERROR: ", message, "\n", file = log_file, append = TRUE)
      stop(message, "\nSee: ", output_dir, call. = FALSE)
    }
  }
  csv <- function(x, name, folder = output_dir) {
    utils::write.csv(x, file.path(folder, name), row.names = FALSE,
                     na = "", fileEncoding = "UTF-8")
  }
  clean <- function(x) trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))
  check_values <- function(bad, details, name, message) {
    bad <- is.na(bad) | bad
    if (any(bad)) {
      csv(details[bad, , drop = FALSE], paste0("ERROR_", name, ".csv"))
      check(FALSE, paste0(message, " See ERROR_", name, ".csv."))
    }
  }

  say("[1/8] Checking packages ...")
  cran <- c("readxl", "data.table", "ggplot2", "ggrepel")
  bioc <- c("limma", "AnnotationDbi", "hgu133plus2.db")
  needed <- function(x) x[!vapply(x, requireNamespace, logical(1), quietly = TRUE)]
  if (length(needed(c(cran, bioc))) && !install_missing) {
    check(FALSE, paste("Missing packages:", paste(needed(c(cran, bioc)), collapse = ", ")))
  }
  if (length(needed(cran))) {
    utils::install.packages(needed(cran), repos = "https://cloud.r-project.org")
  }
  if (length(needed(bioc))) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(needed(bioc), ask = FALSE, update = FALSE)
  }
  check(length(needed(c(cran, bioc))) == 0L, "Some required packages are unavailable.")
  n_threads <- as.integer(n_threads)
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(n_threads)
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)

  # All plots are vector PDF plus 600-dpi PNG by default. No point subsampling.
  save_plot <- function(p, stem, folder, width = 6.5, height = 5.3) {
    if (capabilities("cairo")) {
      ggplot2::ggsave(file.path(folder, paste0(stem, ".pdf")), plot = p,
                     width = width, height = height, units = "in",
                     device = grDevices::cairo_pdf, bg = "white")
    } else {
      ggplot2::ggsave(file.path(folder, paste0(stem, ".pdf")), plot = p,
                     width = width, height = height, units = "in",
                     device = grDevices::pdf, useDingbats = FALSE, bg = "white")
    }
    ggplot2::ggsave(file.path(folder, paste0(stem, ".png")), plot = p,
                   width = width, height = height, units = "in",
                   dpi = png_dpi, device = "png", bg = "white")
  }
  theme_paper <- function() {
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.35),
        axis.text = ggplot2::element_text(colour = "black"),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
        legend.position = "top", legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(10, 12, 10, 10))
  }

  say("[2/8] Reading and validating the revised clinical table ...")
  raw <- as.data.frame(readxl::read_excel(clinical_file, sheet = 1,
                                         col_types = "text", .name_repair = "minimal"))
  names(raw) <- clean(names(raw))
  required <- c("Group", "Accession", "Title", "Source name", "Individual",
                "Brain region", "Gender", "Characteristics")
  check(!anyDuplicated(names(raw)), "Clinical column names are duplicated.")
  check(all(required %in% names(raw)),
        paste("Required clinical columns:", paste(required, collapse = ", ")))
  # Ignore completely blank rows only. Missing AD Individual/Braak/MMSE fields
  # are NOT grounds to discard an otherwise complete age/sex/diagnosis record.
  blank_row <- apply(raw, 1L, function(z) all(is.na(z) | clean(z) == ""))
  raw <- raw[!blank_row, , drop = FALSE]
  check(nrow(raw) > 0L, "Clinical table has no data rows.")
  csv(raw, "clinical_input_snapshot.csv")
  group_key <- tolower(clean(raw$Group))
  group_map <- c(ctrl = "Control", control = "Control", nc = "Control", ad = "AD")
  sex_key <- tolower(clean(raw$Gender))
  sex_map <- c(female = "F", f = "F", male = "M", m = "M")
  age_text <- clean(raw$Characteristics)
  age_pattern <- "^age\\s*\\(yrs\\)\\s*:\\s*([0-9]+(?:\\.[0-9]+)?)\\s*$"
  age_ok <- !is.na(age_text) & grepl(age_pattern, age_text, perl = TRUE,
                                    ignore.case = TRUE)
  age <- rep(NA_real_, nrow(raw))
  age[age_ok] <- as.numeric(sub(age_pattern, "\\1", age_text[age_ok],
                                 perl = TRUE, ignore.case = TRUE))
  meta <- data.frame(
    Sample_ID = clean(raw$Accession), Title = clean(raw$Title),
    Group = unname(group_map[group_key]), Sex = unname(sex_map[sex_key]),
    Age = age, Region = tolower(clean(raw[["Brain region"]])),
    Individual_original = clean(raw$Individual), stringsAsFactors = FALSE)
  bad <- is.na(meta$Sample_ID) | !grepl("^GSM[0-9]+$", meta$Sample_ID) |
    is.na(meta$Group) | is.na(meta$Sex) | !is.finite(meta$Age) |
    meta$Age <= 0 | meta$Age > 120 | is.na(meta$Region) | meta$Region != "hippocampus"
  check_values(bad, cbind(meta, Age_original = age_text), "clinical_fields",
               "An accession, group, sex, age, or brain-region field is invalid.")
  check(!anyDuplicated(meta$Sample_ID), "Clinical accessions are duplicated.")

  parse_title <- function(s) {
    p_control <- "^hippocampus_(female|male)_([0-9]+(?:\\.[0-9]+)?)yrs_indiv([0-9]+(?:-[0-9]+)*)$"
    p_ad <- "^hippocampus_(female|male)_([0-9]+(?:\\.[0-9]+)?)_ad_([0-9]+)$"
    a <- regmatches(s, regexec(p_control, s, ignore.case = TRUE, perl = TRUE))[[1L]]
    b <- regmatches(s, regexec(p_ad, s, ignore.case = TRUE, perl = TRUE))[[1L]]
    if (length(a) == 4L) {
      data.frame(Title_Group = "Control", Title_Sex = unname(sex_map[tolower(a[2L])]),
                 Title_Age = as.numeric(a[3L]), Donor_label = a[4L])
    } else if (length(b) == 4L) {
      data.frame(Title_Group = "AD", Title_Sex = unname(sex_map[tolower(b[2L])]),
                 Title_Age = as.numeric(b[3L]), Donor_label = b[4L])
    } else {
      data.frame(Title_Group = NA_character_, Title_Sex = NA_character_,
                 Title_Age = NA_real_, Donor_label = NA_character_)
    }
  }
  check_values(is.na(meta$Title), meta, "missing_title", "A sample title is missing.")
  title_meta <- do.call(rbind, lapply(meta$Title, parse_title))
  title_bad <- is.na(title_meta$Title_Group) |
    is.na(title_meta$Title_Sex) | !is.finite(title_meta$Title_Age) |
    meta$Group != title_meta$Title_Group | meta$Sex != title_meta$Title_Sex |
    abs(meta$Age - title_meta$Title_Age) > 1e-8
  check_values(title_bad, cbind(meta, title_meta), "title_disagreement",
               "The title disagrees with the clinical group, sex, or age.")
  meta$Donor_label <- title_meta$Donor_label
  # Control and AD donor numbers come from different naming schemes.
  # Keep labels such as 16-07 and 40-04 intact; do not merge across groups.
  meta$Donor_ID <- paste(meta$Group, meta$Donor_label, sep = ":")
  have_individual <- !is.na(meta$Individual_original) & nzchar(meta$Individual_original)
  individual_id <- clean(sub(",.*$", "", meta$Individual_original))
  check_values(have_individual & individual_id != meta$Donor_label, meta,
               "individual_disagreement", "Individual ID and title donor ID disagree.")
  check(!anyDuplicated(meta$Donor_ID),
        "More than one hippocampal sample has the same within-group donor ID.")

  # Prespecified complete hippocampal cohort; catches the earlier wrong-region
  # accessions even if their hand-entered titles say hippocampus.
  expected_control <- paste0("GSM", c(
    300168, 300169, 300170, 300171, 300172, 300174, 300178, 300182,
    300185, 300187, 300190, 300193, 300197, 300205, 300209, 300215,
    300219, 300223, 300227, 300231, 300235, 300239, 300243, 300255,
    300262, 300268, 300272, 300276, 300280, 300286, 300290, 300294,
    300298, 300301, 300305, 300309, 300313, 300317, 300321, 300325,
    300329, 300333, 300339))
  expected_ad <- paste0("GSM", 1176211:1176229)
  expected_ids <- c(expected_control, expected_ad)
  bad_ids <- setdiff(meta$Sample_ID, expected_ids)
  missing_ids <- setdiff(expected_ids, meta$Sample_ID)
  if (length(bad_ids) || length(missing_ids)) {
    csv(data.frame(Sample_ID = c(bad_ids, missing_ids),
                   Problem = c(rep("Unexpected accession", length(bad_ids)),
                               rep("Missing accession", length(missing_ids)))),
         "ERROR_expected_cohort.csv")
    check(FALSE, "Clinical table differs from the audited 62-sample hippocampal cohort.")
  }
  expected_group <- ifelse(meta$Sample_ID %in% expected_ad, "AD", "Control")
  check_values(meta$Group != expected_group, meta, "accession_group",
               "A group label disagrees with the prespecified accession list.")
  meta$Group <- factor(meta$Group, levels = c("Control", "AD"))
  meta$Sex <- factor(meta$Sex, levels = c("F", "M"))
  rownames(meta) <- meta$Sample_ID
  meta$Included_main <- meta$Age >= min_age
  csv(meta, "clinical_validated_all62.csv")

  summarise_cohort <- function(m, analysis) {
    do.call(rbind, lapply(levels(m$Group), function(g) {
      z <- m[m$Group == g, , drop = FALSE]
      data.frame(Analysis = analysis, Group = g, Samples = nrow(z),
                 Donors = length(unique(z$Donor_ID)), Female = sum(z$Sex == "F"),
                 Male = sum(z$Sex == "M"), Age_mean = mean(z$Age),
                 Age_SD = stats::sd(z$Age), Age_min = min(z$Age), Age_max = max(z$Age))
    }))
  }
  main_meta <- meta[meta$Included_main, , drop = FALSE]
  check(all(table(main_meta$Group) >= 3L), "Too few samples in a main-analysis group.")
  main_tag <- paste0("main_age", format(min_age, trim = TRUE, scientific = FALSE))
  cohort_summary <- rbind(summarise_cohort(main_meta, main_tag),
                          summarise_cohort(meta, "all_ages"))
  csv(cohort_summary, "cohort_summary.csv")
  print(cohort_summary, row.names = FALSE)

  say("[3/8] Matching the expression header and taking log2 once ...")
  # This input starts with ID_REF and GSM columns. It is an extracted expression
  # table, not a complete GEO SOFT/series-matrix metadata document.
  header <- names(data.table::fread(matrix_file, nrows = 0L, check.names = FALSE,
                                    data.table = FALSE, showProgress = FALSE))
  check(length(header) > 1L && identical(header[1L], "ID_REF"),
        "Expression file must start with an ID_REF/GSM table header.")
  check(!anyDuplicated(header[-1L]), "Expression sample names are duplicated.")
  missing_expression <- setdiff(meta$Sample_ID, header[-1L])
  check(length(missing_expression) == 0L,
        paste("Clinical samples missing from expression:",
              paste(missing_expression, collapse = ", ")))
  csv(data.frame(Sample_ID = header[-1L],
                 Selected_hippocampus = header[-1L] %in% meta$Sample_ID),
       "expression_sample_selection.csv")
  tab <- data.table::fread(matrix_file, select = c("ID_REF", meta$Sample_ID),
                           colClasses = list(character = "ID_REF"),
                           check.names = FALSE, data.table = FALSE,
                           showProgress = interactive(), nThread = n_threads)
  probe_ids <- clean(tab$ID_REF)
  check(!anyNA(probe_ids) && all(nzchar(probe_ids)) && !anyDuplicated(probe_ids),
        "Expression probe IDs are missing, empty, or duplicated.")
  check(all(vapply(tab[meta$Sample_ID], is.numeric, logical(1))),
        "Some expression columns are not numeric.")
  x <- as.matrix(tab[, meta$Sample_ID, drop = FALSE])
  storage.mode(x) <- "double"
  rownames(x) <- probe_ids
  check(identical(colnames(x), rownames(meta)), "Expression and metadata orders disagree.")
  check(all(is.finite(x)), "Expression contains missing or nonfinite values.")
  check(all(x > 0), "Linear expression must be positive. Check the input scale; no flooring is applied.")
  input_quantiles <- stats::quantile(x, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1))
  csv(data.frame(Quantile = names(input_quantiles), Linear_value = unname(input_quantiles),
                 Log2_value = log2(unname(input_quantiles))), "input_scale_summary.csv")
  say("Selected matrix: ", nrow(x), " probes x ", ncol(x), " samples.")
  x <- log2(x)
  rm(tab)

  say("[4/8] Resolving probe mappings and averaging log2 probes per gene ...")
  db <- get("hgu133plus2.db", envir = asNamespace("hgu133plus2.db"))
  valid_keys <- intersect(probe_ids, AnnotationDbi::keys(db, keytype = "PROBEID"))
  check(length(valid_keys) > 0L, "No GPL570 probe IDs match hgu133plus2.db.")
  ann <- unique(AnnotationDbi::select(db, keys = valid_keys, keytype = "PROBEID",
                                      columns = c("SYMBOL", "ENTREZID")))
  csv(ann, "probe_mapping_raw.csv")
  ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL) &
               !is.na(ann$ENTREZID) & nzchar(ann$ENTREZID), , drop = FALSE]
  check(nrow(ann) > 0L, "No probes have both a symbol and an Entrez ID.")
  n_symbols <- vapply(split(ann$SYMBOL, ann$PROBEID), function(z) length(unique(z)), integer(1))
  n_entrez <- vapply(split(ann$ENTREZID, ann$PROBEID), function(z) length(unique(z)), integer(1))
  audit <- data.frame(PROBEID = probe_ids,
                       n_symbols = as.integer(n_symbols[probe_ids]),
                       n_entrez_ids = as.integer(n_entrez[probe_ids]),
                       stringsAsFactors = FALSE)
  audit$n_symbols[is.na(audit$n_symbols)] <- 0L
  audit$n_entrez_ids[is.na(audit$n_entrez_ids)] <- 0L
  audit$Decision <- "retained"
  audit$Decision[audit$n_symbols == 0L | audit$n_entrez_ids == 0L] <- "unmapped"
  audit$Decision[audit$n_symbols > 1L | audit$n_entrez_ids > 1L] <- "ambiguous_probe"
  audit$Decision[grepl("^AFFX", probe_ids, ignore.case = TRUE)] <- "array_control"
  mapping <- unique(ann[ann$PROBEID %in% audit$PROBEID[audit$Decision == "retained"],
                         c("PROBEID", "SYMBOL", "ENTREZID"), drop = FALSE])
  symbol_n_entrez <- vapply(split(mapping$ENTREZID, mapping$SYMBOL),
                             function(z) length(unique(z)), integer(1))
  ambiguous_symbols <- names(symbol_n_entrez)[symbol_n_entrez > 1L]
  ambiguous_symbol_probes <- mapping$PROBEID[mapping$SYMBOL %in% ambiguous_symbols]
  audit$Decision[audit$PROBEID %in% ambiguous_symbol_probes] <- "ambiguous_gene_symbol"
  mapping <- mapping[!mapping$SYMBOL %in% ambiguous_symbols, , drop = FALSE]
  check(nrow(mapping) > 0L && !anyDuplicated(mapping$PROBEID),
        "Retained probe-to-gene mappings are not unique.")
  csv(audit, "probe_mapping_audit.csv")
  csv(mapping, "probe_mapping_retained.csv")
  gene_all <- limma::avereps(x[mapping$PROBEID, , drop = FALSE], ID = mapping$SYMBOL)
  n_probe <- table(mapping$SYMBOL)
  gene_info <- unique(mapping[, c("SYMBOL", "ENTREZID")])
  gene_info$n_probes <- as.integer(n_probe[gene_info$SYMBOL])
  names(gene_info)[1L] <- "gene_symbol"
  gene_info <- gene_info[match(rownames(gene_all), gene_info$gene_symbol), , drop = FALSE]
  # No absolute-intensity or count-based filtering is appropriate for these ratios.
  # Remove numerically constant genes in the MAIN samples, then use this same
  # tested-gene universe in the all-age sensitivity analysis.
  main_variance <- apply(gene_all[, main_meta$Sample_ID, drop = FALSE], 1L, stats::var)
  tested <- is.finite(main_variance) & main_variance > 1e-12
  gene_info$Variance_main <- main_variance
  gene_info$Tested <- tested
  csv(gene_info, "gene_filter_audit.csv")
  check(sum(tested) >= 10L, "Too few nonconstant mapped genes remain.")
  gene_all <- gene_all[tested, , drop = FALSE]
  gene_info <- gene_info[tested, , drop = FALSE]
  check(!anyDuplicated(rownames(gene_all)), "Gene symbols are duplicated after averaging.")
  say("Retained ", nrow(mapping), " uniquely mapped probes; testing ",
      nrow(gene_all), " genes. No protein-coding filter was applied.")
  rm(x)

  fit_model <- function(m, tag) {
    folder <- file.path(output_dir, tag)
    dir.create(folder, showWarnings = FALSE)
    m$Age_c <- m$Age - mean(m$Age)
    e <- gene_all[, m$Sample_ID, drop = FALSE]
    check(identical(colnames(e), rownames(m)), paste(tag, "sample order mismatch."))
    design <- stats::model.matrix(~ 0 + Group + Age_c + Sex, data = m)
    check(qr(design)$rank == ncol(design), paste(tag, "design is not full rank."))
    check(nrow(design) > ncol(design) + 1L, paste(tag, "insufficient residual degrees of freedom."))
    fit0 <- limma::lmFit(e, design)
    contrast <- limma::makeContrasts(AD_vs_Control = GroupAD - GroupControl, levels = design)
    fit <- limma::eBayes(limma::contrasts.fit(fit0, contrast), trend = FALSE, robust = FALSE)
    result <- limma::topTable(fit, coef = "AD_vs_Control", number = Inf,
                               adjust.method = "BH", sort.by = "P", confint = 0.95)
    result$gene_symbol <- rownames(result)
    result$ENTREZID <- gene_info$ENTREZID[match(result$gene_symbol, gene_info$gene_symbol)]
    result$n_probes <- gene_info$n_probes[match(result$gene_symbol, gene_info$gene_symbol)]
    check(all(is.finite(result$logFC)) && all(is.finite(result$P.Value)) &&
            all(is.finite(result$adj.P.Val)), paste(tag, "nonfinite gene statistics."))
    check(nrow(result) == nrow(gene_all), paste(tag, "unexpected number of tested genes."))
    result$status <- "Not_DEG"
    result$status[result$adj.P.Val < fdr_cutoff & result$logFC >= log2fc_cutoff] <- "UP"
    result$status[result$adj.P.Val < fdr_cutoff & result$logFC <= -log2fc_cutoff] <- "DOWN"
    result <- result[, c("gene_symbol", "ENTREZID", "n_probes", "logFC", "CI.L", "CI.R",
                          "AveExpr", "t", "P.Value", "adj.P.Val", "B", "status")]
    csv(result, "AD_vs_Control_all_genes.csv", folder)
    csv(result[result$status != "Not_DEG", , drop = FALSE], "AD_vs_Control_DEGs.csv", folder)
    csv(result[result$status == "UP", , drop = FALSE], "AD_vs_Control_UP.csv", folder)
    csv(result[result$status == "DOWN", , drop = FALSE], "AD_vs_Control_DOWN.csv", folder)
    csv(m, "metadata_used.csv", folder)
    csv(data.frame(Sample_ID = rownames(design), design, check.names = FALSE), "design_matrix.csv", folder)
    summary <- data.frame(Analysis = tag, Control = sum(m$Group == "Control"),
                           AD = sum(m$Group == "AD"), Tested = nrow(result),
                           DOWN = sum(result$status == "DOWN"),
                           Not_DEG = sum(result$status == "Not_DEG"),
                           UP = sum(result$status == "UP"))
    csv(summary, "DEG_summary.csv", folder)
    say(tag, ": Control = ", summary$Control, ", AD = ", summary$AD,
        "; UP = ", summary$UP, ", DOWN = ", summary$DOWN)
    list(result = result, metadata = m, design = design, fit = fit,
         folder = folder, tag = tag, summary = summary)
  }

  say("[5/8] Main limma model: Group + age + sex ...")
  main <- fit_model(main_meta, main_tag)
  sensitivity <- NULL
  say("[6/8] All-age sensitivity analysis with the same tested genes ...")
  if (run_all_age_sensitivity) sensitivity <- fit_model(meta, "sensitivity_all_ages")
  models <- list(main = main)
  if (!is.null(sensitivity)) models$sensitivity <- sensitivity
  csv(do.call(rbind, lapply(models, function(z) z$summary)), "DEG_summary_all_models.csv")
  if (!is.null(sensitivity)) {
    a <- main$result[, c("gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "status")]
    b <- sensitivity$result[, c("gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "status")]
    comparison <- merge(a, b, by = "gene_symbol", suffixes = c("_main", "_all_ages"), sort = FALSE)
    comparison$Same_direction <- sign(comparison$logFC_main) == sign(comparison$logFC_all_ages)
    comparison$DEG_in_both <- comparison$status_main != "Not_DEG" &
      comparison$status_all_ages != "Not_DEG"
    csv(comparison, "main_vs_all_ages.csv")
  }

  say("[7/8] Saving publication-style volcano plots and QC ...")
  # Common axes across models include ALL tested genes; no fixed -4/+4 clipping.
  x_limit <- max(log2fc_cutoff * 2, max(vapply(models,
                  function(z) max(abs(z$result$logFC)), numeric(1))) * 1.12)
  y_limit <- max(-log10(fdr_cutoff) * 1.5, max(vapply(models,
                  function(z) max(-log10(pmax(z$result$adj.P.Val, .Machine$double.xmin))),
                  numeric(1))) * 1.12)
  make_volcano <- function(model) {
    d <- model$result
    d$neg_log10_fdr <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    d$status <- factor(d$status, levels = c("DOWN", "Not_DEG", "UP"))
    d <- d[order(d$status != "Not_DEG"), , drop = FALSE]
    pick <- function(status) {
      z <- d[d$status == status, , drop = FALSE]
      z <- z[order(z$adj.P.Val, z$P.Value, -abs(z$logFC), z$gene_symbol), , drop = FALSE]
      utils::head(z, n_labels_per_direction)
    }
    labels <- rbind(pick("UP"), pick("DOWN"))
    csv(labels, "volcano_labelled_genes.csv", model$folder)
    counts <- table(d$status)
    legend_labels <- c(DOWN = paste0("Down (", counts["DOWN"], ")"),
                        Not_DEG = paste0("Not DEG (", counts["Not_DEG"], ")"),
                        UP = paste0("Up (", counts["UP"], ")"))
    cohort <- if (identical(model$tag, main_tag)) {
      paste0("Age >= ", min_age, " years")
    } else "All ages"
    p <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = neg_log10_fdr)) +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
                           linewidth = 0.3, colour = "grey50") +
      ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
                           linetype = "dashed", linewidth = 0.3, colour = "grey50") +
      ggplot2::geom_point(ggplot2::aes(colour = status), shape = 16,
                           size = 0.85, alpha = 0.65) +
      ggplot2::scale_colour_manual(
        values = c(DOWN = "#0072B2", Not_DEG = "#B8B8B8", UP = "#D55E00"),
        breaks = c("DOWN", "Not_DEG", "UP"), labels = unname(legend_labels), drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        nrow = 1L, override.aes = list(alpha = 1, size = 2))) +
      ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit)) +
      ggplot2::labs(title = "GSE48350: Hippocampus",
        subtitle = paste0(cohort, "; AD (n = ", sum(model$metadata$Group == "AD"),
                           ") vs Control (n = ", sum(model$metadata$Group == "Control"),
                           "); adjusted for age and sex"),
        x = expression(log[2]~fold~change~(AD/Control)),
        y = expression(-log[10]~FDR)) + theme_paper()
    if (nrow(labels)) {
      p <- p + ggrepel::geom_text_repel(data = labels,
        ggplot2::aes(label = gene_symbol), seed = 42L, size = 2.8,
        box.padding = 0.45, point.padding = 0.15, min.segment.length = 0,
        segment.colour = "grey50", segment.size = 0.3,
        max.overlaps = Inf, max.time = 2, show.legend = FALSE)
    }
    save_plot(p, "Volcano_AD_vs_Control", model$folder, width = 7.0, height = 5.5)
    p
  }
  make_qc <- function(model) {
    m <- model$metadata
    e <- gene_all[, m$Sample_ID, drop = FALSE]
    variance <- apply(e, 1L, stats::var)
    top <- head(order(variance, decreasing = TRUE), min(2000L, nrow(e)))
    # No rank truncation: denominator includes all PCs of the selected genes.
    pc <- stats::prcomp(t(e[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
    percent <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    scores <- data.frame(m, PC1 = pc$x[, 1L], PC2 = pc$x[, 2L])
    csv(scores, "PCA_scores.csv", model$folder)
    csv(data.frame(PC = seq_along(percent), Variance_percent = percent),
         "PCA_variance_explained.csv", model$folder)
    csv(data.frame(gene_symbol = rownames(e)[top]), "PCA_selected_genes.csv", model$folder)
    p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = Group, shape = Sex)) +
      ggplot2::geom_point(size = 2.3, alpha = 0.85) +
      ggplot2::scale_colour_manual(values = c(Control = "#0072B2", AD = "#D55E00")) +
      ggplot2::labs(title = "GSE48350: PCA",
        subtitle = paste0("Unadjusted normalized log2 ratios; top ", length(top), " variable genes"),
        x = sprintf("PC1 (%.1f%%)", percent[1L]), y = sprintf("PC2 (%.1f%%)", percent[2L])) +
      theme_paper()
    save_plot(p, "QC_PCA", model$folder)
    p_age <- ggplot2::ggplot(m, ggplot2::aes(x = Group, y = Age, colour = Group)) +
      ggplot2::geom_boxplot(width = 0.4, outlier.shape = NA, linewidth = 0.35) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.09, seed = 42),
                           size = 1.7, alpha = 0.8) +
      ggplot2::scale_colour_manual(values = c(Control = "#0072B2", AD = "#D55E00")) +
      ggplot2::labs(title = "GSE48350: Age distribution", x = NULL, y = "Age (years)") +
      theme_paper() + ggplot2::theme(legend.position = "none")
    save_plot(p_age, "QC_Age", model$folder, width = 4.4, height = 4.4)
    q <- t(apply(e, 2L, stats::quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    qc <- data.frame(Sample_ID = rownames(q), Group = m$Group,
                       q025 = q[, 1L], q25 = q[, 2L], Median = q[, 3L],
                       q75 = q[, 4L], q975 = q[, 5L])
    csv(qc, "QC_sample_log2_quantiles.csv", model$folder)
    qc$Sample_ID <- factor(qc$Sample_ID, levels = m$Sample_ID)
    p_q <- ggplot2::ggplot(qc, ggplot2::aes(x = Sample_ID, colour = Group)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = q025, ymax = q975), linewidth = 0.3) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = q25, ymax = q75), linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = Median), size = 0.9) +
      ggplot2::scale_colour_manual(values = c(Control = "#0072B2", AD = "#D55E00")) +
      ggplot2::labs(title = "GSE48350: Sample distributions",
        subtitle = "Thin bars: 2.5-97.5%; thick bars: interquartile range; dots: median",
        x = NULL, y = "Normalized log2 ratio") + theme_paper() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    save_plot(p_q, "QC_Expression_distributions", model$folder, width = 9, height = 4.8)
    correlations <- stats::cor(e)
    csv(data.frame(Sample_ID = rownames(correlations), correlations, check.names = FALSE),
         "QC_sample_correlations.csv", model$folder)
    grDevices::pdf(file.path(model$folder, "QC_Limma_variance.pdf"), width = 6, height = 5,
                    useDingbats = FALSE)
    tryCatch(limma::plotSA(model$fit, main = "Residual variance vs mean log2 ratio"),
             finally = grDevices::dev.off())
  }
  main_volcano <- make_volcano(main)
  make_qc(main)
  if (!is.null(sensitivity)) {
    make_volcano(sensitivity)
    make_qc(sensitivity)
  }

  say("[8/8] Saving provenance, settings, and analysis objects ...")
  package_names <- c(cran, bioc)
  versions <- data.frame(Package = package_names,
                          Version = vapply(package_names,
                            function(p) as.character(utils::packageVersion(p)), character(1)))
  csv(versions, "package_versions.csv")
  utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  checksums <- tools::md5sum(c(clinical_file, matrix_file))
  csv(data.frame(File = names(checksums), MD5 = unname(checksums)), "input_checksums.csv")
  settings <- list(min_age = min_age, fdr_cutoff = fdr_cutoff,
                    log2fc_cutoff = log2fc_cutoff, input_scale = "positive linear ratio",
                    transformation = "log2 once; no pseudocount; no extra normalization",
                    probe_summary = "mean log2 values of probes with unique symbol/Entrez mappings",
                    model = "~ 0 + Group + Age_c + Sex", contrast = "GroupAD - GroupControl",
                    eBayes_trend = FALSE, eBayes_robust = FALSE,
                    protein_coding_filter = FALSE, constant_gene_tolerance = 1e-12,
                    sensitivity = run_all_age_sensitivity)
  dput(settings, file = file.path(output_dir, "analysis_settings.R"))
  methods <- c(
    "GSE48350: human hippocampal microarray analysis (GPL570).",
    sprintf("Main cohort: age >= %s; %d AD and %d controls.", min_age,
              sum(main_meta$Group == "AD"), sum(main_meta$Group == "Control")),
    "The authors' normalized positive linear ratios were log2-transformed once.",
    "No additional between-array normalization or count-based filtering was applied.",
    "Array-control and ambiguously mapped probes were excluded using hgu133plus2.db.",
    "Normalized log2 values from retained probes for the same gene symbol were averaged.",
    "Numerically constant genes in the main cohort were excluded before model fitting.",
    "A limma model included diagnosis, continuous age, and sex, with AD minus Control as the contrast.",
    "Ordinary empirical Bayes moderation was used (trend=FALSE; robust=FALSE).",
    sprintf("DEGs required BH FDR < %s and absolute estimated log2 fold change >= %s.",
              fdr_cutoff, log2fc_cutoff),
    "Positive and negative logFC indicate higher and lower expression in AD, respectively.",
    if (run_all_age_sensitivity) paste0(
      "An all-age sensitivity analysis used the same gene set and adjusted for age and sex; ",
      "BH correction was applied separately to each analysis.") else NULL,
    "All mapped gene types were retained; no protein-coding or iLINCS-coverage filter was applied.",
    "Braak stage, APOE genotype and MMSE were not adjustment variables.",
    "PCA used unadjusted normalized log2 ratios; samples were not automatically removed.",
    "RNA integrity, PMI, and sequencing/array-batch covariates were not available in this clinical input.",
    "Package versions, input checksums, metadata and design matrices accompany the outputs."
  )
  writeLines(methods, file.path(output_dir, "analysis_record.txt"))
  if (save_rds) {
    saveRDS(list(expression_log2 = gene_all, gene_annotation = gene_info,
                  metadata_all = meta, main = main, sensitivity = sensitivity,
                  settings = settings, package_versions = versions),
              file.path(output_dir, "GSE48350_analysis_objects.rds"))
  }
  say("Finished. Results: ", output_dir)
  print(main_volcano)
  invisible(list(output_dir = output_dir, main = main,
                  sensitivity = sensitivity, clinical_summary = cohort_summary))
}

# No analysis runs just by sourcing this file. Start it explicitly with:
# GSE48350_results <- run_GSE48350()

