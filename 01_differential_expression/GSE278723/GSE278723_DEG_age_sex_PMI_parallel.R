# GSE278723: complete AD-versus-NC analysis from the three original files.
# Save this script in BASE_DIR and run source(..., encoding = "UTF-8").
# Input files are never overwritten. Each run creates a new results directory.
# Main: region + age + sex + PMI + diagnosis, with donor_id as the block.
# Method: TMM -> voom -> two-pass duplicateCorrelation -> lmFit -> eBayes.
# This uses limma's shared donor correlation, not the dream method.

# -------------------- Settings: normally change only BASE_DIR -----------------
BASE_DIR <- "D:/My Desktop/drug/1.DEG/GSE278723"
COUNTS_FILE <- "GSE278723_Counts.xlsx"
CLINICAL_FILE <- "Supplementary Table 1.xlsx"
GROUP_FILE <- "GSE278723_sample_information.txt"
COUNTS_SHEET <- 1
COUNTS_HEADER_ROW <- 1L
GENE_COLUMN <- "Column1"
CLINICAL_SHEET <- 1
CLINICAL_HEADER_ROW <- 1L
WORKERS <- 4L                 # Upper limit; use 1L for serial execution
FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.58
MAKE_VOLCANO <- TRUE          # Requires ggplot2 and ggrepel
SAVE_ANALYSIS_OBJECTS <- TRUE
EXPECTED_TISSUES <- 110L
EXPECTED_DONORS_PER_GROUP <- 14L
REGIONS <- c("CA1", "CA2", "CA3", "CA4")
set.seed(278723)

# ------------------------------- Helpers ------------------------------------
clean_text <- function(x) {
  z <- trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))
  z[z %in% c("", "NA", "N/A", "NaN", "-", "--")] <- NA_character_
  z
}

donor_key <- function(x) {
  z <- clean_text(x)
  ok <- !is.na(z) & grepl("^[0-9]+(\\.0+)?$", z)
  z[!ok] <- NA_character_
  z[ok] <- sub("\\.0+$", "", z[ok])
  z[ok] <- sub("^0+(?=[0-9])", "", z[ok], perl = TRUE)
  z
}

diagnosis_key <- function(x) {
  z <- gsub("[[:space:]_-]", "", toupper(clean_text(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("NC", "CTRL", "CONTROL")] <- "NC"
  out[z %in% "AD"] <- "AD"
  out[z %in% "CAA"] <- "CAA"
  out[z %in% c("AD+CAA", "CAA+AD", "ADCAA", "CAAAD")] <- "AD+CAA"
  out
}

numeric_value <- function(x) {
  z <- clean_text(x)
  ok <- !is.na(z) & grepl("^[0-9]+(\\.[0-9]+)?$", z)
  ans <- rep(NA_real_, length(z))
  ans[ok] <- as.numeric(z[ok])
  ans
}

read_xlsx_text <- function(path, sheet, header_row) {
  d <- readxl::read_excel(path, sheet = sheet, skip = header_row - 1L,
                          col_types = "text", .name_repair = "minimal",
                          progress = FALSE)
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  names(d) <- trimws(gsub("\u00a0", " ", names(d), fixed = TRUE))
  d
}

require_columns <- function(d, needed, label) {
  if (anyDuplicated(names(d))) stop(label, ": duplicated column names.")
  missing <- setdiff(needed, names(d))
  if (length(missing)) {
    stop(label, ": missing columns: ", paste(missing, collapse = ", "),
         ". Actual columns: ", paste(names(d), collapse = " | "),
         ". Check the sheet and header-row settings at the top of this script.")
  }
}

# The worker is self-contained, so it also works in Windows PSOCK processes.
correlation_chunk <- function(chunk, design, block) {
  limma::duplicateCorrelation(chunk$E, design, block = block,
                             weights = chunk$weights)$atanh.correlations
}

parallel_correlation <- function(v, design, block, cl = NULL) {
  if (is.null(cl)) return(limma::duplicateCorrelation(v, design, block = block))
  n <- nrow(v$E)
  ids <- split(seq_len(n), rep_len(seq_len(min(n, 4L * length(cl))), n))
  chunks <- lapply(ids, function(i) {
    list(E = v$E[i, , drop = FALSE], weights = v$weights[i, , drop = FALSE])
  })
  pieces <- parallel::parLapplyLB(cl, chunks, correlation_chunk,
                                 design = design, block = block)
  z <- rep(NA_real_, n)
  for (j in seq_along(ids)) z[ids[[j]]] <- pieces[[j]]
  if (!any(is.finite(z)) || any(is.infinite(z))) {
    stop("Donor correlations could not be estimated reliably.")
  }
  # Exactly the limma pooling rule: pool all gene-level Fisher-z values first.
  # Averaging the separate chunk consensus correlations would be incorrect.
  rho <- tanh(mean(z, trim = 0.15, na.rm = TRUE))
  if (!is.finite(rho)) stop("The consensus donor correlation is not finite.")
  list(consensus.correlation = rho, atanh.correlations = z)
}

result_table <- function(fit) {
  tt <- limma::topTable(fit, coef = "diagnosisAD", number = Inf,
                        adjust.method = "BH", sort.by = "P")
  ans <- data.frame(gene_symbol = rownames(tt), tt, row.names = NULL,
                    check.names = FALSE)
  ok <- is.finite(ans$logFC) & is.finite(ans$P.Value) & is.finite(ans$adj.P.Val)
  ans$status <- ifelse(ok, "Not_DEG", "Not_estimable")
  ans$status[ok & ans$adj.P.Val < FDR_CUTOFF & ans$logFC >= LOG2FC_CUTOFF] <- "UP"
  ans$status[ok & ans$adj.P.Val < FDR_CUTOFF & ans$logFC <= -LOG2FC_CUTOFF] <- "DOWN"
  ans
}

save_result_tables <- function(res, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  put <- function(x, name) utils::write.csv(x, file.path(directory, name),
                                           row.names = FALSE, na = "NA")
  put(res, "AD_vs_NC_all_genes.csv")
  put(res[res$status %in% c("UP", "DOWN"), ], "AD_vs_NC_DEGs.csv")
  for (direction in c("UP", "DOWN")) {
    sub <- res[res$status == direction, , drop = FALSE]
    put(sub, paste0("AD_vs_NC_", direction, ".csv"))
    writeLines(sub$gene_symbol, file.path(directory, paste0(direction, "_genes.txt")))
  }
}

with_pdf <- function(path, width, height, expr) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

make_volcano <- function(res, out) {
  d <- res[is.finite(res$logFC) & is.finite(res$adj.P.Val), , drop = FALSE]
  d$neg_log10_fdr <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
  d$group <- factor(ifelse(d$status == "UP", "Upregulated",
                          ifelse(d$status == "DOWN", "Downregulated", "Not DEG")),
                     levels = c("Downregulated", "Not DEG", "Upregulated"))
  d <- d[order(d$group != "Not DEG"), , drop = FALSE]
  labels <- do.call(rbind, lapply(c("UP", "DOWN"), function(s) {
    x <- d[d$status == s, , drop = FALSE]
    x <- x[order(x$adj.P.Val, -abs(x$logFC), x$gene_symbol), , drop = FALSE]
    head(x, 5L)
  }))
  cols <- c(Downregulated = "#0072B2", "Not DEG" = "#B8B8B8", Upregulated = "#D55E00")
  sizes <- table(d$group)
  legend_labels <- sprintf("%s (%s)", levels(d$group), format(as.integer(sizes), big.mark = ","))
  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, neg_log10_fdr)) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 0.65, alpha = 0.65, shape = 16) +
    ggplot2::geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
                        linewidth = 0.3, colour = "grey45") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
                        linetype = "dashed", linewidth = 0.3, colour = "grey45") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, labels = legend_labels) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(title = "GSE278723", subtitle = "AD vs NC; age, sex and PMI adjusted",
                  x = expression(log[2]~fold~change~(AD/NC)),
                  y = expression(-log[10]~FDR), colour = NULL) +
    ggplot2::theme_classic(base_size = 8, base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                    plot.subtitle = ggplot2::element_text(size = 7.5),
                    axis.line = ggplot2::element_line(linewidth = 0.35),
                    axis.ticks = ggplot2::element_line(linewidth = 0.3),
                    legend.position = "bottom", legend.text = ggplot2::element_text(size = 7),
                    legend.key.width = grid::unit(0.3, "cm"),
                    plot.margin = ggplot2::margin(8, 10, 6, 8)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2,
                         override.aes = list(size = 1.8, alpha = 1)))
  if (nrow(labels)) {
    p <- p + ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = gene_symbol),
             seed = 278723, size = 2.4, min.segment.length = 0,
             box.padding = 0.35, point.padding = 0.2, segment.size = 0.25,
             segment.colour = "grey45", max.overlaps = Inf, show.legend = FALSE)
  }
  ggplot2::ggsave(file.path(out, "Volcano_AD_vs_NC.pdf"), p,
                  device = "pdf", width = 120, height = 115, units = "mm",
                  bg = "white", useDingbats = FALSE)
  ggplot2::ggsave(file.path(out, "Volcano_AD_vs_NC.png"), p,
                  device = "png", width = 120, height = 115, units = "mm",
                  dpi = 600, bg = "white")
}

# ---------------------------- Complete analysis ------------------------------
run_GSE278723 <- function() {
  if (!dir.exists(BASE_DIR)) stop("BASE_DIR does not exist: ", BASE_DIR)
  if (!is.finite(WORKERS) || WORKERS < 1) stop("WORKERS must be >= 1.")
  if (FDR_CUTOFF <= 0 || FDR_CUTOFF >= 1 || LOG2FC_CUTOFF <= 0) stop("Invalid thresholds.")
  started <- proc.time()[["elapsed"]]
  out <- file.path(BASE_DIR, "adjusted_results", format(Sys.time(), "run_%Y%m%d_%H%M%S"))
  if (dir.exists(out)) stop("Output directory already exists; run again in a moment.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "NA")
  checks <- data.frame(check = character(), status = character(), detail = character())
  check <- function(name, condition, detail) {
    checks <<- rbind(checks, data.frame(check = name,
                     status = if (isTRUE(condition)) "PASS" else "FAIL", detail = detail))
    csv(checks, "validation_results.csv")
    if (!isTRUE(condition)) stop(name, ": ", detail, " [see ", out, "]")
  }
  cl <- NULL
  finished <- FALSE
  stage <- "initialization"
  on.exit({
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    csv(checks, "validation_results.csv")
    writeLines(c(if (finished) "COMPLETE" else paste("STOPPED during", stage),
                 "", capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))
    if (!finished) message("Analysis stopped. Reports already written to: ", out)
  }, add = TRUE)

  stage <- "package installation"
  message("[1/10] Checking packages; the first installation can take time ...")
  cran_repo <- unname(getOption("repos")["CRAN"])
  if (length(cran_repo) != 1L || is.na(cran_repo) || cran_repo == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  cran <- c("readxl", "statmod", if (MAKE_VOLCANO) c("ggplot2", "ggrepel"))
  absent <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) utils::install.packages(absent, repos = "https://cloud.r-project.org")
  bioc <- c("edgeR", "limma")
  absent <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(absent, ask = FALSE, update = FALSE)
  }
  pkgs <- c(cran, bioc)
  check("required_packages", all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)),
        paste(pkgs, collapse = ", "))

  stage <- "clinical metadata"
  message("[2/10] Reading diagnoses and donor-level age, sex and PMI ...")
  files <- file.path(BASE_DIR, c(COUNTS_FILE, CLINICAL_FILE, GROUP_FILE))
  inventory <- data.frame(file = basename(files), exists = file.exists(files),
                          bytes = file.info(files)$size)
  csv(inventory, "input_files.csv")
  check("input_files_exist", all(inventory$exists), paste(basename(files), collapse = "; "))
  csv(data.frame(file = basename(files), md5 = unname(tools::md5sum(files))), "input_MD5.csv")
  info <- utils::read.delim(files[3], colClasses = "character", check.names = FALSE,
                           quote = "", comment.char = "")
  names(info) <- trimws(sub("^\ufeff", "", names(info)))
  require_columns(info, c("Sample ID", "Group"), "Sample information TXT")
  info <- data.frame(donor_id = donor_key(info[["Sample ID"]]),
                      group_original = clean_text(info$Group),
                      diagnosis = diagnosis_key(info$Group))
  csv(info, "group_file_cleaned.csv")
  check("group_file_IDs", !anyNA(info$donor_id) && !anyDuplicated(info$donor_id),
        "Donor IDs must be numeric and unique after normalizing leading zeros.")
  check("group_file_labels", !anyNA(info$diagnosis), "Allowed groups: NC, AD, CAA, AD+CAA.")

  clinical <- read_xlsx_text(files[2], CLINICAL_SHEET, CLINICAL_HEADER_ROW)
  require_columns(clinical, c("Sample ID", "Sex", "Age", "Group", "PMI(h)"), "Clinical XLSX")
  clinical$donor_id <- donor_key(clinical[["Sample ID"]])
  clinical$diagnosis_clinical <- diagnosis_key(clinical$Group)
  clinical$age <- numeric_value(clinical$Age)
  clinical$PMI <- numeric_value(clinical[["PMI(h)"]])
  sex_raw <- toupper(clean_text(clinical$Sex))
  clinical$sex <- unname(c(M = "M", MALE = "M", F = "F", FEMALE = "F")[sex_raw])
  csv(clinical, "clinical_file_cleaned.csv")
  check("clinical_IDs", !anyNA(clinical$donor_id) && !anyDuplicated(clinical$donor_id),
        "Clinical donor IDs must be numeric and unique. Do not match by row order.")
  mismatch <- merge(info, clinical[c("donor_id", "diagnosis_clinical")],
                     by = "donor_id", all = TRUE)
  mismatch$groups_agree <- !is.na(mismatch$diagnosis) &
    !is.na(mismatch$diagnosis_clinical) & mismatch$diagnosis == mismatch$diagnosis_clinical
  csv(mismatch, "diagnosis_crosscheck.csv")

  stage <- "count matrix and sample selection"
  message("[3/10] Reading the integer count matrix and matching sample IDs ...")
  header <- readxl::read_excel(files[1], sheet = COUNTS_SHEET,
                skip = COUNTS_HEADER_ROW - 1L, n_max = 0,
                .name_repair = "minimal", progress = FALSE)
  nms <- trimws(names(header))
  check("counts_header_unique", !anyDuplicated(nms) && !anyNA(nms) && all(nzchar(nms)),
        "Every count-matrix column needs a unique name.")
  check("gene_column", sum(nms == GENE_COLUMN) == 1L,
        paste0("Expected gene column: ", GENE_COLUMN, "; actual first column: ", nms[1]))
  gene_col <- which(nms == GENE_COLUMN)
  count_cols <- setdiff(seq_along(nms), gene_col)
  sample_names <- nms[count_cols]
  types <- rep("numeric", length(nms)); types[gene_col] <- "text"
  dat <- readxl::read_excel(files[1], sheet = COUNTS_SHEET,
               skip = COUNTS_HEADER_ROW - 1L, col_types = types,
               .name_repair = "minimal", progress = FALSE)
  genes <- clean_text(dat[[gene_col]])
  m <- as.matrix(as.data.frame(dat)[, count_cols, drop = FALSE])
  storage.mode(m) <- "double"
  dimnames(m) <- list(genes, sample_names)
  rm(dat)
  csv(data.frame(gene_symbol = genes, missing_id = is.na(genes),
                   duplicate_id = duplicated(genes) | duplicated(genes, fromLast = TRUE)),
      "gene_ID_audit.csv")
  check("gene_IDs", !anyNA(genes) && !anyDuplicated(genes),
        "Gene symbols must be nonmissing and unique; do not sum duplicates blindly.")
  check("finite_counts", all(is.finite(m)), paste("Nonfinite entries:", sum(!is.finite(m))))
  check("nonnegative_counts", all(m >= 0), "Counts must be nonnegative.")
  check("integer_counts", all(abs(m - round(m)) < 1e-8),
        "Use the raw gene counts. TPM/FPKM or log values must not be rounded into counts.")
  # Integer values alone do not establish provenance: use the author's count file.
  rx <- "^([0-9]+)\\.(CA[1-4]|NX|EC)$"
  valid_names <- grepl(rx, toupper(sample_names))
  csv(data.frame(sample_id = sample_names, valid_name = valid_names), "sample_name_audit.csv")
  check("sample_name_format", all(valid_names), "Expected names such as 01.CA1 or 15.NX.")
  meta <- data.frame(sample_id = sample_names,
                     donor_id = donor_key(sub(rx, "\\1", toupper(sample_names))),
                     region = sub(rx, "\\2", toupper(sample_names)))
  meta$diagnosis <- info$diagnosis[match(meta$donor_id, info$donor_id)]
  meta$included <- meta$diagnosis %in% c("NC", "AD") & meta$region %in% REGIONS
  meta$reason <- ifelse(meta$included, "Included",
                 ifelse(is.na(meta$diagnosis), "Unmatched donor ID; not reassigned",
                 ifelse(!meta$diagnosis %in% c("NC", "AD"), "CAA or AD+CAA", "Outside CA1-CA4")))
  csv(meta, "sample_inclusion_audit.csv")
  m <- m[, meta$included, drop = FALSE]
  meta <- meta[meta$included, c("sample_id", "donor_id", "region", "diagnosis")]
  ci <- match(meta$donor_id, clinical$donor_id)
  meta$diagnosis_clinical <- clinical$diagnosis_clinical[ci]
  meta$age <- clinical$age[ci]
  meta$sex <- clinical$sex[ci]
  meta$PMI <- clinical$PMI[ci]
  # Keep pathology scores for auditing only; no pathology-based reclassification.
  for (col in intersect(c("ABCscore", "CAAscore"), names(clinical))) meta[[col]] <- clinical[[col]][ci]
  rownames(meta) <- meta$sample_id
  csv(meta, "metadata_before_validation.csv")
  check("selected_donor_linkage", !anyNA(ci) && !anyNA(meta$diagnosis_clinical) &&
          all(meta$diagnosis == meta$diagnosis_clinical),
        "Every selected donor must match the clinical file and have consistent diagnoses.")
  check("selected_covariates", all(is.finite(meta$age)) && all(meta$age > 0 & meta$age < 130) &&
          !anyNA(meta$sex) && all(is.finite(meta$PMI)) && all(meta$PMI >= 0),
        "Every selected donor needs a numeric age, M/F sex and nonnegative numeric PMI(h). No silent deletion.")
  meta$diagnosis <- factor(meta$diagnosis, levels = c("NC", "AD"))
  meta$region <- factor(meta$region, levels = REGIONS)
  meta$sex <- factor(meta$sex, levels = c("F", "M"))
  meta$donor_id <- factor(meta$donor_id)
  donors <- unique(meta[c("donor_id", "diagnosis", "age", "sex", "PMI")])
  coverage <- with(meta, table(donor_id, region))
  composition <- with(meta, table(region, diagnosis))
  csv(data.frame(region = rownames(composition), as.data.frame.matrix(composition)), "sample_numbers_by_region.csv")
  csv(as.data.frame(coverage), "donor_region_coverage.csv")
  csv(donors, "donor_metadata.csv")
  check("one_sample_per_donor_region", all(coverage <= 1L), "Technical duplicates require separate handling.")
  check("donor_consistency", !anyDuplicated(donors$donor_id), "Each donor must have a single diagnosis and covariate record.")
  check("sample_alignment", identical(colnames(m), rownames(meta)), "Count columns and metadata rows must match exactly.")
  expected <- matrix(c(12L, 14L, 14L, 14L, 14L, 14L, 14L, 14L), nrow = 4L)
  check("expected_cohort", ncol(m) == EXPECTED_TISSUES && nrow(donors) == 28L &&
          all(table(donors$diagnosis) == EXPECTED_DONORS_PER_GROUP) &&
          setequal(as.character(donors$donor_id), as.character(1:28)) &&
          all(composition == expected),
        "Reference: 110 tissues, 14 NC + 14 AD donors; CA1 NC=12/AD=14, other regions NC=14/AD=14.")
  expected_coverage <- matrix(1L, nrow = 28L, ncol = 4L,
                               dimnames = list(as.character(1:28), REGIONS))
  expected_coverage[c("3", "14"), "CA1"] <- 0L
  check("expected_missing_subregions", all(coverage == expected_coverage[rownames(coverage), REGIONS]),
        "Only donor 03 and donor 14 should lack CA1; other selected donor-region combinations should be present.")
  check("nonempty_samples", all(colSums(m) > 0), "Each selected sample must contain counts.")
  meta$age_c <- meta$age - mean(donors$age)
  meta$PMI_c <- meta$PMI - mean(donors$PMI)
  main_formula <- ~ region + age_c + sex + PMI_c + diagnosis
  sensitivity_formula <- ~ region + age_c + sex + diagnosis
  regional_formula <- ~ age_c + sex + PMI_c + diagnosis
  validate_design <- function(design, label) {
    check(paste0(label, "_design"), qr(design)$rank == ncol(design) &&
            nrow(design) > ncol(design) + 2L && "diagnosisAD" %in% colnames(design),
          "Design must be full rank: no constant or perfectly confounded covariates.")
  }
  design <- model.matrix(main_formula, meta)
  design_sens <- model.matrix(sensitivity_formula, meta)
  validate_design(design, "main")
  validate_design(design_sens, "sensitivity")
  for (r in REGIONS) validate_design(model.matrix(regional_formula, meta[meta$region == r, ]), r)
  csv(meta, "metadata_used.csv")
  csv(data.frame(sample_id = rownames(design), design, check.names = FALSE), "design_matrix_main.csv")
  donor_summary <- do.call(rbind, lapply(levels(donors$diagnosis), function(g) {
    d <- donors[donors$diagnosis == g, ]
    data.frame(diagnosis = g, donors = nrow(d), female = sum(d$sex == "F"), male = sum(d$sex == "M"),
                age_mean = mean(d$age), age_SD = sd(d$age), age_min = min(d$age), age_max = max(d$age),
                PMI_mean = mean(d$PMI), PMI_SD = sd(d$PMI), PMI_min = min(d$PMI), PMI_max = max(d$PMI))
  }))
  csv(donor_summary, "donor_covariate_summary.csv")
  print(composition)
  print(donor_summary, row.names = FALSE)

  stage <- "expression filtering and TMM"
  message("[4/10] Filtering low expression and applying TMM normalization ...")
  y <- edgeR::DGEList(counts = m)
  keep <- edgeR::filterByExpr(y, group = interaction(meta$diagnosis, meta$region),
                              min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
  csv(data.frame(gene_symbol = rownames(m), all_zero = rowSums(m) == 0,
                   retained_for_main_analysis = keep), "gene_filter_audit.csv")
  y <- y[keep, , keep.lib.sizes = FALSE]
  check("genes_after_filtering", nrow(y) >= 100L, paste("Genes retained:", nrow(y)))
  y <- edgeR::calcNormFactors(y, method = "TMM")
  qc <- data.frame(meta, total_counts = colSums(m), detected_genes = colSums(m > 0),
                     filtered_library_size = y$samples$lib.size,
                     TMM_factor = y$samples$norm.factors,
                     effective_library_size = y$samples$lib.size * y$samples$norm.factors)
  csv(qc, "QC_samples.csv")
  saveRDS(list(counts = m, metadata = meta, donors = donors, retained_genes = keep,
                 design_main = design, design_no_PMI = design_sens),
          file.path(out, "checked_inputs.rds"), compress = FALSE)

  stage <- "parallel setup"
  message("[5/10] Starting parallel donor-correlation estimation ...")
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 2L
  workers <- max(1L, min(as.integer(WORKERS), as.integer(cores)))
  if (workers > 1L) {
    cl <- tryCatch(parallel::makePSOCKcluster(workers), error = function(e) {
      message("PSOCK startup failed; continuing with one process: ", conditionMessage(e))
      NULL
    })
    if (is.null(cl)) workers <- 1L
  }
  if (!is.null(cl)) {
    invisible(parallel::clusterCall(cl, function(paths) {
      .libPaths(paths)
      NULL
    }, .libPaths()))
  }
  fit_repeated <- function(design_now, label, verify_parallel = FALSE) {
    v0 <- limma::voom(y, design_now, plot = FALSE)
    if (verify_parallel && !is.null(cl)) {
      small <- v0[seq_len(min(32L, nrow(v0))), ]
      ref <- limma::duplicateCorrelation(small, design_now, block = meta$donor_id)
      got <- parallel_correlation(small, design_now, meta$donor_id, cl)
      check("parallel_matches_serial", isTRUE(all.equal(ref$atanh.correlations,
                got$atanh.correlations, tolerance = 1e-8, check.attributes = FALSE)) &&
              isTRUE(all.equal(ref$consensus.correlation, got$consensus.correlation,
                               tolerance = 1e-8, check.attributes = FALSE)),
            "First 32 genes: parallel and serial limma must agree.")
    }
    message("  ", label, ": correlation pass 1 (", workers, " processes, ", nrow(y), " genes)")
    c1 <- parallel_correlation(v0, design_now, meta$donor_id, cl)$consensus.correlation
    v <- with_pdf(file.path(out, paste0("QC_voom_", label, ".pdf")), 6, 5, {
      limma::voom(y, design_now, block = meta$donor_id, correlation = c1, plot = TRUE)
    })
    message("  ", label, ": correlation pass 2")
    c2 <- parallel_correlation(v, design_now, meta$donor_id, cl)$consensus.correlation
    fit <- limma::lmFit(v, design_now, block = meta$donor_id, correlation = c2)
    # Moderation and BH correction are performed on the complete gene set.
    fit <- limma::eBayes(fit)
    list(results = result_table(fit), voom = v, fit = fit, cor1 = c1, cor2 = c2)
  }

  stage <- "main analysis"
  message("[6/10] Main model: region + age + sex + PMI + diagnosis; donor blocked ...")
  main <- fit_repeated(design, "main", verify_parallel = TRUE)
  save_result_tables(main$results, out)

  stage <- "sensitivity analysis"
  message("[7/10] Sensitivity model without PMI, using the same donors and genes ...")
  sens <- fit_repeated(design_sens, "no_PMI")
  save_result_tables(sens$results, file.path(out, "sensitivity_no_PMI"))
  if (!is.null(cl)) { parallel::stopCluster(cl); cl <- NULL }
  si <- match(main$results$gene_symbol, sens$results$gene_symbol)
  comparison <- data.frame(gene_symbol = main$results$gene_symbol,
                  logFC_main = main$results$logFC, FDR_main = main$results$adj.P.Val,
                  status_main = main$results$status,
                  logFC_no_PMI = sens$results$logFC[si], FDR_no_PMI = sens$results$adj.P.Val[si],
                  status_no_PMI = sens$results$status[si])
  comparison$same_direction <- with(comparison, sign(logFC_main) == sign(logFC_no_PMI))
  csv(comparison, "sensitivity_PMI_comparison.csv")
  fc_spearman <- with(comparison, cor(logFC_main, logFC_no_PMI,
                                      method = "spearman", use = "complete.obs"))

  stage <- "regional analyses"
  message("[8/10] Examining the four subregions separately ...")
  regional_results <- setNames(vector("list", length(REGIONS)), REGIONS)
  region_comparison <- main$results[c("gene_symbol", "logFC", "adj.P.Val", "status")]
  names(region_comparison)[-1] <- c("logFC_main", "FDR_main", "status_main")
  regional_counts <- data.frame()
  for (r in REGIONS) {
    take <- meta$region == r
    md <- droplevels(meta[take, , drop = FALSE])
    yr <- edgeR::DGEList(counts = m[, take, drop = FALSE])
    kr <- edgeR::filterByExpr(yr, group = md$diagnosis,
                              min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
    csv(data.frame(gene_symbol = rownames(m), retained = kr), paste0("gene_filter_", r, ".csv"))
    check(paste0(r, "_genes_after_filtering"), sum(kr) >= 100L, paste("Genes retained:", sum(kr)))
    yr <- edgeR::calcNormFactors(yr[kr, , keep.lib.sizes = FALSE], method = "TMM")
    dr <- model.matrix(regional_formula, md)
    vr <- limma::voom(yr, dr, plot = FALSE)
    # One tissue per donor within a subregion: no repeated-measure block needed.
    fr <- limma::eBayes(limma::lmFit(vr, dr))
    rr <- result_table(fr)
    regional_results[[r]] <- rr
    save_result_tables(rr, file.path(out, "by_region", r))
    ii <- match(region_comparison$gene_symbol, rr$gene_symbol)
    region_comparison[[paste0(r, "_logFC")]] <- rr$logFC[ii]
    region_comparison[[paste0(r, "_FDR")]] <- rr$adj.P.Val[ii]
    regional_counts <- rbind(regional_counts, data.frame(region = r,
              NC = sum(md$diagnosis == "NC"), AD = sum(md$diagnosis == "AD"),
              genes_tested = nrow(rr), UP = sum(rr$status == "UP"), DOWN = sum(rr$status == "DOWN")))
  }
  regional_fc <- as.matrix(region_comparison[paste0(REGIONS, "_logFC")])
  region_comparison$n_regions_estimable <- rowSums(is.finite(regional_fc))
  same <- sign(regional_fc) == sign(region_comparison$logFC_main)
  region_comparison$n_regions_same_direction <- rowSums(same & is.finite(regional_fc), na.rm = TRUE)
  region_comparison$all_four_same_direction <- region_comparison$n_regions_estimable == 4L &
    region_comparison$n_regions_same_direction == 4L & region_comparison$logFC_main != 0
  csv(region_comparison, "regional_effect_direction.csv")
  csv(regional_counts, "regional_DEG_summary.csv")

  stage <- "quality control and figures"
  message("[9/10] Saving QC plots and the main volcano plot ...")
  vv <- main$voom$E
  rv <- rowMeans(vv^2) - rowMeans(vv)^2
  idx <- head(order(rv, decreasing = TRUE), 2000L)
  pc <- prcomp(t(vv[idx, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 3L)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  csv(data.frame(sample_id = rownames(pc$x), pc$x, check.names = FALSE), "QC_PCA_scores.csv")
  sample_cor <- cor(vv)
  csv(data.frame(sample_id = rownames(sample_cor), sample_cor, check.names = FALSE), "QC_sample_correlations.csv")
  with_pdf(file.path(out, "QC.pdf"), 10, 7, {
    colours <- c(NC = "#0072B2", AD = "#D55E00")
    shapes <- c(CA1 = 16, CA2 = 17, CA3 = 15, CA4 = 18)
    groups <- interaction(meta$diagnosis, meta$region, sep = "/")
    par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
    boxplot(log10(qc$total_counts) ~ groups, las = 2, xlab = "",
             ylab = "log10(sum of gene counts)", main = "Gene-count totals")
    boxplot(qc$detected_genes ~ groups, las = 2, xlab = "",
             ylab = "Genes with count > 0", main = "Gene detection")
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
    boxplot(age ~ diagnosis, donors, col = colours, xlab = "", ylab = "Age (years)", main = "Donor ages")
    boxplot(PMI ~ diagnosis, donors, col = colours, xlab = "", ylab = "PMI (hours)", main = "Donor PMI")
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
    plot(pc$x[, 1], pc$x[, 2], col = colours[as.character(meta$diagnosis)],
         pch = shapes[as.character(meta$region)],
         xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", pct[2]),
         main = "PCA: 2,000 most variable genes; covariates not removed")
    legend("topright", names(colours), col = colours, pch = 16, bty = "n")
    legend("bottomright", names(shapes), pch = shapes, bty = "n")
    heatmap(sample_cor, symm = TRUE, scale = "none", labRow = FALSE, labCol = FALSE,
              col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(60),
              margins = c(2, 2), main = "Sample correlations: filtered log2 expression")
    plot(comparison$logFC_main, comparison$logFC_no_PMI, pch = 16, cex = 0.35,
         col = grDevices::adjustcolor("grey35", alpha.f = 0.3),
         xlab = "Main model log2FC", ylab = "Model without PMI log2FC",
         main = sprintf("PMI sensitivity: Spearman rho = %.3f", fc_spearman))
    abline(0, 1, col = "#D55E00", lty = 2)
  })
  if (MAKE_VOLCANO) make_volcano(main$results, out)

  stage <- "saving results"
  message("[10/10] Saving analysis objects and the methods record ...")
  correlations <- data.frame(model = c("main", "no_PMI"),
                   pass1 = c(main$cor1, sens$cor1), pass2 = c(main$cor2, sens$cor2))
  csv(correlations, "donor_correlations.csv")
  if (SAVE_ANALYSIS_OBJECTS) {
    saveRDS(list(metadata = meta, donors = donors, filtered_DGEList = y,
            design_main = design, design_no_PMI = design_sens,
            voom_main = main$voom, fit_main = main$fit,
            results_main = main$results, results_no_PMI = sens$results,
            results_by_region = regional_results, correlations = correlations,
            settings = list(FDR = FDR_CUTOFF, abs_log2FC = LOG2FC_CUTOFF, workers = workers)),
            file.path(out, "analysis_objects.rds"), compress = FALSE)
  }
  notes <- c(
    "Dataset: GSE278723; analysis contrast: AD minus NC (positive logFC means higher in AD).",
    "Selection: author-assigned NC/AD only, hippocampal CA1-CA4; CAA/AD+CAA and NX/EC excluded.",
    "Unmatched donor IDs such as 451 are reported and never reassigned to 45.",
    "14 NC and 14 AD independent donors; 110 tissues. NC donors 03 and 14 lack CA1.",
    "Missing subregions were not imputed; other tissues from those donors were retained.",
    "Author-defined NC was retained, including NC with ABCscore L. NC is not claimed pathology-free.",
    "ABCscore and CAAscore are retained for audit, not automatically used for exclusions or adjustment.",
    "Covariates are individual donor age, sex (F reference) and postmortem interval in hours.",
    "Age and PMI were centered at their means across the 28 independent donors.",
    "Main design: ~ region + age_c + sex + PMI_c + diagnosis; donor_id is the repeated-measure block.",
    "Donor_id is not added as a fixed effect because diagnosis is constant within donor.",
    "Low-expression filtering: edgeR::filterByExpr, group=interaction(diagnosis,region),",
    "min.count=10, min.total.count=15, large.n=10, min.prop=0.7; TMM on retained genes.",
    "Main fitting: voom -> duplicateCorrelation -> voom with correlation -> duplicateCorrelation",
    "-> blocked lmFit -> eBayes -> topTable; the consensus donor correlation is shared across genes.",
    "The correlation steps may run in parallel, pooling all Fisher-z values with trim=0.15.",
    "Empirical Bayes moderation and BH correction use the full tested gene set, not separate chunks.",
    "This is limma-voom with duplicateCorrelation, not a gene-specific dream mixed model.",
    "Sensitivity: omit PMI, retain region, age, sex and donor block; same donors and filtered genes.",
    "Regional analyses: separate age/sex/PMI-adjusted models, each with one tissue per donor;",
    "low-expression filtering and TMM are recalculated within each region; BH within each result table.",
    "Regional analyses use overlapping donors and do not count as four independent datasets.",
    "A significant main-model coefficient does not establish significance in every subregion.",
    "RIN, sequencing batch and cell-type composition were not supplied and are not adjusted.",
    "PMI is not a substitute for RIN. No samples were automatically excluded using depth or PCA.",
    paste0("DEG label: BH-adjusted P < ", FDR_CUTOFF, " and absolute estimated log2FC >= ", LOG2FC_CUTOFF, "."),
    "This is an estimated-fold-change filter after a test against zero, not a TREAT threshold test.",
    "All tested genes are retained in the all_genes CSV; filtered-out genes are recorded in filter audits.",
    "The volcano plot uses FDR on the y axis; exactly zero FDR is floored at .Machine$double.xmin for display only.",
    "Positive/negative gene lists are for downstream integration; no arbitrary top-150 truncation is applied.",
    paste("Input genes:", nrow(m), "; all-zero genes:", sum(rowSums(m) == 0), "; main genes tested:", nrow(main$results)),
    paste("Main UP:", sum(main$results$status == "UP"), "; DOWN:", sum(main$results$status == "DOWN")),
    paste("Spearman correlation of logFC with vs without PMI:", signif(fc_spearman, 6)),
    paste("Parallel R processes:", workers),
    paste("Elapsed seconds including installation, I/O and figures:", round(proc.time()[["elapsed"]] - started, 1)),
    "", "Donor correlations:", capture.output(correlations),
    "", "Input MD5 checksums:", capture.output(tools::md5sum(files)),
    "", capture.output(sessionInfo()))
  writeLines(notes, file.path(out, "analysis_summary.txt"))
  writeLines(c("COMPLETE", format(Sys.time())), file.path(out, "COMPLETE.txt"))
  finished <- TRUE
  print(table(main$results$status))
  message("Finished. Results: ", normalizePath(out, winslash = "/"))
  invisible(list(results = main$results, metadata = meta, output_directory = out))
}

GSE278723_results <- run_GSE278723()
