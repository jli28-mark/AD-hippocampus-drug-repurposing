# GSE5281: inspect expression scale, sample linkage and preprocessing evidence.
# This script does not fit differential-expression models, change source files,
# or save a transformed/normalised expression matrix.
# Temporary log2 values are used only for diagnostic summaries/plots.
#
# source("D:/My Desktop/drug/1.DEG/GSE5281/check_GSE5281_input.R", encoding="UTF-8")
# GSE5281_checks <- check_GSE5281_input()

check_GSE5281_input <- function(
    base_dir = "D:/My Desktop/drug/1.DEG/GSE5281",
    matrix_file = file.path(base_dir, "GSE5281_series_matrix.txt"),
    clinical_file = NULL, fetch_geo = TRUE, n_threads = 4L) {

  for (p in c("data.table", "readxl")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      utils::install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  if (!file.exists(matrix_file)) stop("Expression matrix not found: ", matrix_file)
  if (is.null(clinical_file)) {
    candidates <- file.path(base_dir, c("sample(1).xlsx", "sample.xlsx"))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates)) stop("Cannot find sample(1).xlsx or sample.xlsx; set clinical_file explicitly.")
    clinical_file <- candidates[1L]
  }
  if (!file.exists(clinical_file)) stop("Clinical file not found: ", clinical_file)
  parent <- file.path(base_dir, "input_checks")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("check_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  report <- file.path(out, "GSE5281_check_report.txt")
  show <- function(x) {
    s <- utils::capture.output(if (is.character(x)) cat(x, sep = "\n") else print(x))
    cat(paste0(s, "\n"), sep = "")
    cat(paste0(s, "\n"), file = report, append = TRUE, sep = "")
  }
  csv <- function(x, filename) utils::write.csv(x, file.path(out, filename), row.names = FALSE, na = "")
  clean <- function(x) trimws(as.character(x))
  open_text <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  }
  # Handle both extracted ID_REF tables and complete GEO series matrices.
  # For the local file, locate the exact footer so it is not parsed as a probe.
  inspect_structure <- function(path, scan_footer = TRUE) {
    con <- open_text(path)
    on.exit(close(con))
    cursor <- 0L; header <- NA_integer_; footer <- NA_integer_; pre <- character()
    repeat {
      lines <- readLines(con, n = 1000L, warn = FALSE)
      if (!length(lines)) break
      lines <- sub("^\ufeff", "", lines)
      if (is.na(header)) {
        h <- which(grepl('^\\s*"?ID_REF"?\\s*\\t', lines, perl = TRUE))
        if (length(h)) {
          header <- cursor + h[1L]
          if (h[1L] > 1L) pre <- c(pre, lines[seq_len(h[1L] - 1L)])
        } else pre <- c(pre, lines)
      }
      f <- which(grepl('^\\s*"?!series_matrix_table_end"?\\s*$', lines, perl = TRUE))
      if (length(f)) footer <- cursor + f[1L]
      cursor <- cursor + length(lines)
      if (!is.na(footer) || (!scan_footer && !is.na(header))) break
      if (is.na(header) && cursor > 20000L) stop("No ID_REF header found: ", path)
    }
    if (is.na(header)) stop("No ID_REF header found: ", path)
    list(header = header, footer = footer, preamble = pre,
         n_data_lines = if (!is.na(footer)) footer - header - 1L else cursor - header)
  }
  # GEO preamble fields are tab-delimited quoted strings.
  geo_field <- function(lines, key) {
    hits <- lines[startsWith(lines, paste0(key, "\t"))]
    lapply(hits, function(s) {
      d <- utils::read.table(text = s, sep = "\t", quote = '"', header = FALSE,
        comment.char = "", fill = TRUE, colClasses = "character", check.names = FALSE)
      as.character(unlist(d[1L, -1L, drop = FALSE], use.names = FALSE))
    })
  }
  extract_geo <- function(lines, ids) {
    accessions <- geo_field(lines, "!Sample_geo_accession")
    if (!length(accessions)) return(NULL)
    all_ids <- accessions[[1L]]
    index <- match(ids, all_ids)
    join_field <- function(key) {
      rows <- geo_field(lines, key)
      vapply(index, function(j) {
        if (is.na(j) || !length(rows)) return(NA_character_)
        vals <- vapply(rows, function(z) if (length(z) >= j) z[j] else NA_character_, character(1))
        vals <- unique(vals[!is.na(vals) & nzchar(vals)])
        if (length(vals)) paste(vals, collapse = " | ") else NA_character_
      }, character(1))
    }
    data.frame(Sample_ID = ids, Found_in_record = !is.na(index),
      GEO_Title = join_field("!Sample_title"),
      GEO_Processing = join_field("!Sample_data_processing"),
      GEO_Characteristics = join_field("!Sample_characteristics_ch1"),
      GEO_Description = join_field("!Sample_description"), stringsAsFactors = FALSE)
  }

  show("[1/5] Clinical sample checks")
  show(paste("Clinical file:", clinical_file))
  meta <- as.data.frame(readxl::read_excel(clinical_file, col_types = "text"))
  required <- c("Sample_ID", "Sample_Title", "Group", "Sex", "Age", "Region", "Cell_Type")
  if (!all(required %in% names(meta))) stop("Clinical table must contain: ", paste(required, collapse = ", "))
  meta <- meta[!apply(meta, 1L, function(z) all(is.na(z) | clean(z) == "")), , drop = FALSE]
  meta$Sample_ID <- clean(meta$Sample_ID)
  group_map <- c(control = "Control", ctrl = "Control", nc = "Control", ad = "AD")
  sex_map <- c(male = "M", m = "M", female = "F", f = "F")
  meta$Group <- unname(group_map[tolower(clean(meta$Group))])
  meta$Sex <- unname(sex_map[tolower(clean(meta$Sex))])
  meta$Age <- suppressWarnings(as.numeric(clean(meta$Age)))
  if (anyNA(meta[, c("Sample_ID", "Group", "Sex", "Age")]) ||
      anyDuplicated(meta$Sample_ID) || any(!grepl("^GSM[0-9]+$", meta$Sample_ID)) ||
      any(!is.finite(meta$Age) | meta$Age <= 0) ||
      any(is.na(meta$Region) | tolower(clean(meta$Region)) != "hippocampus")) {
    csv(meta, "ERROR_clinical_fields.csv")
    stop("Invalid or duplicated clinical fields; see ", out)
  }
  expected_control <- paste0("GSM", 119628:119640)
  expected_ad <- paste0("GSM", 238799:238808)
  if (!setequal(meta$Sample_ID[meta$Group == "Control"], expected_control) ||
      !setequal(meta$Sample_ID[meta$Group == "AD"], expected_ad)) {
    csv(meta, "ERROR_clinical_cohort.csv")
    stop("Clinical accessions/groups differ from the intended 13-control/10-AD cohort.")
  }
  csv(meta, "clinical_checked.csv")
  cohort <- do.call(rbind, lapply(c("Control", "AD"), function(g) {
    z <- meta[meta$Group == g, , drop = FALSE]
    data.frame(Group = g, N = nrow(z), Female = sum(z$Sex == "F"), Male = sum(z$Sex == "M"),
      Age_mean = mean(z$Age), Age_SD = stats::sd(z$Age), Age_min = min(z$Age), Age_max = max(z$Age))
  }))
  show(cohort); csv(cohort, "cohort_summary.csv")

  show("[2/5] Reading the FULL expression matrix and checking the 23 selected samples")
  structure <- inspect_structure(matrix_file)
  tab <- data.table::fread(matrix_file, skip = structure$header - 1L,
    nrows = structure$n_data_lines, colClasses = list(character = 1L),
    check.names = FALSE, data.table = FALSE, nThread = n_threads, blank.lines.skip = TRUE)
  if (names(tab)[1L] != "ID_REF") stop("First expression column is not ID_REF.")
  if (anyDuplicated(names(tab)[-1L])) stop("Duplicated GSM columns in the expression file.")
  if (!all(grepl("^GSM[0-9]+$", names(tab)[-1L]))) stop("Some expression columns are not GSM samples.")
  if (!all(vapply(tab[-1L], is.numeric, logical(1)))) stop("Some expression columns are not numeric.")
  ids <- clean(tab[[1L]])
  x <- as.matrix(tab[-1L]); storage.mode(x) <- "double"; rownames(x) <- ids
  rm(tab)
  missing_ids <- setdiff(meta$Sample_ID, colnames(x))
  if (length(missing_ids)) stop("Clinical accessions absent from matrix: ", paste(missing_ids, collapse = ", "))
  hip <- x[, meta$Sample_ID, drop = FALSE]
  csv(data.frame(Sample_ID = colnames(x), Selected_hippocampus = colnames(x) %in% meta$Sample_ID),
      "expression_sample_selection.csv")
  basic <- list(full_matrix_probes_samples = dim(x), hippocampal_probes_samples = dim(hip),
    duplicated_probe_ids = sum(duplicated(ids)), missing_or_empty_probe_ids = sum(is.na(ids) | ids == ""),
    missing_values_all = sum(is.na(x)), infinite_values_all = sum(is.infinite(x)),
    zero_values_all = sum(x == 0, na.rm = TRUE), negative_values_all = sum(x < 0, na.rm = TRUE),
    missing_values_HIP = sum(is.na(hip)), infinite_values_HIP = sum(is.infinite(hip)),
    zero_values_HIP = sum(hip == 0, na.rm = TRUE), negative_values_HIP = sum(hip < 0, na.rm = TRUE),
    footer_found = !is.na(structure$footer),
    probe_id_examples = head(ids, 10L))
  show(basic)

  show("[3/5] Raw scale and temporary log2 summaries (no transformation is saved)")
  probs <- c(0, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 1)
  scale_table <- function(e, label) {
    vals <- e[is.finite(e)]
    if (!length(vals)) stop("No finite expression values.")
    q <- stats::quantile(vals, probs = probs, names = FALSE)
    qlog <- if (all(vals > 0)) stats::quantile(log2(vals), probs = probs, names = FALSE) else rep(NA_real_, length(probs))
    data.frame(Matrix = label, Percentile = 100 * probs, Original = q, Temporary_log2 = qlog)
  }
  scales <- rbind(scale_table(x, "all_samples"), scale_table(hip, "HIP_23_samples"))
  show(scales); csv(scales, "scale_quantiles.csv")
  sample_qc <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    v <- x[, j]; finite <- v[is.finite(v)]
    q <- if (length(finite)) stats::quantile(finite, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), names = FALSE) else rep(NA_real_, 7L)
    data.frame(Sample_ID = colnames(x)[j], Selected_HIP = colnames(x)[j] %in% meta$Sample_ID,
      Missing = sum(is.na(v)), Infinite = sum(is.infinite(v)), Zero = sum(v == 0, na.rm = TRUE),
      Negative = sum(v < 0, na.rm = TRUE), Min = q[1L], P1 = q[2L], Q25 = q[3L],
      Median = q[4L], Q75 = q[5L], P99 = q[6L], Max = q[7L])
  }))
  csv(sample_qc, "per_sample_quantiles_all.csv")
  show(sample_qc[match(meta$Sample_ID, sample_qc$Sample_ID), ])
  if (all(is.finite(hip)) && all(hip > 0)) {
    grDevices::pdf(file.path(out, "HIP_raw_and_temporary_log2.pdf"), width = 11, height = 8,
                    useDingbats = FALSE)
    tryCatch({
      graphics::par(mfrow = c(2, 1), mar = c(7, 4, 3, 1))
      cols <- ifelse(meta$Group == "AD", "#D55E00", "#0072B2")
      graphics::boxplot(as.data.frame(hip), outline = FALSE, las = 2, cex.axis = 0.65, col = cols,
        main = "GSE5281 hippocampus: original values", ylab = "Input expression")
      graphics::boxplot(as.data.frame(log2(hip)), outline = FALSE, las = 2, cex.axis = 0.65, col = cols,
        main = "Temporary log2 values for scale inspection", ylab = "Temporary log2 expression")
    }, finally = grDevices::dev.off())
  } else show("Diagnostic log2 plot skipped: selected samples contain nonfinite or nonpositive values.")
  show(c("Interpretation: large values spanning orders of magnitude support linear intensity.",
         "A narrow range is not proof of prior log2 transformation.",
         "Numeric distributions alone cannot identify MAS5, RMA, or quantile normalization."))

  show("[4/5] Local preprocessing records")
  evidence <- structure$preamble[grepl("processing|normaliz|MAS5|RMA|quantile|scal", structure$preamble,
                                       ignore.case = TRUE)]
  writeLines(structure$preamble, file.path(out, "local_GEO_preamble.txt"))
  if (length(evidence)) show(evidence) else show("No preprocessing record in the local preamble; the local file may be an extracted expression table.")
  local_geo <- extract_geo(structure$preamble, meta$Sample_ID)
  if (!is.null(local_geo)) csv(local_geo, "local_GEO_sample_records.csv")

  show("[5/5] Optional GEO reference records and expression spot-check")
  geo_records <- NULL; reference_check <- NULL
  if (fetch_geo) {
    url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/GSE5281_series_matrix.txt.gz"
    reference_file <- file.path(out, "GSE5281_GEO_reference_series_matrix.txt.gz")
    old_timeout <- getOption("timeout")
    options(timeout = max(180, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    tryCatch({
      show(paste("Downloading GEO reference:", url))
      status <- utils::download.file(url, reference_file, mode = "wb", quiet = FALSE)
      if (!identical(status, 0L)) stop("Reference download failed.")
      ref_structure <- inspect_structure(reference_file, scan_footer = FALSE)
      writeLines(ref_structure$preamble, file.path(out, "reference_GEO_preamble.txt"))
      geo_records <- extract_geo(ref_structure$preamble, meta$Sample_ID)
      if (is.null(geo_records)) stop("Reference has no readable sample metadata.")
      csv(geo_records, "GEO_reference_sample_records.csv")
      show(geo_records[, c("Sample_ID", "GEO_Title", "GEO_Processing")])
      # Reading a small reference block avoids loading a second full matrix.
      con <- open_text(reference_file)
      ref <- tryCatch(utils::read.table(con, header = TRUE, sep = "\t", quote = '"',
        skip = ref_structure$header - 1L, nrows = 20L, comment.char = "", check.names = FALSE,
        stringsAsFactors = FALSE), finally = close(con))
      shared_probes <- intersect(as.character(ref[[1L]]), ids)
      shared_samples <- intersect(meta$Sample_ID, names(ref)[-1L])
      if (!length(shared_probes) || !length(shared_samples)) stop("No shared identifiers for reference comparison.")
      r <- as.matrix(ref[match(shared_probes, ref[[1L]]), shared_samples, drop = FALSE])
      local <- x[match(shared_probes, ids), shared_samples, drop = FALSE]
      delta <- abs(local - r)
      tolerance <- 1e-5 * pmax(1, abs(r))
      reference_check <- data.frame(Probes_compared = length(shared_probes),
        Samples_compared = length(shared_samples), Values_compared = length(delta),
        Values_within_tolerance = sum(is.finite(delta) & delta <= tolerance),
        Max_absolute_difference = max(delta, na.rm = TRUE), Relative_tolerance = 1e-5)
      csv(reference_check, "GEO_expression_spotcheck.csv"); show(reference_check)
      show("This is a spot-check, not a full-file identity test. Reference processing records cannot automatically be assigned to a mismatching local matrix.")
    }, error = function(e) show(paste("GEO reference retrieval/comparison incomplete:", conditionMessage(e),
      "The local numerical checks remain available; preprocessing provenance is not fully confirmed.")))
  } else show("GEO download disabled. Preprocessing must be checked against source records separately.")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  show(paste("Finished checking. Output folder:", normalizePath(out, winslash = "/")))
  show("No differential-expression model or normalization was run; source files were not modified.")
  invisible(list(output_dir = normalizePath(out, winslash = "/"), basic_checks = basic,
    cohort_summary = cohort, quantile_comparison = scales, sample_qc = sample_qc,
    local_records = local_geo, GEO_records = geo_records, GEO_spotcheck = reference_check))
}
