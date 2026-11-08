# Extract paired main-model effects/SE, harmonize genes and define a NEW common universe.
# No models are refitted. No significance/effect-size filter or meta-analysis is run.
# Uses the existing HGNC snapshot; never silently downloads/updates it.
# source("D:/My Desktop/drug/2.Gene_harmonization/prepare_four_study_meta_inputs.R", encoding="UTF-8")
# meta_inputs <- prepare_four_study_meta_inputs()

prepare_four_study_meta_inputs <- function(
    base_dir = "D:/My Desktop/drug/2.Gene_harmonization",
    hgnc_file = NULL, n_threads = 4L, tolerance = 1e-7) {
  required_packages <- c("data.table", "limma", "edgeR", "DESeq2")
  for (p in required_packages) if (!requireNamespace(p, quietly = TRUE))
    stop("Required package missing: ", p, ". Install it before running this script.", call. = FALSE)
  stopifnot(dir.exists(base_dir), tolerance > 0)
  studies <- c("GSE278723", "GSE173955", "GSE48350", "GSE5281")
  check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  clean <- function(x) trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  find_input <- function(filename) {
    choices <- list.files(base_dir, full.names = TRUE)
    hit <- choices[key(basename(choices)) == key(filename)]
    check(length(hit) == 1L, paste("Expected exactly one", filename, "in", base_dir,
      "but found", length(hit)))
    hit
  }
  csv_files <- setNames(vapply(studies, function(s)
    find_input(paste0(s, "_all_genes.csv")), character(1)), studies)
  model_names <- c("analysis_objects.rds", "GSE173955_analysis_objects.rds",
    "GSE48350_analysis_objects.rds", "GSE5281_analysis_objects.rds")
  model_files <- setNames(vapply(model_names, find_input, character(1)), studies)
  # Prefer the reference cached by the first harmonization. If absent, use an
  # old run snapshot only when all located snapshots have the same checksum.
  if (is.null(hgnc_file)) {
    cached <- file.path(base_dir, "annotation", "hgnc_complete_set.txt")
    if (file.exists(cached)) hgnc_file <- cached else {
      candidates <- unique(c(
        file.path(base_dir, "hgnc_complete_set.txt"),
        file.path(base_dir, "HGNC_reference_snapshot.txt"),
        list.files(file.path(base_dir, "results"), pattern = "^HGNC_reference_snapshot\\.txt$",
          recursive = TRUE, full.names = TRUE)))
      candidates <- candidates[file.exists(candidates)]
      check(length(candidates) > 0, paste("Existing HGNC snapshot not found. Set hgnc_file explicitly; do not download a new version."))
      check(length(unique(unname(tools::md5sum(candidates)))) == 1L,
        "Several different HGNC snapshots were found. Set hgnc_file to the intended earlier snapshot.")
      hgnc_file <- candidates[1L]
    }
  }
  check(file.exists(hgnc_file), paste("HGNC file not found:", hgnc_file))
  parent <- file.path(base_dir, "results_with_SE")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  old_threads <- data.table::getDTthreads(); data.table::setDTthreads(n_threads)
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  info <- data.frame(Dataset = rep(studies, 2L), Type = rep(c("all_gene_CSV", "model_RDS"), each = 4L),
    File = c(unname(csv_files), unname(model_files)))
  info$MD5 <- unname(tools::md5sum(info$File)); csv(info, "input_files_checksums.csv")
  check(file.copy(hgnc_file, file.path(out, "HGNC_reference_snapshot.txt")), "Cannot save HGNC snapshot.")
  csv(data.frame(File = normalizePath(hgnc_file, winslash = "/"),
    MD5 = unname(tools::md5sum(hgnc_file))), "annotation_record.csv")

  numeric_column <- function(d, choices, required = TRUE) {
    name <- choices[choices %in% names(d)][1L]
    if (is.na(name)) {
      if (required) stop("Missing column: ", paste(choices, collapse = " / "))
      return(rep(NA_real_, nrow(d)))
    }
    suppressWarnings(as.numeric(clean(d[[name]])))
  }
  text_column <- function(d, choices) {
    name <- choices[choices %in% names(d)][1L]
    if (is.na(name)) rep(NA_character_, nrow(d)) else clean(d[[name]])
  }
  same_numbers <- function(x, y, label) {
    check(length(x) == length(y) && identical(is.na(x), is.na(y)),
      paste(label, "has missing-value or length disagreement between CSV and saved model."))
    ok <- !is.na(x)
    check(all(is.finite(x[ok])) && all(is.finite(y[ok])), paste(label, "contains infinite values."))
    check(all(abs(x[ok] - y[ok]) <= tolerance * pmax(1, abs(y[ok]))),
      paste(label, "does not match saved main-model output. Possible wrong model/run/contrast; stopping."))
    if (any(ok)) max(abs(x[ok] - y[ok])) else NA_real_
  }
  preparations <- list(); checks <- list(); models <- list()
  message("[1/5] Extracting paired effects and SE; checking every CSV row against its saved MAIN model ...")
  for (s in studies) {
    raw <- data.table::fread(csv_files[[s]], data.table = FALSE, check.names = FALSE,
      colClasses = "character", na.strings = c("", "NA", "NaN"))
    check(nrow(raw) > 100L && !anyDuplicated(names(raw)), paste(s, "invalid CSV."))
    original_symbol <- text_column(raw, c("gene_symbol", "Gene", "Symbol"))
    effect_csv <- numeric_column(raw, c("log2FoldChange", "logFC", "log2FC"))
    p_csv <- numeric_column(raw, c("pvalue", "P.Value", "P_value"), required = FALSE)
    fdr_csv <- numeric_column(raw, c("padj", "adj.P.Val", "FDR"), required = FALSE)
    z <- readRDS(model_files[[s]])
    if (s == "GSE173955") {
      check(!is.null(z$dds) && !is.null(z$results), "GSE173955 RDS must contain dds and the new DESeq2 results.")
      check(identical(levels(z$dds$Group), c("Control", "AD")), "GSE173955 Group reference is unexpected.")
      fresh <- DESeq2::results(z$dds, contrast = c("Group", "AD", "Control"),
        alpha = if (!is.null(z$settings$FDR)) z$settings$FDR else 0.05,
        pAdjustMethod = "BH", independentFiltering = TRUE, parallel = FALSE)
      model_ids <- rownames(fresh)
      input_ids <- text_column(raw, "GeneID")
      check(!anyNA(input_ids) && !anyDuplicated(model_ids), "Missing input GeneID or nonunique DESeq2 rows.")
      j <- match(input_ids, model_ids)
      check(!anyNA(j) && setequal(input_ids, model_ids), "GSE173955 CSV does not cover the same gene IDs as the saved model.")
      effect <- as.numeric(fresh$log2FoldChange[j]); se <- as.numeric(fresh$lfcSE[j])
      p_model <- as.numeric(fresh$pvalue[j]); fdr_model <- as.numeric(fresh$padj[j])
      same_numbers(numeric_column(raw, c("lfcSE", "SE")), se, "GSE173955 SE")
      se_method <- "DESeq2_Wald_model_SE_unshrunken"
      coefficient <- "Group: AD - Control"
      meta <- z$metadata
    } else {
      if (s == "GSE278723") {
        fit <- z$fit_main; stored <- z$results_main; meta <- z$metadata
        coefficient <- "diagnosisAD"
        check(!is.null(meta) && identical(levels(meta$diagnosis), c("NC", "AD")),
          "GSE278723 diagnosis reference must be NC then AD.")
      } else if (s == "GSE48350") {
        fit <- z$main$fit; stored <- z$main$result; meta <- z$main$metadata
        coefficient <- "AD_vs_Control"
        check(!is.null(meta) && all(meta$Age >= 60) &&
          sum(meta$Group == "Control") == 25L && sum(meta$Group == "AD") == 19L,
          "GSE48350 must be the main age>=60 cohort: 25 Control and 19 AD.")
        check(!is.null(fit$contrasts) && all(c("GroupAD", "GroupControl") %in% rownames(fit$contrasts)) &&
          fit$contrasts["GroupAD", coefficient] == 1 && fit$contrasts["GroupControl", coefficient] == -1,
          "GSE48350 contrast is not the expected AD-Control contrast.")
      } else {
        fit <- z$fit; stored <- z$results; meta <- z$meta
        coefficient <- "GroupAD"
        check(!is.null(meta) && identical(levels(meta$Group), c("Control", "AD")),
          "GSE5281 Group reference must be Control then AD.")
      }
      check(!is.null(fit$coefficients) && coefficient %in% colnames(fit$coefficients) &&
        !is.null(fit$stdev.unscaled) && !is.null(fit$sigma), paste(s, "saved main fit/SE components missing."))
      model_ids <- rownames(fit$coefficients)
      input_ids <- original_symbol
      check(!anyNA(input_ids) && !anyDuplicated(model_ids), paste(s, "missing/nonunique model gene symbols."))
      j <- match(input_ids, model_ids)
      check(!anyNA(j) && setequal(input_ids, model_ids), paste(s, "CSV genes differ from saved main model."))
      effect <- as.numeric(fit$coefficients[j, coefficient])
      # Ordinary model SE, not logFC/moderated_t and not a P-value reconstruction.
      se <- as.numeric(fit$stdev.unscaled[j, coefficient] * fit$sigma[j])
      stored_symbols <- text_column(stored, c("gene_symbol", "Gene"))
      check(!anyNA(stored_symbols) && !anyDuplicated(stored_symbols), paste(s, "stored result identities are ambiguous."))
      sj <- match(input_ids, stored_symbols)
      check(!anyNA(sj), paste(s, "CSV gene absent from saved result table."))
      p_model <- numeric_column(stored, c("P.Value", "pvalue"))[sj]
      fdr_model <- numeric_column(stored, c("adj.P.Val", "padj", "FDR"))[sj]
      se_method <- "limma_ordinary_SE_sigma_times_stdev_unscaled"
    }
    delta <- same_numbers(effect_csv, effect, paste(s, "log2FC"))
    if (any(is.finite(p_csv))) same_numbers(p_csv, p_model, paste(s, "original P values"))
    if (any(is.finite(fdr_csv))) same_numbers(fdr_csv, fdr_model, paste(s, "original FDR"))
    # Keep original statistics; original limma P/FDR use moderated variance whereas
    # the SE here deliberately uses ordinary model variance for inverse weighting.
    d <- data.frame(Dataset = s, Source_row = seq_len(nrow(raw)) + 1L,
      Original_symbol = original_symbol,
      Input_EntrezID = text_column(raw, c("GeneID", "ENTREZID")),
      Model_gene_ID = input_ids, log2FC = effect, SE = se, variance = se^2,
      P_value = p_model, FDR = fdr_model, SE_method = se_method,
      Model_coefficient = coefficient, stringsAsFactors = FALSE)
    valid <- is.finite(d$log2FC) & is.finite(d$SE) & d$SE > 0 & is.finite(d$variance) & d$variance > 0
    # Respect model-level unavailable tests, e.g. DESeq2 Cook's outlier handling.
    # Missing padj alone is NOT an exclusion when effect/SE/raw P are available.
    test_available <- is.finite(d$P_value) & d$P_value >= 0 & d$P_value <= 1
    d$Input_status <- ifelse(!valid, "invalid_effect_or_SE",
      ifelse(!test_available, "model_test_unavailable", ifelse(duplicated(raw), "identical_source_row_duplicate", "valid")))
    csv(d, paste0(s, "_effects_SE_before_harmonization.csv"))
    csv(d[d$Input_status != "valid", ], paste0(s, "_SE_input_exclusions.csv"))
    checks[[s]] <- data.frame(Dataset = s, CSV_rows = nrow(d), Model_rows = length(model_ids),
      Max_abs_log2FC_difference = delta, Coefficient = coefficient, SE_method = se_method,
      Valid_effect_SE = sum(valid), Model_test_unavailable = sum(valid & !test_available),
      Input_eligible = sum(d$Input_status == "valid"), FDR_missing = sum(!is.finite(d$FDR)))
    preparations[[s]] <- d
    models[[s]] <- data.frame(Dataset = s, Metadata_rows = if (is.null(meta)) NA_integer_ else nrow(meta),
      Coefficient = coefficient)
    if (!is.null(meta)) csv(meta, paste0(s, "_model_metadata.csv"))
    message(s, ": CSV/model effects verified; ", sum(d$Input_status == "valid"), " eligible rows.")
    rm(z)
  }
  se_checks <- do.call(rbind, checks); rownames(se_checks) <- NULL
  csv(se_checks, "model_SE_checks.csv"); print(se_checks)

  message("[2/5] Loading the existing HGNC snapshot ...")
  h <- data.table::fread(hgnc_file, colClasses = "character", data.table = FALSE, na.strings = c("", "NA"))
  needed <- c("hgnc_id", "symbol", "status", "locus_group", "locus_type", "entrez_id", "ensembl_gene_id", "prev_symbol", "alias_symbol")
  check(all(needed %in% names(h)), "Existing HGNC snapshot lacks required fields.")
  h <- h[h$status == "Approved", , drop = FALSE]
  check(nrow(h) > 10000L && !anyDuplicated(h$hgnc_id) && !anyDuplicated(h$symbol), "Invalid HGNC snapshot.")
  expand_field <- function(field) {
    z <- strsplit(ifelse(is.na(h[[field]]), "", h[[field]]), "|", fixed = TRUE)
    data.frame(Token = unlist(z, use.names = FALSE), Index = rep(seq_len(nrow(h)), lengths(z)))
  }
  lk <- unique(rbind(expand_field("prev_symbol"), expand_field("alias_symbol")))
  lk <- lk[!is.na(lk$Token) & nzchar(lk$Token), ]
  aliases <- split(lk$Index, lk$Token)
  entrez <- split(seq_len(nrow(h)), h$entrez_id)
  all_input <- do.call(rbind, preparations); rownames(all_input) <- NULL
  identities <- unique(all_input[c("Original_symbol", "Input_EntrezID")])
  identities$Identity_key <- paste(ifelse(is.na(identities$Original_symbol), "", identities$Original_symbol),
    ifelse(is.na(identities$Input_EntrezID), "", identities$Input_EntrezID), sep = "\t")

  message("[3/5] Resolving identifiers and selecting protein-coding genes ...")
  mapping <- do.call(rbind, lapply(seq_len(nrow(identities)), function(i) {
    s <- identities$Original_symbol[i]; e <- identities$Input_EntrezID[i]
    direct <- match(s, h$symbol)
    sym <- if (!is.na(direct)) direct else if (!is.na(s) && nzchar(s)) unique(aliases[[s]]) else integer()
    supplied_id <- !is.na(e) && nzchar(e)
    valid_id <- supplied_id && grepl("^[0-9]+$", e)
    id <- if (valid_id) unique(entrez[[e]]) else integer()
    target <- NA_integer_; method <- "unmapped"; failure <- ""
    if (supplied_id && !valid_id) {
      failure <- "nonunique_or_invalid_source_EntrezID"
    } else if (length(id) > 1L) {
      failure <- "ambiguous_source_EntrezID"
    } else if (length(id) == 1L) {
      # A stable source ID is stronger evidence than a name alias. However, an
      # exact current-symbol disagreement is not silently overwritten.
      if (!is.na(direct) && direct != id) failure <- "source_ID_current_symbol_conflict" else {
        target <- id; method <- if (!is.na(direct)) "source_EntrezID_and_approved_symbol" else "source_EntrezID"
      }
    } else if (length(sym) == 1L) {
      target <- sym
      method <- if (!is.na(direct)) "approved_symbol" else "unique_previous_or_alias"
      if (supplied_id) method <- paste0(method, "_source_ID_unresolved")
    } else if (length(sym) > 1L) failure <- "ambiguous_previous_or_alias" else failure <- "unmapped"
    data.frame(Identity_key = identities$Identity_key[i],
      HGNC_ID = if (!is.na(target)) h$hgnc_id[target] else NA_character_,
      gene_symbol = if (!is.na(target)) h$symbol[target] else NA_character_,
      Mapping_method = method, Mapping_failure = failure,
      Symbol_candidate_HGNC_IDs = if (length(sym)) paste(h$hgnc_id[sym], collapse = "|") else "",
      Source_ID_candidate_HGNC_IDs = if (length(id)) paste(h$hgnc_id[id], collapse = "|") else "",
      locus_group = if (!is.na(target)) h$locus_group[target] else NA_character_,
      locus_type = if (!is.na(target)) h$locus_type[target] else NA_character_,
      ENTREZID = if (!is.na(target)) h$entrez_id[target] else NA_character_,
      Ensembl_gene_ID = if (!is.na(target)) h$ensembl_gene_id[target] else NA_character_,
      stringsAsFactors = FALSE)
  }))
  csv(cbind(identities[c("Original_symbol", "Input_EntrezID")], mapping), "identity_mapping_dictionary.csv")
  identity_keys <- paste(ifelse(is.na(all_input$Original_symbol), "", all_input$Original_symbol),
    ifelse(is.na(all_input$Input_EntrezID), "", all_input$Input_EntrezID), sep = "\t")
  matched <- match(identity_keys, mapping$Identity_key)
  check(!anyNA(matched), "Internal identifier lookup failure.")
  audit <- cbind(all_input, mapping[matched, setdiff(names(mapping), "Identity_key"), drop = FALSE])
  rownames(audit) <- NULL
  audit$Final_status <- audit$Input_status
  eligible <- audit$Input_status == "valid"
  audit$Final_status[eligible & is.na(audit$HGNC_ID)] <- "unresolved_identifier"
  mapped <- eligible & !is.na(audit$HGNC_ID)
  coding <- mapped & !is.na(audit$locus_group) & audit$locus_group == "protein-coding gene"
  coding[is.na(coding)] <- FALSE
  audit$Final_status[mapped & !coding] <- "not_protein_coding"
  audit$Final_status[coding] <- "retained"
  dup_key <- paste(audit$Dataset, audit$HGNC_ID, sep = "|")
  duplicates <- names(which(table(dup_key[coding]) > 1L))
  audit$Final_status[coding & dup_key %in% duplicates] <- "conflicting_duplicate_HGNC_ID"
  csv(audit, "gene_mapping_SE_audit.csv")
  csv(audit[audit$Final_status %in% c("identical_source_row_duplicate", "conflicting_duplicate_HGNC_ID"), ],
    "duplicate_gene_SE_audit.csv")
  csv(audit[audit$Final_status != "retained", ], "all_excluded_records.csv")
  results <- list(); summary_list <- list()
  cols <- c("Dataset", "HGNC_ID", "gene_symbol", "Original_symbol", "Input_EntrezID", "Model_gene_ID",
    "ENTREZID", "Ensembl_gene_ID", "log2FC", "SE", "variance", "P_value", "FDR", "SE_method",
    "Model_coefficient", "Mapping_method", "locus_group", "locus_type", "Source_row")
  for (s in studies) {
    d <- audit[audit$Dataset == s, , drop = FALSE]
    r <- d[d$Final_status == "retained", cols, drop = FALSE]
    r <- r[order(r$HGNC_ID), ]; rownames(r) <- NULL
    check(nrow(r) > 100L && !anyDuplicated(r$HGNC_ID) && !anyDuplicated(r$gene_symbol),
      paste(s, "has too few or nonunique retained genes."))
    csv(r, paste0(s, "_harmonized_meta_input.csv")); results[[s]] <- r
    summary_list[[s]] <- data.frame(Dataset = s, Input_rows = nrow(d),
      Invalid_effect_SE = sum(d$Input_status == "invalid_effect_or_SE"),
      Test_unavailable = sum(d$Input_status == "model_test_unavailable"),
      Exact_duplicate_rows = sum(d$Input_status == "identical_source_row_duplicate"),
      Unresolved_ID_rows = sum(d$Final_status == "unresolved_identifier"),
      Noncoding_rows = sum(d$Final_status == "not_protein_coding"),
      Conflicting_HGNC_IDs = length(unique(d$HGNC_ID[d$Final_status == "conflicting_duplicate_HGNC_ID"])),
      Conflicting_rows = sum(d$Final_status == "conflicting_duplicate_HGNC_ID"),
      Retained_coding_genes_with_SE = nrow(r),
      Retained_FDR_missing = sum(!is.finite(r$FDR)), SE_min = min(r$SE), SE_median = stats::median(r$SE), SE_max = max(r$SE))
  }
  dataset_summary <- do.call(rbind, summary_list); rownames(dataset_summary) <- NULL
  csv(dataset_summary, "dataset_summary_with_SE.csv"); print(dataset_summary)

  message("[4/5] Recomputing four-study coverage using valid paired effects/SE ...")
  union_ids <- sort(unique(unlist(lapply(results, function(r) r$HGNC_ID))))
  coverage <- data.frame(HGNC_ID = union_ids, gene_symbol = h$symbol[match(union_ids, h$hgnc_id)])
  for (s in studies) coverage[[s]] <- union_ids %in% results[[s]]$HGNC_ID
  coverage$N_datasets <- rowSums(coverage[studies])
  csv(coverage, "gene_coverage_with_SE.csv")
  coverage_summary <- data.frame(N_datasets = 1:4,
    N_genes = as.integer(table(factor(coverage$N_datasets, levels = 1:4))))
  csv(coverage_summary, "gene_coverage_summary_with_SE.csv"); print(coverage_summary)
  common <- coverage[coverage$N_datasets == 4L, c("HGNC_ID", "gene_symbol"), drop = FALSE]
  csv(common, "common_genes_four_studies_with_SE.csv")
  common_inputs <- lapply(results, function(r) r[match(common$HGNC_ID, r$HGNC_ID), , drop = FALSE])
  check(all(vapply(common_inputs, function(r) !anyNA(r$HGNC_ID) &&
    identical(r$gene_symbol, common$gene_symbol), logical(1))), "Common gene alignment failure.")
  all_long <- do.call(rbind, results); rownames(all_long) <- NULL
  csv(all_long, "four_studies_all_available_meta_input_long.csv")
  common_long <- do.call(rbind, common_inputs); rownames(common_long) <- NULL
  csv(common_long, "four_studies_common_meta_input_long.csv")
  wide <- common
  for (s in studies) for (v in c("log2FC", "SE", "variance", "P_value", "FDR"))
    wide[[paste(s, v, sep = "_")]] <- common_inputs[[s]][[v]]
  csv(wide, "four_studies_common_effects_SE_wide.csv")
  for (s in studies) csv(common_inputs[[s]], paste0(s, "_common_meta_input.csv"))

  message("[5/5] Saving audit and reproducibility records ...")
  rules <- c(
    "Purpose: prepare paired effects/SE and a new four-study universe. No meta-analysis or differential fitting performed.",
    "All effects are AD minus Control/NC. Every CSV effect is checked against the intended saved main model.",
    "GSE48350 uses main age>=60: 25 Control/19 AD; no all-age sensitivity model is used.",
    "GSE173955 uses the new NCBI-count DESeq2 age/sex model, not the old author's LRT table.",
    "limma ordinary SE = sigma * stdev.unscaled for the actual diagnosis coefficient/contrast.",
    "Original limma P/FDR are moderated-test statistics; they are preserved as study evidence and are not recomputed from ordinary SE.",
    "DESeq2 effects and lfcSE are unshrunken Wald model quantities; no LR-to-SE approximation is used.",
    "No study significance, fold-change or minimum/maximum SE cutoff is applied.",
    "Unavailable model tests (e.g. DESeq2 outlier handling) are excluded; missing adjusted P alone is allowed.",
    "The original HGNC snapshot is reused and copied with its checksum. No new annotation is downloaded.",
    "Approved HGNC records only; protein-coding means locus_group equals protein-coding gene.",
    "Where supplied, a unique source Entrez ID has priority. A conflict with an exact approved symbol is excluded.",
    "If the source ID is unresolved, symbol matching uses exact approved symbol, then a unique combined previous/alias candidate set.",
    "No forced uppercase, fuzzy matching, manual gene splitting, or lowest-P duplicate selection.",
    "Only fully identical CSV source rows are deduplicated. All other repeated within-study coding HGNC IDs are excluded.",
    "Entrez-ID-assisted mapping is now available for the new GSE173955 table and GSE48350; counts need not equal the earlier symbol-only harmonization.",
    "All conversions and identity conflicts are retained in the audit; alias-only candidates still warrant identity review before drug queries.",
    "Study P values/FDR remain unchanged. No FDR is recalculated after protein-coding filtering.",
    "The common gene set is recomputed, not assumed to be the old 14,725-gene universe.",
    "Use four_studies_common_meta_input_long.csv or four_studies_common_effects_SE_wide.csv for the next meta-analysis stage.")
  writeLines(rules, file.path(out, "README_analysis_rules.txt"))
  csv(data.frame(Package = required_packages, Version = vapply(required_packages,
    function(p) as.character(utils::packageVersion(p)), character(1))), "package_versions.csv")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  saveRDS(list(results = results, common_inputs = common_inputs, common_genes = common,
    audit = audit, dataset_summary = dataset_summary, coverage = coverage, model_checks = se_checks,
    hgnc_file = hgnc_file, hgnc_md5 = unname(tools::md5sum(hgnc_file))), file.path(out, "harmonized_meta_inputs.rds"))
  message("Finished. Common genes with valid effects/SE in all four studies: ", nrow(common))
  message("Results: ", normalizePath(out, winslash = "/"))
  invisible(list(output_dir = normalizePath(out, winslash = "/"), dataset_summary = dataset_summary,
    coverage_summary = coverage_summary, common_genes = common, common_inputs = common_inputs,
    results = results, model_checks = se_checks))
}
