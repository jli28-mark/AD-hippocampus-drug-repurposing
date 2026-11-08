# AD37 clinical trial and PubMed evidence retrieval; version 1.1.0 (2026-11-05)
# Searches all 37 original drugs; separately flags the non-PET drugs with >=2 seed sets (Z < -1.5).
# DATA SOURCES: ClinicalTrials.gov API v2 and NCBI PubMed E-utilities.
# THIS IS A SEARCH/REVIEW WORKFLOW, NOT AN AUTOMATED DETERMINATION OF DRUG EFFICACY.
# Drug identity (and aliases), clinical-trial relevance, and direction of evidence require review.
# Source files are read-only. Cached responses prevent repeat downloads unless REFRESH <- TRUE.

BASE_DIR <- "D:/My Desktop/drug/6.AD network"
REFRESH <- FALSE                     # TRUE to obtain current records instead of cached records
MAX_CT_PER_QUERY <- 500L            # safety limit; truncation is reported
MAX_PUBMED_PER_DRUG <- 600L         # safety limit; truncation is reported
PUBMED_BATCH_SIZE <- 100L
NCBI_EMAIL <- ""                     # OPTIONAL: your actual email; never invent an address
NCBI_API_KEY <- Sys.getenv("NCBI_API_KEY", unset = "")  # optional, stored as environment variable

packages <- c("httr2", "jsonlite", "xml2", "readxl")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
if (!all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Missing packages: httr2, jsonlite, xml2, readxl")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
str1 <- function(x, sep = " | ") {
  if (is.null(x) || length(x) == 0L) return("")
  x <- unlist(x, use.names = FALSE)
  if (!length(x)) return("")
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  paste(unique(x), collapse = sep)
}
uni <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}
join <- function(x) paste(uni(x), collapse = " | ")
write_data <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
bind_frames <- function(x) {
  x <- Filter(function(y) is.data.frame(y) && nrow(y) > 0L, x)
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}
get_path <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
req_url <- function(url, query = list()) {
  req <- httr2::request(url)
  if (length(query)) req <- do.call(httr2::req_url_query, c(list(req), query))
  req <- httr2::req_user_agent(req, "AD37EvidenceSearch/1.1 (academic literature review)")
  req <- httr2::req_timeout(req, 65)
  httr2::req_retry(req, max_tries = 4)
}
api_json <- function(url, query = list(), ncbi = FALSE) {
  if (ncbi) Sys.sleep(0.42)  # stay below NCBI 3-requests-per-second limit
  httr2::resp_body_json(httr2::req_perform(req_url(url, query)), simplifyVector = FALSE)
}
api_xml <- function(url, query = list(), ncbi = FALSE) {
  if (ncbi) Sys.sleep(0.42)
  res <- httr2::req_perform(req_url(url, query))
  xml2::read_xml(httr2::resp_body_raw(res))
}

# The manually compiled Excel has exactly DrugBank_Name, Queried_ATC and BBB columns.
# It is descriptive input; BBB+/BBB- is NOT automatically considered experimental evidence.
read_manual_atc_bbb <- function(base_dir, candidates) {
  files <- file.path(base_dir, c("37ATC(1).xlsx", "37ATC.xlsx"))
  files <- files[file.exists(files)]
  if (!length(files)) stop("Missing 37ATC(1).xlsx or 37ATC.xlsx in ", base_dir)
  x <- as.data.frame(readxl::read_excel(files[1], sheet = 1, col_types = "text",
                                          .name_repair = "minimal"), stringsAsFactors = FALSE)
  required <- c("DrugBank_Name", "Queried_ATC", "BBB")
  if (!all(required %in% names(x))) stop("Excel columns must include: ", paste(required, collapse = ", "))
  clean <- function(z) gsub("\\s+", " ", tolower(trimws(as.character(z))), perl = TRUE)
  # Optional future provenance columns are accepted without being required.
  optional <- c("BBB_Source", "BBB_Evidence_type", "BBB_Reference", "BBB_Probability")
  for (col in optional) if (!col %in% names(x)) x[[col]] <- ""
  x <- x[, c(required, optional), drop = FALSE]
  if (nrow(x) != 37L || anyDuplicated(clean(x$DrugBank_Name)))
    stop("Manual Excel must contain 37 distinct drug names")
  ii <- match(clean(candidates$DrugBank_Name), clean(x$DrugBank_Name))
  if (anyNA(ii)) stop("Unmatched manual Excel names: ",
                     paste(candidates$DrugBank_Name[is.na(ii)], collapse = ", "))
  x <- x[ii, , drop = FALSE]
  labels <- toupper(trimws(as.character(x$BBB)))
  labels[is.na(labels) | !nzchar(labels)] <- ""
  if (any(!labels %in% c("BBB+", "BBB-", "")))
    stop("Unexpected BBB labels in Excel: ", paste(unique(labels[!labels %in% c("BBB+", "BBB-", "")]), collapse = ", "))
  ans <- data.frame(
    DrugBank_ID = candidates$DrugBank_ID,
    Drug = candidates$DrugBank_Name,
    ATC_Manual_raw = trimws(as.character(x$Queried_ATC)),
    BBB_Manual_Label = labels,
    BBB_Source = trimws(as.character(x$BBB_Source)),
    BBB_Evidence_type = ifelse(is.na(x$BBB_Evidence_type) | !nzchar(trimws(x$BBB_Evidence_type)),
                               "Unspecified_manual", trimws(as.character(x$BBB_Evidence_type))),
    BBB_Reference = trimws(as.character(x$BBB_Reference)),
    BBB_Probability = trimws(as.character(x$BBB_Probability)),
    BBB_Verification_status = ifelse(is.na(x$BBB_Source) | !nzchar(trimws(x$BBB_Source)),
                                     "NEEDS_SOURCE_REVIEW", "SOURCE_PROVIDED_NEEDS_VERIFICATION"),
    stringsAsFactors = FALSE)
  list(rows = ans, source_file = files[1])
}

# A compact, editable alias list. Add aliases only after confirming exact drug identity.
default_aliases <- function(candidates) {
  extra <- list(
    "DB00746" = c("desferrioxamine", "desferoxamine"),
    "DB00668" = "adrenaline",
    "DB00169" = "vitamin D3",
    "DB02701" = "niacinamide",
    "DB09148" = "florbetaben",
    "DB09149" = "florbetapir",
    "DB09151" = "flutemetamol",
    "DB00071" = "porcine insulin"
  )
  do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) {
    id <- candidates$DrugBank_ID[i]
    x <- uni(c(candidates$DrugBank_Name[i], extra[[id]]))
    data.frame(DrugBank_ID = id, Drug = candidates$DrugBank_Name[i],
               Search_name = x, Enabled = "YES", stringsAsFactors = FALSE)
  }))
}

# ClinicalTrials.gov study parser. All text from registry requires eligibility review.
trial_record <- function(study, id, drug, alias, condition, queried_at) {
  p <- study$protocolSection %||% list()
  idmod <- p$identificationModule %||% list()
  stat <- p$statusModule %||% list()
  cond <- p$conditionsModule %||% list()
  design <- p$designModule %||% list()
  arms <- p$armsInterventionsModule %||% list()
  outcomes <- p$outcomesModule %||% list()
  refs <- p$referencesModule$references %||% list()
  interventions <- arms$interventions %||% list()
  pmids <- uni(vapply(refs, function(z) str1(z$pmid), character(1)))
  int_names <- vapply(interventions, function(z) str1(z$name), character(1))
  int_types <- vapply(interventions, function(z) str1(z$type), character(1))
  primary <- outcomes$primaryOutcomes %||% list()
  nct <- str1(idmod$nctId)
  data.frame(
    DrugBank_ID = id, Drug = drug,
    NCT_ID = nct, Query_alias = alias, Query_condition = condition,
    Study_title = str1(idmod$briefTitle),
    Overall_status = str1(stat$overallStatus),
    Study_type = str1(design$studyType),
    Phase = str1(design$phases),
    Enrollment = str1(design$enrollmentInfo$count),
    Start_date = str1(stat$startDateStruct$date),
    Completion_date = str1(stat$completionDateStruct$date),
    Conditions = str1(cond$conditions),
    Interventions = join(int_names),
    Intervention_types = join(int_types),
    Primary_outcome_names = join(vapply(primary, function(z) str1(z$measure), character(1))),
    Registry_has_results = isTRUE(study$hasResults) || !is.null(study$resultsSection),
    Registry_reference_PMIDs = paste(pmids, collapse = " | "),
    Study_URL = if (nzchar(nct)) paste0("https://clinicaltrials.gov/study/", nct) else "",
    Queried_at_UTC = queried_at,
    stringsAsFactors = FALSE
  )
}

fetch_trials <- function(id, drug, aliases, max_results = MAX_CT_PER_QUERY) {
  base <- "https://clinicaltrials.gov/api/v2/studies"
  conditions <- c("Alzheimer Disease", "Mild Cognitive Impairment")
  records <- list(); errors <- character(); any_truncated <- FALSE
  searched <- 0L; counts <- integer()
  queried_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  for (alias in aliases) for (condition in conditions) {
    q_name <- paste(alias, "/", condition)
    tryCatch({
      token <- ""; retrieved <- 0L; pages <- 0L
      repeat {
        q <- list("query.cond" = condition, "query.intr" = alias,
                  "pageSize" = "100", "countTotal" = "true", "format" = "json")
        if (nzchar(token)) q$pageToken <- token
        z <- api_json(base, q)
        searched <- searched + 1L; pages <- pages + 1L
        if (pages == 1L) counts[q_name] <- as.integer(z$totalCount %||% 0L)
        studies <- z$studies %||% list()
        if (length(studies)) for (s in studies) {
          records[[length(records) + 1L]] <- trial_record(s, id, drug, alias, condition, queried_at)
        }
        retrieved <- retrieved + length(studies)
        token <- str1(z$nextPageToken)
        if (!nzchar(token)) break
        if (retrieved >= max_results) { any_truncated <- TRUE; break }
      }
    }, error = function(e) {
      errors <<- c(errors, paste(q_name, conditionMessage(e), sep = ": "))
    })
  }
  out <- bind_frames(records)
  if (nrow(out)) {
    out <- out[nzchar(out$NCT_ID), , drop = FALSE]
    if (nrow(out)) {
      groups <- split(seq_len(nrow(out)), out$NCT_ID)
      out <- do.call(rbind, lapply(groups, function(ix) {
        row <- out[ix[1], , drop = FALSE]
        row$Query_alias <- join(out$Query_alias[ix])
        row$Query_condition <- join(out$Query_condition[ix])
        row$Registry_reference_PMIDs <- join(unlist(strsplit(
          out$Registry_reference_PMIDs[ix], "\\s*\\|\\s*", perl = TRUE)))
        row
      }))
      rownames(out) <- NULL
    }
  }
  list(rows = out,
       status = if (length(errors)) "PARTIAL_OR_ERROR" else "COMPLETE",
       errors = errors, truncated = any_truncated,
       total_search_hits = sum(counts), queried_at = queried_at)
}

# PubMed queries are strict about drug name in title/abstract and AD/MCI relevance.
# For older studies with only MeSH-indexed drug names, a secondary manual search may be needed.
pubmed_term <- function(aliases) {
  aliases <- gsub('"', '', aliases, fixed = TRUE)
  drug_clause <- paste0('"', aliases, '"[Title/Abstract]', collapse = " OR ")
  disease <- paste(
    '"Alzheimer Disease"[MeSH Terms]',
    'alzheimer*[Title/Abstract]',
    '"Mild Cognitive Impairment"[MeSH Terms]',
    '"mild cognitive impairment"[Title/Abstract]',
    sep = " OR ")
  paste0("(", drug_clause, ") AND (", disease, ")")
}
xml_value <- function(node, path) {
  z <- xml2::xml_find_first(node, path)
  if (inherits(z, "xml_missing")) return("")
  trimws(xml2::xml_text(z))
}
extract_nct <- function(text) {
  hits <- regmatches(text, gregexpr("NCT[0-9]{8}", text, perl = TRUE, ignore.case = TRUE))[[1]]
  if (!length(hits) || identical(hits, character(0))) return("")
  join(toupper(hits))
}
pubmed_article <- function(article, id, drug, query, queried_at) {
  pmid <- xml_value(article, "./MedlineCitation/PMID")
  title <- xml_value(article, "./MedlineCitation/Article/ArticleTitle")
  journal <- xml_value(article, "./MedlineCitation/Article/Journal/Title")
  year <- xml_value(article, "./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year")
  if (!nzchar(year)) {
    medline_date <- xml_value(article, "./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate")
    year <- if (grepl("[12][0-9]{3}", medline_date))
      regmatches(medline_date, regexpr("[12][0-9]{3}", medline_date)) else ""
  }
  if (!nzchar(year)) year <- xml_value(article, "./MedlineCitation/Article/ArticleDate/Year")
  doi <- xml_value(article, ".//ArticleIdList/ArticleId[@IdType='doi']")
  pts <- xml2::xml_find_all(article, "./MedlineCitation/Article/PublicationTypeList/PublicationType")
  absnodes <- xml2::xml_find_all(article, "./MedlineCitation/Article/Abstract/AbstractText")
  abstract <- if (length(absnodes)) paste(xml2::xml_text(absnodes), collapse = " ") else ""
  # Search full XML record to capture NCT numbers in data-bank accession fields, not only abstracts.
  ncts <- extract_nct(xml2::xml_text(article))
  data.frame(
    DrugBank_ID = id, Drug = drug, PMID = pmid, Title = title,
    Year = year, Journal = journal, DOI = doi,
    Publication_types = join(xml2::xml_text(pts)),
    Abstract_excerpt = if (nzchar(abstract)) substr(abstract, 1L, 220L) else "",
    Reported_NCT_IDs = ncts,
    PubMed_URL = if (nzchar(pmid)) paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/") else "",
    Search_query = query, Queried_at_UTC = queried_at,
    stringsAsFactors = FALSE)
}
fetch_pubmed <- function(id, drug, aliases, max_records = MAX_PUBMED_PER_DRUG) {
  prefix <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  query <- pubmed_term(aliases)
  meta <- list(db = "pubmed", term = query, retmode = "json", retmax = as.character(max_records),
               retstart = "0", sort = "pub date", tool = "AD37EvidenceSearch")
  if (nzchar(NCBI_EMAIL)) meta$email <- NCBI_EMAIL
  if (nzchar(NCBI_API_KEY)) meta$api_key <- NCBI_API_KEY
  queried_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  z <- api_json(paste0(prefix, "esearch.fcgi"), meta, ncbi = TRUE)
  es <- z$esearchresult
  if (is.null(es)) stop("PubMed ESearch returned no esearchresult")
  ids <- uni(es$idlist %||% character())
  count <- as.integer(es$count %||% 0L)
  records <- list()
  if (length(ids)) {
    batches <- split(ids, ceiling(seq_along(ids) / PUBMED_BATCH_SIZE))
    for (batch in batches) {
      p <- list(db = "pubmed", id = paste(batch, collapse = ","), retmode = "xml",
                tool = "AD37EvidenceSearch")
      if (nzchar(NCBI_EMAIL)) p$email <- NCBI_EMAIL
      if (nzchar(NCBI_API_KEY)) p$api_key <- NCBI_API_KEY
      doc <- api_xml(paste0(prefix, "efetch.fcgi"), p, ncbi = TRUE)
      articles <- xml2::xml_find_all(doc, ".//PubmedArticle")
      for (a in articles) {
        records[[length(records) + 1L]] <- pubmed_article(a, id, drug, query, queried_at)
      }
    }
  }
  out <- bind_frames(records)
  if (nrow(out)) {
    out <- out[nzchar(out$PMID) & !duplicated(out$PMID), , drop = FALSE]
    rownames(out) <- NULL
  }
  list(rows = out, status = "COMPLETE", total_found = count,
       retrieved = nrow(out), truncated = count > length(ids),
       query = query, queried_at = queried_at)
}

# Cross-reference registry citations and PubMed NCT mentions without assuming efficacy.
make_links <- function(trials, papers) {
  if (!nrow(trials) || !nrow(papers)) return(data.frame())
  rows <- list()
  for (i in seq_len(nrow(trials))) {
    t <- trials[i, ]
    p <- papers[papers$DrugBank_ID == t$DrugBank_ID, , drop = FALSE]
    if (!nrow(p)) next
    ref_ids <- uni(unlist(strsplit(t$Registry_reference_PMIDs, "\\s*\\|\\s*", perl = TRUE)))
    pmid_match <- p$PMID %in% ref_ids
    nct_match <- vapply(p$Reported_NCT_IDs, function(x) {
      t$NCT_ID %in% uni(unlist(strsplit(x, "\\s*\\|\\s*", perl = TRUE)))
    }, logical(1))
    idx <- which(pmid_match | nct_match)
    if (length(idx)) for (j in idx) {
      rows[[length(rows) + 1L]] <- data.frame(
        DrugBank_ID = t$DrugBank_ID, Drug = t$Drug,
        NCT_ID = t$NCT_ID, PMID = p$PMID[j],
        Link_basis = paste(c(if (pmid_match[j]) "Registry_reference_PMID",
                             if (nct_match[j]) "PubMed_reports_NCT_ID"), collapse = " | "),
        Study_URL = t$Study_URL,
        PubMed_URL = p$PubMed_URL[j], stringsAsFactors = FALSE)
    }
  }
  out <- bind_frames(rows)
  if (nrow(out)) out <- out[!duplicated(out[c("DrugBank_ID", "NCT_ID", "PMID")]), ]
  out
}

run_AD37_evidence <- function(base_dir = BASE_DIR) {
  candidates_files <- file.path(base_dir, c("14_nominal_P_candidates.csv", "14_nominal_P_candidates(1).csv"))
  candidates_files <- candidates_files[file.exists(candidates_files)]
  if (!length(candidates_files)) stop("Copy 14_nominal_P_candidates.csv into: ", base_dir)
  candidates <- utils::read.csv(candidates_files[1], stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("DrugBank_ID", "DrugBank_Name") %in% names(candidates))) stop("Missing DrugBank_ID or DrugBank_Name")
  if (nrow(candidates) != 37L || anyDuplicated(candidates$DrugBank_ID)) stop("Expected 37 unique DrugBank IDs")

  results_root <- file.path(base_dir, "AD37_results")
  runs <- if (dir.exists(results_root)) list.dirs(results_root, recursive = FALSE, full.names = TRUE) else character()
  runs <- runs[file.exists(file.path(runs, "06_network_proximity_wide_and_counts.csv"))]
  if (!length(runs)) stop("Completed AD37 network results not found under AD37_results")
  mtime <- vapply(runs, function(r) as.numeric(file.info(file.path(r, "06_network_proximity_wide_and_counts.csv"))$mtime), numeric(1))
  run <- runs[which.max(mtime)]
  network <- utils::read.csv(file.path(run, "06_network_proximity_wide_and_counts.csv"), stringsAsFactors = FALSE)
  if (anyDuplicated(network$DrugBank_ID)) stop("Duplicate DrugBank IDs in network table")
  ind <- match(candidates$DrugBank_ID, network$DrugBank_ID)
  if (anyNA(ind)) stop("37 candidates do not all match previous network results")
  network <- network[ind, , drop = FALSE]
  manual <- read_manual_atc_bbb(base_dir, candidates)
  atc_bbb <- manual$rows

  out <- file.path(run, "12_ClinicalTrials_PubMed")
  cache_dir <- file.path(out, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  write_data(atc_bbb, file.path(out, "00_manual_ATC_BBB_input_snapshot.csv"))
  # One master table, by DrugBank ID. No local network recalculation.
  master <- data.frame(
    DrugBank_ID = candidates$DrugBank_ID, Drug = candidates$DrugBank_Name,
    ATC_Manual_raw = atc_bbb$ATC_Manual_raw,
    BBB_Manual_Label = atc_bbb$BBB_Manual_Label,
    BBB_Source = atc_bbb$BBB_Source,
    BBB_Evidence_type = atc_bbb$BBB_Evidence_type,
    BBB_Reference = atc_bbb$BBB_Reference,
    BBB_Probability = atc_bbb$BBB_Probability,
    BBB_Verification_status = atc_bbb$BBB_Verification_status,
    Meta_AD_Z = network$Meta_AD_Z, AD_Z = network$AD_Z,
    Amyloid_Z = network$Amyloid_Z, Tau_Z = network$Tau_Z,
    N_modules_Z_below_minus1.5 = network$N_modules_Z_below_minus1.5,
    stringsAsFactors = FALSE)
  k <- suppressWarnings(as.integer(master$N_modules_Z_below_minus1.5))
  master$Priority13_nonPET_2plus <- !is.na(k) & k >= 2L &
    !master$DrugBank_ID %in% c("DB09148", "DB09149", "DB09151")
  write_data(master, file.path(out, "00_AD37_master_network_ATC_BBB.csv"))
  aliases_file <- file.path(base_dir, "AD37_drug_search_aliases.csv")
  if (!file.exists(aliases_file)) {
    write_data(default_aliases(candidates), aliases_file)
    message("Wrote editable drug alias table: ", aliases_file)
  }
  aliases <- utils::read.csv(aliases_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("DrugBank_ID", "Search_name", "Enabled")
  if (!all(required %in% names(aliases))) stop("Alias file needs columns: ", paste(required, collapse = ", "))
  aliases$Enabled <- toupper(trimws(as.character(aliases$Enabled)))
  aliases <- aliases[aliases$Enabled %in% c("YES", "TRUE", "1"), , drop = FALSE]
  if (any(!aliases$DrugBank_ID %in% candidates$DrugBank_ID)) stop("Alias table contains an unknown DrugBank ID")
  write_data(aliases, file.path(out, "00_aliases_used.csv"))

  results <- vector("list", nrow(candidates))
  failures <- list()
  for (i in seq_len(nrow(candidates))) {
    id <- candidates$DrugBank_ID[i]; drug <- candidates$DrugBank_Name[i]
    ali <- uni(c(drug, aliases$Search_name[aliases$DrugBank_ID == id]))
    cache_file <- file.path(cache_dir, paste0(id, ".rds"))
    message(sprintf("[%d/37] %s (%s)", i, drug, id))
    if (file.exists(cache_file) && !isTRUE(REFRESH)) {
      z <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (!is.null(z) && identical(z$aliases, ali) &&
          identical(z$clinical$status, "COMPLETE") &&
          identical(z$pubmed$status, "COMPLETE")) {
        results[[i]] <- z
        message("  Reused cached results")
        next
      }
    }
    ct <- tryCatch(fetch_trials(id, drug, ali), error = function(e)
      list(rows = data.frame(), status = "ERROR", errors = conditionMessage(e),
           truncated = NA, total_search_hits = NA_integer_, queried_at = ""))
    pm <- tryCatch(fetch_pubmed(id, drug, ali), error = function(e)
      list(rows = data.frame(), status = "ERROR", total_found = NA_integer_,
           retrieved = NA_integer_, truncated = NA, query = pubmed_term(ali),
           queried_at = "", error = conditionMessage(e)))
    z <- list(id = id, drug = drug, aliases = ali, clinical = ct, pubmed = pm)
    saveRDS(z, cache_file)
    results[[i]] <- z
    if (ct$status != "COMPLETE") failures[[length(failures) + 1L]] <- data.frame(
      DrugBank_ID = id, Drug = drug, Source = "ClinicalTrials.gov",
      Details = str1(ct$errors))
    if (pm$status != "COMPLETE") failures[[length(failures) + 1L]] <- data.frame(
      DrugBank_ID = id, Drug = drug, Source = "PubMed", Details = str1(pm$error))
  }
  trials <- bind_frames(lapply(results, function(x) x$clinical$rows))
  papers <- bind_frames(lapply(results, function(x) x$pubmed$rows))
  links <- make_links(trials, papers)
  failures <- bind_frames(failures)
  if (nrow(trials)) trials <- trials[order(trials$DrugBank_ID, trials$NCT_ID), , drop = FALSE]
  if (nrow(papers)) papers <- papers[order(papers$DrugBank_ID, papers$Year, decreasing = FALSE), , drop = FALSE]
  write_data(trials, file.path(out, "01_ClinicalTrials_all_records.csv"))
  write_data(papers, file.path(out, "02_PubMed_all_records.csv"))
  write_data(links, file.path(out, "03_registry_publication_crosslinks.csv"))
  write_data(failures, file.path(out, "06_query_errors.csv"))

  sumrows <- lapply(seq_len(nrow(candidates)), function(i) {
    z <- results[[i]]; id <- z$id
    t <- if (nrow(trials)) trials[trials$DrugBank_ID == id, , drop = FALSE] else data.frame()
    p <- if (nrow(papers)) papers[papers$DrugBank_ID == id, , drop = FALSE] else data.frame()
    l <- if (nrow(links)) links[links$DrugBank_ID == id, , drop = FALSE] else data.frame()
    row <- network[i, , drop = FALSE]
    k <- suppressWarnings(as.integer(row$N_modules_Z_below_minus1.5))
    flag <- !is.na(k) && k >= 2L && !id %in% c("DB09148", "DB09149", "DB09151")
    data.frame(
      DrugBank_ID = id, Drug = z$drug,
      ATC_Manual_raw = atc_bbb$ATC_Manual_raw[i],
      BBB_Manual_Label = atc_bbb$BBB_Manual_Label[i],
      BBB_Source = atc_bbb$BBB_Source[i],
      BBB_Evidence_type = atc_bbb$BBB_Evidence_type[i],
      BBB_Reference = atc_bbb$BBB_Reference[i],
      BBB_Probability = atc_bbb$BBB_Probability[i],
      BBB_Verification_status = atc_bbb$BBB_Verification_status[i],
      Meta_AD_Z = row$Meta_AD_Z, AD_Z = row$AD_Z,
      Amyloid_Z = row$Amyloid_Z, Tau_Z = row$Tau_Z,
      N_modules_Z_below_minus1.5 = k,
      Priority13_nonPET_2plus = flag,
      NCT_IDs = if (nrow(t)) join(t$NCT_ID) else "",
      Trial_phases_observed = if (nrow(t)) join(t$Phase) else "",
      Trial_statuses_observed = if (nrow(t)) join(t$Overall_status) else "",
      PMIDs = if (nrow(p)) join(p$PMID) else "",
      N_ClinicalTrials_records = nrow(t),
      N_ClinicalTrials_with_posted_results = if (nrow(t)) sum(t$Registry_has_results) else 0L,
      N_PubMed_records_retrieved = nrow(p),
      N_PubMed_database_hits = z$pubmed$total_found,
      N_explicit_NCT_PMID_links = nrow(l),
      ClinicalTrials_search_status = z$clinical$status,
      ClinicalTrials_truncated = z$clinical$truncated,
      PubMed_search_status = z$pubmed$status,
      PubMed_truncated = z$pubmed$truncated,
      ClinicalTrials_queried_at_UTC = z$clinical$queried_at,
      PubMed_queried_at_UTC = z$pubmed$queried_at,
      stringsAsFactors = FALSE)
  })
  summary <- do.call(rbind, sumrows)
  rownames(summary) <- NULL
  write_data(summary, file.path(out, "04_evidence_summary_all37.csv"))
  write_data(summary[summary$Priority13_nonPET_2plus, , drop = FALSE],
             file.path(out, "04a_priority13_evidence_summary.csv"))

  # These fields are deliberately EMPTY: screening / direction of evidence must be verified.
  review <- summary[, c("DrugBank_ID", "Drug", "Priority13_nonPET_2plus",
                         "N_ClinicalTrials_records", "N_PubMed_records_retrieved"), drop = FALSE]
  review$Clinical_trial_screening <- "PENDING_REVIEW"
  review$Relevant_NCT_IDs <- ""
  review$Clinical_trial_results <- "NOT_ASSESSED"
  review$Preclinical_studies_screening <- "PENDING_REVIEW"
  review$Relevant_PMIDs <- ""
  review$Preclinical_evidence_direction <- "NOT_ASSESSED"
  review$Safety_or_exposure_notes <- ""
  review$BBB_Source_verified <- ""
  review$BBB_Evidence_type_verified <- ""
  review$BBB_Reference_verified <- ""
  review$Reviewer <- ""
  review$Review_date <- ""
  review_file <- file.path(out, "05_manual_evidence_review_all37.csv")
  if (!file.exists(review_file)) write_data(review, review_file)
  else message("Preserved existing manual review: ", review_file)

  # Record-level screening files; do not infer outcome direction from titles or PubMed counts.
  if (nrow(trials)) {
    treview <- trials[, c("DrugBank_ID", "Drug", "NCT_ID", "Study_title", "Conditions",
                          "Interventions", "Phase", "Overall_status", "Study_URL"), drop = FALSE]
    treview$AD_MCI_relevant <- "PENDING_REVIEW"
    treview$Candidate_is_study_intervention <- "PENDING_REVIEW"
    treview$Trial_outcome_direction <- "NOT_ASSESSED"
    treview$Review_notes <- ""
  } else treview <- data.frame(DrugBank_ID=character(), Drug=character(), NCT_ID=character(),
                               AD_MCI_relevant=character(), Candidate_is_study_intervention=character(),
                               Trial_outcome_direction=character(), Review_notes=character())
  trial_review_file <- file.path(out, "05a_trial_record_review.csv")
  if (!file.exists(trial_review_file)) write_data(treview, trial_review_file)
  else message("Preserved existing manual trial review: ", trial_review_file)
  if (nrow(papers)) {
    preview <- papers[, c("DrugBank_ID", "Drug", "PMID", "Title", "Year",
                          "Publication_types", "PubMed_URL"), drop = FALSE]
    preview$AD_MCI_relevant <- "PENDING_REVIEW"
    preview$Study_model <- "NOT_ASSESSED"
    preview$Evidence_direction <- "NOT_ASSESSED"
    preview$Review_notes <- ""
  } else preview <- data.frame(DrugBank_ID=character(), Drug=character(), PMID=character(),
                               AD_MCI_relevant=character(), Study_model=character(),
                               Evidence_direction=character(), Review_notes=character())
  pubmed_review_file <- file.path(out, "05b_pubmed_record_review.csv")
  if (!file.exists(pubmed_review_file)) write_data(preview, pubmed_review_file)
  else message("Preserved existing manual PubMed review: ", pubmed_review_file)

  manifest <- data.frame(
    Property = c("Script_version", "Source_candidates", "Source_network", "Source_manual_ATC_BBB",
                 "ClinicalTrials_endpoint", "PubMed_endpoints", "Search_conditions", "PubMed_query_rule",
                 "ClinicalTrials_cap_per_alias_condition", "PubMed_cap_per_drug",
                 "NCBI_rate_limit_handling", "Data_caching", "Retrieved_at_UTC",
                 "Analysis_warning"),
    Value = c("1.1.0", normalizePath(candidates_files[1], winslash = "/"),
              normalizePath(file.path(run, "06_network_proximity_wide_and_counts.csv"), winslash = "/"),
              normalizePath(manual$source_file, winslash = "/"),
              "https://clinicaltrials.gov/api/v2/studies",
              "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi and efetch.fcgi",
              "Alzheimer Disease; Mild Cognitive Impairment",
              "Drug name/aliases in title/abstract AND AD/MCI disease terms",
              as.character(MAX_CT_PER_QUERY), as.character(MAX_PUBMED_PER_DRUG),
              "At most 1 request per 0.42 s; retry transient failures",
              paste("REFRESH =", REFRESH),
              format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
              "Search hits are NOT confirmed efficacy. Manual BBB labels lack source information until reviewed."))
  write_data(manifest, file.path(out, "00_search_settings.csv"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  if (nrow(failures) || any(summary$ClinicalTrials_truncated %in% TRUE) ||
      any(summary$PubMed_truncated %in% TRUE)) {
    writeLines("PARTIAL RETRIEVAL: inspect 06_query_errors.csv and truncation fields in 04_evidence_summary_all37.csv.",
               file.path(out, "INCOMPLETE_OR_TRUNCATED_SEARCH.txt"))
  } else {
    unlink(file.path(out, "INCOMPLETE_OR_TRUNCATED_SEARCH.txt"))
  }
  message("Saved to: ", out)
  print(summary[, c("Drug", "Priority13_nonPET_2plus", "N_ClinicalTrials_records",
                    "N_PubMed_records_retrieved", "ClinicalTrials_search_status",
                    "PubMed_search_status")], row.names = FALSE)
  invisible(out)
}

# Source this file to run. The original 37-drug network results are not changed.
AD37_EVIDENCE_DIR <- run_AD37_evidence(BASE_DIR)
