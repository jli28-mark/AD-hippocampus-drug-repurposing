# endoAD full-library screen of the fixed 109-gene AD expression signature
# Version 1.0.0. Source this file in RStudio; no other R script is required.
# Inputs are read-only. Each run has a new result directory.
# Cached XML metadata, network distances and permutation chunks are reused
# automatically when the corresponding inputs and settings are unchanged.
# Primary statistic: Fang et al., Nature Aging 2021, Eq. 2 (bidirectional closest).
# DOI: 10.1038/s43587-021-00138-z
# This is an adapted gene-set screen, not a replication of the full endoAD study.

BASE_DIR <- "D:/My Desktop/drug/5.endoAD_109_screen"
N_PERM <- 10000L
N_WORKERS <- 2L
RANDOM_SEED <- 20261105L
MIN_BIN_SIZE <- 100L

.e109_version <- "1.0.0"
.e109_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
.e109_col <- function(x, choices, required = FALSE) {
  j <- match(.e109_key(choices), .e109_key(names(x))); j <- j[!is.na(j)]
  if (length(j)) return(names(x)[j[1L]])
  if (required) stop("Missing column: ", paste(choices, collapse = " / "),
    ". Available columns: ", paste(names(x), collapse = ", "), call. = FALSE)
  NULL
}
.e109_trim <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
.e109_present <- function(x) nzchar(.e109_trim(x)) & !toupper(.e109_trim(x)) %in% c("NA", "N/A", "NULL")
.e109_join <- function(x) paste(sort(unique(.e109_trim(x)[.e109_present(x)])), collapse = " | ")
.e109_split <- function(x) {
  z <- trimws(unlist(strsplit(paste(x[!is.na(x)], collapse = "|"), "[|;]")))
  unique(z[.e109_present(z)])
}
.e109_flag <- function(x) {
  z <- tolower(.e109_trim(x))
  ifelse(z %in% c("true", "1", "yes"), TRUE,
    ifelse(z %in% c("false", "0", "no"), FALSE, NA))
}
.e109_df <- function(x, template = data.frame()) {
  if (!length(x)) return(template)
  as.data.frame(data.table::rbindlist(x, fill = TRUE, use.names = TRUE), stringsAsFactors = FALSE)
}
.e109_csv <- function(x, file) {
  if (!ncol(x)) x <- data.frame(Note = character())
  utils::write.csv(x, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
.e109_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, compress = FALSE, version = 3); unname(tools::md5sum(f))
}
.e109_save <- function(x, file, compress = FALSE) {
  tmp <- tempfile("partial_", tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE); saveRDS(x, tmp, compress = compress, version = 3)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("Cannot finalize cache file: ", file)
}
.e109_read <- function(file, sep = "auto", header = TRUE, literal = FALSE) {
  head <- readLines(file, n = 3L, warn = FALSE, encoding = "UTF-8")
  if (any(grepl("<!doctype html|<html|github.*blob.*<", head, ignore.case = TRUE)))
    stop("This appears to be a webpage, not the raw data file: ", file, call. = FALSE)
  x <- data.table::fread(file, sep = sep, header = header, quote = if (literal) "" else "\"",
    colClasses = "character", na.strings = c("", "NA"), check.names = FALSE,
    data.table = FALSE, encoding = "UTF-8", fill = TRUE, showProgress = interactive())
  names(x) <- sub("^\ufeff", "", names(x))
  if (header && anyDuplicated(names(x))) stop("Duplicate column names: ", file)
  x
}
.e109_pairs <- function(file, kind) {
  x <- .e109_read(file, header = FALSE, literal = TRUE)
  if (ncol(x) < 2L || !nrow(x)) stop("Expected at least two columns in ", file)
  a <- .e109_key(.e109_trim(x[1L, 1L])); b <- .e109_key(.e109_trim(x[1L, 2L]))
  header <- if (kind == "drug") a %in% c("drug", "drugid", "drugbankid", "drugbank", "compound") else
    a %in% c("genea", "gene1", "geneid1", "entrez1", "entrezid1", "source", "sourceid", "proteina", "protein1", "node1", "nodea", "from", "a", "v1")
  header <- header || (a == "geneid" && b %in% c("geneid", "geneid2"))
  if (header) { nm <- make.unique(.e109_trim(unlist(x[1L, ], use.names = FALSE))); x <- x[-1L, , drop = FALSE]; names(x) <- nm }
  if (!nrow(x)) stop("No data rows in ", file)
  y <- data.frame(Source_row = seq_len(nrow(x)) + as.integer(header),
    ID1 = .e109_trim(x[[1L]]), ID2 = .e109_trim(x[[2L]]), stringsAsFactors = FALSE)
  if (ncol(x) > 2L) for (j in 3:ncol(x)) y[[paste0("Source_extra_", j)]] <- x[[j]]
  if (kind == "drug" && any(!grepl("^DB[0-9]{5}$", y$ID1)))
    stop("DrugTargetNetwork first column must contain DrugBank IDs. Check the file/header; it was not guessed or swapped.")
  attr(y, "source_header_detected") <- header
  y
}

.e109_hgnc <- function(file) {
  x <- .e109_read(file, "\t")
  take <- function(choices, required = FALSE) {
    z <- .e109_col(x, choices, required)
    if (is.null(z)) rep(NA_character_, nrow(x)) else x[[z]]
  }
  h <- data.frame(HGNC_ID = take("hgnc_id", TRUE), Symbol = take("symbol", TRUE),
    Status = take("status", TRUE), Locus_group = take("locus_group", TRUE),
    Entrez = take(c("entrez_id", "entrezid")), Ensembl = take("ensembl_gene_id"),
    UniProt = take("uniprot_ids"), Previous = take("prev_symbol"), Alias = take("alias_symbol"))
  h <- h[!is.na(h$Status) & h$Status == "Approved", , drop = FALSE]
  if (!nrow(h) || anyNA(h$HGNC_ID) || anyDuplicated(h$HGNC_ID) || anyDuplicated(h$Symbol))
    stop("Invalid approved HGNC records.")
  index <- function(z, split = FALSE) {
    vv <- if (split) lapply(z, .e109_split) else lapply(z, function(a) a[.e109_present(a)])
    ii <- rep(seq_len(nrow(h)), lengths(vv)); zz <- unlist(vv, use.names = FALSE)
    list2env(split(ii, zz), hash = TRUE, parent = emptyenv())
  }
  list(table = h, hgnc = index(h$HGNC_ID), symbol = index(h$Symbol), entrez = index(h$Entrez),
    ensembl = index(h$Ensembl), uniprot = index(h$UniProt, TRUE),
    previous = index(h$Previous, TRUE), alias = index(h$Alias, TRUE))
}
.e109_map <- function(values, ref) {
  values <- unique(.e109_trim(values))
  .e109_df(lapply(values, function(v) {
    getids <- function(name, key = v) get0(key, ref[[name]], inherits = FALSE, ifnotfound = integer())
    k <- integer(); basis <- ""
    if (grepl("^HGNC:[0-9]+$", v)) { k <- getids("hgnc"); basis <- "HGNC_ID" }
    else if (grepl("^[0-9]+$", v)) { k <- getids("entrez"); basis <- "Entrez_ID" }
    else if (grepl("^ENSG[0-9]+(\\.[0-9]+)?$", v)) { k <- getids("ensembl", sub("\\.[0-9]+$", "", v)); basis <- "Ensembl_ID" }
    else if (nzchar(v)) for (nm in c("symbol", "previous", "alias", "uniprot")) {
      k <- getids(nm); basis <- nm; if (length(k)) break
    }
    k <- unique(k); status <- if (length(k) == 1L) "Mapped" else if (length(k)) "Ambiguous" else "Unmapped"
    coding <- length(k) == 1L && identical(ref$table$Locus_group[k], "protein-coding gene")
    data.frame(Input_ID = v, HGNC_ID = if (length(k) == 1L) ref$table$HGNC_ID[k] else NA_character_,
      gene_symbol = if (length(k) == 1L) ref$table$Symbol[k] else NA_character_,
      Mapping_status = status, Mapping_basis = basis, Protein_coding = coding)
  }))
}
.e109_gene_input <- function(x, ref) {
  cc <- .e109_col(x, c("HGNC_ID", "hgnc_accession", "gene_symbol", "GeneSymbol", "Name_GeneSymbol", "Gene", "Symbol", "HGNC_symbol", "approved_symbol", "entrez_id", "ensembl_gene_id"), TRUE)
  vals <- x[[cc]]
  sc <- .e109_col(x, c("gene_symbol", "GeneSymbol", "Name_GeneSymbol", "Gene", "Symbol", "HGNC_symbol", "approved_symbol"))
  if (!is.null(sc)) {
    empty <- !.e109_present(vals); vals[empty] <- x[[sc]][empty]
    both <- which(.e109_present(x[[cc]]) & .e109_present(x[[sc]]))
    a <- .e109_map(x[[cc]][both], ref); b <- .e109_map(x[[sc]][both], ref)
    ia <- match(.e109_trim(x[[cc]][both]), a$Input_ID); ib <- match(.e109_trim(x[[sc]][both]), b$Input_ID)
    conflict <- a$Mapping_status[ia] == "Mapped" & b$Mapping_status[ib] == "Mapped" & a$HGNC_ID[ia] != b$HGNC_ID[ib]
    if (any(conflict, na.rm = TRUE)) stop("Gene ID and gene symbol disagree in an input gene table.")
  }
  .e109_map(vals, ref)
}

.e109_degree_bins <- function(degrees, min_size) {
  groups <- list(); current <- integer()
  for (d in sort(unique(degrees))) {
    current <- c(current, which(degrees == d))
    if (length(current) >= min_size) { groups[[length(groups) + 1L]] <- current; current <- integer() }
  }
  if (length(current)) {
    if (!length(groups)) groups[[1L]] <- current else groups[[length(groups)]] <- c(groups[[length(groups)]], current)
  }
  bins <- integer(length(degrees)); for (i in seq_along(groups)) bins[groups[[i]]] <- i
  bins
}
.e109_distances <- function(g, cache_dir, block_size, max_gb, log) {
  n <- igraph::vcount(g); node_ids <- igraph::V(g)$name
  estimate <- as.double(n)^2 / 1024^3
  if (estimate > max_gb) stop(sprintf("Distance cache needs %.2f GiB before worker copies; exceeds max_distance_gb=%.2f.", estimate, max_gb))
  key <- .e109_hash(list(version = "unweighted_LCC_raw_distances_v1", nodes = node_ids, edges = igraph::as_edgelist(g)))
  folder <- file.path(cache_dir, key); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  whole <- file.path(folder, "distance_matrix_raw.rds")
  if (file.exists(whole)) {
    z <- tryCatch(readRDS(whole), error = function(e) NULL)
    if (!is.null(z) && identical(z$nodes, node_ids) && is.raw(z$D) && identical(dim(z$D), c(n, n))) {
      log("Using cached all-pairs distances."); return(list(D = z$D, key = key, file = whole))
    }
    rm(z); gc(FALSE)
  }
  log(sprintf("Building raw distance matrix: %d nodes, %.2f GiB plus temporary blocks.", n, estimate))
  D <- matrix(as.raw(0), n, n)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
  for (k in seq_along(chunks)) {
    ix <- chunks[[k]]; f <- file.path(folder, sprintf("block_%05d.rds", k))
    z <- if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL
    valid <- !is.null(z) && identical(z$indices, ix) && is.raw(z$values) && identical(dim(z$values), c(length(ix), n))
    if (!valid) {
      dd <- igraph::distances(g, v = ix, to = seq_len(n), mode = "all", weights = NA)
      if (any(!is.finite(dd)) || any(dd < 0 | dd > 254 | dd != floor(dd)))
        stop("Distances must be finite integers 0..254 in the LCC. This network needs a different storage representation.")
      z <- list(indices = ix, values = matrix(as.raw(dd), length(ix), n)); rm(dd)
      .e109_save(z, f)
    }
    D[ix, ] <- z$values; rm(z)
    log(sprintf("Distance blocks completed: %d/%d", k, length(chunks)))
  }
  if (any(diag(D) != as.raw(0))) stop("Distance matrix diagonal is invalid.")
  # Check symmetry in deterministic samples without duplicating the full matrix.
  ix <- unique(round(seq(1, n, length.out = min(n, 50L))))
  if (!identical(D[ix, ix, drop = FALSE], t(D[ix, ix, drop = FALSE]))) stop("Distance matrix is not symmetric.")
  .e109_save(list(nodes = node_ids, D = D), whole)
  list(D = D, key = key, file = whole)
}

# Only direct children of an outer DrugBank drug record define its identity.
.e109_catalog <- function(xml_file, xml_md5, cache_dir, log) {
  key <- .e109_hash(list(parser = "direct_children_SAX_catalog_v1", md5 = xml_md5,
    XML = as.character(utils::packageVersion("XML"))))
  file <- file.path(cache_dir, paste0("drugbank_catalog_", key, ".rds"))
  if (file.exists(file)) {
    z <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!is.null(z) && isTRUE(z$complete) && identical(z$xml_md5, xml_md5)) {
      log("Using cached DrugBank identity/status catalog."); return(z)
    }
  }
  st <- new.env(parent = emptyenv())
  st$n <- 0L; st$root_seen <- FALSE; st$ended <- FALSE
  st$attrs <- character(); st$rows <- list(); st$links <- list()
  children <- function(node, tag)
    XML::getNodeSet(node, paste0("./*[local-name()='", tag, "']"))
  value <- function(node) .e109_trim(XML::xmlValue(node))
  first <- function(nodes) if (length(nodes)) value(nodes[[1L]]) else ""
  process <- function(node) {
    st$n <- st$n + 1L
    ids <- children(node, "drugbank-id")
    id_values <- vapply(ids, value, character(1))
    primary <- vapply(ids, function(x)
      identical(tolower(XML::xmlGetAttr(x, "primary", default = "")), "true"), logical(1))
    p <- unique(id_values[primary & grepl("^DB[0-9]{5,}$", id_values)])
    id_ok <- length(p) == 1L
    primary_id <- if (id_ok) p else NA_character_
    groups <- XML::getNodeSet(node,
      "./*[local-name()='groups']/*[local-name()='group']")
    groups <- tolower(unique(vapply(groups, value, character(1))))
    groups <- groups[nzchar(groups)]
    st$rows[[st$n]] <- data.frame(
      XML_record = st$n, DrugBank_ID = primary_id,
      DrugBank_Name = first(children(node, "name")),
      DrugBank_Groups = paste(sort(groups), collapse = " | "),
      DrugBank_Type = XML::xmlGetAttr(node, "type", default = ""),
      Primary_ID_status = if (id_ok) "ok" else "missing_or_conflicting_primary_ID",
      Direct_DrugBank_IDs = .e109_join(id_values),
      stringsAsFactors = FALSE)
    if (id_ok) {
      take <- grepl("^DB[0-9]{5,}$", id_values)
      st$links[[st$n]] <- unique(data.frame(
        Source_DrugBank_ID = id_values[take],
        DrugBank_ID = primary_id,
        Is_primary = id_values[take] == primary_id,
        stringsAsFactors = FALSE))
    } else st$links[[st$n]] <- data.frame()
    if (st$n %% 1000L == 0L) log("DrugBank records parsed: ", st$n)
    invisible(NULL)
  }
  handlers <- list(
    .startElement = function(name, attrs, ...) {
      if (!st$root_seen) {
        st$root_seen <- TRUE
        if (sub(".*:", "", name) != "drugbank") stop("XML root is not drugbank.")
        st$attrs <- attrs
      }
      invisible(NULL)
    },
    .endDocument = function(...) { st$ended <- TRUE; invisible(NULL) })
  log("Streaming DrugBank XML, one outer drug record at a time ...")
  XML::xmlEventParse(xml_file, handlers = handlers, branches = list(drug = process),
    saxVersion = 2L, useDotNames = TRUE, useTagName = FALSE, addContext = FALSE,
    trim = FALSE, ignoreBlanks = FALSE, replaceEntities = FALSE,
    validate = FALSE, addFinalizer = TRUE)
  if (!st$n || !st$ended) stop("DrugBank XML parsing did not complete.")
  catalog <- .e109_df(st$rows)
  links <- .e109_df(Filter(function(x) nrow(x) > 0L, st$links),
    data.frame(Source_DrugBank_ID = character(), DrugBank_ID = character(), Is_primary = logical()))
  valid <- catalog$Primary_ID_status == "ok"
  if (anyDuplicated(catalog$DrugBank_ID[valid]))
    stop("Duplicate primary DrugBank IDs in the XML; identity cannot be resolved safely.")
  z <- list(complete = TRUE, xml_md5 = xml_md5, catalog = catalog, aliases = links,
    root_attributes = st$attrs, n_records = st$n)
  .e109_save(z, file); z
}

.e109_pool <- function(source_ids, db) {
  source_ids <- sort(unique(source_ids))
  aliases <- split(db$aliases$DrugBank_ID, db$aliases$Source_DrugBank_ID)
  mapping <- .e109_df(lapply(source_ids, function(id) {
    p <- unique(aliases[[id]])
    status <- if (!length(p)) "ID_not_found_in_XML" else if (length(p) > 1L)
      "ambiguous_ID_in_XML" else if (p == id) "exact_primary_ID" else "documented_secondary_ID"
    resolved <- length(p) == 1L
    data.frame(Source_DrugBank_ID = id,
      DrugBank_ID = if (resolved) p else NA_character_,
      Identity_status = status,
      Pool_key = if (resolved) p else paste0("UNRESOLVED:", id),
      stringsAsFactors = FALSE)
  }))
  groups <- split(seq_len(nrow(mapping)), mapping$Pool_key)
  pool <- .e109_df(lapply(names(groups), function(key) {
    a <- mapping[groups[[key]], , drop = FALSE]
    ids <- unique(a$DrugBank_ID[!is.na(a$DrugBank_ID)])
    resolved <- length(ids) == 1L
    row <- if (resolved) match(ids, db$catalog$DrugBank_ID) else NA_integer_
    name <- if (is.na(row)) "" else db$catalog$DrugBank_Name[row]
    labels <- if (is.na(row)) "" else db$catalog$DrugBank_Groups[row]
    type <- if (is.na(row)) "" else db$catalog$DrugBank_Type[row]
    flags <- tolower(.e109_split(labels))
    approved <- if (length(flags)) "approved" %in% flags else NA
    withdrawn <- if (length(flags)) "withdrawn" %in% flags else NA
    eligibility <- if (!resolved) a$Identity_status[1L] else if (!length(flags))
      "DrugBank_groups_missing" else if (withdrawn) "withdrawn_label" else if (!approved)
      "not_approved_label" else "eligible_approved_nonwithdrawn"
    data.frame(Pool_key = key,
      DrugBank_ID = if (resolved) ids else a$Source_DrugBank_ID[1L],
      DrugBank_Name = name,
      NP_Display_name = if (nzchar(name)) name else a$Source_DrugBank_ID[1L],
      DrugBank_Groups = labels, DrugBank_Type = type,
      DrugBank_Approved = approved, DrugBank_Withdrawn = withdrawn,
      Identity_status = .e109_join(a$Identity_status),
      Source_DrugBank_IDs = .e109_join(a$Source_DrugBank_ID),
      N_source_DrugBank_IDs = nrow(a), Eligibility_status = eligibility,
      stringsAsFactors = FALSE)
  }))
  list(mapping = mapping, pool = pool)
}

.e109_distance <- function(D, S, T) {
  if (!length(S) || !length(T)) stop("Empty disease or target set.")
  m <- matrix(as.integer(D[S, T, drop = FALSE]), nrow = length(S), ncol = length(T))
  (sum(matrixStats::rowMins(m, useNames = FALSE)) +
    sum(matrixStats::colMins(m, useNames = FALSE))) / (length(S) + length(T))
}
.e109_draw <- function(counts, pools) {
  unlist(lapply(names(counts), function(b) {
    p <- pools[[b]]; k <- as.integer(counts[[b]])
    if (length(p) < k) stop("Insufficient nodes in a degree bin.")
    p[sample.int(length(p), k, replace = FALSE)]
  }), use.names = FALSE)
}
.e109_null_chunk <- function(task) {
  s <- get(".e109_worker_state", envir = .GlobalEnv)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  K <- length(s$profile_counts)
  values <- matrix(NA_real_, length(task$indices), K)
  disease_degree <- numeric(length(task$indices))
  target_degree_sum <- target_degree_sumsq <- numeric(K)
  for (j in seq_along(task$indices)) {
    iteration <- task$indices[j]
    set.seed(as.integer((as.double(s$seed) + as.double(iteration) * 104729) %% 2147483646 + 1))
    S <- .e109_draw(s$disease_counts, s$pools)
    disease_degree[j] <- mean(s$degrees[S])
    for (p in seq_len(K)) {
      T <- .e109_draw(s$profile_counts[[p]], s$pools)
      values[j, p] <- .e109_distance(s$D, S, T)
      td <- mean(s$degrees[T])
      target_degree_sum[p] <- target_degree_sum[p] + td
      target_degree_sumsq[p] <- target_degree_sumsq[p] + td^2
    }
  }
  list(indices = task$indices, values = values, disease_degree = disease_degree,
    target_degree_sum = target_degree_sum, target_degree_sumsq = target_degree_sumsq,
    profile_ids = s$profile_ids)
}

.e109_permutations <- function(D, S, target_sets, degrees, n_perm, seed,
    n_workers, min_bin_size, chunk_size, graph_key, cache_dir, out, log) {
  n <- nrow(D)
  bins <- .e109_degree_bins(degrees, min_bin_size)
  bin_report <- .e109_df(lapply(sort(unique(bins)), function(b) data.frame(
    Bin = b, Degree_min = min(degrees[bins == b]), Degree_max = max(degrees[bins == b]),
    Network_nodes = sum(bins == b), Signature_nodes = sum(bins[S] == b))))
  .e109_csv(bin_report, file.path(out, "08_degree_bins.csv"))
  target_counts <- lapply(target_sets, function(x) table(bins[x]))
  profile_strings <- vapply(target_counts, function(x)
    paste(names(x), as.integer(x), sep = ":", collapse = "|"), character(1))
  unique_profiles <- sort(unique(profile_strings))
  profile_index <- match(profile_strings, unique_profiles)
  profile_ids <- sprintf("DP%04d", seq_along(unique_profiles))
  profile_counts <- target_counts[match(unique_profiles, profile_strings)]
  observed <- vapply(target_sets, function(T) .e109_distance(D, S, T), numeric(1))
  key <- .e109_hash(list(engine = "full_network_two_set_bidirectional_v1",
    graph = graph_key, S = S, degree_bins = bins, profiles = unique_profiles,
    seed = seed, n_perm = n_perm, chunk_size = chunk_size,
    R = as.character(getRversion()),
    matrixStats = as.character(utils::packageVersion("matrixStats"))))
  folder <- file.path(cache_dir, key)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  tasks <- lapply(split(seq_len(n_perm), ceiling(seq_len(n_perm)/chunk_size)),
    function(i) list(indices = as.integer(i)))
  null <- matrix(NA_real_, n_perm, length(profile_ids),
    dimnames = list(NULL, profile_ids))
  disease_degree <- rep(NA_real_, n_perm)
  target_degree_sum <- target_degree_sumsq <- numeric(length(profile_ids))
  pending <- list()
  valid_chunk <- function(z, task) {
    !is.null(z) && identical(z$indices, task$indices) &&
      identical(z$profile_ids, profile_ids) &&
      identical(dim(z$values), c(length(task$indices), length(profile_ids))) &&
      length(z$disease_degree) == length(task$indices) &&
      length(z$target_degree_sum) == length(profile_ids) &&
      length(z$target_degree_sumsq) == length(profile_ids) &&
      all(is.finite(z$values)) && all(is.finite(z$disease_degree)) &&
      all(is.finite(z$target_degree_sum)) && all(is.finite(z$target_degree_sumsq))
  }
  for (k in seq_along(tasks)) {
    task <- tasks[[k]]
    task$file <- file.path(folder, sprintf("chunk_%05d.rds", k))
    z <- if (file.exists(task$file)) tryCatch(readRDS(task$file), error = function(e) NULL) else NULL
    if (valid_chunk(z, task)) {
      null[task$indices, ] <- z$values
      disease_degree[task$indices] <- z$disease_degree
      target_degree_sum <- target_degree_sum + z$target_degree_sum
      target_degree_sumsq <- target_degree_sumsq + z$target_degree_sumsq
    } else pending[[length(pending) + 1L]] <- task
  }
  state <- list(D = D, disease_counts = table(bins[S]),
    pools = split(seq_len(n), bins), profile_counts = profile_counts,
    profile_ids = profile_ids, degrees = as.numeric(degrees), seed = seed)
  had_state <- exists(".e109_worker_state", envir = .GlobalEnv, inherits = FALSE)
  if (had_state) old_state <- get(".e109_worker_state", envir = .GlobalEnv)
  on.exit({
    if (had_state) assign(".e109_worker_state", old_state, envir = .GlobalEnv)
    else if (exists(".e109_worker_state", envir = .GlobalEnv, inherits = FALSE))
      rm(".e109_worker_state", envir = .GlobalEnv)
  }, add = TRUE)
  assign(".e109_worker_state", state, envir = .GlobalEnv)
  cl <- NULL
  if (length(pending) && n_workers > 1L) {
    log("Starting ", n_workers, " workers; distance matrix per worker: ",
      sprintf("%.2f GiB.", length(D)/1024^3))
    override <- Sys.getenv("ENDO109_RSCRIPT", unset = "")
    cl <- if (nzchar(override)) parallel::makePSOCKcluster(n_workers, rscript = override) else
      parallel::makePSOCKcluster(n_workers)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
    versions <- parallel::clusterCall(cl, function(paths) {
      .libPaths(paths)
      if (!requireNamespace("matrixStats", quietly = TRUE)) stop("Worker cannot load matrixStats.")
      as.character(utils::packageVersion("matrixStats"))
    }, .libPaths())
    if (!all(unlist(versions) == as.character(utils::packageVersion("matrixStats"))))
      stop("Worker matrixStats versions differ.")
    ex <- new.env(parent = emptyenv())
    for (nm in c(".e109_draw", ".e109_distance", ".e109_null_chunk"))
      assign(nm, get(nm, mode = "function"), envir = ex)
    assign(".e109_worker_state", state, envir = ex)
    parallel::clusterExport(cl, ls(ex, all.names = TRUE), envir = ex)
  }
  done <- sum(is.finite(null[, 1L]))
  log("Permutation iterations cached: ", done, "/", n_perm,
    "; distinct target degree profiles: ", length(profile_ids))
  while (length(pending)) {
    k <- min(length(pending), n_workers)
    todo <- pending[seq_len(k)]
    started <- proc.time()[["elapsed"]]
    ans <- if (is.null(cl)) lapply(todo, .e109_null_chunk) else
      parallel::parLapply(cl, todo, .e109_null_chunk)
    for (j in seq_along(ans)) {
      z <- ans[[j]]; task <- todo[[j]]
      if (!valid_chunk(z, task)) stop("Invalid permutation chunk.")
      null[z$indices, ] <- z$values
      disease_degree[z$indices] <- z$disease_degree
      target_degree_sum <- target_degree_sum + z$target_degree_sum
      target_degree_sumsq <- target_degree_sumsq + z$target_degree_sumsq
      .e109_save(z, task$file)
    }
    added <- sum(vapply(todo, function(x) length(x$indices), integer(1)))
    done <- done + added
    elapsed <- proc.time()[["elapsed"]] - started
    log("Permutations completed: ", done, "/", n_perm,
      "; last batch: ", round(elapsed, 1), " s",
      "; approximate remaining: ", round((n_perm - done) * elapsed/max(added, 1)/60, 1), " min.")
    pending <- if (k == length(pending)) list() else pending[-seq_len(k)]
  }
  if (any(!is.finite(null))) stop("Permutation matrix is incomplete.")
  mu <- colMeans(null)
  sd <- vapply(seq_len(ncol(null)), function(j) stats::sd(null[, j]), numeric(1))
  p <- vapply(seq_along(target_sets), function(j)
    (1 + sum(null[, profile_index[j]] <= observed[j] + 1e-12))/(n_perm + 1), numeric(1))
  z <- rep(NA_real_, length(target_sets))
  good_sd <- sd[profile_index] > 0
  z[good_sd] <- (observed[good_sd] - mu[profile_index[good_sd]])/sd[profile_index[good_sd]]
  target_mu <- target_degree_sum/n_perm
  target_sd <- sqrt(pmax(0,
    (target_degree_sumsq - target_degree_sum^2/n_perm)/(n_perm - 1)))
  stats <- data.frame(Pool_key = names(target_sets),
    NP_Distance = unname(observed), NP_Null_mean = mu[profile_index],
    NP_Null_SD = sd[profile_index], NP_Z = z, NP_Empirical_P = p,
    NP_Null_SD_zero = !good_sd,
    NP_Degree_profile_ID = profile_ids[profile_index],
    NP_Observed_target_mean_degree = vapply(target_sets, function(T) mean(degrees[T]), numeric(1)),
    NP_Null_target_mean_degree = target_mu[profile_index],
    NP_Null_target_mean_degree_SD = target_sd[profile_index],
    NP_Observed_signature_mean_degree = mean(degrees[S]),
    NP_Null_signature_mean_degree = mean(disease_degree), row.names = NULL)
  list(stats = stats, null = null, disease_degree = disease_degree,
    profile_definitions = data.frame(Profile_ID = profile_ids, Degree_bin_counts = unique_profiles,
      Null_mean = mu, Null_SD = sd, Null_target_mean_degree = target_mu,
      Null_target_mean_degree_SD = target_sd),
    profile_assignment = data.frame(Pool_key = names(target_sets),
      Profile_ID = profile_ids[profile_index]),
    key = key, folder = folder)
}

.e109_figures <- function(result, out, dpi) {
  tested <- result[result$NP_Status == "ok" & is.finite(result$NP_Z), , drop = FALSE]
  tested <- tested[order(tested$NP_Z, tested$NP_Empirical_P, tested$DrugBank_ID), , drop = FALSE]
  draw <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    if (!nrow(tested)) {
      graphics::plot.new(); graphics::text(.5, .5, "No finite network Z scores")
      graphics::plot.new(); graphics::text(.5, .5, "See the complete result/status tables")
      return(invisible(NULL))
    }
    top <- utils::head(tested, 30L)
    colors <- ifelse(top$NP_BH_supported, "#B2182B",
      ifelse(top$NP_Nominal_candidate, "#EF8A62", "#4393C3"))
    graphics::par(mar = c(4, 11, 3, 1))
    graphics::barplot(rev(top$NP_Z), names.arg = rev(top$NP_Display_name),
      horiz = TRUE, las = 1, col = rev(colors), border = NA, cex.names = .65,
      xlab = "Bidirectional network proximity Z", main = "Top drugs by network Z")
    graphics::abline(v = 0, col = "grey40")
    graphics::abline(v = -1.5, col = "grey40", lty = 2)
    graphics::par(mar = c(4, 4, 3, 1))
    colors <- ifelse(tested$NP_BH_supported, "#B2182B",
      ifelse(tested$NP_Nominal_candidate, "#EF8A62", "#4393C3"))
    graphics::plot(tested$NP_Z, -log10(tested$NP_Empirical_P),
      pch = 16, col = grDevices::adjustcolor(colors, alpha.f = .6),
      xlab = "Bidirectional network proximity Z",
      ylab = expression(-log[10]("empirical P")), main = "All completed drug tests")
    graphics::abline(v = -1.5, h = -log10(.05), lty = 2, col = "grey55")
    graphics::legend("topright",
      legend = c("BH FDR < 0.05, Z < 0", "Nominal P < 0.05, Z < 0", "Other tested"),
      pch = 16, col = c("#B2182B", "#EF8A62", "#4393C3"), bty = "n", cex = .65)
  }
  rows <- list()
  for (kind in c("png", "pdf")) {
    path <- file.path(out, paste0("Network_prioritization.", kind))
    error <- ""; opened <- FALSE
    ok <- tryCatch({
      if (kind == "png") grDevices::png(path, width = 14, height = 8, units = "in", res = dpi)
      else grDevices::pdf(path, width = 14, height = 8)
      opened <- TRUE; draw(); TRUE
    }, error = function(e) { error <<- conditionMessage(e); FALSE },
      finally = if (opened) grDevices::dev.off())
    rows[[length(rows) + 1L]] <- data.frame(File = basename(path), Success = ok, Error = error)
  }
  .e109_df(rows)
}

run_endoAD_109_full_library <- function(
    base_dir = "D:/My Desktop/drug/5.endoAD_109_screen",
    n_perm = 10000L, n_workers = 2L, seed = 20261105L,
    min_bin_size = 100L, expected_signature_n = 109L,
    distance_block_size = 200L, permutation_chunk_size = 100L,
    max_distance_gb = 2, install_missing = TRUE, make_figures = TRUE,
    png_dpi = 300L) {
  if (!dir.exists(base_dir)) stop("Input directory does not exist: ", base_dir)
  base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
  root <- file.path(base_dir, "endoAD_109_results")
  cache <- file.path(base_dir, "endoAD_109_cache")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = root)
  dir.create(out)
  finished <- FALSE
  on.exit(if (!finished) writeLines(
    "INCOMPLETE RUN. Do not interpret this directory as a completed screen. Fix the error and Source the script again; matching cache chunks will be reused.",
    file.path(out, "INCOMPLETE_RUN.txt")), add = TRUE)
  log <- function(...) {
    msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(..., collapse = ""))
    message(msg); cat(msg, "\n", file = file.path(out, "run.log"), append = TRUE)
  }
  log("Results: ", out)
  checks <- list()
  check <- function(name, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(Check = name, Passed = isTRUE(passed), Detail = detail)
    .e109_csv(.e109_df(checks), file.path(out, "analysis_checks.csv"))
    if (!isTRUE(passed)) stop(name, ": ", detail, call. = FALSE)
    invisible(TRUE)
  }
  packages <- c("data.table", "igraph", "XML", "matrixStats")
  lib <- file.path(cache, paste0("R_library_", getRversion()))
  dir.create(lib, showWarnings = FALSE); .libPaths(c(lib, .libPaths()))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && install_missing)
    utils::install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
  check("Required_packages", all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)),
    "Requires data.table, igraph, XML and matrixStats.")
  integer_ok <- function(x, minimum) length(x) == 1L && is.finite(x) &&
    x >= minimum && x == as.integer(x)
  check("Valid_parameters", integer_ok(n_perm, 100L) && integer_ok(n_workers, 1L) &&
      integer_ok(seed, 1L) && seed < 2147483646 && integer_ok(min_bin_size, 2L) &&
      integer_ok(distance_block_size, 1L) && integer_ok(permutation_chunk_size, 1L) &&
      length(max_distance_gb) == 1L && is.finite(max_distance_gb) && max_distance_gb > 0,
    "Permutation count >=100; positive valid worker/seed/bin/block settings.")
  n_perm <- as.integer(n_perm); n_workers <- as.integer(n_workers)
  distance_block_size <- as.integer(distance_block_size)
  permutation_chunk_size <- as.integer(permutation_chunk_size)
  cores <- parallel::detectCores(logical = FALSE)
  if (is.finite(cores)) n_workers <- min(n_workers, max(1L, cores - 1L))
  find_file <- function(pattern, label) {
    paths <- list.files(base_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    check(paste0("Input_", label), length(paths) == 1L,
      paste0("Expected exactly one ", label, " file; found: ",
        if (length(paths)) paste(basename(paths), collapse = ", ") else "none"))
    paths
  }
  files <- c(
    Signature = find_file("^signature_candidates\\.csv$", "signature"),
    HGNC = find_file("^hgnc_complete_set\\.txt$", "HGNC"),
    Interactome = find_file("^HumanInteractome\\.(txt|tsv)$", "interactome"),
    Drug_targets = find_file("^DrugTargetNetwork\\.txt$", "drug_targets"),
    DrugBank_XML = find_file("\\.xml$", "DrugBank_XML"))
  before <- file.info(files)
  log("[1/8] Recording file checksums, including the DrugBank XML ...")
  manifest <- data.frame(Input = names(files), File = unname(files), Bytes = before$size,
    Modified = as.character(before$mtime), MD5 = unname(tools::md5sum(files)))
  check("File_checksums", !anyNA(manifest$MD5), "All five input files have readable checksums.")
  .e109_csv(manifest, file.path(out, "input_manifest.csv"))
  settings <- list(script_version = .e109_version, base_dir = base_dir, n_perm = n_perm,
    n_workers = n_workers, seed = seed, min_bin_size = min_bin_size,
    expected_signature_n = expected_signature_n, distance_block_size = distance_block_size,
    permutation_chunk_size = permutation_chunk_size,
    disease_null_pool = "all nodes of the retained protein-coding interactome LCC",
    target_null_pool = "all nodes of the same interactome LCC",
    primary_distance = "(sum disease-to-nearest-target + sum target-to-nearest-disease)/(n_disease+n_targets)",
    eligibility = "endoAD drug-target source; unambiguous DrugBank ID; approved group present; withdrawn group absent",
    primary_BH_family = "all completed eligible drug tests with valid empirical P",
    source_drug_alias_rule = "only documented direct DrugBank secondary IDs; union targets within the same primary record",
    profile_null_sharing = "same target size and global-degree-bin counts share one null distribution",
    null_sampling = "without replacement within each set; disease and target draws independent and may overlap",
    input_manifest = manifest)
  .e109_save(settings, file.path(out, "analysis_settings.rds"), compress = TRUE)
  log("[2/8] Mapping the original expression signature and interactome ...")
  ref <- .e109_hgnc(files[["HGNC"]])
  signature <- .e109_read(files[["Signature"]])
  .e109_csv(signature, file.path(out, "00_signature_input_snapshot.csv"))
  if (!is.null(expected_signature_n))
    check("Input_signature_rows", nrow(signature) == expected_signature_n,
      paste("Expected", expected_signature_n, "original gene records; found", nrow(signature)))
  signature_map <- .e109_gene_input(signature, ref)
  .e109_csv(signature_map, file.path(out, "01_signature_gene_mapping.csv"))
  check("Signature_identity_mapping",
    nrow(signature_map) == nrow(signature) &&
      all(signature_map$Mapping_status == "Mapped" & signature_map$Protein_coding) &&
      !anyDuplicated(signature_map$HGNC_ID),
    "Each input signature record must resolve to a distinct approved HGNC protein-coding gene. See 01_signature_gene_mapping.csv.")
  network <- .e109_pairs(files[["Interactome"]], "network")
  endpoints <- .e109_map(c(network$ID1, network$ID2), ref)
  .e109_csv(endpoints, file.path(out, "02_network_endpoint_mapping.csv"))
  ia <- match(network$ID1, endpoints$Input_ID); ib <- match(network$ID2, endpoints$Input_ID)
  valid <- endpoints$Mapping_status[ia] == "Mapped" & endpoints$Protein_coding[ia] &
    endpoints$Mapping_status[ib] == "Mapped" & endpoints$Protein_coding[ib]
  valid[is.na(valid)] <- FALSE
  .e109_csv(network[!valid, , drop = FALSE], file.path(out, "02a_excluded_network_edges.csv"))
  a <- endpoints$HGNC_ID[ia[valid]]; b <- endpoints$HGNC_ID[ib[valid]]
  loops <- sum(a == b); keep <- a != b
  edges <- unique(data.frame(from = pmin(a[keep], b[keep]), to = pmax(a[keep], b[keep])))
  check("Nonempty_cleaned_network", nrow(edges) > 0L, "Valid protein-coding interaction edges remain.")
  vertices <- data.frame(name = sort(unique(c(edges$from, edges$to))))
  g0 <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  comp <- igraph::components(g0)
  g <- igraph::induced_subgraph(g0, which(comp$membership == which.max(comp$csize)))
  nodes <- igraph::V(g)$name
  symbols <- ref$table$Symbol[match(nodes, ref$table$HGNC_ID)]
  node_table <- data.frame(Network_index = seq_along(nodes), HGNC_ID = nodes,
    gene_symbol = symbols, Global_degree = as.numeric(igraph::degree(g)))
  .e109_csv(node_table, file.path(out, "02b_network_nodes.csv"))
  signature_map$In_network_LCC <- signature_map$HGNC_ID %in% nodes
  signature_map$Network_index <- match(signature_map$HGNC_ID, nodes)
  .e109_csv(signature_map, file.path(out, "01_signature_gene_mapping.csv"))
  S <- sort(unique(signature_map$Network_index[signature_map$In_network_LCC]))
  check("Signature_network_coverage", length(S) >= 2L,
    paste(length(S), "of", nrow(signature), "signature genes are in the network LCC. No seed-induced LCC restriction is used."))
  .e109_csv(data.frame(HGNC_ID = nodes[S], gene_symbol = symbols[S], Network_index = S),
    file.path(out, "01a_fixed_network_signature.csv"))
  .e109_csv(data.frame(Item = c("Raw_network_rows", "Edges_with_unmapped_or_noncoding_endpoints",
      "Removed_self_loops", "Unique_cleaned_edges", "Network_LCC_nodes", "Network_LCC_edges",
      "Input_signature_genes", "Signature_genes_in_network_LCC"),
    Value = c(nrow(network), sum(!valid), loops, nrow(edges), igraph::vcount(g),
      igraph::ecount(g), nrow(signature), length(S))), file.path(out, "02c_network_summary.csv"))
  log("Fixed network signature: ", length(S), "/", nrow(signature),
    " genes; network: ", length(nodes), " nodes, ", igraph::ecount(g), " edges.")

  log("[3/8] Establishing the complete endoAD source drug pool ...")
  drug_links <- .e109_pairs(files[["Drug_targets"]], "drug")
  db <- .e109_catalog(files[["DrugBank_XML"]],
    manifest$MD5[manifest$Input == "DrugBank_XML"], cache, log)
  .e109_csv(db$catalog, file.path(out, "03_DrugBank_catalog.csv"))
  .e109_csv(db$catalog[db$catalog$Primary_ID_status != "ok", , drop = FALSE],
    file.path(out, "03a_DrugBank_primary_ID_review.csv"))
  .e109_csv(data.frame(Attribute = names(db$root_attributes),
    Value = as.character(db$root_attributes)), file.path(out, "03b_XML_root_attributes.csv"))
  resolved <- .e109_pool(drug_links$ID1, db)
  .e109_csv(resolved$mapping, file.path(out, "04_source_DrugBank_ID_resolution.csv"))
  result <- resolved$pool
  check("Complete_source_pool", nrow(resolved$mapping) == length(unique(drug_links$ID1)),
    "Every source DrugBank ID is represented in the identity-resolution table.")

  log("[4/8] Mapping all source drug-associated proteins ...")
  tm <- .e109_map(drug_links$ID2, ref)
  ii <- match(drug_links$ID2, tm$Input_ID)
  pm <- resolved$mapping$Pool_key[
    match(drug_links$ID1, resolved$mapping$Source_DrugBank_ID)]
  target_map <- data.frame(Source_row = drug_links$Source_row,
    Source_DrugBank_ID = drug_links$ID1, Pool_key = pm,
    Source_protein_ID = drug_links$ID2,
    HGNC_ID = tm$HGNC_ID[ii], gene_symbol = tm$gene_symbol[ii],
    Mapping_status = tm$Mapping_status[ii], Mapping_basis = tm$Mapping_basis[ii],
    Protein_coding = tm$Protein_coding[ii], stringsAsFactors = FALSE)
  target_map$Network_index <- match(target_map$HGNC_ID, nodes)
  target_map$In_network_LCC <- target_map$Mapping_status == "Mapped" &
    target_map$Protein_coding & !is.na(target_map$Network_index)
  target_map$In_network_LCC[is.na(target_map$In_network_LCC)] <- FALSE
  .e109_csv(target_map, file.path(out, "05_all_source_target_mapping.csv"))
  by_drug <- split(seq_len(nrow(target_map)), target_map$Pool_key)
  sets <- vector("list", nrow(result)); names(sets) <- result$Pool_key
  result$NP_N_source_protein_IDs <- result$NP_N_HGNC_PC_targets <- result$NP_N_network_targets <- 0L
  result$NP_N_signature_overlap <- 0L
  result$NP_Network_target_genes <- result$NP_Network_target_HGNC_IDs <- ""
  result$NP_Signature_overlap_genes <- ""
  for (j in seq_len(nrow(result))) {
    t <- target_map[by_drug[[result$Pool_key[j]]], , drop = FALSE]
    sets[[j]] <- sort(unique(t$Network_index[t$In_network_LCC]))
    pc <- t$Mapping_status == "Mapped" & t$Protein_coding
    result$NP_N_source_protein_IDs[j] <- length(unique(t$Source_protein_ID[.e109_present(t$Source_protein_ID)]))
    result$NP_N_HGNC_PC_targets[j] <- length(unique(t$HGNC_ID[pc & !is.na(t$HGNC_ID)]))
    result$NP_N_network_targets[j] <- length(sets[[j]])
    result$NP_Network_target_genes[j] <- .e109_join(symbols[sets[[j]]])
    result$NP_Network_target_HGNC_IDs[j] <- .e109_join(nodes[sets[[j]]])
    overlap <- intersect(sets[[j]], S)
    result$NP_N_signature_overlap[j] <- length(overlap)
    result$NP_Signature_overlap_genes[j] <- .e109_join(symbols[overlap])
  }
  eligible <- result$Eligibility_status == "eligible_approved_nonwithdrawn"
  result$NP_Status <- result$Eligibility_status
  result$NP_Status[eligible & result$NP_N_HGNC_PC_targets == 0L] <- "no_HGNC_protein_coding_targets"
  result$NP_Status[eligible & result$NP_N_HGNC_PC_targets > 0L &
    result$NP_N_network_targets == 0L] <- "no_targets_in_network_LCC"
  result$NP_Status[eligible & result$NP_N_network_targets > 0L] <- "eligible_pending"
  target_keys <- vapply(sets, function(T) if (length(T)) paste(nodes[T], collapse = "|") else "", character(1))
  key_levels <- sort(unique(target_keys[nzchar(target_keys)]))
  result$NP_Target_set_ID <- ifelse(nzchar(target_keys),
    sprintf("TS%05d", match(target_keys, key_levels)), NA_character_)
  .e109_csv(result, file.path(out, "06_all_source_drug_coverage.csv"))
  .e109_csv(result[eligible, , drop = FALSE], file.path(out, "07_approved_nonwithdrawn_pool.csv"))
  test_ix <- which(result$NP_Status == "eligible_pending")
  target_sets <- sets[test_ix]
  log("Source IDs: ", nrow(resolved$mapping), "; consolidated pool entries: ", nrow(result),
    "; eligible approved/nonwithdrawn: ", sum(eligible), "; network-testable: ", length(test_ix), ".")

  result$NP_Distance <- result$NP_Null_mean <- result$NP_Null_SD <- result$NP_Z <-
    result$NP_Empirical_P <- result$NP_Network_FDR <- NA_real_
  result$NP_Null_SD_zero <- NA
  result$NP_Degree_profile_ID <- NA_character_
  result$NP_Observed_target_mean_degree <- result$NP_Null_target_mean_degree <-
    result$NP_Null_target_mean_degree_SD <- result$NP_Observed_signature_mean_degree <-
    result$NP_Null_signature_mean_degree <- NA_real_
  result$NP_Permutations <- NA_integer_
  inference <- NULL; distance_details <- list()
  if (length(test_ix)) {
    log("[5/8] Computing/restoring all-pairs unweighted network distances ...")
    dm <- .e109_distances(g, file.path(cache, "distances"), distance_block_size,
      max_distance_gb, log)
    settings$distance_cache_key <- dm$key; settings$distance_cache_file <- dm$file
    log("[6/8] Two-set degree-bin-matched proximity permutations ...")
    inference <- .e109_permutations(dm$D, S, target_sets, igraph::degree(g),
      n_perm, seed, n_workers, min_bin_size, permutation_chunk_size,
      dm$key, file.path(cache, "permutations"), out, log)
    ix <- match(inference$stats$Pool_key, result$Pool_key)
    for (column in setdiff(names(inference$stats), "Pool_key"))
      result[[column]][ix] <- inference$stats[[column]]
    result$NP_Permutations[ix] <- n_perm
    result$NP_Status[ix] <- "ok"
    .e109_csv(inference$profile_definitions, file.path(out, "09_null_degree_profiles.csv"))
    .e109_csv(inference$profile_assignment, file.path(out, "09a_drug_degree_profile_assignment.csv"))
    settings$permutation_cache_key <- inference$key
    settings$permutation_cache_folder <- inference$folder
    .e109_save(list(null_distances = inference$null,
      profiles = inference$profile_definitions, assignment = inference$profile_assignment,
      disease_mean_degree_by_iteration = inference$disease_degree,
      key = inference$key), file.path(out, "19_null_profile_distributions.rds"), compress = TRUE)
    for (j in seq_along(target_sets)) {
      T <- target_sets[[j]]
      m <- matrix(as.integer(dm$D[S, T, drop = FALSE]), nrow = length(S), ncol = length(T))
      distance_details[[j]] <- data.frame(Pool_key = names(target_sets)[j],
        DrugBank_ID = result$DrugBank_ID[test_ix[j]], Target_HGNC_ID = nodes[T],
        Target_gene = symbols[T],
        Minimum_distance_to_signature = matrixStats::colMins(m, useNames = FALSE),
        In_signature = T %in% S)
    }
  } else log("[5-6/8] No eligible network-testable drugs; complete status tables will be saved.")
  good <- which(result$NP_Status == "ok" & is.finite(result$NP_Empirical_P))
  check("Completed_test_family", length(good) == length(test_ix),
    "Every eligible network-testable drug has a completed empirical P value.")
  check("Valid_empirical_P", !length(good) || all(result$NP_Empirical_P[good] >=
      1/(n_perm + 1) & result$NP_Empirical_P[good] <= 1),
    "Inclusive lower-tail empirical P values with the finite-permutation correction.")
  if (length(good))
    result$NP_Network_FDR[good] <- stats::p.adjust(result$NP_Empirical_P[good], method = "BH")
  result$NP_BH_family_size <- length(good)
  result$NP_Rank_Z <- result$NP_Rank_P <- NA_real_
  result$NP_N_tested_drugs_same_targets <- NA_integer_
  if (length(good)) {
    result$NP_Rank_Z[good] <- rank(result$NP_Z[good], ties.method = "min", na.last = "keep")
    result$NP_Rank_P[good] <- rank(result$NP_Empirical_P[good], ties.method = "min")
    freq <- table(result$NP_Target_set_ID[good])
    result$NP_N_tested_drugs_same_targets[good] <- as.integer(freq[result$NP_Target_set_ID[good]])
  }
  result$NP_Exploratory_Z_candidate <- result$NP_Status == "ok" &
    is.finite(result$NP_Z) & result$NP_Z < -1.5
  result$NP_Nominal_candidate <- result$NP_Status == "ok" &
    is.finite(result$NP_Z) & result$NP_Z < 0 & is.finite(result$NP_Empirical_P) &
    result$NP_Empirical_P < .05
  result$NP_BH_supported <- result$NP_Status == "ok" &
    is.finite(result$NP_Z) & result$NP_Z < 0 & is.finite(result$NP_Network_FDR) &
    result$NP_Network_FDR < .05
  result$NP_Evidence <- result$NP_Status
  result$NP_Evidence[good] <- ifelse(!is.finite(result$NP_Z[good]), "undefined_Z_degenerate_null",
    ifelse(result$NP_Z[good] < 0, "negative_Z_without_nominal_support", "nonnegative_Z"))
  result$NP_Evidence[result$NP_Nominal_candidate] <- "nominal_proximity"
  result$NP_Evidence[result$NP_BH_supported] <- "BH_FDR_supported_proximity"
  tested <- result[good, , drop = FALSE]
  tested <- tested[order(is.na(tested$NP_Z), tested$NP_Z,
    tested$NP_Empirical_P, tested$DrugBank_ID), , drop = FALSE]
  tested$NP_Order_in_Z_table <- seq_len(nrow(tested))
  all_order <- order(result$NP_Status != "ok", is.na(result$NP_Z),
    result$NP_Z, result$NP_Empirical_P, result$DrugBank_ID)
  log("[7/8] Saving all rankings, status records, and candidate flags ...")
  .e109_csv(tested, file.path(out, "11_all_tested_drugs_ranked.csv"))
  .e109_csv(result[all_order, , drop = FALSE], file.path(out, "12_all_drugs_network_results.csv"))
  .e109_csv(tested[tested$NP_Exploratory_Z_candidate, , drop = FALSE],
    file.path(out, "13_exploratory_Z_candidates.csv"))
  .e109_csv(tested[tested$NP_Nominal_candidate, , drop = FALSE],
    file.path(out, "14_nominal_P_candidates.csv"))
  .e109_csv(tested[tested$NP_BH_supported, , drop = FALSE],
    file.path(out, "15_BH_FDR_supported_candidates.csv"))
  .e109_csv(result[result$NP_Status != "ok", , drop = FALSE],
    file.path(out, "16_not_tested_or_not_evaluable.csv"))
  .e109_csv(.e109_df(distance_details, data.frame(Pool_key = character(),
    DrugBank_ID = character(), Target_HGNC_ID = character(), Target_gene = character(),
    Minimum_distance_to_signature = numeric(), In_signature = logical())),
    file.path(out, "17_target_to_signature_distances.csv"))
  .e109_csv(as.data.frame(table(result$NP_Status), stringsAsFactors = FALSE),
    file.path(out, "status_counts.csv"))
  figure_exports <- if (make_figures) .e109_figures(result, out, png_dpi) else
    data.frame(File = character(), Success = logical(), Error = character())
  .e109_csv(figure_exports, file.path(out, "figure_exports.csv"))
  after <- file.info(files)
  check("Inputs_unchanged", identical(before$size, after$size) &&
      identical(as.numeric(before$mtime), as.numeric(after$mtime)),
    "Input file sizes and modification times did not change during the run.")
  summary <- data.frame(Item = c(
    "Input_signature_genes", "Signature_genes_in_network_LCC", "Network_LCC_nodes",
    "Network_LCC_edges", "DrugBank_records_in_XML", "Unique_source_DrugBank_IDs",
    "Consolidated_source_pool_entries", "Approved_nonwithdrawn_pool_drugs",
    "Approved_drugs_with_network_targets", "Network_tested_drugs",
    "Approved_not_evaluable_drugs", "Source_pool_entries_not_tested",
    "Distinct_tested_target_sets", "Distinct_tested_degree_profiles",
    "Tested_drugs_with_undefined_Z", "Negative_Z_drugs", "Exploratory_Z_candidates",
    "Nominal_P_candidates", "BH_FDR_supported_drugs", "Permutation_iterations",
    "BH_family_size", "Failed_figure_exports", "Run_complete"),
    Value = c(nrow(signature), length(S), igraph::vcount(g), igraph::ecount(g),
      db$n_records, nrow(resolved$mapping), nrow(result), sum(eligible),
      sum(eligible & result$NP_N_network_targets > 0L), length(good),
      sum(eligible & result$NP_Status != "ok"), sum(result$NP_Status != "ok"),
      length(unique(result$NP_Target_set_ID[good])),
      if (is.null(inference)) 0L else ncol(inference$null),
      sum(result$NP_Status == "ok" & !is.finite(result$NP_Z)),
      sum(result$NP_Status == "ok" & is.finite(result$NP_Z) & result$NP_Z < 0),
      sum(result$NP_Exploratory_Z_candidate), sum(result$NP_Nominal_candidate),
      sum(result$NP_BH_supported), if (length(good)) n_perm else 0L, length(good),
      sum(!figure_exports$Success), 1L))
  .e109_csv(summary, file.path(out, "analysis_summary.csv"))
  .e109_save(settings, file.path(out, "analysis_settings.rds"), compress = TRUE)
  .e109_save(list(settings = settings, summary = summary, results = result,
    ranked = tested, graph = g, network_nodes = node_table,
    signature_mapping = signature_map, fixed_signature_indices = S,
    drug_ID_mapping = resolved$mapping,
    source_target_mapping = target_map, target_sets = sets),
    file.path(out, "analysis_objects.rds"), compress = TRUE)
  .e109_csv(data.frame(Package = packages,
    Version = vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))),
    file.path(out, "package_versions.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(c(
    "METHODS AND INTERPRETATION",
    paste("Script version:", .e109_version),
    "This adapted endoAD analysis screens the fixed expression-signature gene set.",
    "No seed-induced LCC significance requirement or seed-LCC restriction was applied.",
    "The original signature was not reselected by network or drug-screen results.",
    "The same local HGNC annotation was used for signature, interactome endpoints, and source drug-associated proteins.",
    "Only unique approved HGNC protein-coding genes were used.",
    "The unweighted undirected interactome was deduplicated, self-loops removed, and its largest connected component retained.",
    "The source drug universe was exactly the supplied endoAD DrugTargetNetwork table.",
    "This is not an exhaustive screen of every drug in DrugBank.",
    "Source DrugBank IDs were matched only to direct primary or unambiguous documented secondary DrugBank IDs.",
    "Source IDs mapping to the same primary DrugBank record were consolidated and their source target links unioned.",
    "Different primary DrugBank drugs remained separate even when their mapped targets were identical.",
    "Primary eligibility required the DrugBank approved group and absence of the withdrawn group.",
    "Approved is a DrugBank group label; an FDA-specific approval jurisdiction was not inferred.",
    "Status-excluded, unresolved and target-unassessable entries remained in the complete source-pool table.",
    "Drug-associated proteins came from the supplied table; mechanism action, binding affinities and species annotations were not re-inferred from its two columns.",
    "No transcriptomic, BBB, target-overlap, Z, empirical-P, or FDR filter was used to preselect the tested approved pool.",
    "The primary statistic followed Fang et al., Nature Aging 2021, Eq. 2; DOI 10.1038/s43587-021-00138-z.",
    "d_bi = [sum_s min_t shortest_path(s,t) + sum_t min_s shortest_path(t,s)] / (n_signature + n_targets).",
    "Both sets were randomized on the fixed interactome, preserving each set's cardinality and global-degree-bin counts.",
    "The random pool for both sets was the entire retained protein-coding network LCC.",
    "Sampling was without replacement within a set; independently sampled disease and target sets could overlap.",
    paste("Degree bins accumulated complete degree classes to at least", min_bin_size,
      "nodes, with a small final bin merged into its predecessor. Matching is by bins, not by exact degree."),
    "Degree-range and observed/null mean-degree diagnostics were saved.",
    paste("Each target degree profile received", n_perm, "permutations; master seed", seed, "."),
    "Drugs with identical target degree-bin count profiles shared the same valid marginal null distribution as a computational rule.",
    "Shared target sets or null distributions are not independent pharmacological evidence.",
    "The permutation draw for a given iteration was deterministic and independent of worker scheduling.",
    "Z = (observed distance - null mean) / null SD; negative Z indicates greater network proximity.",
    "P = [1 + count(null distance <= observed distance, with tolerance 1e-12)] / (1 + permutations).",
    "P was empirical and one-sided; no normal approximation from Z was used.",
    "If null SD was zero, Z remained NA but the completed finite empirical P was retained in the BH family.",
    "BH correction included every completed eligible drug test with a valid empirical P, including drugs with positive Z.",
    "Exploratory Z flag: Z < -1.5. Nominal flag: Z < 0 and P < 0.05. Corrected flag: Z < 0 and BH FDR < 0.05.",
    "All completed tests were retained in the ranked table regardless of these flags.",
    "The distance is unsigned and does not establish transcriptomic reversal, therapeutic action direction, brain exposure, or clinical benefit.",
    "Source references: https://github.com/ChengF-Lab/endoAD ; https://doi.org/10.1038/s43587-021-00138-z"
  ), file.path(out, "analysis_methods.txt"))
  if (exists(".e109_this_file", inherits = TRUE)) {
    src <- get(".e109_this_file", inherits = TRUE)
    if (length(src) == 1L && !is.na(src) && file.exists(src))
      file.copy(src, file.path(out, "analysis_script.R"), overwrite = FALSE)
  }
  log("[8/8] Screen complete. All source-pool records and all completed tests were retained.")
  finished <- TRUE
  print(summary, row.names = FALSE)
  if (nrow(tested)) {
    cat("\nTop 20 by network Z (full ranking saved):\n")
    print(utils::head(tested[, c("DrugBank_ID", "NP_Display_name", "NP_N_network_targets",
      "NP_Distance", "NP_Z", "NP_Empirical_P", "NP_Network_FDR")], 20L), row.names = FALSE)
  }
  cat("\nFinished. Results:\n", out, "\n", sep = "")
  invisible(out)
}

.e109_this_file <- tryCatch({
  f <- sys.frame(1)$ofile
  if (is.null(f)) NA_character_ else normalizePath(f, winslash = "/", mustWork = FALSE)
}, error = function(e) NA_character_)

if (!identical(Sys.getenv("ENDO109_SKIP_AUTORUN"), "1")) {
  ENDOAD_109_DIR <- run_endoAD_109_full_library(
    base_dir = BASE_DIR, n_perm = N_PERM, n_workers = N_WORKERS,
    seed = RANDOM_SEED, min_bin_size = MIN_BIN_SIZE)
}
