# AD37_ATC_seed_network_B3DB_one_click.R
# Version: 1.0.0 (2026-11-05)
# Purpose: ATC classification + 3 Fang 2021 seed gene sets + 10,000
# degree-bin-matched bidirectional network-proximity permutations + B3DB matching.
# This is an ADAPTED seed-set screen, NOT an exact replication of GPSnet / Fang.
# It evaluates 37 previously preselected drugs: conditional and exploratory.
# Sources:
# Fang et al., Nature Aging (2021), DOI: 10.1038/s43587-021-00138-z
# Meng et al., Scientific Data (2021), DOI: 10.1038/s41597-021-01069-5
# B3DB: https://github.com/theochem/B3DB
# ATC: https://atcddd.fhi.no/
# EndoAD: https://github.com/ChengF-Lab/endoAD
#
# REQUIRED in BASE_DIR:
#   37ATC.xlsx
#   Table S13.xlsx
#   B3DB_classification.tsv
#   14_nominal_P_candidates.csv  [or 14_nominal_P_candidates(1).csv]
# For the exact same PPI as the previous run, the script automatically looks for
# analysis_objects.rds in BASE_DIR or the latest completed run under PREVIOUS_DIR.
# Otherwise it reconstructs PPI using original HumanInteractome.txt/.tsv and
# hgnc_complete_set.txt in BASE_DIR or PREVIOUS_DIR (same files as prior run).
# Outputs are in a NEW run directory. Inputs never overwritten.

BASE_DIR <- "D:/My Desktop/drug/6.AD network"
PREVIOUS_DIR <- "D:/My Desktop/drug/5.endoAD_109_screen"
N_PERM <- 10000L
RANDOM_SEED <- 20261105L
N_WORKERS <- 2L                    # change to 1 if insufficient RAM
MIN_BIN_SIZE <- 100L
CHUNK_SIZE <- 100L
DISTANCE_BLOCK_SIZE <- 200L
MAX_DISTANCE_GB <- 2.0            # raw byte distance matrix alone (not total RAM)
MODULE_MODE <- "mapped_seeds"      # alternative: "seed_induced_LCC"
INSTALL_MISSING <- TRUE
SCRIPT_VERSION <- "1.0.0"

run_AD37 <- function() {
  stopifnot(dir.exists(BASE_DIR), N_PERM >= 100L, CHUNK_SIZE >= 1L,
            MIN_BIN_SIZE >= 2L, N_WORKERS >= 1L,
            MODULE_MODE %in% c("mapped_seeds", "seed_induced_LCC"))
  base <- normalizePath(BASE_DIR, winslash = "/", mustWork = TRUE)
  root <- file.path(base, "AD37_results")
  cache <- file.path(base, "AD37_cache")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = root)
  dir.create(out)
  finished <- FALSE
  on.exit(if (!finished) writeLines(
    "INCOMPLETE RUN: check run.log, fix the problem, rerun. Cached chunks can be reused.",
    file.path(out, "INCOMPLETE_RUN.txt")), add = TRUE)
  logit <- function(...) {
    msg <- paste0(format(Sys.time(), "%H:%M:%S"), " ", paste0(..., collapse = ""))
    message(msg)
    cat(msg, "\n", file = file.path(out, "run.log"), append = TRUE)
  }
  write_csv <- function(x, file) {
    if (!is.data.frame(x)) x <- as.data.frame(x)
    data.table::fwrite(x, file.path(out, file), na = "", bom = TRUE)
  }
  clean <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
  split_pipe <- function(x) {
    x <- clean(x)
    if (!nzchar(x)) return(character())
    z <- unique(clean(unlist(strsplit(x, "[|;]", perl = TRUE))))
    z[nzchar(z)]
  }
  name_key <- function(x) {
    x <- tolower(clean(x))
    x <- gsub("[[:space:]]+", " ", x)
    x <- gsub("[[:space:]_.-]+", "", x)
    x
  }
  hash_obj <- function(x) {
    f <- tempfile(); on.exit(unlink(f), add = TRUE)
    saveRDS(x, f, compress = FALSE, version = 3)
    unname(tools::md5sum(f))
  }
  atomic_rds <- function(x, dest) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("partial_", tmpdir = dirname(dest))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = FALSE, version = 3)
    if (file.exists(dest)) unlink(dest)
    if (!file.rename(tmp, dest)) stop("Could not save ", dest)
  }
  find_one <- function(dir, pattern, label, required = TRUE) {
    if (!dir.exists(dir)) {
      if (required) stop("Folder missing: ", dir)
      return(NA_character_)
    }
    fs <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(fs) != 1L) {
      if (required) stop("Expected exactly one ", label, " in ", dir,
                         "; found ", length(fs), ": ", paste(basename(fs), collapse = ", "))
      return(NA_character_)
    }
    normalizePath(fs, winslash = "/", mustWork = TRUE)
  }
  choose_source <- function(pattern, label) {
    for (dir in unique(c(base, PREVIOUS_DIR))) {
      ff <- find_one(dir, pattern, label, required = FALSE)
      if (!is.na(ff)) return(ff)
    }
    stop("Cannot locate ", label, " in BASE_DIR or PREVIOUS_DIR.")
  }
  # Install only required packages. This script uses base graphics, not ggplot2.
  pkgs <- c("data.table", "readxl", "igraph", "matrixStats")
  lib <- file.path(cache, paste0("R_lib_", getRversion()))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && INSTALL_MISSING)
    utils::install.packages(missing, repos = "https://cloud.r-project.org", lib = lib)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

  logit("Output directory: ", out)
  atc_file <- find_one(base, "^37ATC\\.xlsx$", "37ATC.xlsx")
  seed_file <- find_one(base, "^Table S13\\.xlsx$", "Table S13.xlsx")
  bbb_file <- find_one(base, "^B3DB_classification\\.tsv$", "B3DB_classification.tsv")
  cand_file <- find_one(base, "^14_nominal_P_candidates(\\([0-9]+\\))?\\.csv$", "37-drug candidate CSV")
  drug <- data.table::fread(cand_file, data.table = FALSE, encoding = "UTF-8")
  must_have <- c("DrugBank_ID", "DrugBank_Name", "NP_Z", "NP_Network_target_HGNC_IDs",
                 "NP_N_network_targets", "NP_Empirical_P")
  if (!all(must_have %in% names(drug)))
    stop("Candidate CSV missing required columns: ", paste(setdiff(must_have, names(drug)), collapse = ", "))
  if (nrow(drug) != 37L || anyDuplicated(drug$DrugBank_ID))
    stop("Expected 37 distinct DrugBank IDs, found ", nrow(drug), " rows; duplicates? ",
         anyDuplicated(drug$DrugBank_ID) > 0)
  drug$DrugBank_ID <- clean(drug$DrugBank_ID)
  drug$DrugBank_Name <- clean(drug$DrugBank_Name)

  # ------- (1) ATC: all valid codes, no arbitrary primary code assignment -------
  logit("[1/5] ATC classification ...")
  atc <- as.data.frame(readxl::read_excel(atc_file, sheet = 1), stringsAsFactors = FALSE)
  if (!all(c("DrugBank_Name", "Queried_ATC") %in% names(atc)))
    stop("37ATC.xlsx requires columns DrugBank_Name and Queried_ATC")
  atc$DrugBank_Name <- clean(atc$DrugBank_Name)
  atc$Queried_ATC <- clean(atc$Queried_ATC)
  if (anyDuplicated(name_key(atc$DrugBank_Name))) stop("Duplicate names in ATC input")
  atc_pos <- match(name_key(drug$DrugBank_Name), name_key(atc$DrugBank_Name))
  if (anyNA(atc_pos)) stop("ATC names do not match: ",
                         paste(drug$DrugBank_Name[is.na(atc_pos)], collapse = ", "))
  if (nrow(atc) != nrow(drug)) stop("ATC row count differs from 37 drugs")
  atc_letters <- c(A="Alimentary tract and metabolism", B="Blood and blood forming organs",
                   C="Cardiovascular system", D="Dermatologicals", G="Genito-urinary system and sex hormones",
                   H="Systemic hormonal preparations (excl. sex hormones and insulins)",
                   J="Antiinfectives for systemic use", L="Antineoplastic and immunomodulating agents",
                   M="Musculo-skeletal system", N="Nervous system",
                   P="Antiparasitic products, insecticides and repellents",
                   R="Respiratory system", S="Sensory organs", V="Various")
  parse_atc <- function(s) {
    s <- toupper(gsub("[[:space:]]+", "", clean(s)))
    parts <- unique(unlist(strsplit(s, "[/;,|]", perl = TRUE)))
    parts <- parts[nzchar(parts)]
    if (length(parts) && any(!grepl("^[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}$", parts)))
      stop("Invalid ATC code(s): ", s)
    parts
  }
  atc_codes <- lapply(atc$Queried_ATC[atc_pos], parse_atc)
  names(atc_codes) <- drug$DrugBank_ID
  atc_long <- do.call(rbind, lapply(seq_len(nrow(drug)), function(i) {
    codes <- atc_codes[[i]]
    if (!length(codes)) codes <- NA_character_
    data.frame(DrugBank_ID = drug$DrugBank_ID[i], Drug = drug$DrugBank_Name[i],
               ATC_code = codes, ATC_level1 = substr(codes, 1L, 1L), stringsAsFactors = FALSE)
  }))
  atc_long$ATC_level1_name <- unname(atc_letters[atc_long$ATC_level1])
  atc_drug <- data.frame(DrugBank_ID = drug$DrugBank_ID, Drug = drug$DrugBank_Name,
                         All_ATC_codes = vapply(atc_codes, paste, collapse = " | ", character(1)),
                         ATC_code_count = lengths(atc_codes), stringsAsFactors = FALSE)
  atc_drug$ATC_level1_codes <- vapply(seq_len(nrow(drug)), function(i)
    paste(unique(na.omit(atc_long$ATC_level1[atc_long$DrugBank_ID == drug$DrugBank_ID[i]])),
          collapse = " | "), character(1))
  atc_drug$Multiple_ATC_review <- atc_drug$ATC_code_count > 1L
  atc_dist <- unique(atc_long[!is.na(atc_long$ATC_level1), c("DrugBank_ID", "ATC_level1", "ATC_level1_name")])
  atc_dist <- aggregate(DrugBank_ID ~ ATC_level1 + ATC_level1_name,
                        atc_dist, FUN = length)
  names(atc_dist)[3] <- "N_unique_drugs"
  atc_dist <- atc_dist[order(-atc_dist$N_unique_drugs, atc_dist$ATC_level1), ]
  write_csv(atc_long, "01_ATC_all_codes_long.csv")
  write_csv(atc_drug, "02_ATC_by_drug.csv")
  write_csv(atc_drug[atc_drug$Multiple_ATC_review, ], "02a_ATC_multiple_codes_review.csv")
  write_csv(atc_dist, "02b_ATC_level1_distribution.csv")
  pdf(file.path(out, "ATC_level1_distribution.pdf"), width = 10, height = 6)
  par(mar = c(4, 17, 3, 2))
  barplot(rev(atc_dist$N_unique_drugs), horiz = TRUE, las = 1,
          names.arg = rev(paste(atc_dist$ATC_level1, atc_dist$ATC_level1_name)),
          xlab = "Number of unique drugs", main = "ATC level 1 (multi-category drugs counted in each class)")
  dev.off()

  # ------- (2) Network: reuse exact earlier graph or reconstruct from raw inputs -------
  logit("[2/5] Loading original interactome and mapping 3 Fang seed sets ...")
  prev_runs <- if (dir.exists(file.path(PREVIOUS_DIR, "endoAD_109_results")))
    list.files(file.path(PREVIOUS_DIR, "endoAD_109_results"), pattern = "^run_",
               full.names = TRUE, recursive = FALSE) else character()
  prev_runs <- prev_runs[file.exists(file.path(prev_runs, "analysis_objects.rds")) &
                           file.exists(file.path(prev_runs, "analysis_summary.csv")) &
                           !file.exists(file.path(prev_runs, "INCOMPLETE_RUN.txt"))]
  if (length(prev_runs)) {
    fi <- file.info(file.path(prev_runs, "analysis_objects.rds"))
    prev_runs <- prev_runs[order(fi$mtime, decreasing = TRUE)]
  }
  obj_files <- c(file.path(base, "analysis_objects.rds"),
                 file.path(base, "endoAD_109_results", "analysis_objects.rds"),
                 file.path(prev_runs, "analysis_objects.rds"))
  obj_files <- obj_files[file.exists(obj_files)]
  prev_obj_file <- if (length(obj_files)) obj_files[1L] else NA_character_
  graph_origin <- ""
  hgnc_file <- NA_character_
  if (!is.na(prev_obj_file)) {
    logit("Using exact graph from prior analysis: ", prev_obj_file)
    obj <- readRDS(prev_obj_file)
    if (is.null(obj$graph)) stop("analysis_objects.rds missing $graph")
    g <- obj$graph
    if (!igraph::is_igraph(g)) stop("Saved graph object is not an igraph")
    if (igraph::is_directed(g)) stop("Saved graph is directed; expected undirected")
    if (is.null(igraph::V(g)$name)) stop("Saved graph has no node identifiers")
    graph_origin <- prev_obj_file
    rm(obj); gc(FALSE)
  } else {
    hgnc_file <- choose_source("^hgnc_complete_set\\.txt$", "HGNC")
    network_file <- choose_source("^HumanInteractome\\.(txt|tsv)$", "HumanInteractome")
    logit("No analysis_objects.rds; reconstructing graph using raw interactome and HGNC")
    h <- data.table::fread(hgnc_file, sep = "\t", data.table = FALSE,
                           colClasses = "character", encoding = "UTF-8")
    need_h <- c("hgnc_id", "symbol", "status", "locus_group", "entrez_id")
    if (!all(need_h %in% names(h))) stop("HGNC input needs: ", paste(need_h, collapse=", "))
    h <- h[h$status == "Approved", , drop = FALSE]
    h <- h[h$locus_group == "protein-coding gene", , drop = FALSE]
    idkey <- function(x) {
      x <- clean(x); x[!nzchar(x)] <- NA_character_; x
    }
    idx <- function(z) split(h$hgnc_id[!is.na(idkey(z))], idkey(z)[!is.na(idkey(z))])
    by_id <- idx(h$hgnc_id)
    by_entrez <- idx(h$entrez_id)
    by_symbol <- idx(h$symbol)
    by_ensembl <- if ("ensembl_gene_id" %in% names(h)) idx(h$ensembl_gene_id) else list()
    add_split <- function(col) {
      if (!(col %in% names(h))) return(list())
      # Keep row alignment explicitly (the alias fields have different lengths).
      keys <- lapply(h[[col]], function(s) clean(unlist(strsplit(clean(s), "[|;]"))))
      ii <- rep(seq_len(nrow(h)), lengths(keys))
      kk <- unlist(keys, use.names = FALSE)
      use <- nzchar(kk) & !is.na(kk)
      split(h$hgnc_id[ii[use]], kk[use])
    }
    by_prev <- add_split("prev_symbol")
    by_alias <- add_split("alias_symbol")
    map_identifier <- function(x) {
      x <- clean(x)
      if (!nzchar(x)) return(NA_character_)
      vv <- if (grepl("^HGNC:[0-9]+$", x)) by_id[[x]]
      else if (grepl("^[0-9]+$", x)) by_entrez[[x]]
      else if (grepl("^ENSG[0-9]+(\\.[0-9]+)?$", x)) by_ensembl[[sub("\\.[0-9]+$", "", x)]]
      else {
        for (tab in list(by_symbol, by_prev, by_alias)) {
          vv <- tab[[x]]
          if (length(vv)) break
        }
        vv
      }
      vv <- unique(vv)
      if (length(vv) != 1L) NA_character_ else vv
    }
    edges_raw <- data.table::fread(network_file, header = FALSE, sep = "\t",
                                   data.table = FALSE, colClasses = "character", fill = TRUE)
    if (ncol(edges_raw) < 2L) stop("HumanInteractome needs >=2 columns")
    endpoints <- unique(c(clean(edges_raw[[1]]), clean(edges_raw[[2]])))
    endpoint_map <- vapply(endpoints, map_identifier, character(1))
    aa <- unname(endpoint_map[match(clean(edges_raw[[1]]), endpoints)])
    bb <- unname(endpoint_map[match(clean(edges_raw[[2]]), endpoints)])
    valid <- !is.na(aa) & !is.na(bb) & aa != bb
    ed <- unique(data.frame(from = pmin(aa[valid], bb[valid]),
                            to = pmax(aa[valid], bb[valid])))
    if (!nrow(ed)) stop("No valid protein-coding network edges after HGNC mapping")
    g0 <- igraph::graph_from_data_frame(ed, directed = FALSE)
    cmp <- igraph::components(g0)
    g <- igraph::induced_subgraph(g0, which(cmp$membership == which.max(cmp$csize)))
    graph_origin <- network_file
    rm(g0, ed, edges_raw, endpoints, aa, bb, endpoint_map); gc(FALSE)
  }
  nodes <- igraph::V(g)$name
  if (!length(nodes) || anyDuplicated(nodes)) stop("Invalid graph node names")
  degrees <- as.numeric(igraph::degree(g))
  n <- length(nodes)
  logit("HumanInteractome LCC: ", n, " nodes; ", igraph::ecount(g), " edges")

  # Seed mapping uses original NCBI Entrez IDs when possible. A disagreement
  # with a mapped symbol is recorded for manual inspection, never hidden.
  if (is.na(hgnc_file)) hgnc_file <- choose_source("^hgnc_complete_set\\.txt$", "HGNC")
  h <- data.table::fread(hgnc_file, sep = "\t", data.table = FALSE,
                         colClasses = "character", encoding = "UTF-8")
  need_h <- c("hgnc_id", "symbol", "status", "locus_group", "entrez_id")
  if (!all(need_h %in% names(h))) stop("HGNC file missing required columns")
  h <- h[h$status == "Approved", , drop = FALSE]
  make_idx <- function(k) {
    k <- clean(k); use <- nzchar(k)
    split(seq_len(nrow(h))[use], k[use])
  }
  ent_idx <- make_idx(h$entrez_id)
  sym_idx <- make_idx(h$symbol)
  prev_idx <- if ("prev_symbol" %in% names(h)) {
    ks <- lapply(h$prev_symbol, function(s) clean(unlist(strsplit(clean(s), "[|;]"))))
    kk <- unlist(ks, use.names = FALSE); ii <- rep(seq_len(nrow(h)), lengths(ks))
    split(ii[nzchar(kk)], kk[nzchar(kk)])
  } else list()
  alias_idx <- if ("alias_symbol" %in% names(h)) {
    ks <- lapply(h$alias_symbol, function(s) clean(unlist(strsplit(clean(s), "[|;]"))))
    kk <- unlist(ks, use.names = FALSE); ii <- rep(seq_len(nrow(h)), lengths(ks))
    split(ii[nzchar(kk)], kk[nzchar(kk)])
  } else list()
  map_seed <- function(sym, entrez) {
    e <- unique(ent_idx[[entrez]])
    s <- unique(sym_idx[[sym]])
    if (!length(s)) s <- unique(prev_idx[[sym]])
    if (!length(s)) s <- unique(alias_idx[[sym]])
    ei <- if (length(e) == 1L) e else integer()
    si <- if (length(s) == 1L) s else integer()
    conflict <- length(ei) && length(si) && ei != si
    if (conflict) return(c(HGNC_ID=NA_character_, Approved_symbol=NA_character_,
                           Mapping_status="Entrez_symbol_conflict"))
    chosen <- if (length(ei)) ei else si
    if (!length(chosen)) return(c(HGNC_ID=NA_character_, Approved_symbol=NA_character_,
                                 Mapping_status="Unmapped_or_ambiguous"))
    if (h$locus_group[chosen] != "protein-coding gene")
      return(c(HGNC_ID=h$hgnc_id[chosen], Approved_symbol=h$symbol[chosen],
               Mapping_status="Non_protein_coding"))
    c(HGNC_ID=h$hgnc_id[chosen], Approved_symbol=h$symbol[chosen], Mapping_status="Mapped")
  }
  seed_x <- as.data.frame(readxl::read_excel(seed_file, sheet = 1,
                               col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE)
  if (ncol(seed_x) < 22L) stop("Table S13 has unexpected layout (need >=22 columns)")
  spec <- list(AD = c(1L, 2L), Amyloid = c(8L, 9L), Tau = c(21L, 22L))
  seed_raw <- do.call(rbind, lapply(names(spec), function(m) {
    cc <- spec[[m]]
    sy <- clean(seed_x[-c(1,2), cc[1]])
    en <- clean(seed_x[-c(1,2), cc[2]])
    good <- nzchar(sy) | nzchar(en)
    data.frame(Module = m, Input_symbol = sy[good], Input_Entrez = en[good],
               stringsAsFactors = FALSE)
  }))
  mapped <- t(vapply(seq_len(nrow(seed_raw)), function(i)
    map_seed(seed_raw$Input_symbol[i], seed_raw$Input_Entrez[i]), character(3)))
  seed_map <- cbind(seed_raw, as.data.frame(mapped, stringsAsFactors = FALSE))
  seed_map$In_global_LCC <- seed_map$Mapping_status == "Mapped" &
                             !is.na(seed_map$HGNC_ID) & seed_map$HGNC_ID %in% nodes
  seed_map$In_seed_induced_LCC <- FALSE
  module_sets <- list(); module_summary <- list()
  for (m in names(spec)) {
    ii <- which(seed_map$Module == m)
    glob <- sort(unique(match(seed_map$HGNC_ID[ii[seed_map$In_global_LCC[ii]]], nodes)))
    glob <- glob[!is.na(glob)]
    if (length(glob) < 2L) stop("Module ", m, " has fewer than 2 mapped protein nodes")
    sg <- igraph::induced_subgraph(g, glob)
    cc <- igraph::components(sg)
    sel <- which(cc$membership == which.max(cc$csize))
    seed_lcc <- sort(glob[sel])
    seed_map$In_seed_induced_LCC[ii] <- seed_map$In_global_LCC[ii] &
      seed_map$HGNC_ID[ii] %in% nodes[seed_lcc]
    S <- if (MODULE_MODE == "mapped_seeds") glob else seed_lcc
    if (length(S) < 2L) stop("Module ", m, " has <2 nodes under MODULE_MODE=", MODULE_MODE,
                             "; change to mapped_seeds")
    module_sets[[m]] <- S
    module_summary[[m]] <- data.frame(Module = m, Input_seed_rows = length(ii),
      Unique_input_Entrez = length(unique(seed_map$Input_Entrez[ii])),
      HGNC_mapped_protein_coding = sum(seed_map$Mapping_status[ii] == "Mapped"),
      Seeds_in_global_LCC = length(glob), Seed_induced_LCC_nodes = length(seed_lcc),
      Seed_induced_LCC_edges = igraph::ecount(igraph::induced_subgraph(g, seed_lcc)),
      Nodes_used_for_proximity = length(S), Mode = MODULE_MODE,
      stringsAsFactors = FALSE)
  }
  module_summary <- do.call(rbind, module_summary)
  write_csv(seed_map, "03_seed_gene_mapping_and_LCC.csv")
  write_csv(module_summary, "04_module_counts_and_definition.csv")
  for (m in names(module_sets)) {
    write_csv(data.frame(Module = m, Network_index = module_sets[[m]],
                         HGNC_ID = nodes[module_sets[[m]]]),
              paste0("04a_module_nodes_", m, ".csv"))
  }
  logit("Module nodes used: ", paste(paste(names(module_sets),
        lengths(module_sets), sep="="), collapse=", "))

  target_sets <- lapply(drug$NP_Network_target_HGNC_IDs, function(s) {
    ids <- split_pipe(s)
    ix <- match(ids, nodes)
    ix <- sort(unique(ix[!is.na(ix)]))
    ix
  })
  names(target_sets) <- drug$DrugBank_ID
  target_n <- lengths(target_sets)
  if (any(target_n == 0L)) stop("Some candidate drugs have zero mapped network targets")
  if (any(target_n != as.integer(drug$NP_N_network_targets))) {
    bad <- which(target_n != as.integer(drug$NP_N_network_targets))
    stop("Network mismatch: candidate targets differ from original analysis for: ",
         paste(drug$DrugBank_Name[bad], collapse=", "),
         ". Use matching previous analysis_objects.rds / original network files.")
  }
  write_csv(data.frame(DrugBank_ID = drug$DrugBank_ID, Drug = drug$DrugBank_Name,
                       N_network_targets = target_n), "04b_candidate_network_target_check.csv")

  # ------- (3) Shortest-path cache / constrained degree-bin permutations -------
  logit("[3/5] Restoring/building all-pairs distances and performing 10,000 permutations ...")
  graph_key <- hash_obj(list(nodes = nodes, edges = igraph::as_edgelist(g)))
  dcache_dir <- file.path(cache, "distance", graph_key)
  dir.create(dcache_dir, recursive=TRUE, showWarnings=FALSE)
  dist_file <- file.path(dcache_dir, "distance_matrix_raw.rds")
  candidate_cache_files <- c(dist_file)
  prior_run_dir <- if (!is.na(prev_obj_file)) dirname(prev_obj_file) else NA_character_
  if (!is.na(prior_run_dir)) {
    st <- file.path(prior_run_dir, "analysis_settings.rds")
    if (file.exists(st)) {
      x <- tryCatch(readRDS(st), error=function(e) NULL)
      if (!is.null(x) && !is.null(x$distance_cache_file))
        candidate_cache_files <- c(candidate_cache_files, x$distance_cache_file)
    }
  }
  old_dir <- file.path(PREVIOUS_DIR, "endoAD_109_cache", "distances")
  if (dir.exists(old_dir))
    candidate_cache_files <- c(candidate_cache_files,
       list.files(old_dir, pattern="^distance_matrix_raw\\.rds$", recursive=TRUE,
                  full.names=TRUE))
  D <- NULL
  for (cf in unique(candidate_cache_files[file.exists(candidate_cache_files)])) {
    zz <- tryCatch(readRDS(cf), error=function(e) NULL)
    if (is.list(zz) && identical(zz$nodes, nodes) && is.raw(zz$D) &&
        identical(dim(zz$D), c(n,n))) {
      D <- zz$D
      logit("Using existing validated distance matrix: ", cf)
      break
    }
    rm(zz); gc(FALSE)
  }
  if (is.null(D)) {
    size_gb <- as.double(n)^2 / 1024^3
    if (size_gb > MAX_DISTANCE_GB)
      stop(sprintf("Raw distance matrix would require %.2f GiB; MAX_DISTANCE_GB=%.2f.",
                   size_gb, MAX_DISTANCE_GB))
    logit(sprintf("Building new raw distance matrix: %.2f GiB", size_gb))
    D <- matrix(as.raw(0), nrow=n, ncol=n)
    chunks <- split(seq_len(n), ceiling(seq_len(n)/DISTANCE_BLOCK_SIZE))
    for (k in seq_along(chunks)) {
      ix <- chunks[[k]]
      bf <- file.path(dcache_dir, sprintf("block_%05d.rds", k))
      bz <- if (file.exists(bf)) tryCatch(readRDS(bf), error=function(e) NULL) else NULL
      if (is.null(bz) || !identical(bz$ix, ix) || !is.raw(bz$v) ||
          !identical(dim(bz$v), c(length(ix), n))) {
        dd <- igraph::distances(g, v=ix, to=seq_len(n), weights=NA, mode="all")
        if (any(!is.finite(dd) | dd < 0 | dd > 254 | dd != floor(dd)))
          stop("Path lengths outside 0..254, raw representation invalid")
        bz <- list(ix=ix, v=matrix(as.raw(dd), nrow=length(ix), ncol=n))
        atomic_rds(bz, bf)
        rm(dd)
      }
      D[ix,] <- bz$v
      rm(bz)
      if (k %% 10L == 0L || k == length(chunks)) logit("Distance blocks: ", k, "/", length(chunks))
    }
    atomic_rds(list(nodes=nodes, D=D), dist_file)
  }
  if (!all(diag(D) == as.raw(0))) stop("Invalid distance matrix diagonal")

  degree_bins <- function(d, minsize) {
    groups <- list(); current <- integer()
    for (val in sort(unique(d))) {
      current <- c(current, which(d == val))
      if (length(current) >= minsize) {
        groups[[length(groups)+1L]] <- current
        current <- integer()
      }
    }
    if (length(current)) {
      if (!length(groups)) groups[[1L]] <- current
      else groups[[length(groups)]] <- c(groups[[length(groups)]], current)
    }
    bins <- integer(length(d))
    for (k in seq_along(groups)) bins[groups[[k]]] <- k
    bins
  }
  bins <- degree_bins(degrees, MIN_BIN_SIZE)
  bins_report <- do.call(rbind, lapply(sort(unique(bins)), function(b)
    data.frame(Bin=b, Degree_min=min(degrees[bins==b]),
               Degree_max=max(degrees[bins==b]), N_network_nodes=sum(bins==b))))
  write_csv(bins_report, "05_degree_bin_definitions.csv")
  pools <- split(seq_len(n), bins)
  dist_measure <- function(D, S, T) {
    a <- matrix(as.integer(D[S, T, drop=FALSE]), nrow=length(S), ncol=length(T))
    (sum(matrixStats::rowMins(a, useNames=FALSE)) +
     sum(matrixStats::colMins(a, useNames=FALSE))) / (length(S)+length(T))
  }
  draw_degree <- function(ct, pools) {
    unlist(lapply(names(ct), function(b) {
      p <- pools[[b]]; k <- as.integer(ct[[b]])
      if (length(p) < k) stop("A degree bin has too few nodes")
      p[sample.int(length(p), size=k, replace=FALSE)]
    }), use.names=FALSE)
  }
  null_chunk <- function(task) {
    ss <- get(".ad37_state", envir=.GlobalEnv)
    K <- length(ss$target_profiles)
    val <- matrix(NA_real_, nrow=length(task), ncol=K)
    for (j in seq_along(task)) {
      iteration <- task[j]
      set.seed(as.integer(((as.double(ss$seed) +
         as.double(iteration)*104729 + as.double(ss$module_no)*19739) %% 2147483646) + 1))
      S0 <- draw_degree(ss$disease_count, ss$pools)
      for (p in seq_len(K)) {
        T0 <- draw_degree(ss$target_profiles[[p]], ss$pools)
        val[j,p] <- dist_measure(ss$D, S0, T0)
      }
    }
    list(indices=as.integer(task), values=val, profile_keys=ss$profile_keys)
  }
  # Avoid serializing the entire enclosing analysis (graph and full input tables)
  # with worker functions when using PSOCK clusters.
  environment(draw_degree) <- baseenv()
  environment(dist_measure) <- baseenv()
  .worker_fn_env <- new.env(parent=baseenv())
  assign("draw_degree", draw_degree, envir=.worker_fn_env)
  assign("dist_measure", dist_measure, envir=.worker_fn_env)
  environment(null_chunk) <- .worker_fn_env

  # Serial/PSOCK worker state is exported once per module, not once per task.
  run_module <- function(module, S, module_no) {
    targ_counts <- lapply(target_sets, function(T) table(bins[T]))
    profile_strings <- vapply(targ_counts, function(z)
      paste(names(z), as.integer(z), sep=":", collapse="|"), character(1))
    uniq <- sort(unique(profile_strings))
    prof_ix <- match(profile_strings, uniq)
    profiles <- targ_counts[match(uniq, profile_strings)]
    obs <- vapply(target_sets, function(T) dist_measure(D, S, T), numeric(1))
    key <- hash_obj(list(version=SCRIPT_VERSION, graph=graph_key, module=module,
       S=sort(S), bins=bins, profiles=uniq, nperm=N_PERM,
       seed=RANDOM_SEED, chunk=CHUNK_SIZE, mode=MODULE_MODE,
       R=as.character(getRversion()), matrixStats=as.character(packageVersion("matrixStats"))))
    folder <- file.path(cache, "null_permutations", key)
    dir.create(folder, recursive=TRUE, showWarnings=FALSE)
    tasks <- split(seq_len(N_PERM), ceiling(seq_len(N_PERM) / CHUNK_SIZE))
    val <- matrix(NA_real_, nrow=N_PERM, ncol=length(uniq))
    pending <- list()
    for (j in seq_along(tasks)) {
      f <- file.path(folder, sprintf("chunk_%04d.rds", j))
      z <- if (file.exists(f)) tryCatch(readRDS(f), error=function(e) NULL) else NULL
      ok <- !is.null(z) && identical(z$indices, as.integer(tasks[[j]])) &&
            identical(z$profile_keys, uniq) &&
            identical(dim(z$values), c(length(tasks[[j]]), length(uniq))) &&
            all(is.finite(z$values))
      if (ok) val[tasks[[j]],] <- z$values else pending[[length(pending)+1L]] <- list(ix=j, task=tasks[[j]], file=f)
    }
    state <- list(D=D, pools=pools, disease_count=table(bins[S]),
                  target_profiles=profiles, profile_keys=uniq,
                  seed=RANDOM_SEED, module_no=module_no)
    assign(".ad37_state", state, envir=.GlobalEnv)
    on.exit(if (exists(".ad37_state", envir=.GlobalEnv, inherits=FALSE))
               rm(".ad37_state", envir=.GlobalEnv), add=TRUE)
    cl <- NULL
    if (length(pending) && N_WORKERS > 1L) {
      cores <- parallel::detectCores(logical=FALSE)
      nw <- if (is.finite(cores)) min(N_WORKERS, max(1L, cores-1L)) else N_WORKERS
      if (nw > 1L) {
        cl <- parallel::makePSOCKcluster(nw)
        on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add=TRUE)
        worker_setup <- function(lib) {
          .libPaths(lib)
          if (!requireNamespace("matrixStats", quietly=TRUE)) stop("matrixStats absent on worker")
          NULL
        }
        environment(worker_setup) <- baseenv()
        parallel::clusterCall(cl, worker_setup, .libPaths())
        ex <- new.env(parent=emptyenv())
        for (nm in c("draw_degree", "dist_measure", "null_chunk"))
          assign(nm, get(nm, mode="function"), envir=ex)
        assign(".ad37_state", state, envir=ex)
        parallel::clusterExport(cl, ls(ex, all.names=TRUE), envir=ex)
        logit("  ", module, ": ", nw, " parallel workers; ", length(uniq), " distinct target degree profiles")
      }
    }
    done <- N_PERM - sum(vapply(pending, function(x) length(x$task), integer(1)))
    logit("  ", module, ": permutations cached ", done, "/", N_PERM)
    while (length(pending)) {
      count <- if (is.null(cl)) 1L else min(length(cl), length(pending))
      todo <- pending[seq_len(count)]
      task_indices <- lapply(todo, `[[`, "task")
      zlist <- if (is.null(cl)) lapply(task_indices, null_chunk) else
        parallel::parLapply(cl, task_indices, null_chunk)
      for (k in seq_along(todo)) {
        z <- zlist[[k]]
        if (!identical(z$indices, as.integer(todo[[k]]$task)) ||
            !identical(z$profile_keys, uniq) || any(!is.finite(z$values)))
          stop("Invalid permutation batch for ", module)
        val[z$indices,] <- z$values
        atomic_rds(z, todo[[k]]$file)
        done <- done + length(z$indices)
      }
      if (done %% 1000L == 0L || done == N_PERM) logit("  ", module, ": ", done, "/", N_PERM)
      pending <- if (count == length(pending)) list() else pending[-seq_len(count)]
    }
    if (!is.null(cl)) { parallel::stopCluster(cl); cl <- NULL }
    if (any(!is.finite(val))) stop("Incomplete permutations for ", module)
    mean0 <- colMeans(val)
    sd0 <- apply(val, 2L, stats::sd)
    nn <- length(target_sets)
    pvalue <- numeric(nn); zscore <- rep(NA_real_, nn)
    for (i in seq_len(nn)) {
      jj <- prof_ix[i]
      pvalue[i] <- (1 + sum(val[,jj] <= obs[i] + 1e-12))/(N_PERM+1)
      if (sd0[jj] > 0) zscore[i] <- (obs[i]-mean0[jj])/sd0[jj]
    }
    fdr <- stats::p.adjust(pvalue, method="BH")
    data.frame(DrugBank_ID=drug$DrugBank_ID, Drug=drug$DrugBank_Name,
               Module=module, Module_mode=MODULE_MODE, N_seed_nodes=length(S),
               N_drug_targets=target_n, NP_Distance=obs, NP_Null_mean=mean0[prof_ix],
               NP_Null_SD=sd0[prof_ix], NP_Z=zscore, NP_Empirical_P=pvalue,
               NP_BH_FDR_within_37=fdr, NP_Exploratory_Z_below_minus1.5=is.finite(zscore) & zscore < -1.5,
               NP_Nominal=is.finite(zscore) & zscore < 0 & pvalue < 0.05,
               N_permutations=N_PERM, Null_profile=prof_ix,
               stringsAsFactors=FALSE)
  }
  by_module <- lapply(seq_along(module_sets), function(j) {
    m <- names(module_sets)[j]
    logit("Network proximity module ", j, "/", length(module_sets), ": ", m)
    z <- run_module(m, module_sets[[j]], j)
    write_csv(z, paste0("05_network_proximity_", m, ".csv"))
    z
  })
  prox <- do.call(rbind, by_module)
  write_csv(prox, "05_network_proximity_all_modules_long.csv")
  wide <- data.frame(DrugBank_ID=drug$DrugBank_ID, Drug=drug$DrugBank_Name,
                     Meta_AD_Z=as.numeric(drug$NP_Z),
                     Meta_AD_empirical_P=as.numeric(drug$NP_Empirical_P),
                     stringsAsFactors=FALSE)
  for (m in names(module_sets)) {
    zz <- prox[prox$Module==m,]
    pos <- match(wide$DrugBank_ID, zz$DrugBank_ID)
    wide[[paste0(m,"_Z")]] <- zz$NP_Z[pos]
    wide[[paste0(m,"_P")]] <- zz$NP_Empirical_P[pos]
    wide[[paste0(m,"_FDR_within_37")]] <- zz$NP_BH_FDR_within_37[pos]
  }
  score_columns <- paste0(names(module_sets), "_Z")
  score_m <- as.matrix(wide[,score_columns,drop=FALSE])
  wide$N_modules_Z_below_minus1.5 <- rowSums(score_m < -1.5, na.rm=TRUE)
  wide$N_modules_negative_Z <- rowSums(score_m < 0, na.rm=TRUE)
  write_csv(wide, "06_network_proximity_wide_and_counts.csv")
  heatmap_pdf <- file.path(out,"Network_proximity_heatmap.pdf")
  plot_heatmap <- function() {
    mm <- as.matrix(wide[,score_columns,drop=FALSE]); storage.mode(mm) <- "double"
    rownames(mm) <- wide$Drug
    ord <- order(-wide$N_modules_Z_below_minus1.5, wide$Meta_AD_Z)
    mm <- mm[ord,,drop=FALSE]
    # Display missing Z as NA; no replacement or imputation.
    graphics::heatmap(mm, Rowv=NA, Colv=NA, scale="none", margins=c(9,15),
      cexRow=0.65, cexCol=1.0, main=paste("AD seed-set network proximity Z (",MODULE_MODE,")",sep=""),
      col=grDevices::colorRampPalette(c("#B2182B","#F7F7F7","#2166AC"))(101))
  }
  pdf(heatmap_pdf, width=10, height=12)
  tryCatch(plot_heatmap(), error=function(e) logit("Heatmap PDF failed: ",conditionMessage(e)))
  dev.off()
  png(file.path(out,"Network_proximity_heatmap.png"), width=3000, height=3600, res=300)
  tryCatch(plot_heatmap(), error=function(e) logit("Heatmap PNG failed: ",conditionMessage(e)))
  dev.off()

  # ------- (4) B3DB: conservative name match, without guessing salt/isomers -------
  logit("[4/5] Matching candidates against original B3DB classification TSV ...")
  b <- data.table::fread(bbb_file, sep="\t", data.table=FALSE,
                         encoding="UTF-8", fill=TRUE, colClasses="character")
  b_names <- names(b)
  if (!"compound_name" %in% b_names || !"BBB+/BBB-" %in% b_names)
    stop("B3DB input needs compound_name and BBB+/BBB- columns. Found: ",
         paste(b_names, collapse=", "))
  b$compound_name <- clean(b$compound_name)
  b$key <- name_key(b$compound_name)
  get_b <- function(row, col) {
    if (!col %in% names(b)) return(NA_character_)
    as.character(b[[col]][row])
  }
  b_rows <- list(); b_sum <- list()
  for (i in seq_len(nrow(drug))) {
    nm <- drug$DrugBank_Name[i]
    key <- name_key(nm)
    idx <- which(nzchar(key) & b$key == key)
    match_type <- if (length(idx)) "normalized_name_exact" else "not_found_by_name"
    if (!length(idx)) {
      b_sum[[i]] <- data.frame(DrugBank_ID=drug$DrugBank_ID[i], Drug=nm,
        B3DB_status="Not_matched", B3DB_records=0L, B3DB_label=NA_character_,
        B3DB_logBB=NA_character_, B3DB_match_method=match_type,
        B3DB_review_required=TRUE, stringsAsFactors=FALSE)
      next
    }
    for (j in idx) b_rows[[length(b_rows)+1L]] <- data.frame(
      DrugBank_ID=drug$DrugBank_ID[i], Drug=nm, B3DB_compound_name=b$compound_name[j],
      B3DB_label=get_b(j,"BBB+/BBB-"), logBB=get_b(j,"logBB"),
      CID=get_b(j,"CID"), SMILES=get_b(j,"SMILES"), Inchi=get_b(j,"Inchi"),
      Reference=get_b(j,"reference"), Group=get_b(j,"group"),
      B3DB_match_method=match_type, B3DB_row=j, stringsAsFactors=FALSE)
    labels <- unique(clean(b[["BBB+/BBB-"]][idx])); labels <- labels[nzchar(labels)]
    is_unambig <- length(idx) == 1L && length(labels) == 1L && labels %in% c("BBB+", "BBB-")
    b_sum[[i]] <- data.frame(DrugBank_ID=drug$DrugBank_ID[i], Drug=nm,
      B3DB_status=if (is_unambig) "One_name_match_needs_identity_review" else "Multiple_or_inconsistent_matches_review",
      B3DB_records=length(idx), B3DB_label=if (is_unambig) labels else NA_character_,
      B3DB_logBB=if (is_unambig) get_b(idx,"logBB") else NA_character_,
      B3DB_match_method=match_type, B3DB_review_required=TRUE,
      stringsAsFactors=FALSE)
  }
  bbb_match <- if (length(b_rows)) do.call(rbind,b_rows) else
     data.frame(DrugBank_ID=character(), Drug=character(), B3DB_compound_name=character(),
       B3DB_label=character(), logBB=character(), CID=character(), SMILES=character(),
       Inchi=character(), Reference=character(), Group=character(),
       B3DB_match_method=character(), B3DB_row=integer())
  bbb_summary <- do.call(rbind,b_sum)
  write_csv(bbb_match,"07_B3DB_name_matches_all_records.csv")
  write_csv(bbb_summary,"08_B3DB_drug_summary.csv")
  write_csv(bbb_summary[bbb_summary$B3DB_status == "Not_matched",],"08a_B3DB_unmatched_for_review.csv")
  write_csv(bbb_summary[bbb_summary$B3DB_status != "Not_matched",],"08b_B3DB_matched_identity_review.csv")
  # B3DB is not a prediction algorithm; no-match is NOT equivalent to BBB-negative.
  # Name matches alone do not prove molecular identity (salt, isomer, tracer).

  # ------- (5) Integrated evidence + transparent provenance -------
  logit("[5/5] Writing integration, provenance and reproducibility reports ...")
  result <- merge(wide, atc_drug, by="DrugBank_ID", all.x=TRUE, sort=FALSE)
  result <- merge(result, bbb_summary[,setdiff(names(bbb_summary),"Drug")],
                  by="DrugBank_ID", all.x=TRUE, sort=FALSE)
  result <- result[match(drug$DrugBank_ID, result$DrugBank_ID),]
  write_csv(result,"09_integrated_ATC_seed_modules_B3DB.csv")
  source_files <- c(Candidates=cand_file, ATC=atc_file, Seeds=seed_file, B3DB=bbb_file,
                    HGNC=hgnc_file, PPI_graph_source=graph_origin)
  md5 <- tools::md5sum(source_files)
  write_csv(data.frame(Input=names(source_files), File=unname(source_files),
        Size_bytes=file.info(source_files)$size, MD5=unname(md5), stringsAsFactors=FALSE),
        "input_provenance_MD5.csv")
  write_csv(data.frame(Package=pkgs, Version=vapply(pkgs, function(x)
    as.character(utils::packageVersion(x)), character(1))), "package_versions.csv")
  writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"))
  settings <- c(
    paste("Script version:",SCRIPT_VERSION), paste("Base directory:",base),
    paste("Network source:",graph_origin), paste("HGNC file:",hgnc_file),
    paste("Module mode:",MODULE_MODE),
    "mapped_seeds = all unique approved HGNC protein-coding seed genes present in the original interactome global LCC.",
    "seed_induced_LCC = largest connected component of the subgraph induced by the mapped seed nodes only.",
    "Neither construction is asserted to reproduce GPSnet's original validated disease module exactly.",
    paste("Permutation N:",N_PERM), paste("Random seed:",RANDOM_SEED),
    paste("Degree-bin minimum nodes:",MIN_BIN_SIZE), paste("Workers:",N_WORKERS),
    "Unweighted, undirected interactome; distances are bidirectional closest distances.",
    "d(S,T) = (sum_s min_t dist(s,t) + sum_t min_s dist(s,t))/(length(S)+length(T)).",
    "Z = (observed distance - null mean)/null standard deviation.",
    "Empirical one-sided P = (1 + number of random distances <= observed)/(N_perm + 1).",
    "For each permutation both disease and drug sets are independently sampled with degree-bin and size matching.",
    "Within-module BH FDR across 37 drugs is CONDITIONAL on preselection; it does not replace the original 1,440-test FDR.",
    "ATC level-1 counts are per unique drug/category; drugs with multiple categories appear in multiple category counts.",
    "ATC mappings are supplied by 37ATC.xlsx; NOT independently verified against WHO ATC/DDD releases.",
    "B3DB matches are conservative normalized-name matches only; structure/CID identity needs independent review.",
    "B3DB not matched means unknown, NOT BBB-negative; B3DB label is not direct proof of brain exposure.",
    "No clinical trial, PubMed, SwissADME or drug efficacy assessment is performed.",
    "37 drug inputs were selected from the 1,440-drug original screen and are not independent replication.",
    "Fang et al. 2021, DOI 10.1038/s43587-021-00138-z; B3DB Meng et al. 2021, DOI 10.1038/s41597-021-01069-5."
  )
  writeLines(settings,file.path(out,"analysis_methods_and_settings.txt"))
  summary <- data.frame(Item=c("Candidate_drugs", "ATC_codes_total", "ATC_multi_code_drugs",
         "Network_nodes", "Network_edges", "Module_mode", "AD_module_nodes",
         "Amyloid_module_nodes", "Tau_module_nodes", "Permutations_per_module_profile",
         "AD_module_Z_below_minus1.5", "Amyloid_module_Z_below_minus1.5",
         "Tau_module_Z_below_minus1.5", "B3DB_one_record_name_matches", "B3DB_unmatched",
         "Run_complete"),
       Value=c(nrow(drug),nrow(atc_long),sum(atc_drug$Multiple_ATC_review),
               n, igraph::ecount(g), MODULE_MODE, lengths(module_sets)[c("AD","Amyloid","Tau")],
               N_PERM,
               vapply(c("AD","Amyloid","Tau"), function(m)
                 sum(prox$Module == m & prox$NP_Exploratory_Z_below_minus1.5, na.rm=TRUE), integer(1)),
               sum(bbb_summary$B3DB_status=="One_name_match_needs_identity_review"),
               sum(bbb_summary$B3DB_status=="Not_matched"), "YES"))
  write_csv(summary,"00_analysis_summary.csv")
  logit("Completed. Results: ",out)
  finished <- TRUE
  print(summary, row.names=FALSE)
  invisible(out)
}

# Source in RStudio or run: Rscript "AD37_ATC_seed_network_B3DB_one_click.R"
AD37_RESULT_DIR <- run_AD37()
