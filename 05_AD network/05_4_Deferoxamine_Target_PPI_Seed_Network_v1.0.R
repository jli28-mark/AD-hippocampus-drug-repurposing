# ============================================================
# Deferoxamine (DB00746): drug-associated proteins -> actual
# PPI shortest paths -> AD / Amyloid / Tau seed genes
# Version 1.0.0 | 2026-11-06
#
# This is a DISPLAY/AUDIT of existing network proximity results.
# It does not repeat permutation tests or prove drug-target binding.
# Disease modules may overlap; each protein occurs only once.
# ============================================================

BASE_DIR <- "D:/My Desktop/drug/6.AD network"
PREVIOUS_DIR <- "D:/My Desktop/drug/5.endoAD_109_screen"
DRUGBANK_ID <- "DB00746"
DRUG_NAME <- "Deferoxamine"

# Leave blank for automatic discovery, or specify the exact directory/file.
AD37_RUN_DIR <- ""
PRIOR_ANALYSIS_RDS <- ""
DRUG_TARGET_FILE <- ""

# Display controls. All seed-to-nearest-target distances are saved regardless.
N_SEEDS_PER_MODULE <- 7L
MAX_DISPLAY_NODES <- 85L
SHOW_BRIDGE_LABELS <- FALSE
RANDOM_SEED <- 20261106L
INSTALL_MISSING <- TRUE
SCRIPT_VERSION <- "1.0.0"

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
req_packages <- c("igraph", "data.table", "ggplot2", "ggrepel")
missing <- req_packages[!vapply(req_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) && INSTALL_MISSING) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
missing <- req_packages[!vapply(req_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

read_csv <- function(f) data.table::fread(
  f, data.table = FALSE, encoding = "UTF-8", check.names = FALSE,
  na.strings = c("", "NA")
)
write_csv <- function(d, file) data.table::fwrite(d, file, na = "", bom = TRUE)
trim <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
fmt <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "NA"
verify_cols <- function(d, cols, label) {
  absent <- setdiff(cols, names(d))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
}
completed_run <- function(dir, requirements) {
  dir.exists(dir) &&
    !file.exists(file.path(dir, "INCOMPLETE_RUN.txt")) &&
    all(file.exists(file.path(dir, requirements)))
}
latest_completed <- function(root, requirements) {
  if (!dir.exists(root)) return(NA_character_)
  d <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  d <- d[grepl("^run_", basename(d))]
  d <- d[vapply(d, completed_run, logical(1), requirements = requirements)]
  if (!length(d)) return(NA_character_)
  d[order(file.info(d)$mtime, decreasing = TRUE)][1L]
}
find_input <- function(filename, override, dirs) {
  if (nzchar(trim(override))) {
    if (!file.exists(override)) stop("Specified input not found: ", override)
    return(normalizePath(override, winslash = "/", mustWork = TRUE))
  }
  options <- file.path(dirs, filename)
  options <- options[file.exists(options)]
  if (!length(options)) stop(
    "Cannot locate ", filename, ". Expected in one of: ",
    paste(dirs, collapse = " | "), ". Set its override path at the top."
  )
  normalizePath(options[1L], winslash = "/", mustWork = TRUE)
}
collapse_one <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = " | ")

# ------------------------------------------------------------
# Source files: previous RDS + raw DrugTargetNetwork + seed maps
# ------------------------------------------------------------
run_deferoxamine_network <- function() {
  stopifnot(dir.exists(BASE_DIR), N_SEEDS_PER_MODULE >= 1L,
            MAX_DISPLAY_NODES >= 10L)
  base <- normalizePath(BASE_DIR, winslash = "/", mustWork = TRUE)
  module_required <- c("03_seed_gene_mapping_and_LCC.csv",
                       "04_module_counts_and_definition.csv",
                       paste0("04a_module_nodes_", c("AD", "Amyloid", "Tau"), ".csv"))
  run_dir <- if (nzchar(trim(AD37_RUN_DIR))) AD37_RUN_DIR else
    latest_completed(file.path(base, "AD37_results"), module_required)
  if (length(run_dir) != 1L || is.na(run_dir) ||
      !completed_run(run_dir, module_required))
    stop("Completed AD37 run with seed mapping files not found. Set AD37_RUN_DIR explicitly.")

  prev_run <- latest_completed(
    file.path(PREVIOUS_DIR, "endoAD_109_results"), "analysis_objects.rds"
  )
  rds_file <- if (nzchar(trim(PRIOR_ANALYSIS_RDS))) {
    PRIOR_ANALYSIS_RDS
  } else if (file.exists(file.path(base, "analysis_objects.rds"))) {
    file.path(base, "analysis_objects.rds")
  } else if (!is.na(prev_run)) {
    file.path(prev_run, "analysis_objects.rds")
  } else ""
  if (!nzchar(rds_file) || !file.exists(rds_file))
    stop("Cannot find original analysis_objects.rds. Set PRIOR_ANALYSIS_RDS.")
  rds_file <- normalizePath(rds_file, winslash = "/", mustWork = TRUE)
  raw_file <- find_input("DrugTargetNetwork.txt", DRUG_TARGET_FILE,
                         c(base, PREVIOUS_DIR))

  output_root <- file.path(run_dir, "14_Deferoxamine_Target_PPI")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                  tmpdir = output_root)
  dir.create(out)
  log <- function(...) {
    s <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(..., collapse = ""))
    message(s)
    cat(s, "\n", file = file.path(out, "run.log"), append = TRUE)
  }
  finished <- FALSE
  on.exit(if (!finished) writeLines("INCOMPLETE: consult run.log", file.path(out, "INCOMPLETE_RUN.txt")), add = TRUE)
  log("Deferoxamine shortest-path network; version ", SCRIPT_VERSION)
  log("Original analysis object: ", rds_file)
  log("Seed mapping run: ", run_dir)
  log("Raw drug-target source: ", raw_file)

  obj <- readRDS(rds_file)
  if (is.null(obj$graph) || !igraph::is_igraph(obj$graph))
    stop("analysis_objects.rds must contain an igraph at $graph")
  g <- obj$graph
  if (igraph::is_directed(g)) stop("Expected undirected PPI network")
  node_ids <- igraph::V(g)$name
  if (is.null(node_ids) || anyDuplicated(node_ids) || anyNA(node_ids))
    stop("PPI network vertex names must contain unique HGNC IDs")
  idx <- setNames(seq_along(node_ids), node_ids)

  for (slot in c("results", "source_target_mapping", "target_sets")) {
    if (is.null(obj[[slot]])) stop("analysis_objects.rds missing $", slot)
  }
  results <- obj$results
  verify_cols(results, c("DrugBank_ID", "Pool_key"), "Original drug results")
  drug_row <- results[trim(results$DrugBank_ID) == DRUGBANK_ID, , drop = FALSE]
  if (nrow(drug_row) != 1L) stop("Expected exactly one ", DRUGBANK_ID,
                                " in original analysis results; got ", nrow(drug_row))
  poolkey <- trim(drug_row$Pool_key[1L])
  T_saved <- obj$target_sets[[poolkey]]
  if (is.null(T_saved) || !length(T_saved)) stop("Drug has no saved PPI targets in original RDS")
  T_saved <- sort(unique(as.integer(T_saved)))
  if (anyNA(T_saved) || any(T_saved < 1L | T_saved > length(node_ids)))
    stop("Saved target indices are not valid PPI vertex indices")

  target_map <- obj$source_target_mapping
  verify_cols(target_map,
              c("Pool_key", "Source_DrugBank_ID", "Source_protein_ID", "HGNC_ID",
                "gene_symbol", "In_network_LCC"), "Original target mapping")
  mapped <- target_map[trim(target_map$Pool_key) == poolkey, , drop = FALSE]
  if (!nrow(mapped)) stop("No source protein links for pool key: ", poolkey)
  trows <- mapped[!is.na(mapped$In_network_LCC) & mapped$In_network_LCC &
                    !is.na(mapped$HGNC_ID), , drop = FALSE]
  T_ids <- sort(unique(trim(trows$HGNC_ID)))
  T_ids <- T_ids[nzchar(T_ids)]
  if (!setequal(T_ids, node_ids[T_saved]))
    stop("Saved drug target set differs from source_target_mapping; inputs may be mixed across runs")
  T <- sort(unique(unname(as.integer(idx[T_ids]))))
  log("Verified ", length(T), " drug-associated proteins in original PPI")

  # Independently check original raw DrugTargetNetwork links, including aliases.
  raw <- data.table::fread(raw_file, sep = "auto", header = FALSE, quote = "",
                           fill = TRUE, data.table = FALSE,
                           colClasses = "character", encoding = "UTF-8")
  if (ncol(raw) < 2L) stop("Raw DrugTargetNetwork.txt has fewer than 2 columns")
  raw_a <- trim(raw[[1L]]); raw_b <- trim(raw[[2L]])
  raw_pairs <- paste(raw_a[grepl("^DB[0-9]{5}$", raw_a)],
                     raw_b[grepl("^DB[0-9]{5}$", raw_a)], sep = "|||", collapse = NULL)
  saved_pairs <- paste(trim(mapped$Source_DrugBank_ID),
                       trim(mapped$Source_protein_ID), sep = "|||")
  absent_raw <- setdiff(unique(saved_pairs), unique(raw_pairs))
  if (length(absent_raw)) stop(
    "The raw DrugTargetNetwork.txt differs from the links recorded in analysis_objects.rds. ",
    "First unmatched pair: ", absent_raw[1L], ". Use the exact original source file."
  )
  write_csv(mapped, file.path(out, "01_original_drug_target_links_all.csv"))
  log("Raw DrugTargetNetwork cross-check passed (", length(unique(saved_pairs)), " unique source pairs)")

  # Symbol mapping: use HGNC IDs, never gene-symbol string joins.
  symbols <- setNames(rep(NA_character_, length(node_ids)), node_ids)
  if (is.data.frame(obj$network_nodes) &&
      all(c("HGNC_ID", "gene_symbol") %in% names(obj$network_nodes))) {
    nn <- obj$network_nodes
    positions <- match(node_ids, trim(nn$HGNC_ID))
    symbols[] <- trim(nn$gene_symbol[positions])
  }
  for (j in seq_len(nrow(trows))) {
    k <- trim(trows$HGNC_ID[j]); s <- trim(trows$gene_symbol[j])
    if (nzchar(k) && nzchar(s) && k %in% node_ids &&
        (is.na(symbols[[k]]) || !nzchar(symbols[[k]]))) symbols[[k]] <- s
  }

  seed_map <- read_csv(file.path(run_dir, "03_seed_gene_mapping_and_LCC.csv"))
  verify_cols(seed_map, c("Module", "HGNC_ID", "Approved_symbol",
                          "Mapping_status", "In_global_LCC"), "Seed mapping")
  seed_sum <- read_csv(file.path(run_dir, "04_module_counts_and_definition.csv"))
  verify_cols(seed_sum, c("Module", "Mode", "Nodes_used_for_proximity"),
              "Seed module summary")
  modules <- c("AD", "Amyloid", "Tau")
  seed_sets <- setNames(vector("list", length(modules)), modules)
  for (m in modules) {
    f <- file.path(run_dir, paste0("04a_module_nodes_", m, ".csv"))
    mm <- read_csv(f)
    verify_cols(mm, "HGNC_ID", f)
    hids <- unique(trim(mm$HGNC_ID)); hids <- hids[nzchar(hids)]
    if (!length(hids) || any(!hids %in% node_ids))
      stop("Seed set ", m, " has absent HGNC IDs in original graph")
    seed_sets[[m]] <- sort(as.integer(idx[hids]))
    row <- seed_sum[seed_sum$Module == m, , drop = FALSE]
    if (nrow(row) != 1L || length(hids) != as.integer(row$Nodes_used_for_proximity))
      stop("Seed set ", m, " disagrees with module summary")
    if (trim(row$Mode[1]) != "mapped_seeds") {
      log("WARNING: Seed module mode is ", row$Mode[1],
          "; the plot follows the saved mode, not necessarily mapped_seeds")
    }
    x <- seed_map[trim(seed_map$Module) == m & trim(seed_map$HGNC_ID) %in% hids,
                  , drop = FALSE]
    for (j in seq_len(nrow(x))) {
      k <- trim(x$HGNC_ID[j]); s <- trim(x$Approved_symbol[j])
      if (nzchar(k) && nzchar(s) &&
          (is.na(symbols[[k]]) || !nzchar(symbols[[k]]))) symbols[[k]] <- s
    }
    log("Validated ", m, " seed genes: ", length(seed_sets[[m]]))
  }
  symbols[is.na(symbols) | !nzchar(symbols)] <- names(symbols)[is.na(symbols) | !nzchar(symbols)]
  get_sym <- function(ids) unname(symbols[node_ids[ids]])
  is_in <- function(v, m) v %in% seed_sets[[m]]

  # Compute full bidirectional closest distances, exactly as used in the
  # original test. No degree-matched permutations are rerun here.
  all_seed_audit <- list()
  all_target_audit <- list()
  selected_paths <- list()
  optional_paths <- list()
  distances <- list()
  addpath <- function(a, b, m, direction, reason, distance) {
    pp <- as.integer(igraph::shortest_paths(
      g, from = a, to = b, mode = "all", weights = NA, output = "vpath"
    )$vpath[[1L]])
    if (!length(pp) || pp[1L] != a || tail(pp, 1L) != b || length(pp) - 1L != distance)
      stop("Cannot reconstruct validated shortest path for ", m, ": ",
           node_ids[a], " -> ", node_ids[b])
    list(Module = m, Direction = direction, Reason = reason,
         From = a, To = b, Distance = distance, Path = pp)
  }
  choose_index <- function(values, candidates) {
    # Deterministic tie-breaking: lowest HGNC ID among equal distances.
    best <- which(values == min(values))
    best[order(node_ids[candidates[best]])][1L]
  }
  for (m in modules) {
    S <- seed_sets[[m]]
    dt <- igraph::distances(g, v = T, to = S, mode = "all", weights = NA)
    if (!all(is.finite(dt))) stop("Disconnected nodes for module ", m)
    dist_target <- apply(dt, 1L, min)
    dist_seed <- apply(dt, 2L, min)
    observed <- (sum(dist_target) + sum(dist_seed)) / (length(T) + length(S))
    distances[[m]] <- observed
    log(m, " observed bidirectional closest distance = ", fmt(observed))

    # Each associated target -> nearest seed: retain one real shortest path.
    for (i in seq_along(T)) {
      k <- choose_index(dt[i, ], S)
      all_target_audit[[length(all_target_audit) + 1L]] <- data.frame(
        Module = m, DrugBank_ID = DRUGBANK_ID, Target_HGNC_ID = node_ids[T[i]],
        Target_symbol = get_sym(T[i]), Nearest_seed_HGNC_ID = node_ids[S[k]],
        Nearest_seed_symbol = get_sym(S[k]), Minimum_hops = dt[i, k],
        N_tied_nearest_seeds = sum(dt[i, ] == dt[i, k]),
        stringsAsFactors = FALSE)
      selected_paths[[length(selected_paths) + 1L]] <- addpath(
        T[i], S[k], m, "Drug_associated_protein_to_seed",
        "Every mapped drug-associated protein, nearest seed", dt[i, k])
    }
    # Each seed -> nearest target: calculate ALL distances for audit.
    seed_best_target <- integer(length(S))
    for (i in seq_along(S)) {
      k <- choose_index(dt[, i], T)
      seed_best_target[i] <- k
      all_seed_audit[[length(all_seed_audit) + 1L]] <- data.frame(
        Module = m, Seed_HGNC_ID = node_ids[S[i]], Seed_symbol = get_sym(S[i]),
        Nearest_target_HGNC_ID = node_ids[T[k]],
        Nearest_target_symbol = get_sym(T[k]), Minimum_hops = dt[k, i],
        N_tied_nearest_targets = sum(dt[, i] == dt[k, i]),
        stringsAsFactors = FALSE)
    }

    # Representative seed paths: favor module-distinctive genes for display.
    # This choice affects only the visualization, not the network distance.
    distinctive <- switch(m,
      AD = setdiff(S, union(seed_sets$Amyloid, seed_sets$Tau)),
      Amyloid = setdiff(S, seed_sets$Tau),
      Tau = setdiff(S, seed_sets$Amyloid))
    rank_seeds <- function(ix) {
      if (!length(ix)) return(integer())
      ix[order(dist_seed[ix], node_ids[S[ix]])]
    }
    primary <- rank_seeds(which(S %in% distinctive))
    secondary <- rank_seeds(setdiff(seq_along(S), primary))
    choice <- head(c(primary, secondary), N_SEEDS_PER_MODULE)
    for (i in choice) {
      t_idx <- seed_best_target[i]
      optional_paths[[length(optional_paths) + 1L]] <- addpath(
        S[i], T[t_idx], m, "Seed_to_nearest_drug_associated_protein",
        "Representative seed-to-target shortest path (display only)", dt[t_idx, i])
    }
  }
  targets_audit <- do.call(rbind, all_target_audit)
  seeds_audit <- do.call(rbind, all_seed_audit)
  write_csv(targets_audit, file.path(out, "02_target_to_nearest_seed_ALL.csv"))
  write_csv(seeds_audit, file.path(out, "03_seed_to_nearest_target_ALL.csv"))

  # Verify observed distances match the saved AD37 network calculation.
  prior_prox <- file.path(run_dir, "05_network_proximity_all_modules_long.csv")
  validation <- data.frame(Module = modules,
                           N_seed_nodes = vapply(seed_sets, length, integer(1)),
                           N_associated_proteins = length(T),
                           Distance_recalculated = unlist(distances[modules]),
                           Distance_original = NA_real_,
                           Absolute_difference = NA_real_)
  if (file.exists(prior_prox)) {
    prior <- read_csv(prior_prox)
    verify_cols(prior, c("DrugBank_ID", "Module", "NP_Distance"), "AD37 proximity results")
    for (i in seq_along(modules)) {
      w <- prior[prior$DrugBank_ID == DRUGBANK_ID & prior$Module == modules[i], , drop = FALSE]
      if (nrow(w) != 1L) stop("Missing/duplicate prior distance for ", modules[i])
      validation$Distance_original[i] <- as.numeric(w$NP_Distance[1L])
      validation$Absolute_difference[i] <- abs(validation$Distance_recalculated[i] -
                                                validation$Distance_original[i])
    }
    if (any(!is.finite(validation$Absolute_difference)) ||
        any(validation$Absolute_difference > 1e-8)) stop(
          "Recalculated network distance does not match original AD37 result. ",
          "Check RDS / seed mapping provenance."
        )
    log("All three recalculated observed distances match original AD37 output")
  } else log("WARNING: original long proximity table absent; independent distance cross-check skipped")
  write_csv(validation, file.path(out, "04_distance_recalculation_validation.csv"))

  # Display: include ALL nearest-seed paths for each drug-associated protein;
  # then optional seed-to-target paths up to the configured display-node cap.
  shown <- selected_paths
  present <- unique(unlist(lapply(shown, function(x) x$Path), use.names = FALSE))
  added <- 0L; skipped <- 0L
  for (p in optional_paths) {
    new <- union(present, p$Path)
    if (length(new) <= MAX_DISPLAY_NODES) {
      shown[[length(shown) + 1L]] <- p
      present <- new
      added <- added + 1L
    } else skipped <- skipped + 1L
  }
  if (length(present) > MAX_DISPLAY_NODES) log(
    "WARNING: even mandatory target-to-seed paths require ", length(present),
    " nodes. Raise MAX_DISPLAY_NODES or simplify using Cytoscape."
  )
  log("Displayed paths: ", length(shown), "; supplementary seed paths added ",
      added, "; skipped ", skipped, "; PPI nodes ", length(present))
  path_table <- do.call(rbind, lapply(shown, function(p) data.frame(
    Module = p$Module, Direction = p$Direction, Reason = p$Reason,
    From_HGNC_ID = node_ids[p$From], From_gene = get_sym(p$From),
    To_HGNC_ID = node_ids[p$To], To_gene = get_sym(p$To),
    Minimum_hops = p$Distance, Path_HGNC_IDs = paste(node_ids[p$Path], collapse = " | "),
    Path_gene_symbols = paste(get_sym(p$Path), collapse = " -> "),
    stringsAsFactors = FALSE)))
  write_csv(path_table, file.path(out, "05_displayed_shortest_paths.csv"))

  edge_records <- list()
  for (p in shown) {
    pp <- p$Path
    if (length(pp) <= 1L) next
    edge_records[[length(edge_records) + 1L]] <- data.frame(
      a = node_ids[head(pp, -1L)], b = node_ids[tail(pp, -1L)],
      Module = p$Module, Direction = p$Direction,
      stringsAsFactors = FALSE)
  }
  if (length(edge_records)) {
    e0 <- do.call(rbind, edge_records)
    e0$key <- paste(pmin(e0$a, e0$b), pmax(e0$a, e0$b), sep = "|||" )
    split_edges <- split(e0, e0$key)
    ppi_edges <- do.call(rbind, lapply(split_edges, function(z) data.frame(
      from = pmin(z$a[1], z$b[1]), to = pmax(z$a[1], z$b[1]),
      Link_type = "PPI_shortest_path", Modules = collapse_one(z$Module),
      N_displayed_path_uses = nrow(z),
      stringsAsFactors = FALSE)))
  } else ppi_edges <- data.frame(from = character(), to = character(),
                                Link_type = character(), Modules = character(),
                                N_displayed_path_uses = integer())
  drug_node <- paste0("DRUG_", DRUGBANK_ID)
  drug_edges <- data.frame(from = drug_node, to = node_ids[T],
                           Link_type = "DrugTargetNetwork_association_UNVERIFIED",
                           Modules = "", N_displayed_path_uses = NA_integer_,
                           stringsAsFactors = FALSE)
  edges <- rbind(ppi_edges, drug_edges)
  present_ids <- sort(unique(c(node_ids[present], node_ids[T])))
  vertex_ids <- c(drug_node, present_ids)
  member <- function(id, set) id %in% node_ids[seed_sets[[set]]]
  classify <- function(id) {
    if (id == drug_node) return("Drug")
    if (id %in% node_ids[T]) {
      if (any(vapply(modules, function(m) member(id, m), logical(1))))
        return("Drug-associated seed (unverified)")
      return("Drug-associated protein (unverified)")
    }
    inA <- member(id, "AD"); inB <- member(id, "Amyloid"); inC <- member(id, "Tau")
    if (inB && inC) return("Amyloid + Tau seed")
    if (inB) return("Amyloid seed")
    if (inC) return("Tau seed")
    if (inA) return("AD-only seed")
    "PPI connector"
  }
  group <- vapply(vertex_ids, classify, character(1))
  node_df <- data.frame(
    name = vertex_ids,
    Gene_symbol = ifelse(vertex_ids == drug_node, DRUG_NAME, unname(symbols[vertex_ids])),
    Group = group,
    Is_drug = vertex_ids == drug_node,
    Is_drug_associated = vertex_ids %in% node_ids[T],
    Is_AD_seed = vertex_ids %in% node_ids[seed_sets$AD],
    Is_Amyloid_seed = vertex_ids %in% node_ids[seed_sets$Amyloid],
    Is_Tau_seed = vertex_ids %in% node_ids[seed_sets$Tau],
    Drug_association_direct_binding_verified = FALSE,
    stringsAsFactors = FALSE)
  node_df$Gene_symbol[node_df$name == drug_node] <- DRUG_NAME
  node_df$Gene_symbol[is.na(node_df$Gene_symbol) | !nzchar(node_df$Gene_symbol)] <-
    node_df$name[is.na(node_df$Gene_symbol) | !nzchar(node_df$Gene_symbol)]
  edge_df <- edges
  edge_df$From_symbol <- node_df$Gene_symbol[match(edge_df$from, node_df$name)]
  edge_df$To_symbol <- node_df$Gene_symbol[match(edge_df$to, node_df$name)]
  write_csv(node_df, file.path(out, "06_Cytoscape_nodes.csv"))
  write_csv(edge_df, file.path(out, "07_Cytoscape_edges.csv"))
  vis_g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = node_df)
  igraph::write_graph(vis_g, file.path(out, "08_Deferoxamine_shortest_paths.graphml"), format = "graphml")

  # Figure: reproducible undirected PPI network with distinct dashed
  # drug-associated protein links (NOT proven molecular binding).
  set.seed(RANDOM_SEED)
  coords <- igraph::layout_with_fr(vis_g, niter = 1800,
                                   weights = ifelse(igraph::E(vis_g)$Link_type ==
                                    "PPI_shortest_path", 1, 1.4))
  xy <- as.data.frame(coords)
  names(xy) <- c("x", "y")
  xy$name <- igraph::V(vis_g)$name
  plot_nodes <- merge(node_df, xy, by = "name", all.x = TRUE, sort = FALSE)
  plot_edges <- merge(edges, xy, by.x = "from", by.y = "name", all.x = TRUE, sort = FALSE)
  names(plot_edges)[names(plot_edges) %in% c("x", "y")] <- c("x", "y")
  to_xy <- xy[match(plot_edges$to, xy$name), c("x", "y")]
  plot_edges$xend <- to_xy$x
  plot_edges$yend <- to_xy$y
  if (anyNA(plot_edges$xend) || anyNA(plot_edges$yend)) stop("Graph plotting endpoints missing")
  plot_nodes$point_size <- ifelse(plot_nodes$Is_drug, 8,
                                 ifelse(plot_nodes$Is_drug_associated, 5.4,
                                        ifelse(plot_nodes$Group == "PPI connector", 2.8, 4.2)))
  plot_nodes$label <- ifelse(plot_nodes$Group == "PPI connector" & !SHOW_BRIDGE_LABELS,
                             "", plot_nodes$Gene_symbol)
  pal <- c(
    "Drug" = "#60439B",
    "Drug-associated protein (unverified)" = "#EFAC54",
    "Drug-associated seed (unverified)" = "#D77235",
    "AD-only seed" = "#4F7AB9",
    "Amyloid seed" = "#23998D",
    "Tau seed" = "#C66897",
    "Amyloid + Tau seed" = "#8A61AA",
    "PPI connector" = "#B2B8BE"
  )
  # use source file Z values only as a label; never infer significance from the subset graph
  ztable <- file.path(run_dir, "06_network_proximity_wide_and_counts.csv")
  subtitle <- paste0("Observed bidirectional distance | AD: ", fmt(distances$AD),
    "   Amyloid: ", fmt(distances$Amyloid), "   Tau: ", fmt(distances$Tau))
  if (file.exists(ztable)) {
    z <- read_csv(ztable)
    r <- z[z$DrugBank_ID == DRUGBANK_ID, , drop = FALSE]
    if (nrow(r) == 1L && all(paste0(modules, "_Z") %in% names(r))) {
      subtitle <- paste0("Network proximity Z | AD: ", fmt(as.numeric(r$AD_Z)),
                         "   Amyloid: ", fmt(as.numeric(r$Amyloid_Z)),
                         "   Tau: ", fmt(as.numeric(r$Tau_Z)))
    }
  }
  px <- plot_edges[plot_edges$Link_type == "PPI_shortest_path", , drop = FALSE]
  dx <- plot_edges[plot_edges$Link_type != "PPI_shortest_path", , drop = FALSE]
  ordinary <- plot_nodes[!plot_nodes$Is_drug, , drop = FALSE]
  compound <- plot_nodes[plot_nodes$Is_drug, , drop = FALSE]
  labelled <- plot_nodes[nzchar(plot_nodes$label) & !plot_nodes$Is_drug, , drop = FALSE]
  plot_title <- paste0(DRUG_NAME, " (", DRUGBANK_ID, ") | associated proteins and shortest PPI paths")
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = px, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                          colour="#A3ADB6", linewidth=0.65, alpha=0.75) +
    ggplot2::geom_segment(data = dx, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                          colour="#8166A7", linetype="dashed", linewidth=0.9, alpha=0.8) +
    ggplot2::geom_point(data = ordinary,
                        ggplot2::aes(x=x, y=y, fill=Group, size=point_size),
                        shape=21, colour="#283540", stroke=0.55) +
    ggplot2::geom_point(data = compound, ggplot2::aes(x=x, y=y),
                        shape=23, size=8.5, fill=pal[["Drug"]], colour="#342359", stroke=0.7) +
    ggrepel::geom_text_repel(data = labelled,
                             ggplot2::aes(x=x, y=y, label=label),
                             seed=RANDOM_SEED, size=3.5, fontface="plain",
                             colour="#23313D", box.padding=0.3, point.padding=0.18,
                             max.overlaps=Inf, min.segment.length=0,
                             segment.alpha=0.55, force=1.15, max.time=3) +
    ggrepel::geom_label_repel(data = compound,
                              ggplot2::aes(x=x, y=y, label=Gene_symbol),
                              seed=RANDOM_SEED, size=4.0, colour="#432C70",
                              fill="white", label.size=0.2,
                              box.padding=0.5, point.padding=0.35,
                              min.segment.length=0, max.overlaps=Inf,
                              force=1.0, max.time=3) +
    ggplot2::scale_fill_manual(values=pal, drop=TRUE, name="Node category") +
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal(clip="off") +
    ggplot2::labs(
      title=plot_title,
      subtitle=subtitle,
      caption=paste(
        "Solid = observed HumanInteractome PPI edge on a selected shortest path; dashed = original DrugTargetNetwork association.",
        "Drug associations are not automatically validated direct-binding targets.",
        "Only representative shortest paths are displayed; full gene sets were used to compute the scores.",
        sep="\n")
    ) +
    ggplot2::theme_void(base_size=12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size=15, face="bold", colour="#23313D", margin=ggplot2::margin(b=5)),
      plot.subtitle = ggplot2::element_text(size=11, colour="#536171", margin=ggplot2::margin(b=12)),
      plot.caption = ggplot2::element_text(size=8.7, colour="#627181", hjust=0, margin=ggplot2::margin(t=15)),
      legend.position="bottom", legend.title=ggplot2::element_text(size=10, face="bold"),
      legend.text=ggplot2::element_text(size=9),
      legend.box="vertical", legend.key.size=grid::unit(0.42, "cm"),
      plot.margin=ggplot2::margin(t=16, r=20, b=15, l=20)
    ) + ggplot2::guides(fill=ggplot2::guide_legend(nrow=2, byrow=TRUE,
                                                  override.aes=list(size=4)))
  pdf_path <- file.path(out, "Figure_Deferoxamine_PPI_network.pdf")
  png_path <- file.path(out, "Figure_Deferoxamine_PPI_network_400dpi.png")
  ggplot2::ggsave(pdf_path, plot=p, width=15.5, height=10.5, units="in", bg="white", device="pdf")
  ggplot2::ggsave(png_path, plot=p, width=15.5, height=10.5, units="in", dpi=400, bg="white")
  write_csv(plot_nodes, file.path(out, "09_plot_node_coordinates.csv"))

  manifest <- data.frame(
    File_role=c("Original_109_analysis_object", "Original_DrugTargetNetwork",
                "AD37_seed_mapping", "AD37_module_summary"),
    File=c(rds_file, raw_file,
           file.path(run_dir,"03_seed_gene_mapping_and_LCC.csv"),
           file.path(run_dir,"04_module_counts_and_definition.csv")),
    stringsAsFactors=FALSE)
  manifest$MD5 <- unname(tools::md5sum(manifest$File))
  write_csv(manifest, file.path(out, "10_input_manifest.csv"))
  summary <- data.frame(
    Metric=c("Drug", "DrugBank_ID", "Script_version", "Original_PPI_nodes",
             "Original_PPI_edges", "Drug_associated_proteins_in_PPI",
             "AD_seed_nodes", "Amyloid_seed_nodes", "Tau_seed_nodes",
             "Displayed_PPI_nodes", "Displayed_PPI_edges", "Displayed_paths",
             "Optional_paths_added", "Optional_paths_omitted_by_node_limit"),
    Value=as.character(c(DRUG_NAME, DRUGBANK_ID, SCRIPT_VERSION,
      igraph::vcount(g), igraph::ecount(g), length(T),
      vapply(seed_sets, length, integer(1)), length(present_ids), nrow(ppi_edges),
      length(shown), added, skipped)), stringsAsFactors=FALSE)
  write_csv(summary, file.path(out,"11_analysis_summary.csv"))
  writeLines(c(
    paste0("Script version: ", SCRIPT_VERSION),
    "This is a targeted visual audit, not a new network drug screen.",
    "Input interactome is the original undirected unweighted network from analysis_objects.rds.",
    "Drug-associated proteins use the original source_target_mapping and are independently checked against raw DrugTargetNetwork.txt.",
    "The three saved seed sets are used as-is from the matching completed AD37 run.",
    "For each disease seed, the minimum PPI hop distance to any drug-associated protein is calculated.",
    "For each drug-associated protein, the minimum PPI hop distance to any disease seed is calculated.",
    "Observed distance = (sum(seed->nearest target hops) + sum(target->nearest seed hops)) / (N seeds + N targets).",
    "The observed distances are compared against saved AD37 distances and the script stops if they disagree.",
    "One shortest path is chosen for each shown pair; equally short alternative paths are not shown.",
    "Every target-to-nearest-seed path is included; representative seed-to-nearest-target paths are added within the display node limit.",
    "Full nearest-distance records are saved separately, including seeds not drawn in the network graphic.",
    "AD/Amyloid/Tau seed sets overlap and are not independent; shared proteins are shown only once.",
    "Drug-target lines show source database associations, not necessarily direct binding, mechanism or clinical efficacy.",
    "No disease gene expression direction, drug perturbation or causal regulatory effect can be inferred from these paths.",
    "The plot is a representative subnetwork, not a visualization of every term contributing to the full proximity Z-score."
  ), file.path(out, "12_methods_and_interpretation.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "13_sessionInfo.txt"))
  finished <- TRUE
  log("COMPLETE. Figure PDF: ", pdf_path)
  log("Cytoscape GraphML: ", file.path(out, "08_Deferoxamine_shortest_paths.graphml"))
  cat("\nOUTPUT FOLDER: ", out, "\n", sep="")
  print(validation, row.names=FALSE)
  invisible(out)
}

# Run automatically when Source is clicked in RStudio.
DEFEROXAMINE_NETWORK_OUTPUT <- run_deferoxamine_network()
