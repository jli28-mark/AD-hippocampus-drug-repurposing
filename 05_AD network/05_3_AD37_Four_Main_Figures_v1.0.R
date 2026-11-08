# ============================================================================
# AD37: Four manuscript figure panels from previously completed analyses
# Version: 1.0.0 | 2026-11-06
# Input: the completed AD37_results run, 14_nominal_P_candidates.csv,
#        37ATC(1).xlsx and 04_evidence_summary_all37.csv.
# Output: 4 high-resolution PNG + 4 vector PDF + data and review template.
# IMPORTANT: no network recalculation, no API calls, no pharmacological claims.
# Run in RStudio: source("D:/My Desktop/drug/6.AD network/AD37_Four_Main_Figures_v1.0.R", encoding="UTF-8")
# ============================================================================

BASE_DIR <- "D:/My Desktop/drug/6.AD network"
MANUAL_EXCEL <- ""  # Optional full path; empty = auto-detect 37ATC(1).xlsx
EVIDENCE_CSV <- ""  # Optional full path; empty = auto-detect completed run
SHOW_RAW_PUBMED_AS_BUBBLE_SIZE <- TRUE
PET_IDS <- c("DB09148", "DB09149", "DB09151")
SEED_Z_THRESHOLD <- -1.5

# ---- Package installation (CRAN; only missing packages) ----------------------
needed <- c("ggplot2", "readxl")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!all(vapply(needed, requireNamespace, logical(1), quietly = TRUE)))
  stop("Packages missing: ggplot2 / readxl")

# ---- Helpers ----------------------------------------------------------------
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
read_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8")
}
first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  paths <- paths[file.exists(paths)]
  if (length(paths)) paths[1L] else ""
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
int <- function(x) suppressWarnings(as.integer(as.character(x)))
str1 <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | toupper(z) %in% c("NA", "NULL", "N/A")] <- ""
  z
}
key_name <- function(x) gsub("[^a-z0-9]", "", tolower(str1(x)))
truth <- function(x) toupper(str1(x)) %in% c("TRUE", "YES", "1")
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
pad <- function(df, key, val = NA) { if (!key %in% names(df)) df[[key]] <- rep(val, nrow(df)); df }
sanitize_bbb <- function(z) {
  z <- toupper(str1(z)); z[!z %in% c("BBB+", "BBB-")] <- "Unknown"
  z
}
parse_atc_level1 <- function(value) {
  a <- str1(value)
  vapply(a, function(s) {
    if (!nzchar(s)) return("Unassigned")
    codes <- trimws(unlist(strsplit(toupper(s), "[|;,/]+")))
    codes <- codes[nzchar(codes)]
    valid <- codes[grepl("^[A-Z][0-9]{2}", codes)]
    if (!length(valid)) return("Unassigned")
    paste(sort(unique(substr(valid, 1, 1))), collapse = "/")
  }, character(1))
}

save_panel <- function(plot, stem, width, height, out) {
  pdf_file <- file.path(out, paste0(stem, ".pdf"))
  png_file <- file.path(out, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height,
                  units = "in", device = "pdf", useDingbats = FALSE, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(png_file, plot = plot, width = width, height = height,
                  units = "in", dpi = 400, bg = "white", limitsize = FALSE)
  log_msg("Saved: ", basename(pdf_file), " + ", basename(png_file))
}

theme_manuscript <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.24), color = "#15263A"),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.9), color = "#50627A"),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(0.75), color = "#586980", hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.margin = ggplot2::margin(15, 24, 12, 14),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = "#354A60"),
      axis.text = ggplot2::element_text(color = "#23364A"),
      legend.title = ggplot2::element_text(face = "bold"))
}

# ---- Load most recent completed network/evidence run -------------------------
run_ad37_four_figures <- function(base = BASE_DIR) {
  base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  runs_root <- file.path(base, "AD37_results")
  runs <- if (dir.exists(runs_root)) list.dirs(runs_root, recursive = FALSE, full.names = TRUE) else character()
  network_path <- file.path(runs, "06_network_proximity_wide_and_counts.csv")
  runs <- runs[file.exists(network_path)]
  if (!length(runs)) stop("No completed run containing 06_network_proximity_wide_and_counts.csv under AD37_results")
  ev_files <- file.path(runs, "12_ClinicalTrials_PubMed", "04_evidence_summary_all37.csv")
  ev_runs <- runs[file.exists(ev_files)]
  if (length(ev_runs) && !nzchar(EVIDENCE_CSV)) {
    ev_times <- file.info(file.path(ev_runs, "12_ClinicalTrials_PubMed", "04_evidence_summary_all37.csv"))$mtime
    run <- ev_runs[which.max(ev_times)]
  } else {
    network_times <- file.info(file.path(runs, "06_network_proximity_wide_and_counts.csv"))$mtime
    run <- runs[which.max(network_times)]
  }
  log_msg("Using network run: ", run)

  if (nzchar(EVIDENCE_CSV)) {
    if (!file.exists(EVIDENCE_CSV)) stop("EVIDENCE_CSV does not exist: ", EVIDENCE_CSV)
    evidence_path <- EVIDENCE_CSV
  } else {
    evidence_path <- first_existing(c(
      file.path(run, "12_ClinicalTrials_PubMed", "04_evidence_summary_all37.csv"),
      file.path(base, "04_evidence_summary_all37.csv")
    ))
  }

  candidate_path <- first_existing(file.path(base, c("14_nominal_P_candidates.csv",
                                                       "14_nominal_P_candidates(1).csv")))
  if (!nzchar(candidate_path)) stop("Place 14_nominal_P_candidates.csv in: ", base)
  candidates <- read_csv(candidate_path)
  required_candidates <- c("DrugBank_ID", "DrugBank_Name", "NP_Z")
  if (!all(required_candidates %in% names(candidates)))
    stop("Candidate table needs: ", paste(required_candidates, collapse = ", "))
  if (anyDuplicated(candidates$DrugBank_ID)) stop("DrugBank_ID duplicates in candidate CSV")

  network <- read_csv(file.path(run, "06_network_proximity_wide_and_counts.csv"))
  required_network <- c("DrugBank_ID", "Meta_AD_Z", "AD_Z", "Amyloid_Z", "Tau_Z",
                        "N_modules_Z_below_minus1.5")
  if (!all(required_network %in% names(network)))
    stop("Missing network columns: ", paste(setdiff(required_network, names(network)), collapse = ", "))
  if (anyDuplicated(network$DrugBank_ID)) stop("DrugBank_ID duplicates in network results")
  j <- match(candidates$DrugBank_ID, network$DrugBank_ID)
  if (anyNA(j)) stop("Candidate IDs absent from network results: ", paste(candidates$DrugBank_ID[is.na(j)], collapse = ", "))
  net <- network[j, , drop = FALSE]

  dat <- data.frame(
    DrugBank_ID = candidates$DrugBank_ID,
    Drug = candidates$DrugBank_Name,
    Meta_AD_Z = num(net$Meta_AD_Z),
    AD_Z = num(net$AD_Z),
    Amyloid_Z = num(net$Amyloid_Z),
    Tau_Z = num(net$Tau_Z),
    N_modules = int(net$N_modules_Z_below_minus1.5),
    stringsAsFactors = FALSE
  )
  if (nrow(dat) != 37L) warning("Expected 37 candidates, found ", nrow(dat))
  if (anyNA(dat$N_modules) || anyNA(dat$Meta_AD_Z)) stop("Missing module counts or Meta-AD Z in network CSV")
  recomputed <- as.integer(dat$AD_Z < SEED_Z_THRESHOLD) +
    as.integer(dat$Amyloid_Z < SEED_Z_THRESHOLD) +
    as.integer(dat$Tau_Z < SEED_Z_THRESHOLD)
  if (anyNA(recomputed) || !identical(recomputed, dat$N_modules))
    stop("Module hit counts disagree with AD/Amyloid/Tau Z values; check input/threshold")
  dat$PET_tracer <- dat$DrugBank_ID %in% PET_IDS
  dat$Shortlisted16 <- dat$N_modules >= 2L
  dat$Priority13 <- dat$Shortlisted16 & !dat$PET_tracer

  # Prefer the most recent manually edited Excel values for ATC and BBB.
  xlsx_path <- if (nzchar(MANUAL_EXCEL)) MANUAL_EXCEL else first_existing(
    file.path(base, c("37ATC(1).xlsx", "37ATC.xlsx")))
  dat$ATC_raw <- ""
  dat$BBB_label <- "Unknown"
  if (nzchar(xlsx_path)) {
    if (!file.exists(xlsx_path)) stop("MANUAL_EXCEL does not exist: ", xlsx_path)
    ex <- as.data.frame(readxl::read_excel(xlsx_path), stringsAsFactors = FALSE)
    if (!all(c("DrugBank_Name", "Queried_ATC", "BBB") %in% names(ex)))
      stop("Manual Excel must have DrugBank_Name, Queried_ATC, BBB")
    if (anyDuplicated(key_name(ex$DrugBank_Name))) stop("Duplicate drug names in Excel")
    jj <- match(key_name(dat$Drug), key_name(ex$DrugBank_Name))
    if (anyNA(jj)) stop("Excel misses drugs: ", paste(dat$Drug[is.na(jj)], collapse = ", "))
    dat$ATC_raw <- str1(ex$Queried_ATC[jj])
    dat$BBB_label <- sanitize_bbb(ex$BBB[jj])
    log_msg("Using manually compiled ATC/BBB Excel: ", xlsx_path)
  }
  dat$ATC_level1 <- parse_atc_level1(dat$ATC_raw)

  # Drug-level API retrieval counts: search hits are NOT confirmed AD efficacy.
  dat$NCT_hits <- NA_integer_
  dat$PubMed_hits <- NA_integer_
  dat$PubMed_truncated <- FALSE
  dat$CT_search_status <- "NOT_RETRIEVED"
  dat$PubMed_search_status <- "NOT_RETRIEVED"
  if (nzchar(evidence_path)) {
    ev <- read_csv(evidence_path)
    if (anyDuplicated(ev$DrugBank_ID)) stop("DrugBank_ID duplicates in evidence summary")
    jj <- match(dat$DrugBank_ID, ev$DrugBank_ID)
    if (anyNA(jj)) stop("Evidence summary lacks drugs: ", paste(dat$Drug[is.na(jj)], collapse = ", "))
    ev <- ev[jj, , drop = FALSE]
    for (cname in c("N_ClinicalTrials_records", "N_PubMed_database_hits", "PubMed_truncated",
                    "ClinicalTrials_search_status", "PubMed_search_status")) {
      ev <- pad(ev, cname, NA)
    }
    dat$NCT_hits <- int(ev$N_ClinicalTrials_records)
    dat$PubMed_hits <- int(ev$N_PubMed_database_hits)
    dat$PubMed_truncated <- truth(ev$PubMed_truncated)
    dat$CT_search_status <- str1(ev$ClinicalTrials_search_status)
    dat$PubMed_search_status <- str1(ev$PubMed_search_status)
    # Only a recorded COMPLETE search supports treating zero hits as no hits.
    dat$NCT_hits[!grepl("^COMPLETE", toupper(dat$CT_search_status))] <- NA_integer_
    dat$PubMed_hits[!grepl("^COMPLETE", toupper(dat$PubMed_search_status))] <- NA_integer_
    log_msg("Using ClinicalTrials/PubMed summary: ", evidence_path)
  } else {
    warning("No 04_evidence_summary_all37.csv located. Figures 3/4 will show retrieval as unavailable.")
  }

  n_tested <- 1440L  # fallback from this project's completed analysis
  if ("NP_BH_family_size" %in% names(candidates)) {
    candidates_n <- unique(int(candidates$NP_BH_family_size))
    candidates_n <- candidates_n[is.finite(candidates_n) & candidates_n > 0]
    if (length(candidates_n) == 1L) n_tested <- candidates_n
  }
  n_all <- nrow(dat)
  n_16 <- sum(dat$Shortlisted16)
  n_pet <- sum(dat$Shortlisted16 & dat$PET_tracer)
  n_13 <- sum(dat$Priority13)
  n_3 <- sum(dat$Priority13 & dat$N_modules == 3L)
  if (!n_13) stop("No non-PET follow-up candidates; cannot create four figures")
  lead_pool <- dat[dat$Priority13 & dat$N_modules == 3L, , drop = FALSE]
  if (!nrow(lead_pool)) lead_pool <- dat[dat$Priority13, , drop = FALSE]
  lead <- lead_pool$Drug[which.min(lead_pool$Meta_AD_Z)]
  lead_id <- lead_pool$DrugBank_ID[which.min(lead_pool$Meta_AD_Z)]
  dat$Lead_by_network_rule <- dat$DrugBank_ID == lead_id
  log_msg("Observed: ", n_tested, " tested -> ", n_all, " nominal -> ", n_16,
          " >=2/3 -> ", n_13, " non-PET. Network-ranking lead: ", lead)

  out <- file.path(run, "13_Manuscript_Four_Figures")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv(dat, file.path(out, "00_figure_source_37_drugs.csv"))
  write_csv(dat[dat$Shortlisted16, , drop = FALSE], file.path(out, "00a_shortlisted_16.csv"))
  write_csv(dat[dat$Priority13, , drop = FALSE], file.path(out, "00b_nonPET_13.csv"))

  # ---- Figure 1: auditable candidate-screening flow -------------------------
  main <- data.frame(
    y = c(5.1, 4.0, 2.9, 1.8, 0.7),
    n = c(n_tested, n_all, n_16, n_13, 1L),
    label = c("approved/non-withdrawn drugs tested",
              "nominal proximity: Z < 0 and empirical P < 0.05",
              "Z < -1.5 in at least 2 of 3 seed sets",
              "non-PET exploratory candidates",
              paste0("network-ranking lead: ", lead)),
    stage = c("Universe", "Initial", "Cross-set", "Non-tracer", "Lead"),
    stringsAsFactors = FALSE
  )
  main$display <- paste0(format(main$n, big.mark = ","), "\n", main$label)
  branches <- data.frame(
    y = c(4.0, 2.9),
    display = c(paste0(n_all - n_16, " below 2-of-3 threshold"),
                paste0(n_pet, " PET imaging tracers")),
    stringsAsFactors = FALSE
  )
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = data.frame(y = c(4.55, 3.45, 2.35, 1.25)),
      ggplot2::aes(x = 0, xend = 0, y = y + 0.12, yend = y - 0.12),
      linewidth = 0.75, color = "#7890A8",
      arrow = grid::arrow(type = "closed", length = grid::unit(0.13, "in"))) +
    ggplot2::geom_rect(data = main,
      ggplot2::aes(xmin = -2.35, xmax = 2.35, ymin = y - 0.42, ymax = y + 0.42, fill = stage),
      color = "#FFFFFF", linewidth = 1.2) +
    ggplot2::geom_text(data = main, ggplot2::aes(x = 0, y = y, label = display),
                       size = 4.2, lineheight = 1.2, fontface = "bold", color = "#142C43") +
    ggplot2::geom_segment(data = branches,
      ggplot2::aes(x = 2.38, xend = 3.10, y = y, yend = y),
      color = "#A2ABB4", linewidth = 0.65,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.11, "in"))) +
    ggplot2::geom_rect(data = branches,
      ggplot2::aes(xmin = 3.12, xmax = 6.75, ymin = y - 0.36, ymax = y + 0.36),
      fill = "#EFF2F6", color = "#CBD4DC", linewidth = 0.65) +
    ggplot2::geom_text(data = branches,
      ggplot2::aes(x = 4.94, y = y, label = display), size = 3.65, color = "#506277") +
    ggplot2::scale_fill_manual(values = c("Universe"="#DDEBF4", "Initial"="#CAE3F2",
       "Cross-set"="#B9E1DA", "Non-tracer"="#A7D4C9", "Lead"="#F2D1A5"), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-2.75, 7.10), ylim = c(0.15, 5.65), expand = FALSE) +
    ggplot2::labs(title = "Exploratory drug prioritization workflow",
      subtitle = "From a fixed AD transcriptomic signature to a network-ranking lead",
      caption = paste0("Two-of-three is an exploratory ranking criterion; the AD, amyloid and tau seed sets overlap. ",
                       "No drug passed BH-FDR < 0.05 in the original library-wide screen.")) +
    ggplot2::theme_void(base_size = 12, base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, color = "#15263A"),
                   plot.subtitle = ggplot2::element_text(size = 11, color = "#50627A"),
                   plot.caption = ggplot2::element_text(size = 9, color = "#637183", hjust = 0),
                   plot.margin = ggplot2::margin(20, 30, 18, 20))
  save_panel(p1, "Figure1_Drug_screening_flow", 12.4, 8.3, out)

  # ---- Figure 2: network proximity heatmap for all 16 -----------------------
  sixteen <- dat[dat$Shortlisted16, , drop = FALSE]
  sixteen <- sixteen[order(-sixteen$N_modules, sixteen$Meta_AD_Z, sixteen$Drug), , drop = FALSE]
  sixteen$Drug_plot <- factor(sixteen$Drug, levels = rev(sixteen$Drug))
  zvars <- c("Meta_AD_Z", "AD_Z", "Amyloid_Z", "Tau_Z")
  hm <- do.call(rbind, lapply(seq_along(zvars), function(k) {
    data.frame(Drug_plot = sixteen$Drug_plot, X = k,
               Z = sixteen[[zvars[k]]], stringsAsFactors = FALSE)
  }))
  hm$z_label <- sprintf("%.2f", hm$Z)
  hm$Z_fill <- pmax(-5.5, pmin(0.3, hm$Z))
  sixteen$BBB_label <- sanitize_bbb(sixteen$BBB_label)
  sixteen$flag <- ifelse(sixteen$PET_tracer, "PET tracer",
                         ifelse(sixteen$Lead_by_network_rule, "Lead", ""))
  p2 <- ggplot2::ggplot(hm, ggplot2::aes(x = X, y = Drug_plot)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Z_fill), color = "white", linewidth = 1.1,
                       width = 0.93, height = 0.89) +
    ggplot2::geom_text(ggplot2::aes(label = z_label), size = 3.5, color = "#122B45") +
    ggplot2::geom_text(data = sixteen,
      ggplot2::aes(x = 5.18, y = Drug_plot, label = paste0(N_modules, "/3")),
      inherit.aes = FALSE, size = 3.7, fontface = "bold", color = "#21394D") +
    ggplot2::geom_text(data = sixteen,
      ggplot2::aes(x = 6.23, y = Drug_plot, label = BBB_label, color = BBB_label),
      inherit.aes = FALSE, size = 3.6, fontface = "bold") +
    ggplot2::geom_text(data = sixteen,
      ggplot2::aes(x = 7.20, y = Drug_plot, label = ATC_level1),
      inherit.aes = FALSE, size = 3.3, color = "#455C74") +
    ggplot2::geom_text(data = sixteen,
      ggplot2::aes(x = 8.25, y = Drug_plot, label = flag),
      inherit.aes = FALSE, size = 3.1, fontface = "bold", color = "#8D4F32") +
    ggplot2::scale_fill_gradient2(low = "#6799C8", mid = "#F3F6FA", high = "#D78C82",
      midpoint = -1.5, limits = c(-5.5, 0.3), name = "Network Z\n(negative = closer)") +
    ggplot2::scale_color_manual(values = c("BBB+" = "#21806C", "BBB-" = "#AE695A",
                                          "Unknown" = "#89929D"), guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5.18, 6.23, 7.20, 8.25),
      labels = c("Meta-AD", "AD seed", "Amyloid", "Tau", "# seed\nsets", "BBB*", "ATC 1", "Flag"),
      limits = c(0.48, 8.77), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = 0.55)) +
    ggplot2::labs(title = "Network proximity profiles of 16 shortlisted drugs",
      subtitle = "Drug-level scores in the original Meta-AD signature and three overlapping seed gene sets",
      x = NULL, y = NULL,
      caption = "PET tracer = excluded from therapeutic follow-up; Lead = first by the stated network-ranking rule.\n*BBB labels are manually compiled and not yet confirmed by drug-level PK evidence. ATC 1 = first-level code(s).") +
    theme_manuscript(12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 10, angle = 12, hjust = 0.55),
                   axis.text.y = ggplot2::element_text(size = 10.5),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right")
  save_panel(p2, "Figure2_16_drug_network_heatmap", 14.7, 9.3, out)

  # ---- Figure 3: network-score/bibliography retrieval bubble chart ----------
  thirteen <- dat[dat$Priority13, , drop = FALSE]
  thirteen <- thirteen[order(thirteen$Meta_AD_Z, thirteen$Drug), , drop = FALSE]
  thirteen$Drug_plot <- factor(thirteen$Drug, levels = rev(thirteen$Drug))
  thirteen$PubMed_size <- ifelse(is.na(thirteen$PubMed_hits), 1, pmax(1, thirteen$PubMed_hits))
  thirteen$NCT_label <- ifelse(is.na(thirteen$NCT_hits), "NCT ?",
                                paste0("NCT ", thirteen$NCT_hits))
  thirteen$counts_label <- paste0(thirteen$N_modules, "/3   |   ", thirteen$NCT_label)
  x_right <- max(-0.17, max(thirteen$Meta_AD_Z, na.rm = TRUE) + 0.9)
  x_text <- -0.95
  p3 <- ggplot2::ggplot(thirteen, ggplot2::aes(y = Drug_plot, x = Meta_AD_Z)) +
    ggplot2::geom_segment(ggplot2::aes(x = -1.0, xend = Meta_AD_Z, y = Drug_plot, yend = Drug_plot),
                          color = "#D6DEE8", linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(size = PubMed_size, color = BBB_label), alpha = 0.85) +
    ggplot2::geom_point(data = thirteen[thirteen$Lead_by_network_rule, , drop = FALSE],
      size = 12.5, shape = 21, fill = NA, color = "#BB6540", stroke = 1.0,
      inherit.aes = TRUE) +
    ggplot2::geom_text(ggplot2::aes(x = x_text, label = counts_label),
                       hjust = 0, size = 3.5, color = "#344C64") +
    ggplot2::scale_color_manual(values = c("BBB+" = "#21806C", "BBB-" = "#BF7968",
                                          "Unknown" = "#8895A3"), name = "Manual BBB label*") +
    ggplot2::scale_size_continuous(range = c(3.2, 10.5),
      breaks = c(1, 10, 50, 150), labels = c("0 / 1", "10", "50", "150"),
      name = "Raw PubMed hits**") +
    ggplot2::scale_x_continuous(limits = c(min(-3.15, min(thirteen$Meta_AD_Z) - 0.25), x_right),
       breaks = seq(-3.0, -1.0, by = 0.5)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = 0.65)) +
    ggplot2::labs(title = "Thirteen non-tracer candidates: network-priority landscape",
      subtitle = "Drugs ordered by proximity to the Meta-AD gene signature",
      x = "Meta-AD network proximity Z (more negative = closer)", y = NULL,
      caption = paste0("Lead (outlined): ", lead, ". The right column shows number of seed sets < -1.5 and raw NCT search hits.\n",
           "*Manual BBB is not a verified measure of human brain exposure. **Raw PubMed search hits are not counts of supportive studies; ",
           "missing searches are displayed with a minimal marker.")) +
    theme_manuscript(12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#F0F3F6"),
      panel.grid.major.x = ggplot2::element_line(color = "#EBEEF2"),
      axis.text.y = ggplot2::element_text(size = 10.2), legend.position = "right")
  save_panel(p3, "Figure3_13_drug_prioritization_bubbles", 14.0, 8.4, out)

  # ---- Manual review template (NEVER replaced once created) -----------------
  review_path <- file.path(out, "Figure4_manual_evidence_review.csv")
  if (!file.exists(review_path)) {
    review_template <- data.frame(
      DrugBank_ID = thirteen$DrugBank_ID, Drug = thirteen$Drug,
      Clinical_AD_result = "Not reviewed",
      Preclinical_AD_result = "Not reviewed",
      Direct_target_status = "Not reviewed",
      Notes = "",
      stringsAsFactors = FALSE
    )
    write_csv(review_template, review_path)
    log_msg("Created manual review template: ", review_path)
  }
  review <- read_csv(review_path)
  mandatory <- c("DrugBank_ID", "Clinical_AD_result", "Preclinical_AD_result", "Direct_target_status")
  if (!all(mandatory %in% names(review))) stop("Manual evidence review file missing required columns")
  if (anyDuplicated(review$DrugBank_ID)) stop("Duplicate DrugBank_ID in Figure4 manual review template")
  ji <- match(thirteen$DrugBank_ID, review$DrugBank_ID)
  if (anyNA(ji)) stop("Figure4 manual review template misses: ", paste(thirteen$Drug[is.na(ji)], collapse = ", "))
  review <- review[ji, , drop = FALSE]
  norm_review <- function(x, allowed) {
    z <- str1(x)
    z[!nzchar(z)] <- "Not reviewed"
    bad <- !z %in% allowed
    if (any(bad)) stop("Unexpected manual evidence status: ", paste(unique(z[bad]), collapse = ", "),
                       "; allowed: ", paste(allowed, collapse = ", "))
    z
  }
  outcome_levels <- c("Not reviewed", "Positive", "Negative", "Mixed", "Inconclusive", "No relevant study")
  target_levels <- c("Not reviewed", "Verified", "Partly verified", "Not verified")
  review$Clinical_AD_result <- norm_review(review$Clinical_AD_result, outcome_levels)
  review$Preclinical_AD_result <- norm_review(review$Preclinical_AD_result, outcome_levels)
  review$Direct_target_status <- norm_review(review$Direct_target_status, target_levels)

  # ---- Figure 4: evidence-STATUS matrix (avoid treating search hits as efficacy)
  columns <- c("AD seed", "Amyloid", "Tau", "BBB*", "ClinicalTrials\nretrieval",
               "PubMed\nretrieval", "AD clinical\nconclusion", "AD preclinical\nconclusion", "Direct target\nQC")
  cells <- list()
  append_cell <- function(id, drug, col, status, label) {
    cells[[length(cells) + 1L]] <<- data.frame(Drug = drug, DrugBank_ID = id,
      Column = col, Status = status, Label = label, stringsAsFactors = FALSE)
  }
  for (i in seq_len(nrow(thirteen))) {
    r <- thirteen[i, , drop = FALSE]
    drug <- r$Drug; id <- r$DrugBank_ID
    zs <- c(r$AD_Z, r$Amyloid_Z, r$Tau_Z)
    for (k in 1:3) {
      z <- as.numeric(zs[k]); yes <- is.finite(z) && z < SEED_Z_THRESHOLD
      append_cell(id, drug, columns[k], if (yes) "Network near" else "Network other",
                  if (yes) "Z < -1.5" else "Z >= -1.5")
    }
    append_cell(id, drug, columns[4], r$BBB_label,
                ifelse(r$BBB_label == "Unknown", "Unknown", r$BBB_label))
    nct <- r$NCT_hits
    append_cell(id, drug, columns[5], if (is.na(nct)) "Not retrieved" else if (nct > 0) "Search hits" else "No hits",
                if (is.na(nct)) "?" else paste0(nct, " hit", ifelse(nct == 1, "", "s")))
    np <- r$PubMed_hits
    append_cell(id, drug, columns[6], if (is.na(np)) "Not retrieved" else if (np > 0) "Search hits" else "No hits",
                if (is.na(np)) "?" else paste0(np, if (r$PubMed_truncated) "+" else ""))
    append_cell(id, drug, columns[7], review$Clinical_AD_result[i], review$Clinical_AD_result[i])
    append_cell(id, drug, columns[8], review$Preclinical_AD_result[i], review$Preclinical_AD_result[i])
    append_cell(id, drug, columns[9], review$Direct_target_status[i], review$Direct_target_status[i])
  }
  mat <- do.call(rbind, cells)
  mat$Column <- factor(mat$Column, levels = columns)
  mat$Drug <- factor(mat$Drug, levels = rev(thirteen$Drug))
  pal <- c("Network near"="#B3D8E8", "Network other"="#EFF3F6",
           "BBB+"="#B5DFC8", "BBB-"="#F1C9B1", "Unknown"="#EEEEF2",
           "Search hits"="#C8D9EF", "No hits"="#EFF2F6", "Not retrieved"="#ECECEF",
           "Not reviewed"="#F2F3F5", "Positive"="#9CD5BC", "Negative"="#E8AEA6",
           "Mixed"="#CFBCE8", "Inconclusive"="#F5E0AB", "No relevant study"="#ECEFF1",
           "Verified"="#9CD5BC", "Partly verified"="#F5E0AB", "Not verified"="#E8AEA6")
  missing_pal <- setdiff(unique(mat$Status), names(pal))
  if (length(missing_pal)) stop("No colors for statuses: ", paste(missing_pal, collapse = ", "))
  mat$Text_color <- ifelse(mat$Status %in% c("Network near", "Search hits", "BBB+"), "#173E55", "#344456")
  p4 <- ggplot2::ggplot(mat, ggplot2::aes(x = Column, y = Drug)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Status), width = 0.95, height = 0.90,
                       color = "#FFFFFF", linewidth = 1.1) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 2.95,
                       color = "#20354A", lineheight = 1.0) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE, expand = ggplot2::expansion(add = 0.55)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = 0.54)) +
    ggplot2::labs(title = "Evidence retrieval and review status of 13 non-tracer candidates",
      subtitle = "Network proximity is calculated; clinical and preclinical outcome classifications require record-level review",
      x = NULL, y = NULL,
      caption = paste0("ClinicalTrials.gov / PubMed cells report SEARCH HITS, not treatment efficacy. ",
         "The three seed sets overlap.\n",
         "*BBB is a manual label whose evidence provenance must be verified. ",
         "Clinical outcome / preclinical outcome / target QC are 'Not reviewed' until you edit Figure4_manual_evidence_review.csv.")) +
    theme_manuscript(11.5) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9.6, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 10.2),
      axis.ticks = ggplot2::element_blank())
  save_panel(p4, "Figure4_13_drug_evidence_status_matrix", 18.5, 8.8, out)
  write_csv(mat, file.path(out, "Figure4_cell_source.csv"))

  # ---- Export reproducibility notes -----------------------------------------
  summary <- data.frame(
    Metric = c("Screened_drugs", "Nominal_network_candidates", "Cross_module_2plus",
               "PET_tracers_removed", "Nontracer_followup", "Nontracer_3of3",
               "Network_ranking_lead", "Lead_DrugBank_ID", "Latest_network_run",
               "Evidence_summary_used", "Manual_ATC_BBB_used",
               "Seed_Z_threshold", "Publication_ready_efficacy_figure4"),
    Value = c(n_tested, n_all, n_16, n_pet, n_13, n_3, lead, lead_id, run,
              if (nzchar(evidence_path)) evidence_path else "NONE",
              if (nzchar(xlsx_path)) xlsx_path else "NONE",
              SEED_Z_THRESHOLD,
              if (any(review$Clinical_AD_result != "Not reviewed" |
                      review$Preclinical_AD_result != "Not reviewed"))
                "PARTIAL manual review; verify record-level citations" else
                "NO: evidence-STATUS visualization only"),
    stringsAsFactors = FALSE)
  write_csv(summary, file.path(out, "00_figure_method_and_input_summary.csv"))
  writeLines(c(
    "AD37 four figures v1.0.0, created from saved analytical outputs, no PPI recomputation.",
    "Figure 1: 1,440 primary drug tests (or the NP_BH_family_size in input), 37 nominal candidates,",
    "          16 with >=2/3 module Z < -1.5, 3 PET tracers removed, 13 non-PET candidates.",
    "          Lead selected as minimum original Meta-AD Z among 3/3 non-PET candidates.",
    "          Selection logic is exploratory and was defined after results were seen.",
    "Figure 2: four actual network Z-score columns, 16 candidates, manual BBB and ATC annotations.",
    "Figure 3: 13 non-tracer drug network score vs raw PubMed search retrieval; counts are not efficacy.",
    "Figure 4: explicit evidence-STATUS matrix. Clinical/preclinical columns are not filled by hit counts.",
    "BBB+/- in the manual Excel are not automatically validated as human brain exposure.",
    "No drug passed BH-FDR < 0.05 in the primary network screen.",
    "All three seed modules overlap. A hit in multiple modules is not independent validation.",
    "ClinicalTrials.gov and PubMed records should be manually reviewed for intervention and AD relevance.",
    "To update Figure 4: edit Figure4_manual_evidence_review.csv (do not alter DrugBank_ID), rerun script.",
    "Allowed AD clinical/preclinical: Not reviewed, Positive, Negative, Mixed, Inconclusive, No relevant study.",
    "Allowed direct target QC: Not reviewed, Verified, Partly verified, Not verified.",
    "Do not label the fourth figure 'efficacy evidence' until review is completed."
  ), file.path(out, "README_figure_interpretation.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  log_msg("ALL FOUR PANELS EXPORTED: ", out)
  invisible(out)
}

AD37_FIGURES_DIR <- run_ad37_four_figures(BASE_DIR)
