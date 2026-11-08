# =============================================================================
# GSE281803: STZ + DFO vs STZ + saline
# Heatmap of DFO-contrast leading-edge genes for APP catabolic process
# Uses previously generated GSEA + DESeq2 VST outputs; no rerun of DESeq2.
# Date: 2026-11-07
# =============================================================================

BASE_DIR <- "D:/My Desktop/drug/10.GSE281803"
RUN_DIR  <- file.path(BASE_DIR, "results_GSE281803_DFO_PPI_v1")
CORE_DIR <- file.path(RUN_DIR, "01_QC_DESeq2")
GSEA_DIR <- file.path(RUN_DIR, "04_Pathway_GSEA")
OUT_DIR  <- file.path(RUN_DIR, "06_DFO_APP_13_LeadingEdge")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

required_packages <- c("data.table", "ggplot2", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Could not load required package: ", pkg)
  }
}

file_gsea <- file.path(GSEA_DIR, "04_FOCUSED_APP_iron_HIF_and_neural_GSEA.csv")
file_vst  <- file.path(CORE_DIR, "07_vst_expression.csv")
file_meta <- file.path(CORE_DIR, "01_sample_metadata_verified.csv")
file_dfo  <- file.path(CORE_DIR, "05_STZDFO_vs_STZSal_all_genes.csv")
file_stz  <- file.path(CORE_DIR, "04_STZSal_vs_Sham_all_genes.csv")
all_inputs <- c(file_gsea, file_vst, file_meta, file_dfo, file_stz)
missing_inputs <- all_inputs[!file.exists(all_inputs)]
if (length(missing_inputs)) stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))

# 1. Select the ACTUAL DFO-treatment leading-edge genes, NOT the STZ-model list.
pathway_id <- "GO_BP::GOBP_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS"
contrast_id <- "STZDFO_vs_STZSal"
gsea <- data.table::fread(file_gsea, data.table = FALSE, check.names = FALSE)
stopifnot(all(c("pathway", "Contrast", "NES", "padj", "LeadingEdge") %in% names(gsea)))
g <- gsea[gsea$pathway == pathway_id & gsea$Contrast == contrast_id, , drop = FALSE]
if (nrow(g) != 1L) stop("Expected exactly one DFO-contrast APP catabolic GSEA entry.")
if (is.na(g$LeadingEdge[1]) || !nzchar(trimws(g$LeadingEdge[1]))) {
  stop("No leading-edge genes found in the selected GSEA row.")
}
genes <- trimws(strsplit(g$LeadingEdge[1], "|", fixed = TRUE)[[1]])
genes <- unique(genes[nzchar(genes)])
cat("DFO leading-edge genes (n=", length(genes), "):\n",
    paste(genes, collapse = ", "), "\n", sep = "")
cat("Pathway NES = ", round(g$NES[1], 3), "; FDR = ",
    signif(g$padj[1], 4), "\n", sep = "")
if (length(genes) != 13L) {
  warning("Leading-edge count differs from the previously observed 13; use current data as source.")
}

# 2. Read/verify the nine rat samples.
meta <- data.table::fread(file_meta, data.table = FALSE, check.names = FALSE)
stopifnot(all(c("Count_sample", "Group", "GSM") %in% names(meta)))
expected_groups <- c("Sham", "STZSal", "STZDFO")
if (nrow(meta) != 9L || !setequal(unique(meta$Group), expected_groups) ||
    any(table(factor(meta$Group, levels = expected_groups)) != 3L)) {
  stop("Expected exactly three Sham, three STZSal and three STZDFO samples.")
}
if (anyDuplicated(meta$Count_sample)) stop("Duplicate count-matrix sample names.")
meta <- meta[order(match(meta$Group, expected_groups)), , drop = FALSE]
meta$GroupLabel <- factor(meta$Group,
                          levels = expected_groups,
                          labels = c("Sham (n=3)", "STZ + saline (n=3)",
                                     "STZ + DFO (n=3)"))
meta$RatLabel <- sub(".*_rat_", "Rat ", meta$Count_sample)
if (anyDuplicated(meta$RatLabel)) stop("Rat labels must be unique.")

vst <- data.table::fread(file_vst, data.table = FALSE, check.names = FALSE)
if (!"Rat_gene" %in% names(vst)) stop("VST file must contain Rat_gene.")
if (anyDuplicated(vst$Rat_gene)) stop("Duplicate gene symbols in VST data.")
if (!all(meta$Count_sample %in% names(vst))) {
  stop("VST sample columns do not match the verified sample metadata.")
}
missing_genes <- setdiff(genes, vst$Rat_gene)
if (length(missing_genes)) stop("Leading-edge genes not found in VST:\n",
                                paste(missing_genes, collapse = ", "))
expr <- as.matrix(vst[match(genes, vst$Rat_gene), meta$Count_sample, drop = FALSE])
rownames(expr) <- genes
storage.mode(expr) <- "double"
if (any(!is.finite(expr))) stop("Non-finite VST values found.")

# 3. Within-gene Z-score calculated across ALL 9 animals, not per group.
z <- t(scale(t(expr), center = TRUE, scale = TRUE))
if (any(!is.finite(z))) stop("A leading-edge gene has zero VST variance.")
long <- data.frame(
  Rat_gene = rep(rownames(z), times = ncol(z)),
  Count_sample = rep(colnames(z), each = nrow(z)),
  VST = as.vector(expr),
  Row_Z = as.vector(z),
  stringsAsFactors = FALSE
)
long$Group <- meta$Group[match(long$Count_sample, meta$Count_sample)]
long$GroupLabel <- meta$GroupLabel[match(long$Count_sample, meta$Count_sample)]
long$RatLabel <- meta$RatLabel[match(long$Count_sample, meta$Count_sample)]
long$Rat_gene <- factor(long$Rat_gene, levels = rev(genes))
long$RatLabel <- factor(long$RatLabel, levels = meta$RatLabel)
long$GroupLabel <- factor(long$GroupLabel, levels = levels(meta$GroupLabel))
data.table::fwrite(long, file.path(OUT_DIR, "01_DFO_13_genes_9_rats_VST_rowZ.csv"))

# 4. Export BOTH contrasts' real DESeq2 estimates for interpretation.
read_effect <- function(file, prefix) {
  x <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  cols <- c("Rat_gene", "log2FoldChange", "pvalue", "padj")
  if (!all(cols %in% names(x))) stop("Missing DESeq2 columns in: ", file)
  y <- x[match(genes, x$Rat_gene), cols, drop = FALSE]
  names(y) <- c("Rat_gene", paste0(prefix, c("_log2FC", "_P", "_FDR")))
  if (anyNA(y$Rat_gene)) stop("One or more leading-edge genes are missing from DESeq2 table")
  y
}
stz_effect <- read_effect(file_stz, "STZ")
dfo_effect <- read_effect(file_dfo, "DFO")
effects <- cbind(stz_effect, dfo_effect[, -1, drop = FALSE])
effects$Opposite_log2FC_direction <- with(effects,
  ifelse(is.na(STZ_log2FC) | is.na(DFO_log2FC), NA,
         STZ_log2FC * DFO_log2FC < 0))
data.table::fwrite(effects, file.path(OUT_DIR, "02_DFO_13_genes_DESeq2_statistics.csv"))

# 5. Publication-oriented plot with fixed group/sample order.
p <- ggplot2::ggplot(long, ggplot2::aes(x = RatLabel, y = Rat_gene, fill = Row_Z)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5, width = 0.97, height = 0.97) +
  ggplot2::facet_grid(. ~ GroupLabel, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradient2(
    low = "#457AA7", mid = "#FFFFFF", high = "#C55C50", midpoint = 0,
    limits = c(-2.5, 2.5), oob = scales::squish,
    name = "VST row Z-score"
  ) +
  ggplot2::labs(
    title = "DFO-contrast leading-edge genes in APP catabolic process",
    subtitle = sprintf("STZ + DFO vs STZ + saline: NES = %.2f; FDR = %.3f (exploratory)",
                       g$NES[1], g$padj[1]),
    x = NULL, y = NULL,
    caption = paste0("Genes selected from the DFO vs STZ GSEA contrast. Each column is one rat. ",
                     "Row Z-scores are descriptive and do not indicate statistical significance.")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.12, "cm"),
    axis.text.y = ggplot2::element_text(size = 10.5, colour = "#242424"),
    axis.text.x = ggplot2::element_text(size = 9.5, colour = "#242424",
                                       angle = 45, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "#EFF1F3", colour = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10),
    plot.caption = ggplot2::element_text(size = 8.5, hjust = 0),
    legend.position = "right"
  )

stem <- file.path(OUT_DIR, "Fig_DFO_APP_catabolic_13_leading_edge")
ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 10.5, height = 6.8,
                units = "in", bg = "white")
ggplot2::ggsave(paste0(stem, ".png"), p, width = 10.5, height = 6.8,
                units = "in", dpi = 400, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(paste0(stem, ".svg"), p, width = 10.5, height = 6.8,
                  units = "in", bg = "white", device = svglite::svglite)
}
writeLines(c(
  paste0("GSEA: ", pathway_id),
  paste0("Contrast: ", contrast_id),
  paste0("NES: ", g$NES[1]),
  paste0("FDR: ", g$padj[1]),
  paste0("Leading-edge count: ", length(genes)),
  paste0("Genes: ", paste(genes, collapse = ", ")),
  "This gene set did not reach FDR < 0.05 in the DFO comparison.",
  "The leading edge is derived from ranked RNA-seq statistics; it does not identify direct drug targets.",
  "Samples are displayed without clustering, with gene-wise scaling across all nine animals."
), file.path(OUT_DIR, "03_analysis_notes.txt"))

cat("\nCompleted. Heatmap and tables saved to:\n", OUT_DIR, "\n", sep = "")
