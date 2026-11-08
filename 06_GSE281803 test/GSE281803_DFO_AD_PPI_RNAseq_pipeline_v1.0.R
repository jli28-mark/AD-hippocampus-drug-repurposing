# =============================================================================
# GSE281803: rat hippocampal RNA-seq reanalysis centered on Deferoxamine
# Drug-associated genes (8), PPI-neighboring AD seed genes, reversal, GSEA
# Version 1.0.0 | 2026-11-07
# This is secondary analysis of published animal RNA-seq; NOT direct-binding
# validation, NOT an independent experiment, NOT evidence of clinical efficacy.
# =============================================================================

# ----------------------------- USER SETTINGS ---------------------------------
BASE_DIR <- "D:/My Desktop/drug/10.GSE281803"
DATA_DIR <- file.path(BASE_DIR, "data")
COUNT_FILE <- file.path(DATA_DIR, "GSE281803_9FineRawcounts.txt")
META_FILE  <- file.path(DATA_DIR, "GSE281803_series_matrix.txt.xlsx")

AD37_RESULTS_ROOT <- "D:/My Desktop/drug/6.AD network/AD37_results"
# If desired, set the exact AD37 run. Empty means newest run with pathway audit.
AD37_RUN_OVERRIDE <- ""
# The PPI audit contains ALL 137 seed genes and nearest target distances.
PPI_AUDIT_NAME <- "03_seed_to_nearest_target_ALL.csv"
MAX_PPI_HOPS <- 2L   # a priori, seed genes within <=2 PPI steps of DFO associations

ADJUST_FOR_BATCH <- FALSE  # set TRUE ONLY after verifying what M3B1/M3B2 mean
COUNT_MIN <- 10L           # gene must have >=10 reads in at least MIN_SAMPLES samples
MIN_SAMPLES <- 3L
FDR_CUTOFF <- 0.05
MIN_ABS_LFC_FOR_DIRECTION <- 0.25 # descriptive ONLY, not a significance rule

RUN_PPI_ANALYSIS <- TRUE
RUN_GSEA <- TRUE
RUN_GO_BP <- TRUE
AUTO_INSTALL_PACKAGES <- TRUE
SEED <- 20261107L

DFO_HUMAN_GENES <- c("APP", "HIF1A", "GMNN", "RGS4", "RGS7", "RGS8", "RGS16", "RGS19")
# -----------------------------------------------------------------------------

run_GSE281803 <- function() {
  set.seed(SEED)
  out <- file.path(BASE_DIR, "results_GSE281803_DFO_PPI_v1")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  core_dir <- file.path(out, "01_QC_DESeq2")
  target_dir <- file.path(out, "02_Eight_associated_genes")
  network_dir <- file.path(out, "03_PPI_neighboring_AD_genes")
  gsea_dir <- file.path(out, "04_Pathway_GSEA")
  fig_dir <- file.path(out, "Figures")
  for (d in c(core_dir,target_dir,network_dir,gsea_dir,fig_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(out, "run.log")
  logit <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
    message(msg)
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  save_csv <- function(x, path) {
    if (!is.data.frame(x)) x <- as.data.frame(x)
    data.table::fwrite(x, path, bom = TRUE, na = "")
  }
  save_plot <- function(p, name, width, height) {
    ggplot2::ggsave(file.path(fig_dir, paste0(name,".pdf")), plot = p,
                    width = width, height = height, units = "in", device = "pdf",
                    bg = "white", limitsize = FALSE)
    ggplot2::ggsave(file.path(fig_dir, paste0(name,"_400dpi.png")), plot = p,
                    width = width, height = height, units = "in", dpi = 400,
                    bg = "white", limitsize = FALSE)
  }
  logit("Starting GSE281803 reanalysis")

  # R package installation: Bioconductor packages from BiocManager; others CRAN.
  core_cran <- c("data.table", "readxl", "ggplot2", "babelgene")
  core_bioc <- c("DESeq2")
  optional_cran <- if (RUN_GSEA) "msigdbr" else character()
  optional_bioc <- if (RUN_GSEA) "fgsea" else character()
  install_cran <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
    if (length(missing) && AUTO_INSTALL_PACKAGES)
      tryCatch(install.packages(missing,repos="https://cloud.r-project.org"),
               error=function(e) logit("CRAN install warning: ",conditionMessage(e)))
  }
  install_bioc <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
    if (length(missing) && AUTO_INSTALL_PACKAGES)
      tryCatch({
        if (!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager",repos="https://cloud.r-project.org")
        BiocManager::install(missing,ask=FALSE,update=FALSE)
      },error=function(e) logit("Bioconductor install warning: ",conditionMessage(e)))
  }
  install_cran(c(core_cran,optional_cran))
  install_bioc(c(core_bioc,optional_bioc))
  required <- c(core_cran,core_bioc)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Missing CORE R packages: ", paste(missing,collapse=", "))
  gsea_packages_present <- !RUN_GSEA ||
     all(vapply(c(optional_cran,optional_bioc),requireNamespace,logical(1),quietly=TRUE))
  if (RUN_GSEA && !gsea_packages_present)
    logit("WARNING: optional GSEA packages missing. Core DESeq2/8-gene analysis will continue.")
  logit("Required core packages loaded: ", paste(required, collapse=", "))

  # ------------------------ 1. INPUT QC / METADATA ---------------------------
  if (!file.exists(COUNT_FILE)) stop("Counts file missing: ", COUNT_FILE)
  if (!file.exists(META_FILE)) stop("Metadata file missing: ", META_FILE)
  logit("[1/7] Reading 9-sample count matrix and GEO series matrix")

  raw <- data.table::fread(COUNT_FILE, data.table=FALSE, check.names=FALSE,
                           sep="\t", encoding="UTF-8")
  if (ncol(raw)!=10L || names(raw)[1]!="FeatureID")
    stop("Expected FeatureID plus 9 sample columns. Actual headers: ",
         paste(names(raw),collapse=" / "))
  gene <- trimws(as.character(raw[[1L]]))
  if (anyNA(gene) || any(!nzchar(gene))) stop("Empty gene identifiers detected")
  sample_cols <- names(raw)[-1L]
  dat <- as.matrix(raw[,-1L,drop=FALSE])
  suppressWarnings(storage.mode(dat) <- "double")
  if (anyNA(dat) || any(!is.finite(dat)) || any(dat<0) || any(dat!=round(dat)))
    stop("Counts contain NA, negative, non-integer or non-finite values")
  if (any(dat > .Machine$integer.max)) stop("Read counts exceed R integer capacity")
  storage.mode(dat) <- "integer"
  rownames(dat) <- gene
  colnames(dat) <- sample_cols

  if (anyDuplicated(gene)) {
    dups <- data.frame(FeatureID=unique(gene[duplicated(gene)]))
    save_csv(dups, file.path(core_dir,"00_duplicate_gene_symbols.csv"))
    logit("WARNING: duplicate gene symbols: ", nrow(dups),
          "; summing raw counts per exact symbol (audit file saved)")
    dat <- rowsum(dat, group=gene, reorder=FALSE)
    storage.mode(dat) <- "integer"
  }
  allzero <- rowSums(dat)==0L
  n_allzero <- sum(allzero)
  dat <- dat[!allzero,,drop=FALSE]
  logit("Gene symbols after duplicate aggregation and all-zero removal: ",nrow(dat))

  md <- as.data.frame(readxl::read_excel(META_FILE, sheet=1,
                     col_types="text", .name_repair="minimal"), stringsAsFactors=FALSE)
  if (ncol(md)<4L || nrow(md)!=9L)
    stop("Metadata must have 9 sample rows and >=4 columns (title, GSM, genotype, treatment)")
  # Use positional columns because the two characteristics columns share a name.
  md <- data.frame(
    Sample_title = trimws(as.character(md[[1L]])),
    GSM = trimws(as.character(md[[2L]])),
    Genotype_raw = trimws(as.character(md[[3L]])),
    Treatment_raw = trimws(as.character(md[[4L]])), stringsAsFactors=FALSE)
  md$Sample_core <- sub("_S[0-9]+$", "", md$Sample_title)
  if (anyDuplicated(md$Sample_core) || anyDuplicated(md$GSM))
    stop("Duplicate sample core IDs or GSM records in metadata")
  ix <- match(sample_cols, md$Sample_core)
  if (anyNA(ix) || anyDuplicated(ix))
    stop("Count columns do not uniquely map to GEO titles: ",
         paste(sample_cols[is.na(ix)],collapse=", "))
  md <- md[ix,,drop=FALSE]
  md$Count_sample <- sample_cols
  md$Group <- ifelse(grepl("_Sham_",md$Sample_core),"Sham",
              ifelse(grepl("_STZSal_",md$Sample_core),"STZSal",
              ifelse(grepl("_STZDFO_",md$Sample_core),"STZDFO",NA_character_)))
  if (anyNA(md$Group)) stop("Unrecognized treatment labels in GEO sample titles")
  if (any(md$Group=="Sham" & !grepl("wild type",md$Genotype_raw,ignore.case=TRUE)) ||
      any(md$Group!="Sham" & !grepl("STZ",md$Genotype_raw,ignore.case=TRUE)) ||
      any(md$Group=="STZDFO" & !grepl("DFO",md$Treatment_raw,ignore.case=TRUE)) ||
      any(md$Group!="STZDFO" & !grepl("saline",md$Treatment_raw,ignore.case=TRUE)))
    stop("GEO genotype/treatment descriptions disagree with inferred groups")
  md$Group <- factor(md$Group, levels=c("Sham","STZSal","STZDFO"))
  if (!identical(as.integer(table(md$Group)),c(3L,3L,3L)))
    stop("Expected 3/3/3 replicates; observed ", paste(table(md$Group),collapse=","))
  md$Batch_candidate <- sub(".*_(M3B[0-9]+)_rat_.*", "\\1", md$Sample_core)
  if (any(!grepl("^M3B[0-9]+$",md$Batch_candidate)))
    md$Batch_candidate <- "UNKNOWN"
  md$Batch_candidate <- factor(md$Batch_candidate)
  rownames(md) <- md$Count_sample
  save_csv(md, file.path(core_dir,"01_sample_metadata_verified.csv"))

  raw_library <- colSums(dat)
  qc <- data.frame(Sample=colnames(dat), GSM=md$GSM, Group=as.character(md$Group),
       Batch_candidate=as.character(md$Batch_candidate),
       Library_counts=as.numeric(raw_library),
       Genes_detected_gt0=colSums(dat>0),
       Genes_with_10plus_counts=colSums(dat>=10L), stringsAsFactors=FALSE)
  if (any(raw_library<=0)) stop("At least one sample has zero library counts")
  save_csv(qc,file.path(core_dir,"02_sample_count_QC.csv"))
  keep <- rowSums(dat>=COUNT_MIN)>=MIN_SAMPLES
  filter_info <- data.frame(Input_unique_genes=nrow(dat)+n_allzero,
          Removed_all_zero=n_allzero, After_nonzero=nrow(dat),
          Kept_for_DESeq2=sum(keep),
          Removed_low_expression=sum(!keep),
          Count_min=COUNT_MIN, Min_samples=MIN_SAMPLES)
  save_csv(filter_info,file.path(core_dir,"03_gene_filter_summary.csv"))
  if (sum(keep)<1000L) stop("Too few genes after filtering; check counts input")
  dat_keep <- dat[keep,,drop=FALSE]
  logit("Genes passing >=",COUNT_MIN," reads in >=",MIN_SAMPLES,
        " samples: ",nrow(dat_keep))

  # ------------------------ 2. DESeq2 MODEL / QC -----------------------------
  logit("[2/7] Fitting DESeq2 negative-binomial model")
  design_formula <- if (ADJUST_FOR_BATCH) stats::as.formula("~ Batch_candidate + Group")
                    else stats::as.formula("~ Group")
  if (ADJUST_FOR_BATCH) {
    if (length(unique(md$Batch_candidate))<2L ||
        qr(stats::model.matrix(design_formula, data=md))$rank <
        ncol(stats::model.matrix(design_formula,data=md)))
      stop("Batch-adjusted design is not identifiable. Check batch annotations")
    logit("WARNING: batch adjustment ENABLED based on M3B tags; verify their interpretation")
  } else logit("M3B1/M3B2 tags shown for QC only, NOT included in differential-expression design")

  dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat_keep,
                 colData=md, design=design_formula)
  dds <- DESeq2::DESeq(dds, quiet=FALSE)
  logit("DESeq2 coefficients: ", paste(DESeq2::resultsNames(dds), collapse=", "))
  res_stz <- DESeq2::results(dds, contrast=c("Group","STZSal","Sham"),
                           alpha=FDR_CUTOFF)
  res_dfo <- DESeq2::results(dds, contrast=c("Group","STZDFO","STZSal"),
                           alpha=FDR_CUTOFF)
  res_to_df <- function(res) {
    y <- as.data.frame(res)
    y$Rat_gene <- rownames(y)
    y <- y[,c("Rat_gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    rownames(y) <- NULL
    y
  }
  stz <- res_to_df(res_stz)
  dfo <- res_to_df(res_dfo)
  save_csv(stz,file.path(core_dir,"04_STZSal_vs_Sham_all_genes.csv"))
  save_csv(dfo,file.path(core_dir,"05_STZDFO_vs_STZSal_all_genes.csv"))
  saveRDS(dds,file.path(core_dir,"DESeq2_model.rds"))
  norm_counts <- DESeq2::counts(dds,normalized=TRUE)
  vst_obj <- DESeq2::vst(dds,blind=TRUE)
  vst_m <- SummarizedExperiment::assay(vst_obj)
  norm_log2 <- log2(norm_counts+1)
  save_csv(data.frame(Rat_gene=rownames(norm_counts),norm_counts,check.names=FALSE),
           file.path(core_dir,"06_normalized_counts.csv"))
  save_csv(data.frame(Rat_gene=rownames(vst_m),vst_m,check.names=FALSE),
           file.path(core_dir,"07_vst_expression.csv"))
  save_csv(data.frame(Sample=colnames(dds),Size_factor=DESeq2::sizeFactors(dds)),
           file.path(core_dir,"08_size_factors.csv"))

  vsd_pca <- stats::prcomp(t(vst_m),center=TRUE,scale.=FALSE)
  pca_var <- 100*vsd_pca$sdev^2/sum(vsd_pca$sdev^2)
  pd <- data.frame(Sample=rownames(vsd_pca$x),
       PC1=vsd_pca$x[,1], PC2=vsd_pca$x[,2],
       Group=md$Group, Batch=md$Batch_candidate,
       GSM=md$GSM,stringsAsFactors=FALSE)
  save_csv(pd,file.path(core_dir,"09_PCA_coordinates.csv"))
  pca_plot <- ggplot2::ggplot(pd,ggplot2::aes(x=PC1,y=PC2,color=Group,shape=Batch))+
      ggplot2::geom_point(size=3.2)+
      ggplot2::geom_text(ggplot2::aes(label=GSM),vjust=-1.05,size=2.7,show.legend=FALSE)+
      ggplot2::labs(x=sprintf("PC1 (%.1f%%)",pca_var[1]),
                    y=sprintf("PC2 (%.1f%%)",pca_var[2]),
                    title="GSE281803: PCA of variance-stabilized counts")+
      ggplot2::theme_classic(base_size=11)+ggplot2::theme(legend.position="bottom")
  save_plot(pca_plot,"Fig_QC_PCA",8.2,6.0)
  corr <- stats::cor(vst_m,method="pearson")
  grDevices::pdf(file.path(fig_dir,"Fig_QC_sample_correlation.pdf"),width=9,height=8)
  graphics::par(mar=c(10,10,3,2))
  graphics::image(seq_len(ncol(corr)),seq_len(ncol(corr)),corr,
          col=grDevices::colorRampPalette(c("#EEF2F6","#4D8DAC","#14364D"))(100),
          axes=FALSE, xlab="",ylab="",main="VST sample Pearson correlation",
          zlim=c(-1,1))
  graphics::axis(1,at=seq_len(ncol(corr)),labels=md$GSM,las=2,cex.axis=0.75)
  graphics::axis(2,at=seq_len(ncol(corr)),labels=md$GSM,las=2,cex.axis=0.75)
  grDevices::dev.off()
  save_csv(data.frame(Sample=rownames(corr),corr,check.names=FALSE),
           file.path(core_dir,"10_sample_correlation.csv"))

  # ------------------------ 3. EIGHT GENES -----------------------------------
  logit("[3/7] Mapping and extracting 8 a priori DFO-associated human proteins")
  get_orthologs <- function(human_genes) {
    human_genes <- unique(human_genes[!is.na(human_genes)&nzchar(human_genes)])
    if (!length(human_genes)) return(data.frame())
    ans <- as.data.frame(babelgene::orthologs(genes=human_genes,
         species="rat",human=TRUE,min_support=3,top=FALSE),stringsAsFactors=FALSE)
    required_cols <- c("human_symbol","symbol","support_n")
    if (!all(required_cols %in% colnames(ans)))
      stop("babelgene returned unexpected columns; cannot safely map orthologs")
    if (!nrow(ans)) return(data.frame(Human_gene=character(),Rat_gene=character(),
                       Ortholog_support=integer(),Mapping_status=character()))
    ans$human_symbol <- as.character(ans$human_symbol)
    ans$symbol <- as.character(ans$symbol)
    ans <- ans[!is.na(ans$human_symbol)&nzchar(ans$human_symbol)&
               !is.na(ans$symbol)&nzchar(ans$symbol),,drop=FALSE]
    ans <- unique(ans[,c("human_symbol","symbol","support_n")])
    # Strong conservative rule: reject gene pairs with non-unique human OR rat partner.
    hm_n <- table(ans$human_symbol); rt_n <- table(ans$symbol)
    ans$Mapping_status <- ifelse(hm_n[ans$human_symbol]==1L & rt_n[ans$symbol]==1L,
                                 "One_to_one", "Ambiguous_orthology")
    names(ans) <- c("Human_gene","Rat_gene","Ortholog_support","Mapping_status")
    ans
  }
  ortho8 <- get_orthologs(DFO_HUMAN_GENES)
  map8 <- merge(data.frame(Human_gene=DFO_HUMAN_GENES,stringsAsFactors=FALSE),
                ortho8,by="Human_gene",all.x=TRUE,sort=FALSE)
  map8 <- map8[match(DFO_HUMAN_GENES,map8$Human_gene),,drop=FALSE]
  map8$Mapping_status[is.na(map8$Mapping_status)] <- "No_mapping"
  map8$Detected_in_raw_counts <- !is.na(map8$Rat_gene) & map8$Rat_gene %in% rownames(dat)
  map8$Included_in_DESeq2 <- !is.na(map8$Rat_gene) &
         map8$Rat_gene %in% rownames(dat_keep) & map8$Mapping_status=="One_to_one"
  save_csv(map8,file.path(target_dir,"01_eight_human_rat_ortholog_QC.csv"))
  # Extract all eight, retaining NA for genes absent or excluded.
  sel <- map8$Mapping_status=="One_to_one" & map8$Included_in_DESeq2
  ri <- match(map8$Rat_gene,stz$Rat_gene)
  rj <- match(map8$Rat_gene,dfo$Rat_gene)
  target_tab <- data.frame(map8,
       STZ_log2FC=stz$log2FoldChange[ri], STZ_P=stz$pvalue[ri], STZ_FDR=stz$padj[ri],
       DFO_log2FC=dfo$log2FoldChange[rj], DFO_P=dfo$pvalue[rj], DFO_FDR=dfo$padj[rj],
       stringsAsFactors=FALSE)
  for (col in c("STZ_log2FC","STZ_P","STZ_FDR","DFO_log2FC","DFO_P","DFO_FDR"))
    target_tab[[col]][!sel] <- NA_real_
  # Separate exploratory BH across 8 preselected associated genes (drug contrast).
  target_tab$DFO_panel8_BH <- stats::p.adjust(target_tab$DFO_P,method="BH",n=8L)
  target_tab$STZ_panel8_BH <- stats::p.adjust(target_tab$STZ_P,method="BH",n=8L)
  target_tab$All9_raw_sum <- vapply(map8$Rat_gene, function(g) {
      if (is.na(g) || !(g %in% rownames(dat))) return(NA_real_)
      as.double(sum(dat[g,]))
  }, numeric(1))
  target_tab$N_samples_detected <- vapply(map8$Rat_gene,function(g) {
      if (is.na(g) || !(g %in% rownames(dat))) return(NA_integer_)
      as.integer(sum(dat[g,]>0L))
  },integer(1))
  save_csv(target_tab,file.path(target_dir,"02_eight_genes_two_contrasts_and_FDR.csv"))

  exp_rows <- list()
  for (i in seq_len(nrow(map8))) {
    g <- map8$Rat_gene[i]
    for (j in seq_len(nrow(md))) {
      exp_rows[[length(exp_rows)+1L]] <- data.frame(
         Human_gene=map8$Human_gene[i],Rat_gene=ifelse(is.na(g),"",g),
         Sample=md$Count_sample[j], GSM=md$GSM[j],Group=as.character(md$Group[j]),
         Raw_count=if (!is.na(g) && g %in% rownames(dat)) as.numeric(dat[g,j]) else NA_real_,
         Normalized_count=if (!is.na(g) && g %in% rownames(norm_counts))
            as.numeric(norm_counts[g,j]) else NA_real_,
         VST=if (!is.na(g) && g %in% rownames(vst_m))
            as.numeric(vst_m[g,j]) else NA_real_,stringsAsFactors=FALSE)
    }
  }
  eight_exp <- do.call(rbind,exp_rows)
  save_csv(eight_exp,file.path(target_dir,"03_eight_genes_9sample_expression.csv"))
  plot_df <- eight_exp[is.finite(eight_exp$Normalized_count),,drop=FALSE]
  if (nrow(plot_df)) {
    plot_df$Group <- factor(plot_df$Group,levels=c("Sham","STZSal","STZDFO"))
    g1 <- ggplot2::ggplot(plot_df,ggplot2::aes(x=Group,y=log2(Normalized_count+1),color=Group))+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.07,height=0,seed=SEED),
                          size=2.1,alpha=.9)+
      ggplot2::stat_summary(fun=median,geom="point",shape=95,size=7,color="black")+
      ggplot2::facet_wrap(~Human_gene,ncol=4,scales="free_y")+
      ggplot2::labs(x=NULL,y="log2(normalized count + 1)",
                    title="Expression of eight DFO-associated genes",
                    subtitle="Each point represents one rat; descriptive expression, not protein binding")+
      ggplot2::theme_classic(base_size=10)+
      ggplot2::theme(legend.position="none",axis.text.x=ggplot2::element_text(angle=32,hjust=1))
    save_plot(g1,"Fig_A_Eight_genes_expression",11,6)
  }
  eff8 <- rbind(data.frame(Gene=target_tab$Human_gene,Contrast="STZSal vs Sham",
                           log2FC=target_tab$STZ_log2FC,FDR=target_tab$STZ_FDR),
                data.frame(Gene=target_tab$Human_gene,Contrast="STZDFO vs STZSal",
                           log2FC=target_tab$DFO_log2FC,FDR=target_tab$DFO_FDR))
  eff8$Gene <- factor(eff8$Gene,levels=rev(DFO_HUMAN_GENES))
  eff8$Sig <- ifelse(!is.na(eff8$FDR) & eff8$FDR<FDR_CUTOFF,"*","")
  heat <- ggplot2::ggplot(eff8,ggplot2::aes(Contrast,Gene,fill=log2FC))+
      ggplot2::geom_tile(color="white",linewidth=.8)+
      ggplot2::geom_text(ggplot2::aes(label=Sig),color="black",size=6)+
      ggplot2::scale_fill_gradient2(low="#3A78A7",mid="#F7F7F7",high="#C05A45",
                                    midpoint=0,na.value="#D9DDE1")+
      ggplot2::labs(x=NULL,y=NULL,fill="log2FC",title="Eight drug-associated genes: differential expression",
          subtitle="Asterisk: genome-wide BH-FDR < 0.05; grey: not tested / filtered")+
      ggplot2::theme_minimal(base_size=11)+ggplot2::theme(panel.grid=ggplot2::element_blank())
  save_plot(heat,"Fig_A2_Eight_genes_effect_heatmap",8.5,5.0)

  # ---------------------- 4. DESCRIPTIVE REVERSAL ----------------------------
  logit("[4/7] Disease-associated / DFO expression direction analysis")
  rev <- merge(stz,dfo,by="Rat_gene",suffixes=c("_STZ","_DFO"),sort=FALSE)
  # group-mean VST: descriptive distance to sham, NOT a valid independent test
  grp <- as.character(md$Group)
  means <- sapply(c("Sham","STZSal","STZDFO"), function(g)
     rowMeans(vst_m[,grp==g,drop=FALSE]))
  means <- as.data.frame(means,check.names=FALSE)
  means$Rat_gene <- rownames(means)
  names(means)[1:3] <- c("Mean_VST_Sham","Mean_VST_STZSal","Mean_VST_STZDFO")
  rev <- merge(rev,means,by="Rat_gene",all.x=TRUE,sort=FALSE)
  rev$Opposite_direction <- with(rev,
       is.finite(log2FoldChange_STZ) & is.finite(log2FoldChange_DFO) &
       abs(log2FoldChange_STZ)>=MIN_ABS_LFC_FOR_DIRECTION &
       abs(log2FoldChange_DFO)>=MIN_ABS_LFC_FOR_DIRECTION &
       log2FoldChange_STZ*log2FoldChange_DFO<0)
  rev$Both_nominal_P_below_005 <- with(rev,
       !is.na(pvalue_STZ) & !is.na(pvalue_DFO) &
       pvalue_STZ<.05 & pvalue_DFO<.05)
  rev$Both_genomewide_FDR_below_005 <- with(rev,
       !is.na(padj_STZ) & !is.na(padj_DFO) &
       padj_STZ<FDR_CUTOFF & padj_DFO<FDR_CUTOFF)
  rev$Distance_STZ_from_Sham <- abs(rev$Mean_VST_STZSal-rev$Mean_VST_Sham)
  rev$Distance_DFO_from_Sham <- abs(rev$Mean_VST_STZDFO-rev$Mean_VST_Sham)
  rev$DFO_closer_to_Sham_descriptively <- with(rev,
       Distance_DFO_from_Sham < Distance_STZ_from_Sham)
  rev$Direction_class <- ifelse(rev$Opposite_direction,"Opposite (descriptive)",
             ifelse(is.finite(rev$log2FoldChange_STZ) & is.finite(rev$log2FoldChange_DFO) &
                    abs(rev$log2FoldChange_STZ)>=MIN_ABS_LFC_FOR_DIRECTION &
                    abs(rev$log2FoldChange_DFO)>=MIN_ABS_LFC_FOR_DIRECTION,
                    "Same (descriptive)","Effect too small / unavailable"))
  save_csv(rev,file.path(core_dir,"11_all_genes_model_DFO_direction_and_VST_distance.csv"))
  eight_rev <- merge(target_tab[,c("Human_gene","Rat_gene","Mapping_status")],
       rev[,c("Rat_gene","Direction_class","Both_nominal_P_below_005",
               "Both_genomewide_FDR_below_005","DFO_closer_to_Sham_descriptively",
               "Distance_STZ_from_Sham","Distance_DFO_from_Sham")],
       by="Rat_gene",all.x=TRUE,sort=FALSE)
  save_csv(eight_rev,file.path(target_dir,"04_eight_genes_descriptive_reversal.csv"))
  logit("Direction analysis remains exploratory: contrasts share the STZSal reference group")

  # ---------------------- 5. OPTIONAL PPI NEIGHBORS --------------------------
  logit("[5/7] Mapping network-connected human AD seed genes to rat orthologs")
  ppi_status <- "SKIPPED"
  if (RUN_PPI_ANALYSIS) {
    tryCatch({
      audit_file <- ""
      if (nzchar(trimws(AD37_RUN_OVERRIDE))) {
        audit_file <- file.path(AD37_RUN_OVERRIDE,"14_Deferoxamine_Target_PPI",PPI_AUDIT_NAME)
      } else if (dir.exists(AD37_RESULTS_ROOT)) {
        runs <- list.dirs(AD37_RESULTS_ROOT, recursive=FALSE, full.names=TRUE)
        audit_files <- file.path(runs,"14_Deferoxamine_Target_PPI",PPI_AUDIT_NAME)
        audit_files <- audit_files[file.exists(audit_files)]
        if (length(audit_files))
          audit_file <- audit_files[which.max(file.info(audit_files)$mtime)]
      }
      if (!nzchar(audit_file) || !file.exists(audit_file)) {
        ppi_status <- "NO_PPI_AUDIT_FILE"
        logit("PPI seed-to-target audit not found; PPI stage skipped, core results retained")
      } else {
        audit <- data.table::fread(audit_file,data.table=FALSE)
        needed <- c("Module","Seed_symbol","Nearest_target_symbol","Minimum_hops")
        if (!all(needed %in% names(audit)))
          stop("PPI audit missing columns: ",paste(setdiff(needed,names(audit)),collapse=","))
        audit$Minimum_hops <- as.numeric(audit$Minimum_hops)
        audit$Seed_symbol <- trimws(as.character(audit$Seed_symbol))
        audit <- audit[audit$Module %in% c("AD","Amyloid","Tau") &
           !is.na(audit$Minimum_hops) & nzchar(audit$Seed_symbol),,drop=FALSE]
        if (!nrow(audit)) stop("No valid PPI seed distance rows")
        audit <- audit[order(audit$Seed_symbol,audit$Minimum_hops),,drop=FALSE]
        seed_unique <- audit[!duplicated(audit$Seed_symbol),
                            c("Seed_symbol","Minimum_hops","Nearest_target_symbol"),drop=FALSE]
        flags <- stats::aggregate(Module~Seed_symbol,audit,
                    FUN=function(z) paste(sort(unique(z)),collapse=" | "))
        seed_unique <- merge(seed_unique,flags,by="Seed_symbol",sort=FALSE)
        names(seed_unique)[names(seed_unique)=="Module"] <- "Seed_membership"
        seed_unique$Within_hop_threshold <- seed_unique$Minimum_hops<=MAX_PPI_HOPS
        save_csv(seed_unique,file.path(network_dir,"01_all_AD_seed_network_distance_audit.csv"))
        chosen <- seed_unique[seed_unique$Within_hop_threshold,,drop=FALSE]
        # Predefined set uses network topology only: NO expression-based selection.
        if (!nrow(chosen)) {
          ppi_status <- "NO_SEEDS_WITHIN_HOP_THRESHOLD"
          logit("No AD seed within ",MAX_PPI_HOPS," PPI hops; skip network expression stage")
        } else {
          mapped <- get_orthologs(chosen$Seed_symbol)
          names(chosen)[names(chosen)=="Seed_symbol"] <- "Human_gene"
          m <- merge(chosen,mapped,by="Human_gene",all.x=TRUE,sort=FALSE)
          m$Mapping_status[is.na(m$Mapping_status)] <- "No_mapping"
          m$Expressed_in_DESeq2 <- !is.na(m$Rat_gene) &
                 m$Rat_gene %in% stz$Rat_gene & m$Mapping_status=="One_to_one"
          z <- match(m$Rat_gene,rev$Rat_gene)
          m$STZ_log2FC <- rev$log2FoldChange_STZ[z]
          m$STZ_FDR <- rev$padj_STZ[z]
          m$DFO_log2FC <- rev$log2FoldChange_DFO[z]
          m$DFO_FDR <- rev$padj_DFO[z]
          m$Direction_class <- rev$Direction_class[z]
          m$Closer_to_Sham <- rev$DFO_closer_to_Sham_descriptively[z]
          m$STZ_log2FC[!m$Expressed_in_DESeq2] <- NA_real_
          m$DFO_log2FC[!m$Expressed_in_DESeq2] <- NA_real_
          save_csv(m,file.path(network_dir,"02_nearby_AD_seed_rat_ortholog_expression.csv"))
          ppi_status <- "COMPLETED"
          logit("PPI neighboring AD seed genes (<=",MAX_PPI_HOPS,"): ",nrow(chosen),
                "; one-to-one mapped/tested: ",sum(m$Expressed_in_DESeq2))
          show <- m[m$Expressed_in_DESeq2 & is.finite(m$DFO_log2FC),,drop=FALSE]
          show <- show[order(show$Minimum_hops,show$Human_gene),,drop=FALSE]
          show <- head(show,35L)  # readable subset only; full table preserved
          if (nrow(show)) {
            zz <- rbind(data.frame(Gene=show$Human_gene,Contrast="STZ vs Sham",
                                   log2FC=show$STZ_log2FC,FDR=show$STZ_FDR),
                        data.frame(Gene=show$Human_gene,Contrast="DFO vs STZ",
                                   log2FC=show$DFO_log2FC,FDR=show$DFO_FDR))
            zz$Gene <- factor(zz$Gene,levels=rev(show$Human_gene))
            zz$Star <- ifelse(!is.na(zz$FDR) & zz$FDR<FDR_CUTOFF,"*","")
            p <- ggplot2::ggplot(zz,ggplot2::aes(Contrast,Gene,fill=log2FC))+
                 ggplot2::geom_tile(color="white",linewidth=.5)+
                 ggplot2::geom_text(ggplot2::aes(label=Star),size=3.7)+
                 ggplot2::scale_fill_gradient2(low="#3A78A7",mid="white",high="#C05A45",
                         midpoint=0,na.value="grey88")+
                 ggplot2::labs(title="Predefined AD seed genes near DFO-associated proteins",
                         subtitle=paste0("Shortest PPI distance <=",MAX_PPI_HOPS,
                         "; subset for readability; * genome-wide BH-FDR <0.05"),
                         x=NULL,y=NULL,fill="log2FC")+
                 ggplot2::theme_minimal(base_size=10)+
                 ggplot2::theme(panel.grid=ggplot2::element_blank())
            save_plot(p,"Fig_B_PPI_adjacent_AD_seed_effects",8.4,max(5.5,nrow(show)*.22+2))
          }
        }
      }
    },error=function(e) {
      ppi_status <<- paste0("ERROR: ",conditionMessage(e))
      logit("WARNING: optional PPI stage failed: ",conditionMessage(e))
    })
  }

  # ---------------------- 6. PATHWAY GSEA ------------------------------------
  logit("[6/7] GSEA using full ranked transcriptome and rat MSigDB ortholog sets")
  gsea_status <- "SKIPPED"
  if (RUN_GSEA && gsea_packages_present) {
    tryCatch({
      collections <- list(Reactome=list(collection="C2",subcollection="CP:REACTOME"))
      if (RUN_GO_BP)
        collections$GO_BP <- list(collection="C5",subcollection="GO:BP")
      path_list <- list(); provenance <- list()
      for (nm in names(collections)) {
        cc <- collections[[nm]]
        tmp <- as.data.frame(msigdbr::msigdbr(db_species="HS",species="Rattus norvegicus",
                    collection=cc$collection,subcollection=cc$subcollection))
        if (!all(c("gs_name","gene_symbol") %in% colnames(tmp)))
          stop("msigdbr returned unexpected column names")
        tmp <- tmp[!is.na(tmp$gene_symbol)&nzchar(tmp$gene_symbol),,drop=FALSE]
        pths <- split(tmp$gene_symbol,tmp$gs_name)
        pths <- lapply(pths,unique)
        names(pths) <- paste(nm,names(pths),sep="::")
        path_list <- c(path_list,pths)
        provenance[[nm]] <- data.frame(Collection=nm,
                   Gene_sets=length(pths),
                   MSigDB_version=if ("db_version" %in% names(tmp))
                          paste(unique(tmp$db_version),collapse="|") else "unknown")
      }
      save_csv(do.call(rbind,provenance),file.path(gsea_dir,"01_gene_set_provenance.csv"))
      save_csv(data.frame(Pathway=rep(names(path_list),lengths(path_list)),
           Rat_gene=unlist(path_list,use.names=FALSE)),
           file.path(gsea_dir,"02_gene_set_membership.csv"))
      run_fgsea <- function(tbl,label) {
        stats0 <- tbl$stat
        names(stats0) <- tbl$Rat_gene
        stats0 <- stats0[is.finite(stats0)&!is.na(names(stats0))&nzchar(names(stats0))]
        stats0 <- sort(stats0,decreasing=TRUE)
        if (anyDuplicated(names(stats0))) stop("Duplicate gene symbols in GSEA rank vector")
        if (length(stats0)<1000) stop("Too few ranked genes for GSEA")
        fg <- fgsea::fgseaMultilevel(pathways=path_list,stats=stats0,
                minSize=10,maxSize=500,eps=0, nPermSimple=1000)
        fg <- as.data.frame(fg)
        if ("leadingEdge" %in% names(fg))
          fg$LeadingEdge <- vapply(fg$leadingEdge,function(x) paste(x,collapse=" | "),character(1))
        fg$leadingEdge <- NULL
        fg$Contrast <- label
        fg
      }
      fg_stz <- run_fgsea(stz,"STZSal_vs_Sham")
      fg_dfo <- run_fgsea(dfo,"STZDFO_vs_STZSal")
      fg_all <- rbind(fg_stz,fg_dfo)
      # fgsea padj is calculated within each contrast across ALL tested gene sets.
      save_csv(fg_all,file.path(gsea_dir,"03_ALL_pathways_fgsea_both_contrasts.csv"))
      relevant <- grepl("APP_|AMYLOID|SECRETASE|BACE|IRON|FERR|HIF|HYPOXI|INFLAMMAT|SYNAP|NEURON",
                       fg_all$pathway,ignore.case=TRUE)
      focused <- fg_all[relevant,,drop=FALSE]
      save_csv(focused,file.path(gsea_dir,"04_FOCUSED_APP_iron_HIF_and_neural_GSEA.csv"))
      gsea_status <- "COMPLETED"
      logit("Full GSEA sets tested: ",nrow(fg_all),"; thematic subset: ",nrow(focused))
      if (nrow(focused)) {
        focus_plot <- focused[is.finite(focused$NES) & !is.na(focused$padj),,drop=FALSE]
        focus_plot <- focus_plot[order(focus_plot$padj,-abs(focus_plot$NES)),,drop=FALSE]
        focus_plot <- head(focus_plot,24L)
        if (nrow(focus_plot)) {
          focus_plot$Pathway_label <- sub("^[^:]+::", "", focus_plot$pathway)
          focus_plot$Pathway_label <- gsub("^(REACTOME_|GOBP_)","",focus_plot$Pathway_label)
          focus_plot$Pathway_label <- factor(focus_plot$Pathway_label,
                        levels=rev(unique(focus_plot$Pathway_label)))
          p <- ggplot2::ggplot(focus_plot,ggplot2::aes(x=NES,y=Pathway_label,
                    size=-log10(pmax(padj,1e-300)),color=Contrast))+
               ggplot2::geom_point(alpha=.85)+
               ggplot2::geom_vline(xintercept=0,linetype="dashed",color="grey55")+
               ggplot2::labs(title="Predefined pathway themes in full-set GSEA",
                       subtitle="Shown pathways selected by keywords; adjusted P from all tested sets",
                       y=NULL,x="Normalized enrichment score",size="-log10(FDR)")+
               ggplot2::theme_bw(base_size=10)+ggplot2::theme(legend.position="bottom")
          save_plot(p,"Fig_C_Pathway_GSEA_themes",11,max(6,7+.14*nrow(focus_plot)))
        }
      }
    },error=function(e) {
      gsea_status <<- paste0("ERROR: ",conditionMessage(e))
      logit("WARNING: optional GSEA stage failed: ",conditionMessage(e))
    })
  }

  # ---------------------- 7. MANIFEST ---------------------------------------
  logit("[7/7] Writing provenance, summary and software versions")
  summary <- data.frame(
    Item=c("RNAseq_samples","Sham_n","STZSal_n","STZDFO_n",
        "Gene_symbols_in_count_matrix_after_zeros", "Genes_retained_for_DESeq2",
        "STZ_genomewide_FDR_lt_005","DFO_genomewide_FDR_lt_005",
        "Eight_genes_one_to_one_and_tested","PPI_stage","GSEA_stage",
        "Batch_adjustment","Count_filter_rule","Direction_effect_threshold"),
    Value=c(ncol(dat),3,3,3,nrow(dat),nrow(dat_keep),
         sum(stz$padj<FDR_CUTOFF,na.rm=TRUE),
         sum(dfo$padj<FDR_CUTOFF,na.rm=TRUE),
         sum(target_tab$Included_in_DESeq2),ppi_status,gsea_status,
         as.character(ADJUST_FOR_BATCH),
         paste0("counts>=",COUNT_MIN," in >=",MIN_SAMPLES," samples"),
         MIN_ABS_LFC_FOR_DIRECTION),stringsAsFactors=FALSE)
  save_csv(summary,file.path(out,"00_analysis_summary.csv"))
  installed_used <- c(required,c(optional_cran,optional_bioc)[vapply(c(optional_cran,optional_bioc),
                                requireNamespace,logical(1),quietly=TRUE)])
  package_version <- data.frame(Package=installed_used,
       Version=vapply(installed_used,function(x) as.character(utils::packageVersion(x)),character(1)))
  save_csv(package_version,file.path(out,"package_versions.csv"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  writeLines(c(
    "GSE281803 DFO rat hippocampus RNA-seq reanalysis, v1.0.0",
    paste("Counts:",normalizePath(COUNT_FILE,winslash="/")),
    paste("GEO sample metadata:",normalizePath(META_FILE,winslash="/")),
    "3 experimental groups, 3 rats per group; all 9 columns matched by core sample ID.",
    paste("Design:",deparse(design_formula)),
    paste("Filter: counts >=",COUNT_MIN,"in >=",MIN_SAMPLES,"samples"),
    "DESeq2 comparisons: STZSal/Sham and STZDFO/STZSal; padj over whole gene set.",
    "Eight DFO-associated genes defined from the original drug-protein network.",
    "Babelgene ortholog mapping: min_support = 3; only unambiguous one-to-one pairs evaluated.",
    "Panel-of-eight BH is exploratory and reported separately from genome-wide BH-FDR.",
    "Shared STZSal group induces statistical dependence between the two contrasts.",
    "Opposite direction and closer-to-Sham measures are descriptive, NOT efficacy tests.",
    "Network-connected AD genes come from original PPI shortest-distance audit, if available.",
    "Reactome/GO GSEA uses all testable pathways; theme filter applied only after global FDR.",
    "Transcript changes neither establish direct drug-protein binding nor APP processing at protein level.",
    "Published dataset reanalysis is not independent biological validation.",
    "Citation: Fine et al. Front Neurosci. DOI 10.3389/fnins.2024.1528374",
    "Official methods: https://bioconductor.org/packages/DESeq2/",
    "Orthologs: https://cran.r-project.org/package=babelgene",
    "MSigDB: https://cran.r-project.org/package=msigdbr",
    "GSEA: https://bioconductor.org/packages/fgsea/"
  ),file.path(out,"analysis_methods_and_cautions.txt"))
  logit("RUN COMPLETED; results saved at: ",out)
  print(summary,row.names=FALSE)
  invisible(out)
}

GSE281803_RESULT_DIR <- run_GSE281803()
