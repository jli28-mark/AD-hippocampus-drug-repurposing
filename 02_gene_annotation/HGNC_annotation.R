
base_dir <- "D:/My Desktop/drug/1.DEG"

# datasets <- c("GSE278723", "GSE173955")
datasets <- c("GSE278723", "GSE173955", "GSE48350", "GSE5281")

annotation_dir <- file.path(base_dir, "shared_HGNC_annotation")
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

out_dir <- file.path(
  base_dir, "protein_coding_results",
  paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. 下载一次，之后始终复用同一份 HGNC 注释 ----------

hgnc_url <- paste0(
  "https://storage.googleapis.com/",
  "public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)

hgnc_file <- file.path(annotation_dir, "hgnc_complete_set.txt")

if (!file.exists(hgnc_file)) {
  message("首次运行：下载 HGNC 注释……")
  options(timeout = max(300, getOption("timeout")))

  tmp_file <- paste0(hgnc_file, ".download")

  status <- download.file(
    hgnc_url, tmp_file,
    method = "libcurl", mode = "wb"
  )

  if (status != 0L || !file.rename(tmp_file, hgnc_file)) {
    stop("HGNC 注释下载失败，请检查网络后重试。")
  }
}

hgnc <- read.delim(
  hgnc_file,
  colClasses = "character",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = "",
  fileEncoding = "UTF-8"
)

required <- c(
  "hgnc_id", "symbol", "status", "locus_group", "locus_type",
  "prev_symbol", "alias_symbol", "entrez_id", "ensembl_gene_id"
)

if (!all(required %in% names(hgnc))) {
  stop("HGNC 文件缺少必要字段，请核查注释文件。")
}

hgnc <- hgnc[
  !is.na(hgnc$status) & hgnc$status == "Approved",
  , drop = FALSE
]

stopifnot(
  !anyDuplicated(hgnc$symbol),
  !anyDuplicated(hgnc$hgnc_id)
)


writeLines(
  c(
    paste("Annotation source:", hgnc_url),
    paste("Local annotation:", hgnc_file),
    paste("MD5:", unname(tools::md5sum(hgnc_file))),
    paste("Analysis time:", format(Sys.time(), tz = "UTC")),
    "Selection: Approved HGNC records; locus_group = protein-coding gene",
    "Existing differential-expression statistics were not recalculated."
  ),
  file.path(out_dir, "annotation_record.txt")
)

# ---------- 2. 建立旧名称、别名与标准名称的对应关系 ----------

make_alias_table <- function(field) {
  x <- hgnc[[field]]
  x[is.na(x)] <- ""

  parts <- strsplit(x, "|", fixed = TRUE)

  data.frame(
    key = unlist(parts, use.names = FALSE),
    index = rep(seq_len(nrow(hgnc)), lengths(parts)),
    stringsAsFactors = FALSE
  )
}

alias_map <- unique(rbind(
  make_alias_table("prev_symbol"),
  make_alias_table("alias_symbol")
))

alias_map$key <- trimws(alias_map$key)
alias_map <- unique(
  alias_map[nzchar(alias_map$key), , drop = FALSE]
)


ambiguous_keys <- unique(alias_map$key[
  duplicated(alias_map$key) |
    duplicated(alias_map$key, fromLast = TRUE)
])

alias_unique <- alias_map[
  !(alias_map$key %in% ambiguous_keys),
  , drop = FALSE
]

# ---------- 3. 注释和筛选函数 ----------

extract_protein_coding <- function(dat, dataset, gene_col = NULL) {

 
  if (is.null(gene_col)) {
    clean_names <- tolower(
      gsub("[^[:alnum:]]", "", names(dat))
    )

    priority <- c("genesymbol", "symbol", "gene")
    found <- priority[priority %in% clean_names]

    if (!length(found)) {
      stop(
        dataset, "：没有识别到基因符号列。\n",
        "请将确认包含基因符号的列改名为 Gene 后重新运行。\n",
        "当前列名：", paste(names(dat), collapse = ", ")
      )
    }

    gene_col <- names(dat)[match(found[1], clean_names)]
  }

  stopifnot(gene_col %in% names(dat))

  message(dataset, "：使用基因符号列 ", gene_col)

  symbol_input <- trimws(as.character(dat[[gene_col]]))
  symbol_input[is.na(symbol_input)] <- ""

  missing_symbol <- symbol_input %in% c("", "---", "NA")

  # 多基因注释单独保留，避免一个探针被随意分配给某个基因
  multiple_symbols <- grepl("[/;,|]", symbol_input)

  valid <- !missing_symbol & !multiple_symbols

  # 优先匹配当前标准名称
  index <- match(symbol_input, hgnc$symbol)
  index[!valid] <- NA_integer_

  map_status <- rep("unmapped", nrow(dat))
  map_status[missing_symbol] <- "missing_symbol"
  map_status[multiple_symbols & !missing_symbol] <- "multiple_symbols"
  map_status[!is.na(index)] <- "approved_symbol"

  # 当前名称未匹配时，再匹配唯一的旧名称或别名
  unresolved <- valid & is.na(index)

  j <- match(symbol_input[unresolved], alias_unique$key)
  index[unresolved] <- alias_unique$index[j]

  map_status[unresolved] <- ifelse(
    !is.na(index[unresolved]),
    "previous_or_alias_symbol",
    ifelse(
      symbol_input[unresolved] %in% ambiguous_keys,
      "ambiguous_alias",
      "unmapped"
    )
  )

  if (all(is.na(index))) {
    stop(dataset, "：没有任何基因匹配，请核查选择的基因符号列。")
  }

  annotation <- data.frame(
    HGNC_input_symbol = symbol_input,
    HGNC_mapping_status = map_status,
    HGNC_ID = hgnc$hgnc_id[index],
    HGNC_symbol = hgnc$symbol[index],
    HGNC_locus_group = hgnc$locus_group[index],
    HGNC_locus_type = hgnc$locus_type[index],
    HGNC_Entrez_ID = hgnc$entrez_id[index],
    HGNC_Ensembl_ID = hgnc$ensembl_gene_id[index],
    stringsAsFactors = FALSE
  )

  # 原始结果列全部保留，注释追加在后面
  annotated <- cbind(dat, annotation)

  keep <- !is.na(index) &
    annotation$HGNC_locus_group %in% "protein-coding gene"

  # 重复基因暂不合并，不按最小 P/FDR 随意选取记录
  repeated <- keep & (
    duplicated(annotation$HGNC_ID) |
      duplicated(annotation$HGNC_ID, fromLast = TRUE)
  )

  dataset_dir <- file.path(out_dir, dataset)
  dir.create(dataset_dir, showWarnings = FALSE)

  save_csv <- function(x, filename) {
    write.csv(
      x, file.path(dataset_dir, filename),
      row.names = FALSE,
      na = "",
      fileEncoding = "UTF-8"
    )
  }

  save_csv(annotated, "all_records_annotated.csv")

  save_csv(
    annotated[keep, , drop = FALSE],
    "protein_coding_all_records.csv"
  )

  save_csv(
    annotated[!keep, , drop = FALSE],
    "excluded_or_unresolved_records.csv"
  )

  save_csv(
    annotated[repeated, , drop = FALSE],
    "repeated_protein_coding_records.csv"
  )

  data.frame(
    dataset = dataset,
    input_records = nrow(dat),
    matched_records = sum(!is.na(index)),
    protein_coding_records = sum(keep),
    unique_protein_coding_genes =
      length(unique(annotation$HGNC_ID[keep])),
    mapped_non_protein_coding_records =
      sum(!is.na(index) & !keep),
    unresolved_records = sum(is.na(index)),
    repeated_protein_coding_genes =
      length(unique(annotation$HGNC_ID[repeated]))
  )
}

# ---------- 4. 依次选择各数据集的完整结果表 ----------

summaries <- vector("list", length(datasets))
input_files <- character(length(datasets))

for (k in seq_along(datasets)) {

  ds <- datasets[k]

  message(
    "\n请选择 ", ds,
    " 的完整差异分析结果：CSV、TSV 或 TXT"
  )

  input_file <- file.choose()
  input_files[k] <- input_file

  # 所有原始列按字符读取，避免在注释步骤改写统计值
  if (grepl("\\.csv$", input_file, ignore.case = TRUE)) {

    dat <- read.csv(
      input_file,
      colClasses = "character",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM"
    )

  } else if (grepl("\\.(txt|tsv)$", input_file, ignore.case = TRUE)) {

    dat <- read.delim(
      input_file,
      colClasses = "character",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "\"",
      comment.char = "",
      fileEncoding = "UTF-8-BOM"
    )

  } else {
    stop(
      "请选择完整结果的 CSV/TSV/TXT 文件。",
      "GSE173955 请使用 author_statistics_extracted.csv。"
    )
  }

  summaries[[k]] <- extract_protein_coding(dat, ds)
}

summary_table <- do.call(rbind, summaries)

write.csv(
  summary_table,
  file.path(out_dir, "four_datasets_annotation_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(dataset = datasets, input_file = input_files),
  file.path(out_dir, "input_files.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

print(summary_table, row.names = FALSE)
message("\n完成。结果保存在：", out_dir)