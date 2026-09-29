# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# inferCNV 独立分析入口
# 输入：主 Seurat 流程产生的注释后 RDS，以及与物种/基因标识匹配的 gene_order.tsv。
# 输出：inferCNV 官方结果、输入审计表、对象和中文说明。

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：把逗号分隔的正常参考细胞类型转换为字符向量。
split_values <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# 功能：从 Seurat v5 对象提取合并后的原始 counts layer。
seurat_counts <- function(obj, assay) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay)
  SeuratObject::LayerData(obj, assay = assay, layer = "counts")
}

# 功能：构建 inferCNV 输入并运行官方 infercnv::run。
run_infercnv <- function(
  input_rds,
  gene_order_file,
  output_dir,
  assay,
  annotation_col,
  reference_types,
  cutoff = 0.1,
  denoise = TRUE,
  HMM = FALSE,
  cluster_by_groups = FALSE,
  num_threads = 8L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(input_rds)
  counts <- seurat_counts(obj, assay)
  metadata <- obj[[]][colnames(counts), , drop = FALSE]

  annotations <- data.frame(
    cell = colnames(counts),
    cell_type = as.character(metadata[[annotation_col]]),
    stringsAsFactors = FALSE
  )
  annotations <- annotations[!is.na(annotations$cell_type) & annotations$cell_type != "", ]
  counts <- counts[, annotations$cell, drop = FALSE]

  missing_reference <- setdiff(reference_types, unique(annotations$cell_type))
  if (length(missing_reference)) {
    stop("以下 reference 不在注释列中：", paste(missing_reference, collapse = ", "))
  }

  annotation_file <- file.path(output_dir, "infercnv_cell_annotations.tsv")
  write.table(
    annotations,
    annotation_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  cell_counts <- as.data.frame(table(annotations$cell_type), stringsAsFactors = FALSE)
  colnames(cell_counts) <- c("cell_type", "n_cells")
  write.csv(cell_counts, file.path(output_dir, "infercnv_input_cell_counts.csv"), row.names = FALSE)

  # gene_order.tsv 通常为 gene、chromosome、start、end 四列且无表头；必须与 counts 行名一致。
  infer_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = reference_types
  )

  # cutoff 常用 0.1（10x）或 1（Smart-seq2）；HMM 是否启用应按研究目的显式决定。
  infer_obj <- infercnv::run(
    infer_obj,
    cutoff = cutoff,
    out_dir = output_dir,
    cluster_by_groups = cluster_by_groups,
    denoise = denoise,
    HMM = HMM,
    num_threads = num_threads,
    plot_steps = FALSE
  )
  saveRDS(infer_obj, file.path(object_dir, "infercnv_object.rds"))

  writeLines(
    c(
      "# inferCNV 结果说明",
      "",
      paste0("- 输入 RDS：", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
      paste0("- assay：", assay, "；注释列：", annotation_col, "。"),
      paste0("- 正常 reference：", paste(reference_types, collapse = ", "), "。"),
      paste0("- cutoff = ", cutoff, "；denoise = ", denoise, "；HMM = ", HMM, "。"),
      "- 输入使用原始 counts；gene_order.tsv 必须与物种、基因 ID 和基因组版本一致。",
      "- 正常 reference 的身份是分析假设，应结合 marker 和样本来源检查。",
      "- inferCNV 是基于转录表达的相对 CNV 推断，不能替代 DNA 层面的拷贝数检测。"
    ),
    file.path(output_dir, "README_results_CN.md")
  )
  invisible(infer_obj)
}

input_rds <- arg_value("--input")
gene_order_file <- arg_value("--gene-order")
output_dir <- arg_value("--output")
object_dir <- arg_value("--object-dir", file.path("data/processed", "infercnv"))
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
if (is.null(input_rds) || is.null(gene_order_file) || is.null(output_dir)) {
  stop("用法：Rscript infercnv_workflow.R --input annotated.rds --gene-order gene_order.tsv --output results/infercnv --references 'T cell,B cell'")
}

run_infercnv(
  input_rds = input_rds,
  gene_order_file = gene_order_file,
  output_dir = output_dir,
  assay = arg_value("--assay", "RNA"),
  annotation_col = arg_value("--annotation-col", "cell_type"),
  reference_types = split_values(arg_value("--references", "T cell,B cell,Endothelial cell")),
  cutoff = as.numeric(arg_value("--cutoff", "0.1")),
  denoise = as.logical(arg_value("--denoise", "TRUE")),
  HMM = as.logical(arg_value("--HMM", "FALSE")),
  cluster_by_groups = as.logical(arg_value("--cluster-by-groups", "FALSE")),
  num_threads = as.integer(arg_value("--threads", "8"))
)

}
if (sys.nframe() == 0L) main()
