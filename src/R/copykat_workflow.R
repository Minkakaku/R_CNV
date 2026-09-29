# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# CopyKAT 独立分析入口
# 输入：主 Seurat 流程产生的注释后 RDS；脚本从 RNA assay 提取 raw UMI counts。
# 多样本对象按样本分别运行，避免官方文档指出的跨样本批次效应。

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：把逗号分隔的正常参考细胞类型转换为字符向量。
split_values <- function(x) {
  if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# 功能：同时保存矢量 PDF 和 600 dpi PNG。
save_plot <- function(plot, stem, width = 8, height = 7) {
  plot <- plot + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
}

# 功能：提取 Seurat v5 合并后的原始 counts layer。
seurat_counts <- function(obj, assay) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay)
  SeuratObject::LayerData(obj, assay = assay, layer = "counts")
}

# 功能：对一个生物学样本运行 CopyKAT 并标准化 prediction 表。
run_one_sample <- function(
  obj,
  sample_name,
  output_dir,
  assay,
  annotation_col,
  reference_types,
  genome,
  id_type,
  ngene_chr,
  win_size,
  ks_cut,
  distance,
  n_cores
) {
  sample_dir <- file.path(output_dir, make.names(sample_name))
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  counts <- seurat_counts(obj, assay)
  counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]

  # 指定已知正常细胞可稳定基线；不指定时传空字符串让 CopyKAT 自动寻找 diploid cells。
  normal_cells <- ""
  if (length(reference_types)) {
    metadata <- obj[[]]
    normal_cells <- rownames(metadata)[metadata[[annotation_col]] %in% reference_types]
    if (!length(normal_cells)) stop(sample_name, " 中没有配置的正常 reference cells。")
  }

  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(sample_dir)
  set.seed(1001)
  result <- copykat::copykat(
    rawmat = as.matrix(counts),
    id.type = id_type,
    ngene.chr = ngene_chr,
    win.size = win_size,
    KS.cut = ks_cut,
    sam.name = make.names(sample_name),
    distance = distance,
    norm.cell.names = normal_cells,
    n.cores = n_cores,
    genome = genome
  )

  prediction <- as.data.frame(result$prediction, stringsAsFactors = FALSE)
  cell_column <- intersect(c("cell.names", "cell", "barcode"), colnames(prediction))[1]
  if (is.na(cell_column)) {
    prediction$cell <- rownames(prediction)
    cell_column <- "cell"
  }
  pred_column <- intersect(c("copykat.pred", "prediction"), colnames(prediction))[1]
  if (is.na(pred_column)) stop("无法识别 CopyKAT prediction 列。")
  prediction$cell <- prediction[[cell_column]]
  prediction$sample <- sample_name
  prediction$copykat_prediction <- prediction[[pred_column]]
  write.csv(prediction, file.path(sample_dir, "copykat_predictions.csv"), row.names = FALSE)
  list(result = result, prediction = prediction)
}

# 功能：按样本运行 CopyKAT、合并标签并保存带预测的 Seurat 对象和 UMAP。
run_copykat <- function(
  input_rds,
  output_dir,
  assay = "RNA",
  sample_col = "orig.ident",
  annotation_col = "cell_type",
  reference_types = character(0),
  genome = "hg20",
  id_type = "S",
  ngene_chr = 5L,
  win_size = 25L,
  ks_cut = 0.1,
  distance = "euclidean",
  n_cores = 8L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(input_rds)
  metadata <- obj[[]]
  sample_ids <- if (sample_col %in% colnames(metadata)) as.character(metadata[[sample_col]]) else rep("sample", nrow(metadata))
  names(sample_ids) <- rownames(metadata)

  # CopyKAT 官方建议每个样本单独计算；合并样本可能把批次差异误判成 CNV。
  runs <- lapply(unique(sample_ids), function(sample_name) {
    cells <- names(sample_ids)[sample_ids == sample_name]
    run_one_sample(
      obj = subset(obj, cells = cells),
      sample_name = sample_name,
      output_dir = output_dir,
      assay = assay,
      annotation_col = annotation_col,
      reference_types = reference_types,
      genome = genome,
      id_type = id_type,
      ngene_chr = ngene_chr,
      win_size = win_size,
      ks_cut = ks_cut,
      distance = distance,
      n_cores = n_cores
    )
  })
  prediction <- do.call(rbind, lapply(runs, `[[`, "prediction"))
  write.csv(prediction, file.path(output_dir, "copykat_predictions_all_samples.csv"), row.names = FALSE)

  obj$copykat_prediction <- prediction$copykat_prediction[match(colnames(obj), prediction$cell)]
  saveRDS(obj, file.path(object_dir, "copykat_annotated.rds"))

  if ("umap" %in% Seurat::Reductions(obj)) {
    plot <- Seurat::DimPlot(
      obj,
      reduction = "umap",
      group.by = "copykat_prediction",
      raster = FALSE,
      cols = c(aneuploid = "#E64B35", diploid = "#4DBBD5", "not.defined" = "#BDBDBD")
    ) + ggplot2::coord_fixed()
    save_plot(plot, file.path(output_dir, "UMAP_copykat_prediction"))
  }

  writeLines(
    c(
      "# CopyKAT 结果说明",
      "",
      paste0("- 输入 RDS：", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
      paste0("- 样本列：", sample_col, "；每个样本分别运行后再合并标签。"),
      paste0("- genome：", genome, "；gene ID type：", id_type, "。"),
      paste0("- 正常 reference：", if (length(reference_types)) paste(reference_types, collapse = ", ") else "算法自动寻找", "。"),
      paste0("- ngene.chr = ", ngene_chr, "；win.size = ", win_size, "；KS.cut = ", ks_cut, "。"),
      "- 输入必须是 raw UMI counts；不要使用 log-normalized 或 scale.data。",
      "- aneuploid/diploid 是 RNA 推断标签；小鼠预测尤其需要目视检查。",
      "- 结论应结合正常细胞、肿瘤 marker 和 DNA CNV 证据，不能只凭 CopyKAT 定义恶性细胞。"
    ),
    file.path(output_dir, "README_results_CN.md")
  )
  invisible(list(object = obj, runs = runs, prediction = prediction))
}

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
object_dir <- arg_value("--object-dir", file.path("data/processed", "copykat"))
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
if (is.null(input_rds) || is.null(output_dir)) {
  stop("用法：Rscript copykat_workflow.R --input annotated.rds --output results/copykat --sample-col orig.ident")
}

run_copykat(
  input_rds = input_rds,
  output_dir = output_dir,
  assay = arg_value("--assay", "RNA"),
  sample_col = arg_value("--sample-col", "orig.ident"),
  annotation_col = arg_value("--annotation-col", "cell_type"),
  reference_types = split_values(arg_value("--references", "")),
  genome = arg_value("--genome", "hg20"),
  id_type = arg_value("--id-type", "S"),
  ngene_chr = as.integer(arg_value("--ngene-chr", "5")),
  win_size = as.integer(arg_value("--win-size", "25")),
  ks_cut = as.numeric(arg_value("--ks-cut", "0.1")),
  distance = arg_value("--distance", "euclidean"),
  n_cores = as.integer(arg_value("--cores", "8"))
)

}
if (sys.nframe() == 0L) main()
