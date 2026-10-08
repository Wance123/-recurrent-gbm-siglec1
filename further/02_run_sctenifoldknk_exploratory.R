suppressPackageStartupMessages({library(Matrix); library(scTenifoldKnk)})
root <- file.path(Sys.getenv('JOB34_ROOT',unset='.'),'further')
input_dir <- file.path(root, "input")
result_dir <- file.path(root, "results", "sctenifoldknk")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

read_condition <- function(condition) {
  prefix <- tolower(condition)
  x <- readMM(file.path(input_dir, paste0(prefix, "_genes_by_cells.mtx")))
  rownames(x) <- readLines(file.path(input_dir, paste0(prefix, "_genes.tsv")), encoding = "UTF-8")
  colnames(x) <- paste0(prefix, "_cell_", seq_len(ncol(x)))
  x
}

targets <- c("SIGLEC1", "STAT1", "IRF7", "IFNAR1", "FCGR1A", "CCL3")
rows <- list()
set.seed(34)
for (condition in c("Primary", "Recurrent")) {
  counts <- read_condition(condition)
  for (target in targets) {
    if (!target %in% rownames(counts)) next
    fit <- scTenifoldKnk(countMatrix = counts, qc = FALSE, gKO = target,
      nc_nNet = 5, nc_nCells = min(300, ncol(counts)), nc_nComp = 3,
      nc_q = 0.9, td_K = 3, td_maxIter = 200, ma_nDim = 2, nCores = 4)
    saveRDS(fit, file.path(result_dir, paste0(condition, "_", target, ".rds")))
    dr <- as.data.frame(fit$diffRegulation)
    dr <- dr[order(dr$p.adj, -dr$distance), ]
    write.csv(dr, file.path(result_dir, paste0(condition, "_", target, "_diffRegulation.csv")), row.names = FALSE)
    rows[[length(rows) + 1]] <- data.frame(condition = condition, target = target,
      cells = ncol(counts), genes = nrow(counts),
      downstream_FDR_0_05 = sum(dr$p.adj < 0.05 & dr$gene != target, na.rm = TRUE))
  }
}
write.csv(do.call(rbind, rows), file.path(result_dir, "run_summary.csv"), row.names = FALSE)
