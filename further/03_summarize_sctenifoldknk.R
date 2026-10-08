result_dir <- file.path(file.path(Sys.getenv('JOB34_ROOT',unset='.'),'further/results/sctenifoldknk'))

rds_files <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE)
summaries <- list()
gene_sets <- list()

for (rds_file in rds_files) {
  stem <- sub("\\.rds$", "", basename(rds_file))
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  condition <- parts[1]
  target <- paste(parts[-1], collapse = "_")
  fit <- readRDS(rds_file)
  dr <- as.data.frame(fit$diffRegulation)
  stopifnot(all(c("gene", "distance", "Z", "FC", "p.value", "p.adj") %in% names(dr)))
  dr <- dr[order(dr$p.adj, -dr$distance), ]
  write.csv(dr, file.path(result_dir, paste0(stem, "_diffRegulation.csv")), row.names = FALSE)

  downstream <- dr$gene[dr$p.adj < 0.05 & dr$gene != target]
  gene_sets[[paste(condition, target, sep = "_")]] <- downstream
  summaries[[length(summaries) + 1]] <- data.frame(
    condition = condition,
    target = target,
    tested_genes = nrow(dr),
    downstream_FDR_0_05 = length(downstream),
    min_downstream_FDR = if (any(dr$gene != target)) min(dr$p.adj[dr$gene != target], na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summaries)
summary_df <- summary_df[order(summary_df$target, summary_df$condition), ]
write.csv(summary_df, file.path(result_dir, "run_summary_corrected.csv"), row.names = FALSE)

targets <- sort(unique(summary_df$target))
comparisons <- lapply(targets, function(target) {
  p <- gene_sets[[paste("Primary", target, sep = "_")]]
  r <- gene_sets[[paste("Recurrent", target, sep = "_")]]
  u <- union(p, r)
  data.frame(
    target = target,
    primary_n = length(p),
    recurrent_n = length(r),
    overlap_n = length(intersect(p, r)),
    jaccard = if (length(u)) length(intersect(p, r)) / length(u) else NA_real_,
    primary_specific_n = length(setdiff(p, r)),
    recurrent_specific_n = length(setdiff(r, p)),
    recurrent_specific_genes = paste(setdiff(r, p), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
comparison_df <- do.call(rbind, comparisons)
write.csv(comparison_df, file.path(result_dir, "primary_recurrent_comparison.csv"), row.names = FALSE)

top_rows <- do.call(rbind, lapply(rds_files, function(rds_file) {
  stem <- sub("\\.rds$", "", basename(rds_file))
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  condition <- parts[1]
  target <- paste(parts[-1], collapse = "_")
  dr <- as.data.frame(readRDS(rds_file)$diffRegulation)
  dr <- dr[dr$gene != target, ]
  dr <- dr[order(dr$p.adj, -dr$distance), ]
  dr <- head(dr, 20)
  cbind(condition = condition, target = target, dr, stringsAsFactors = FALSE)
}))
write.csv(top_rows, file.path(result_dir, "top20_downstream_by_run.csv"), row.names = FALSE)

print(summary_df)
print(comparison_df[, setdiff(names(comparison_df), "recurrent_specific_genes")])
