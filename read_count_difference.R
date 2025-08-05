#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("jsonlite")

rm(list=ls())

print("Loading arguments")
args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]


# setwd("/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss")
# case_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds"
# ctrl_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds"
# outdir <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m2/default/D1.somefile.txt"

print("Loading seurat objects")
case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Read counts
print("Converting to expr matrices")
expr_orig <- ctrl_obj[["RNA"]]$counts
expr_anon <- case_obj[["RNA"]]$counts

print("Getting readcounts/rowsums")
rs_orig <- data.frame(rs_orig = rowSums(expr_orig)) %>% rownames_to_column(var = "gene")
rs_anon <- data.frame(rs_anon = rowSums(expr_anon)) %>% rownames_to_column(var = "gene")

print("Joining readcount matrices")
out_rs_both <- inner_join(rs_orig,rs_anon,"gene") %>% 
  mutate(rcd = rs_orig - rs_anon)
rm(rs_orig,rs_anon)

## Cellcounts
print("Getting cellcounts")
out_case_count <- ncol(case_obj)
out_ctrl_count <- ncol(ctrl_obj)

print("Getting common cells")
common_cells <- intersect(colnames(case_obj),colnames(ctrl_obj))
out_comm_count <- length(common_cells)
unique_in_ctrl <- colnames(ctrl_obj)[!colnames(ctrl_obj) %in% common_cells]
unique_in_anon <- colnames(case_obj)[!colnames(case_obj) %in% common_cells]

## Common cells
print("Subset to common cells")
ctrl_obj <- subset(ctrl_obj,cells=common_cells)
case_obj <- subset(case_obj,cells=common_cells)

# Convert Seurat objects to expression matrices
print("Converting subset seurat objects to expr matrix")
expr_orig <- ctrl_obj[["RNA"]]$counts
expr_anon <- case_obj[["RNA"]]$counts

## Mean absolute error
print("Calculating MAE")
expr_absdiff <- abs(expr_orig - expr_anon)
out_MAE <- mean(expr_absdiff)
row_MAE <- data.frame(MAE = rowMeans(expr_absdiff)) %>% 
  rownames_to_column(var = "gene")
row_MAE_1pct <- quantile(row_MAE$MAE,probs = 0.99)
row_MAE <-  row_MAE %>% 
  arrange(desc(MAE)) %>% 
  slice_head(n = 20)

cell_MAE <- data.frame(MAE = colMeans(expr_absdiff)) %>% 
  rownames_to_column(var = "cell") 
cell_MAE_1pct <- quantile(cell_MAE$MAE,probs = 0.99)
cell_MAE <- cell_MAE %>% 
  arrange(desc(MAE)) %>% 
  slice_head(n = 20)

# Format output
json_obj <- list(
  module = "R2_ReadCountDifference",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
  metrics = list(
    total_cells = list(
      ctrl = out_ctrl_count,
      anon = out_case_count,
      diff = out_case_count-out_ctrl_count,
      in_common = out_comm_count,
      unique_in_ctrl = unique_in_ctrl,
      unique_in_anon = unique_in_anon
    ),
    mean_absolute_error = list(
      overall = out_MAE,
      gene_99th_percentile = row_MAE_1pct,
      top_20_gene_values = row_MAE$MAE,
      top_20_gene_names = row_MAE$gene,
      cell_99th_percentile = cell_MAE_1pct,
      top_20_cell_values = cell_MAE$MAE,
      top_20_cell_names = cell_MAE$cell
    ),
    count_differences_genewise = list(
      mean = mean(out_rs_both$rcd),
      median = median(out_rs_both$rcd),
      min = min(out_rs_both$rcd),
      max = max(out_rs_both$rcd),
      q25 = quantile(out_rs_both$rcd,probs = 0.25),
      q75 = quantile(out_rs_both$rcd,probs = 0.75)
    )
  )
)

# Save output
write_json(json_obj, outdir, pretty = TRUE, auto_unbox = TRUE)






