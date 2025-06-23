#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")

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
out_case_count <- length(colnames(case_obj))
out_ctrl_count <- length(colnames(ctrl_obj))

print("Getting common cells")
common_cells <- intersect(colnames(case_obj),colnames(ctrl_obj))
out_comm_count <- length(common_cells)

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
row_MAE <- data.frame(MAE = rowMeans(expr_absdiff)) %>% rownames_to_column(var = "gene")

print("Creating out dataframe")
out_df <- data.frame(
  names = c("Control cell count","Anonymous cell count","Common cell count","Mean absolute error"),
  values = c(out_ctrl_count,out_case_count,out_comm_count,out_MAE)
)

print("Writing outfile")
write.table(out_df,file = outdir)



