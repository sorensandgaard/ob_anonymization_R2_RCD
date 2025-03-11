#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("aricode")

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]


# setwd("/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss")
# case_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds"
# ctrl_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds"
# outdir <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m2/default/D1.somefile.txt"

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Read counts
expr_orig <- ctrl_obj[["RNA"]]$counts
expr_anon <- case_obj[["RNA"]]$counts
rs_orig <- data.frame(rs_orig = rowSums(expr_orig)) %>% rownames_to_column(var = "gene")
rs_anon <- data.frame(rs_anon = rowSums(expr_anon)) %>% rownames_to_column(var = "gene")
out_rs_both <- inner_join(rs_orig,rs_anon,"gene") %>% 
  mutate(rcd = rs_orig - rs_anon)
rm(rs_orig,rs_anon)

## Cellcounts
out_case_count <- length(colnames(case_obj))
out_ctrl_count <- length(colnames(ctrl_obj))
common_cells <- intersect(colnames(case_obj),colnames(ctrl_obj))
out_comm_count <- length(common_cells)

## Common cells
ctrl_obj <- subset(ctrl_obj,cells=common_cells)
case_obj <- subset(case_obj,cells=common_cells)

# Convert Seurat objects to expression matrices
expr_orig <- ctrl_obj[["RNA"]]$counts
expr_anon <- case_obj[["RNA"]]$counts

## Mean absolute error
expr_absdiff <- abs(expr_orig - expr_anon)
out_MAE <- mean(expr_absdiff)
row_MAE <- data.frame(MAE = rowMeans(expr_absdiff)) %>% rownames_to_column(var = "gene")

## Clustering with Aricode
case_obj <- NormalizeData(case_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)
ctrl_obj <- NormalizeData(ctrl_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.5)


ctrl_clusters <- FetchData(case_obj,vars = c("seurat_clusters")) %>% 
  rownames_to_column(var = "cell")
case_clusters <- FetchData(ctrl_obj,vars = c("seurat_clusters")) %>% 
  rownames_to_column(var = "cell")
cluster_mat <- full_join(x = ctrl_clusters,
                         y = case_clusters,
                         by = "cell")

out_ARI <- ARI(cluster_mat$seurat_clusters.x,cluster_mat$seurat_clusters.y)

out_df <- data.frame(
  names = c("Control cell count","Anonymous cell count","Common cell count","Mean absolute error","Adjust rand index"),
  values = c(out_ctrl_count,out_case_count,out_comm_count,out_MAE,out_ARI)
)

write.table(out_df,file = outdir)



