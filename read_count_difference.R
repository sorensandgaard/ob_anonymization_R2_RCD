#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("EMDomics")

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]


# setwd("/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss")
# case_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds"
# ctrl_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds"
# outdir <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m1/default/D1.somefile.txt"

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Do something - write earth movers distance

# Normalize
case_obj <- NormalizeData(case_obj)
ctrl_obj <- NormalizeData(ctrl_obj)

# Convert Seurat objects to expression matrices
expr_orig <- case_obj[["RNA"]]$data
expr_anon <- case_obj[["RNA"]]$data

# Check genelists are identical and save genelist and genecount
genelist_orig <- rownames(expr_orig)
genelist_anon <- rownames(expr_anon)

try(if(!identical(genelist_orig,genelist_anon)) stop("Genelists are not identical"))
genelist <- genelist_orig
rm(genelist_orig,genelist_anon)

# Keep only genes that are non-zero in at least one matrix
non_zero_anon <- rowSums(expr_anon > 0) > 0
non_zero_orig <- rowSums(expr_orig > 0) > 0
common_genes <- non_zero_anon & non_zero_orig

expr_anon <- as.matrix(expr_anon[common_genes, ])
expr_orig <- as.matrix(expr_orig[common_genes, ])

# Combine the expression matrices
combined_matrix <- cbind(expr_anon, expr_orig)
colnames(combined_matrix) <- paste("sample", 1:ncol(combined_matrix), sep="")

# Create labels vector
labels <- c(rep("Group1", ncol(expr_anon)), rep("Group2", ncol(expr_orig)))
names(labels) <- colnames(combined_matrix)

# Calculate EMD
emd_results <- calculate_emd(combined_matrix, labels,nperm = 2,parallel = F)

emd_results <- emd_results$emd

write.csv(emd_results,outdir)
