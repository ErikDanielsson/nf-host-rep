#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)
interaction_csv_fn <- args[1]
interaction_nex_fn <- args[2]

interactions <- read.csv(interaction_csv_fn, sep = "\t", row.names = 1)

write.nexus.data(as.matrix(interactions), interaction_nex_fn, format = "standard")
