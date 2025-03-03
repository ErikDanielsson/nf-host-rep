#!/usr/bin/env Rscript

library(tidyverse)
library(ape)
library(evolnets)
library(treepplr)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

symbiont_tree_fn <- args[1]
host_tree_fn <- args[2]
interactions_csv_fn <- args[3]
phyjson_path <- args[4]

symbiont_tree <- read_tree_from_revbayes(symbiont_tree_fn)

host_tree <- read.tree(host_tree_fn)
m <- as.matrix(read.csv(interactions_csv_fn, row.names = 1))

###### Prepare input for treepplr ######
#
# symbiont_tree: TreeLabeled, ntips: Int, nhosts: Int,
# interactions: Int[], host_distances: Real[],
# dMean: Real, tune: Real

ntips <- Ntip(symbiont_tree)
nhosts <- Ntip(host_tree)
interactions <- m
host_distances <- cophenetic.phylo(host_tree)
dMean <- sum(host_distances) / factorial(nhosts)
tune <- 0.9

phyjson_tree <- treepplr::tp_phylo_2_phyjson(symbiont_tree)
phyjson_tree_rec <- tp_phyjson_list(phyjson_tree)
json_obj <- fromJSON(RJSONIO::toJSON(phyjson_tree_rec))
json_obj$ntips <- ntips
json_obj$nhosts <- nhosts
json_obj$ntips <- ntips
json_obj$nhosts <- nhosts
json_obj$interactions <- c(t(interactions))
json_obj$host_distances <- c(t(host_distances))
json_obj$dMean <- dMean
json_obj$tune <- tune
write_json(json_obj, phyjson_path, pretty = TRUE, auto_unbox = TRUE)
