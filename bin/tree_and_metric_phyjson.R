#! /usr/bin/env Rscript
library(tidyverse)
library(ape)
library(evolnets)
library(treepplr)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

symbiont_tree_fn <- args[1]
host_tree_fn <- args[2]
phyjson_path <- args[3]
interaction_params_path <- args[4]

symbiont_tree <- read_tree_from_revbayes(symbiont_tree_fn)
host_tree <- read.tree(host_tree_fn)

ntips <- Ntip(symbiont_tree)
nhosts <- Ntip(host_tree)
host_distances <- cophenetic.phylo(host_tree)
# The mean distance is the sum of all distances divided by the number of ordered pairs
meanDistance <- sum(host_distances) / (nhosts * (nhosts - 1))
tune <- 1.0

phyjson_tree <- treepplr::tp_phylo_2_phyjson(symbiont_tree)
phyjson_tree_rec <- tp_phyjson_list(phyjson_tree)
json_obj <- fromJSON(RJSONIO::toJSON(phyjson_tree_rec))
json_obj$ntips <- ntips
json_obj$nhosts <- nhosts
json_obj$ntips <- ntips
json_obj$nhosts <- nhosts
json_obj$host_distances <- c(t(host_distances))
json_obj$dMean <- meanDistance
json_obj$tune <- tune
write_json(json_obj, phyjson_path, pretty = TRUE, auto_unbox = TRUE)
