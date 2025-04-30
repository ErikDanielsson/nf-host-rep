#!/usr/bin/env Rscript

#### Small simulated dataset to compare TreePPL and RevBayes ####

library(ape)

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
genid <- args[1]
ntips_symbiont <- as.integer(args[2])
ntips_host <- as.integer(args[3])

# Set a random seed to make sure that we get the same tree every time
set.seed(genid)

# Create the symbiont tree
tree <- rcoal(ntips_symbiont, rooted = TRUE)
tree$tip.label <- paste0("Symbiont", 1:ntips_symbiont)

sim_height <- node.depth.edgelength(tree)[1]
height <- 10.0 * log(ntips_symbiont)
tree$edge.length <- tree$edge.length * height / sim_height

is.binary(tree)
is.ultrametric(tree)
is.rooted(tree)

# Create the host tree
host_tree <- rcoal(ntips_host, rooted = TRUE)
host_tree$tip.label <- paste0("Host", 1:ntips_host)

# Write tree to temp file
tmp_symbiont_tree_fn <- paste0("tmp", genid, ".tre")
write.tree(tree, tmp_symbiont_tree_fn)

# Add tiny subroot branch to tree
subroot_length <- height * 1e-12
tree_string <- readLines(tmp_symbiont_tree_fn)
tree_tiny_stem_string <- paste0(
    sub(");$", "):", tree_string), subroot_length, ";"
)
host_tree_fn <- paste0("host_tree.", genid, ".tre")
symbiont_tree_fn <- paste0("symbiont_tree.", genid, ".tre")
write.tree(host_tree, host_tree_fn)
writeLines(tree_tiny_stem_string, symbiont_tree_fn)
