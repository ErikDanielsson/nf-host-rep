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

height <- node.depth.edgelength(tree)[1]
scaling_factor <- 2.0 / height
tree$edge.length <- tree$edge.length * scaling_factor

is.binary(tree)
is.ultrametric(tree)
is.rooted(tree)

# Create the host tree
host_tree <- rcoal(ntips_host, rooted = TRUE)
host_tree$tip.label <- paste0("Host", 1:ntips_host)

host_tree_fn <- paste0("host_tree.", genid, ".tre")
symbiont_tree_fn <- paste0("symbiont_tree.", genid, ".tre")
write.tree(host_tree, host_tree_fn)
write.tree(tree, symbiont_tree_fn)

# Add subroot branch to tree
tree_string <- readLines(symbiont_tree_fn)
tree_tiny_stem_string <- sub(");$", "):0.01;", tree_string)
tiny_tree_fn <- paste0("tiny_stem_tree.", genid, ".tre")
writeLines(tree_tiny_stem_string, tiny_tree_fn)
tree_long_stem_string <- sub(");$", "):2.0;", tree_string)
long_tree_fn <- paste0("long_stem_tree.", genid, ".tre")
writeLines(tree_long_stem_string, long_tree_fn)
