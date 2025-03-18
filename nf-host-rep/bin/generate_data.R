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
tree$tip.label <- paste0("S", 1:ntips_symbiont)

height <- node.depth.edgelength(tree)[1]
scaling_factor <- 2.0 / height
tree$edge.length <- tree$edge.length * scaling_factor

is.binary(tree)
is.ultrametric(tree)
is.rooted(tree)

# Create the host tree
host_tree <- rcoal(ntips_host, rooted = TRUE)
host_tree$tip.label <- paste0("H", 1:ntips_host)

# Create interaction matrix
# NOTE: here I am not generating a interaction matrix that conveys the
# phylogenetic signal in the host trees, I am just sampling random entries
# according to the distribution defined below.
# Perhaps it makes more sense to sample conditional on the host tree if we want
# to guarantee that the interactions has a phylogenetic signal

# First create a basic matrix where each row and col has atleast one 2
n <- max(ntips_symbiont, ntips_host)
m <- min(ntips_host, ntips_symbiont)
interactions <- diag(2, m, m)
if (n > m) {
    a <- diag(2, n - m, m)
    interactions <- rbind(a, interactions)
    if (ntips_symbiont < ntips_host) {
        interactions <- t(interactions)
    }
}
# Randomly permute the columns of the matrix
interactions <- interactions[, sample(ncol(interactions))]

# Add in some random number for an extra dose of randomness
state_dist <- c(0.75, 0.0, 0.25)
random_entries <- matrix(
    data = sample(
        0:2, ntips_symbiont * ntips_host,
        replace = TRUE, prob = state_dist
    ), nrow = ntips_symbiont, ncol = ntips_host, byrow = TRUE
)
interactions <- pmax(interactions, random_entries)
print(interactions)

rownames(interactions) <- tree$tip.label
colnames(interactions) <- host_tree$tip.label

# Write output files
interaction_csv_fn <- paste0("interactions.", genid, ".csv")
interaction_nex_fn <- paste0("interactions.", genid, ".nex")
host_tree_fn <- paste0("host_tree.", genid, ".tre")
symbiont_tree_fn <- paste0("symbiont_tree.", genid, ".tre")
write.csv(interactions, interaction_csv_fn, row.names = TRUE)
write.nexus.data(interactions, interaction_nex_fn, format = "standard")
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
