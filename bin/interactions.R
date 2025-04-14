#!/usr/bin/env Rscript

# Create interaction matrix
# NOTE: here I am not generating a interaction matrix that conveys the
# phylogenetic signal in the host trees, I am just sampling random entries
# according to the distribution defined below.
# Perhaps it makes more sense to sample conditional on the host tree if we want
# to guarantee that the interactions has a phylogenetic signal

library(ape)

args <- commandArgs(trailingOnly = TRUE)
genid <- args[1]
ntips_symbiont <- as.integer(args[2])
ntips_host <- as.integer(args[3])

set.seed(genid)

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

rownames(interactions) <- paste0("Symbiont", seq_len(nrow(interactions)))
colnames(interactions) <- paste0("Host", seq_len(ncol(interactions)))

# Write output files
interaction_csv_fn <- paste0("interactions.", genid, ".csv")
interaction_nex_fn <- paste0("interactions.", genid, ".nex")
write.table(interactions, interaction_csv_fn, row.names = TRUE, quote = FALSE, sep = "\t")
write.nexus.data(interactions, interaction_nex_fn, format = "standard")
