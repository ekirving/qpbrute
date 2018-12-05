#!/usr/bin/env Rscript
require(admixturegraph)
library(MASS)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
graph_file <- args[1]
output_prefix <- args[2]
dstat_file <- args[3]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/qpbrute')
graph_file <- 'graphs/pygmyhog-1d9676e.graph'
output_prefix <- 'pygmyhog'
dstat_file <- 'dstats/pygmyhog.csv'

# load the Dstat data
dstats <- read.csv(dstat_file)
dstats <- dstats[dstats$D > 0,] # drop negative results (because they are symetrical)

convert_qpgraph <- function(graph_file) {
    # Convert qpGraph model to R format

    leaves <- c("SCWB", "NCWB", "LIB", "ISEA_SVSV", "EUWB", "AF")
    inner_nodes <- c('R', 'a1', 'a1a', 'a1b', 'a2', 'a2a', 'a2b', 'a3', 'a3a', 'a3b',
                     'n1', 'n2', 'n3', 'n4', 'n5', 'n6')

    edges <- parent_edges(c(edge('a3b', 'n6'), edge('LIB', 'a1'), edge('NCWB', 'n4'),
                             edge('n4', 'n1'), edge('a2a', 'n5'), edge('a3a', 'a2'),
                             edge('a1a', 'n2'), edge('SCWB', 'a2'), edge('AF', 'R'),
                             edge('a1b', 'n5'), edge('EUWB', 'a3'), edge('n5', 'n3'),
                             edge('n6', 'n2'), edge('n1', 'n6'), edge('n2', 'R'),
                             edge('ISEA_SVSV', 'n3'), edge('n3', 'n1'), edge('a2b', 'n4'),

                            admixture_edge('a3', 'a3a', 'a3b'),
                            admixture_edge('a2', 'a2a', 'a2b'),
                            admixture_edge('a1', 'a1a', 'a1b')))

    agraph(leaves, inner_nodes, edges)
}

# perform the conversion
graph <- convert_qpgraph(graph_file)

# plot the graph
pdf(file=paste0('bayes/', output_prefix, '-graph.pdf'))
plot(graph)
dev.off()

# plot the D-stats
pdf(file=paste0('bayes/', output_prefix, '-dstat.pdf'))
plot(f4stats(dstats))
dev.off()

# fit the graph
graph_fit <- fit_graph(dstats, graph)

# plot the fit
pdf(file=paste0('bayes/', output_prefix, '-fit.pdf'))
plot(graph_fit)
dev.off()

# run the MCMC on the params
mcmc <- make_mcmc_model(graph, dstats)
initial <- rep(0.5, length(mcmc$parameter_names))
chain <- run_metropolis_hasting(mcmc, initial, iterations = 10000, verbose = TRUE)

# save the chain (just in case we want it later)
write.matrix(chain, file=paste0('bayes/', output_prefix, '-chain.mtx'), sep = "\t")

# burn in and thin the chain
thinned <- thinning(burn_in(chain, 4000), 100)

# save the thin chain
write.matrix(thinned, file=paste0('bayes/', output_prefix, '-thinned.mtx'), sep = "\t")
