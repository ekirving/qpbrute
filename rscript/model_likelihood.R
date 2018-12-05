#!/usr/bin/env Rscript
library(admixturegraph)
library(MASS)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
graph_code <- args[2]
csv_file <- args[3]

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# graph_code <- '1d9676e'
# csv_file <- 'dstats/pygmyhog.csv'

# load the Dstat data
dstats <- read.csv(csv_file)
dstats <- dstats[dstats$D >= 0,] # drop negative results (because they are symetrical)

# Convert qpGraph model to R format
convert_qpgraph <- function(graph_file) {
    leaves <- c()
    nodes <- c()
    edges <- c()

    # load the graph into memory
    fin <- file(graph_file, open="r")
    lines <-readLines(fin)
    close(fin)

    # iterate over each line
    for (i in 1:length(lines)) {
        line <- strsplit(lines[i], split="\t")[[1]]

        # process the different types of entries in the graph file
        if (line[1] == 'root') {
            nodes <- c(nodes, line[2])
        } else if (line[1] == 'label') {
            leaves <- c(leaves, line[2])
        } else if (line[1] == 'edge') {
            edges <- c(edges, edge(line[4], line[3]))
            nodes <- c(nodes, line[4])
        } else if (line[1] == 'admix') {
            edges <- c(edges, admixture_edge(line[2], line[3], line[4]))
            nodes <- c(nodes, line[2])
        }
    }

    # get all the inner nodes
    inner_nodes <- setdiff(nodes, leaves)

    # make the graph object
    agraph(leaves, inner_nodes, parent_edges(edges))
}

# perform the conversion
graph <- convert_qpgraph(paste0("graphs/", prefix, "-", graph_code, ".graph"))

# plot the graph
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-graph.pdf'))
plot(graph)
dev.off()

# plot the D-stats
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-dstat.pdf'))
plot(f4stats(dstats))
dev.off()

# fit the graph
graph_fit <- fit_graph(dstats, graph)

# plot the fit
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-fit.pdf'))
plot(graph_fit)
dev.off()

# run the MCMC on the params
mcmc <- make_mcmc_model(graph, dstats)
initial <- rep(0.5, length(mcmc$parameter_names))
chain <- run_metropolis_hasting(mcmc, initial, iterations = 10000, verbose = TRUE)

# TODO check ESS to make sure the MCMC has converged

# save the chain (just in case we want it later)
write.matrix(chain, file=paste0('bayes/', prefix, "-", graph_code, '-chain.mtx'), sep = "\t")

# burn in and thin the chain
thinned <- thinning(burn_in(chain, 4000), 100)

# save the thin chain
write.matrix(thinned, file=paste0('bayes/', prefix, "-", graph_code, '-thinned.mtx'), sep = "\t")
