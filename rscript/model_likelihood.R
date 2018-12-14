#!/usr/bin/env Rscript
library(admixturegraph)
library(coda)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
graph_code <- args[2]
csv_file <- args[3]
num_temps <- strtoi(args[4])
num_iters <- strtoi(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# graph_code <- '501575a'
# csv_file <- 'dstats/pygmyhog.csv'
# num_temps <- 3
# num_iters <- 1e5

# load the Dstat data
dstats <- read.csv(csv_file)
dstats.pos <- dstats[dstats$D >= 0,] # drop negative results (because they are symetrical)

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
            edges <- c(edges, admixture_edge(line[2], line[3], line[4], line[2]))
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
plot(f4stats(dstats.pos))
dev.off()

# fit the graph
graph_fit <- fit_graph(dstats.pos, graph)

# plot the fit
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-fit.pdf'))
plot(graph_fit)
dev.off()

# setup the MCMC model
mcmc <- make_mcmc_model(graph, dstats)

# choose some random starting values for the params
initial <- runif(length(mcmc$parameter_names))

chain.csv <- paste0('bayes/', prefix, '-', graph_code, '-chain.csv')

# the MCMC can take a long time, so let's break it into chunks
chunk_size <- 1e3
num_chunks <- num_iters / chunk_size

for (i in 1:num_chunks) {
    # see https://www.rdocumentation.org/packages/admixturegraph/versions/1.0.2/topics/run_metropolis_hasting
    chunk <- run_metropolis_hasting(mcmc, initial, no_temperatures = num_temps,
                                    iterations = chunk_size, verbose = FALSE)
    # save the current chunk
    write.table(chunk, file=chain.csv, sep = ",", row.names = FALSE, col.names = (i == 1), append = (i != 1))

    # update the starting position
    initial <- as.numeric(as.vector(tail(subset(chunk, select=-c(prior, likelihood, posterior)), 1)))
}

# load the whole chain
chain <- read.csv(file=chain.csv)

# check ESS to make sure the MCMC has converged
ess <- t(effectiveSize(subset(chain, select=-c(prior, likelihood, posterior))))
write.csv(ess, file=paste0('bayes/', prefix, '-', graph_code, '-ess.csv'), row.names = FALSE)

# burn in and thin the chain
thinned <- thinning(burn_in(chain, 4000), 100)

# save the thin chain
write.csv(thinned, file=paste0('bayes/', prefix, "-", graph_code, '-thinned.csv'))
