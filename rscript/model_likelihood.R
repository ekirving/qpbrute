#!/usr/bin/env Rscript
library(admixturegraph)
library(coda)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
graph_code <- args[2]
dstats_file <- args[3]
num_temps <- strtoi(args[4])
num_iters <- strtoi(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# graph_code <- '1d9676e'
# dstats_file <- 'dstats/pygmyhog.csv'
# num_temps <- 3
# num_iters <- 1e6

# load the Dstat data
dstats <- read.csv(dstats_file)
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

# number of independent MCMC chains
num_chains <- 2

# make the run repeatable
seed <- floor(runif(1, min=0, max=1e5))
set.seed(seed)

# burn in by 10%
burn <- num_iters * 0.1

# thin to 1%
thin <- 100

print("Starting MCMC...")
print(paste0("Graph code: ", graph_code))
print(paste0("Num chains: ", num_chains))
print(paste0("Num temperatures: ", num_temps))
print(paste0("Num iterations: ", num_iters))
print(paste0("Burn in: ", burn))
print(paste0("Thin: ", thin))
print(paste0("Random seed: ", seed))

chains <- c()

for (i in 1:num_chains) {
    print(paste0("Starting chain: ", i))

    # choose some random starting values for the params
    initial <- runif(length(mcmc$parameter_names))

    # see https://www.rdocumentation.org/packages/admixturegraph/versions/1.0.2/topics/run_metropolis_hasting
    chain <- run_metropolis_hasting(mcmc, initial, no_temperatures = num_temps,
                                    iterations = num_iters, verbose = TRUE)

    print(paste0("Finished chain: ", i))

    # save the full chain
    write.csv(chain, file=paste0('bayes/', prefix, '-', graph_code, '-chain-', i, '.csv'), row.names = FALSE)

    # burn in and thin the chain
    thinned <- thinning(burn_in(chain, burn), thin)
    mcmc.thin <- mcmc(thinned, start=burn, thin=thin)

    # save the thin chain
    write.csv(thinned, file=paste0('bayes/', prefix, "-", graph_code, '-thinned-', i, '.csv'), row.names = FALSE)

    # compute the ESS for the thinned chain
    ess.thin <- t(effectiveSize(subset(thinned, select=-c(prior, likelihood, posterior))))
    write.csv(ess.thin, file=paste0('bayes/', prefix, '-', graph_code, '-ess-', i, '.csv'), row.names = FALSE)

    # plot the trace
    pdf(file=paste0('bayes/', prefix, "-", graph_code, '-trace-', i, '.pdf'))
    plot(mcmc.thin)
    dev.off()

    # add the chain to the list
    chains[[i]] <- mcmc.thin
}

chains.all = mcmc.list(chains)

# plot the combined traces
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-trace-0.pdf'))
plot(chains.all)
dev.off()

# plot the Gelman and Rubin's convergence diagnostic
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-gelman.pdf'))
gelman.plot(chains.all)
dev.off()

# print some summary details
print("Summary...")
summary(chains.all)

print("Acceptance Rate...")
print(1 - rejectionRate(mcmc.thin))

print("Effective Sample Size...")
effectiveSize(chains.all)

# NB. values substantially above 1 indicate lack of convergence.
print("Gelman and Rubin's convergence diagnostic...")
print(gelman.diag(chains.all))
