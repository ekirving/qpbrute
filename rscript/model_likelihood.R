#!/usr/bin/env Rscript
library(admixturegraph)
library(coda)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
graph_code <- args[2]
dstats_file <- args[3]
num_chains <- strtoi(args[4])
num_temps <- strtoi(args[5])
num_iters <- strtoi(args[6])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# graph_code <- '1d9676e'
# dstats_file <- 'dstats/pygmyhog.csv'
# num_chains <- 2
# num_temps <- 10
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
off <- dev.off()

# plot the D-stats
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-dstat.pdf'))
plot(f4stats(dstats.pos))
off <- dev.off()

# fit the graph
graph_fit <- fit_graph(dstats.pos, graph)

# plot the fit
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-fit.pdf'))
plot(graph_fit)
off <- dev.off()

# setup the MCMC model
mcmc <- make_mcmc_model(graph, dstats)

# burn in by 10%
burn <- num_iters * 0.1

# thin to 1%
thin <- 100

# make the run repeatable
seed <- floor(runif(1, min=0, max=1e5))
set.seed(seed)

cat("Starting MCMC...\n")
cat("Graph: ", graph_code, "\n")
cat("Chains: ", num_chains, "\n")
cat("Temps: ", num_temps, "\n")
cat("Iters: ", num_iters, "\n")
cat("Burn in: ", burn, "\n")
cat("Thin: ", thin, "\n")
cat("Seed: ", seed, "\n", "\n")

chains <- c()

for (i in 1:num_chains) {
    cat("Starting chain: ", i, "\n")

    # choose some random starting values for the params
    initial <- runif(length(mcmc$parameter_names))

    # see https://www.rdocumentation.org/packages/admixturegraph/versions/1.0.2/topics/run_metropolis_hasting
    chain <- run_metropolis_hasting(mcmc, initial, no_temperatures = num_temps,
                                    iterations = num_iters, verbose = TRUE)
    cat("\n")

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
    off <- dev.off()

    # add the chain to the list
    chains[[i]] <- mcmc.thin
}

chains.all = mcmc.list(chains)

# plot the combined traces
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-trace-0.pdf'))
plot(chains.all)
off <- dev.off()

# plot the Gelman and Rubin's convergence diagnostic
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-gelman.pdf'))
gelman.plot(chains.all)
off <- dev.off()

# print some summary details
summary(chains.all)

# NB. ideal is 0.234
cat("Acceptance Rate...\n")
cat(1 - rejectionRate(mcmc.thin), "\n\n")

cat("Effective Sample Size...\n")
effectiveSize(chains.all)
cat("\n")

# NB. values substantially above 1 indicate lack of convergence.
cat("Gelman and Rubin's convergence diagnostic...\n")
cat(gelman.diag(chains.all))
