#!/usr/bin/env Rscript
library(admixturegraph, quietly=T)
library(coda, quietly=T)
library(fitR, quietly=T)  # see http://sbfnk.github.io/mfiidd/introduction.html

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
# graph_code <- 'ghost2'
# dstats_file <- 'dstats/pygmyhog.csv'
# num_chains <- 2
# num_temps <- 5
# num_iters <- 1e6

# default to 10% burn in and thinning
burn <- num_iters / 10
thin <- 10

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

# make the run repeatable by setting a random seed
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

    fullchain.file = paste0('bayes/', prefix, '-', graph_code, '-chain-', i, '.csv')

    # don't rerun completed chains
    if (file.exists(fullchain.file)) {
        cat("Loading chain: ", i, "\n")

        chain <- read.csv(fullchain.file)

    } else {
        cat("Starting chain: ", i, "\n")

        # choose some random starting values for the params
        initial <- runif(length(mcmc$parameter_names))

        # see https://www.rdocumentation.org/packages/admixturegraph/versions/1.0.2/topics/run_metropolis_hasting
        chain <- run_metropolis_hasting(mcmc, initial, no_temperatures = num_temps,
                                        iterations = num_iters, verbose = TRUE)

        # save the full chain
        write.csv(chain, file=fullchain.file, row.names = FALSE)
    }

    # convert to MCMC object
    mcmc.chain <- mcmc(chain)

    cat("\n")

    # check the acceptance rate (ideal is 0.234)
    cat("Acceptance Rate =", 1 - rejectionRate(mcmc.chain)[1], "\n\n")

    # plot ESS vs. burn-in
    cat("Plotting ESS vs. burn-in.", "\n\n")
    pdf(file=paste0('bayes/', prefix, "-", graph_code, '-ess-burn-', i, '.pdf'), width=21, height=14)
    plotESSBurn(mcmc.chain, step.size=burn/2)
    off <- dev.off()

    cat("Thinning chain: ", i, "\n")

    # burn in and thin the chain
    chain.thin <- thinning(burn_in(mcmc.chain, k=burn), k=thin)
    mcmc.thin <- mcmc(chain.thin, start=burn, thin=thin)

    # save the thin chain
    write.csv(mcmc.thin, file=paste0('bayes/', prefix, "-", graph_code, '-thinned-', i, '.csv'), row.names = FALSE)

    # print the summary stats
    print(summary(mcmc.thin))

    cat("Effective Sample Size.\n")
    print(effectiveSize(mcmc.thin))
    cat("\n")

    cat("Plotting autocorrelation.\n\n")
    pdf(file=paste0('bayes/', prefix, "-", graph_code, '-thin-autocorr-', i, '.pdf'))
    autocorr.plot(mcmc.thin)
    off <- dev.off()

    cat("Plotting the trace.\n\n")
    pdf(file=paste0('bayes/', prefix, "-", graph_code, '-thin-trace-', i, '.pdf'))
    plot(mcmc.thin)
    off <- dev.off()

    # add the thinnned chain to the list
    chains[[i]] <- mcmc.thin
}

cat("\n\n", "--------------", "\n\n")
cat("Analysing combined chains.", "\n\n")

chains.all = mcmc.list(chains)

# print the summary stats
print(summary(chains.all))

cat("Effective Sample Size.\n")
print(effectiveSize(chains.all))
cat("\n")

cat("Plotting autocorrelation.", "\n\n")
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-thin-autocorr-0.pdf'))
autocorr.plot(chains.all)
off <- dev.off()

# plot the combined traces
cat("Plotting combined traces.", "\n\n")
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-thin-trace-0.pdf'))
plot(chains.all)
off <- dev.off()

cat("Plotting the Gelman and Rubin's convergence diagnostic.", "\n\n")
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-thin-gelman.pdf'))
gelman.plot(chains.all)
off <- dev.off()

# NB. values substantially above 1 indicate lack of convergence.
cat("Gelman and Rubin's convergence diagnostic.", "\n")
print(gelman.diag(chains.all))
