#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2018
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(admixturegraph))
quiet(library(coda))
quiet(library(fitR))
quiet(library(data.table))
quiet(library(rjson))
quiet(library(stringr))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
graph_code <- args[2]
dstats_file <- args[3]
num_chains <- strtoi(args[4])
num_temps <- strtoi(args[5])
num_iters <- strtoi(args[6])
num_burn <- strtoi(args[7])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# graph_code <- '1d9676e'
# dstats_file <- 'dstats/pygmyhog.csv'
# num_chains <- 2
# num_temps <- 5
# num_iters <- 2e6
# num_burn <- 1.1e6

# maximum PSRF
max_psrf <- 1.2

# load any custom burn in values
burn_file <- paste0(prefix, '-burnin.csv')
if (file.exists(burn_file)) {
    burn <- read.csv(burn_file, row.names = 1, col.names = c('', 'burn'), header = F)
    num_burn <- ifelse(is.na(burn[graph_code,]), num_burn, burn[graph_code,])
}

# load the Dstat data
dstats <- read.csv(dstats_file)

# Convert qpGraph model to R format
convert_qpgraph <- function(graph_file) {
    leaves <- c()
    nodes <- c()
    edges <- c()

    # load the graph into memory
    fin <- file(graph_file, open="r")
    lines <- readLines(fin)
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
graph <- convert_qpgraph(paste0(prefix, "/graphs/", prefix, "-", graph_code, ".graph"))

# plot the graph
pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-graph.pdf'))
plot(graph, show_admixture_labels = TRUE)
dev.off()

# plot the D-stats
pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-dstat.pdf'))
plot(f4stats(dstats))
dev.off()

# number of inner nodes
n <- length(graph$inner_nodes)

# fit the graph
graph_fit <- fit_graph(
  dstats, 
  graph, 
  # allow some extra time to find the best fit
  optimisation_options = neldermead::optimset(
    method='fminbnd', 
    MaxFunEvals=400*n, 
    MaxIter=400*n
  )
)

cat("Summary of graph fit...\n")
cat(summary(graph_fit))
cat("\n")

# plot the fit
pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-fit.pdf'))
plot(graph_fit)
dev.off()

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
cat("Burn in: ", num_burn, "\n")
cat("Seed: ", seed, "\n", "\n")

chains <- c()

run_chain <- function(i, num_iters) {

    chain.prev <- NULL
    chain.new <- NULL

    fullchain.file <- paste0(prefix, "/bayes/", prefix, '-', graph_code, '-chain-', i, '.csv')

    # check if the chain exists
    if (file.exists(fullchain.file) && file.info(fullchain.file)$size > 0) {

        cat("Loading chain: ", i, "\n")
        chain.prev <- fread(fullchain.file, header = T, sep = ',')
        row.count <- nrow(chain.prev)

        if (row.count == num_iters) {
            cat("Chain has correct number of iterations.", "\n")
            num_iters <- 0

        } else if (row.count < num_iters) {
            # update the number of iterations
            num_iters <- num_iters - row.count

            cat("Restarting chain with", num_iters, "more iterations.", "\n")

            # set the Markov chain to restart from the last position
            initial <- as.numeric(setDF(tail(chain.prev, 1))[mcmc$parameter_names])

            # starting with a value of 1 breaks the MCMC
            initial <- replace(initial,initial==1, 0.9999999)

        } else {
            cat("Truncating chain to length", num_iters, "iterations.", "\n")
            chain.prev <- head(chain.prev, num_iters)
            num_iters <- 0
        }

    } else {
        cat("Starting new chain: ", i, "\n")

        # choose some random starting values for the params
        initial <- runif(length(mcmc$parameter_names))
    }

    if (num_iters > 0) {
        # run the MCMC chain for the required number of iterations
        # see https://www.rdocumentation.org/packages/admixturegraph/versions/1.0.2/topics/run_metropolis_hasting
        chain.new <- run_metropolis_hasting(mcmc, initial, no_temperatures = num_temps,
                                            iterations = num_iters, verbose = TRUE)
    }

    cat("\n\n")

    # merge any old and new chains
    chain <- rbind(chain.prev, chain.new)

    if (!is.null(chain.new)) {
        # save the full chain
        write.csv(chain, file=fullchain.file, row.names = FALSE)
    }

    # convert to MCMC object and return
    mcmc(chain)
}

for (i in 1:num_chains) {

    # run the MCMC chain
    mcmc.chain <- run_chain(i, num_iters)

    # check the acceptance rate (ideal is 0.234)
    cat("Acceptance Rate =", 1 - rejectionRate(mcmc.chain)[1], "\n\n")

    # burn in the chain
    mcmc.burn <- mcmc(burn_in(mcmc.chain, k=num_burn), start=num_burn)

    # calculate the ESS for all params
    ess <- effectiveSize(mcmc.burn)

    if (min(ess) < 100) {
        cat(paste0("WARNING: ESS below threshold. min(ess) = ", min(ess),
                       " ./", prefix, "/bayes/", prefix, '-', graph_code, '-chain-', i, '.csv'))
    }

    # NB. we do not thin the chain because there is no need
    # see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00131.x

    # print the summary stats
    print(summary(mcmc.burn))

    cat("Effective Sample Size.\n")
    print(ess)
    cat("\n")

    # plot ESS vs. burn-in
    cat("Plotting ESS vs. burn-in.", "\n\n")
    pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-ess-burn-', i, '.pdf'), width=21, height=14)
    plotESSBurn(mcmc.chain, step.size=round(num_burn/2, 0))
    dev.off()

    cat("Plotting thinned autocorrelation.\n\n")
    pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-burn-thin-autocorr-', i, '.pdf'))
    autocorr.plot(thinning(mcmc.burn, k=1000), lag.max=50)
    dev.off()

    cat("Plotting the trace.\n\n")
    png(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-burn-trace-', i, '-pt%d.png'), width=7, height=7, units='in', res=300)
    plot(mcmc.burn)
    dev.off()

    # add the burned-in chain to the list (and drop posterior)
    chains[[i]] <- mcmc.burn[ , -which(colnames(mcmc.burn) %in% c("posterior"))]
}

cat("\n\n", "--------------", "\n\n")
cat("Analysing combined chains.", "\n\n")

chains.all <- mcmc.list(chains)

# calculate the summary stats
chains.summary <- summary(chains.all)
print(chains.summary)
write.table(chains.summary$quantiles, file = paste0(prefix, "/bayes/", prefix, "-", graph_code, '-chainAll-summary.tsv'), sep = "\t")

# calculate the ESS for all params across all the replicates
ess <- effectiveSize(chains.all)
cat(toJSON(ess, indent = 2), file = paste0(prefix, "/bayes/", prefix, "-", graph_code, '-chainAll-ess.json'))

cat("Effective Sample Size.\n")
print(ess)
cat("\n")

# plot the combined traces
cat("Plotting combined traces.", "\n\n")
png(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-burn-trace-0-pt%d.png'), width=7, height=7, units='in', res=300)
plot(chains.all)
dev.off()

gelman <- tryCatch({
    # NB. values substantially above 1 indicate lack of convergence.
    gelman.diag(chains.all, multivariate = TRUE, autoburnin = FALSE)
}, error = function(e) {
    # not all models can compute a value multivariate PSRF
    gelman.diag(chains.all, multivariate = FALSE, autoburnin = FALSE)
})

gelman.list <- gelman$psrf[, 1]
gelman.list['mpsrf'] <- if (is.null(gelman$mpsrf)) NA else gelman$mpsrf
cat(toJSON(gelman.list, indent = 2), file = paste0(prefix, "/bayes/", prefix, "-", graph_code, '-chainAll-psrf.json'))

cat("Gelman and Rubin's convergence diagnostic.", "\n")
print(gelman)

if (!is.na(gelman$mpsrf)) {
  if (gelman$mpsrf > max_psrf) {
    cat(paste0("WARNING: MPSRF of model is above threshold = ", round(gelman$mpsrf, 3)), "\n")
  }
} else if (gelman$psrf['lnL', 1] > max_psrf) {
    cat(paste0("WARNING: PSRF of likelihood is above threshold = ", round(gelman$psrf['lnL', 1], 3)), "\n")
}

# do any params have infinite variance (very bad!)
inf.params <- names(gelman$psrf[, 1][gelman$psrf[, 1] == 'Inf'])

if (length(inf.params) > 0) {
    # drop the offending columns so we can print the Gelman plot
    chains.all <- chains.all[, !(colnames(chains.all[[1]]) %in% inf.params)]
}

cat("Plotting the Gelman and Rubin's convergence diagnostic.", "\n\n")
pdf(file=paste0(prefix, "/bayes/", prefix, "-", graph_code, '-burn-gelman.pdf'), width = 7, height = 7)
gelman.plot(chains.all)
dev.off()

cat("FINISHED!")