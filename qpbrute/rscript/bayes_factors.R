#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2018
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(admixturegraph))
quiet(library(gtools))
quiet(library(stringr))
quiet(library(ggplot2))
quiet(library(reshape2))
quiet(library(viridis))
quiet(library(scales))
quiet(library(raster))
quiet(library(data.table))
quiet(library(rjson))
quiet(library(dplyr))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
num_burn <- strtoi(args[2])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# num_burn <- 1e5

# load any custom burn in values
burn_over <- paste0(prefix, "/", prefix, '-burnin.csv')
if (file.exists(burn_over)) {
    burn <- read.csv(burn_over, row.names = 1, col.names = c('', 'burn'), header = F)
}

# load all the MCMC chains
mcmc.regex <- paste0(prefix, "-(.+)-chain-(\\d+).csv")
mcmc.files <- list.files(path=paste0(prefix, "/bayes"), pattern=mcmc.regex, full.names=TRUE)
names(mcmc.files) <- str_match(mcmc.files, mcmc.regex)[,2]
graphs <- unique(names(mcmc.files))

# load all the chains, and burn them in
chains.all <- lapply(mcmc.files, function(x) {
    graph <- str_match(x, mcmc.regex)[,2]
    offset <- ifelse(!exists('burn') || is.na(burn[graph,]), num_burn, burn[graph,])
    cat("Loading chain: ", x, "\t", "Burn in:", offset, "\n")
    burn_in(fread(x, header = T, sep = ','), k=offset)
})

chains <- list()
cat("Merging replicate chains.", "\n")
for (graph in graphs) {
    chains[[graph]] <- rbindlist(chains.all[names(chains.all) == graph])
}
remove(chains.all)

cat("Computing the likelihoods for all models.", "\n")
ll <- data.frame()
for (graph_code in graphs) {
  # get the likelihood
  likelihood <- model_likelihood_n(chains[[graph_code]][,'likelihood'], 100)
  
  # get the PSRF scores
  psrf <- data.frame(fromJSON(file = paste0(prefix, "/bayes/", prefix, "-", graph_code, '-chainAll-psrf.json')))
  names(psrf) <- paste0(ifelse(names(psrf) != "mpsrf", "psrf_", ""), names(psrf))
  
  # merge into a dataframe
  ll <- dplyr::bind_rows(ll, cbind(likelihood, psrf))
}
rownames(ll) <- graphs

# reverse sort by likelihood
ll <- ll[order(-ll$mean),]
chains <- chains[row.names(ll)]

# save the model likelihoods
write.csv(ll, file=paste0(prefix, "/bayes/", prefix, "-likelihoods.csv"))

num_chains <- length(chains)

if (num_chains >= 2) {

    # make a matrix to hold all the pairwise combinations
    perms <- permutations(num_chains, 2, names(chains))
    mtx <- matrix(nrow=num_chains, ncol=num_chains)
    rownames(mtx) <- names(chains)
    colnames(mtx) <- names(chains)

    cat("Comparing the Bayes factors for all model pairs.", "\n")
    for(i in 1:nrow(perms)) {
        x <- perms[i,1]
        y <- perms[i,2]
        K <- model_bayes_factor_n(chains[[x]][, "likelihood"],
                                  chains[[y]][, "likelihood"], 100)
        mtx[x,y] <- K[,'mean']
    }

    # save the Bayes factors
    write.csv(mtx, file=paste0(prefix, "/bayes/", prefix, "-bayes-factors.csv"))
    
    # convert the matrix into the 3-column format required by geom_tile
    melted_mtx <- reshape2::melt(mtx)

    # we consider K > 150 to be very strong statistical support, so cap at that threshold
    # see https://en.wikipedia.org/wiki/Bayes_factor#Interpretation
    melted_mtx$value <- clamp(melted_mtx$value, lower=-150, upper=150)

    # set the sort oder for the heatmap by reordering the factors
    melted_mtx$Var1 <- factor(melted_mtx$Var1, levels=rev(names(chains)))
    melted_mtx$Var2 <- factor(melted_mtx$Var2, levels=names(chains))

    plt <- ggplot(data = melted_mtx, aes(x=Var2, y=Var1, fill=value)) +
      geom_tile() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text.align = 1) +
      scale_fill_viridis(name = "K", na.value = 'gainsboro', option='viridis')
      # rescale the color palette so the zero threshold is obvious
      # values=rescale(c(-1, 0-.Machine$double.eps, 0, 0+.Machine$double.eps,1)))
  
    # plot the heatmap
    pdf(file=paste0(prefix, "/bayes/", prefix, "-heatmap.pdf"), width=9, height=7)
    print(plt)
    dev.off()
}

cat("FINISHED!")