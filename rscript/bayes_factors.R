#!/usr/bin/env Rscript
require(admixturegraph)
library(gtools)
library(stringr)
library(ggplot2)
library(reshape2)
library(viridis)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[2]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/qpbrute')
prefix <- 'pygmyhog'

# regex pattern for the MCMC chains
regx <- paste0(prefix, "-(.+)-thinned.mtx")

# load all the thinned chains
files <- list.files(path='bayes', pattern=regx, full.names=TRUE)
chains <- lapply(files, function(x) x=data.matrix(read.csv(x, sep = "\t")))
graphs <- str_match(files,regx)[,2]

# make a matrix to hold all the pairwise combinations
num_chains <- length(chains)
perms <- permutations(num_chains, 2)
mtx <- matrix(nrow=num_chains, ncol=num_chains)
rownames(mtx) <- graphs
colnames(mtx) <- graphs

# compare the Bayes factors for all model pairs
for(i in 1:nrow(perms)) {
    x <- perms[i,1]
    y <- perms[i,2]
    bayes <- model_bayes_factor_n(chains[[x]][, "likelihood"],
                                  chains[[y]][, "likelihood"], 100)
    mtx[x,y] <- bayes[,'mean']
}

melted_mtx <- melt(mtx)

# plot the heatmap
pdf(file=paste0('bayes/', prefix, "-heatmap.pdf"), width=9, height=7)
ggplot(data = melted_mtx, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_fill_viridis(name = "Bayes factor", na.value = 'gainsboro')
dev.off()

# model_likelihood_n(thinned[, "likelihood"], 100)
