#!/usr/bin/env Rscript
require(admixturegraph)
library(gtools)
library(stringr)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)
library(raster)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'

# regex pattern for the MCMC chains
regx <- paste0(prefix, "-(.+)-thinned.mtx")

# load all the thinned chains
files <- list.files(path='bayes', pattern=regx, full.names=TRUE)
chains <- lapply(files, function(x) data.matrix(read.csv(x, sep = "\t")))
graphs <- str_match(files,regx)[,2]

# make a matrix to hold all the pairwise combinations
num_chains <- length(chains)
perms <- permutations(num_chains, 2)
mtx <- matrix(nrow=num_chains, ncol=num_chains)
rownames(mtx) <- graphs
colnames(mtx) <- graphs

# compute the likelihoods for all models
ll <- data.frame()
for(i in 1:length(chains)) {
    ll <- rbind(ll, model_likelihood_n(chains[[i]][, "likelihood"], 100))
}
rownames(ll) <- graphs
write.csv(ll, file=paste0('bayes/', prefix, "-likelihood.csv"))

# compare the Bayes factors for all model pairs
for(i in 1:nrow(perms)) {
    x <- perms[i,1]
    y <- perms[i,2]
    bayes <- model_bayes_factor_n(chains[[x]][, "likelihood"],
                                  chains[[y]][, "likelihood"], 100)
    mtx[x,y] <- bayes[,'mean']
}

# we consider K > 150 to be very strong statistical support
# see https://en.wikipedia.org/wiki/Bayes_factor#Interpretation
melted_mtx <- melt(mtx)
melted_mtx$value <- clamp(melted_mtx$value, lower=-150, upper=150)

# plot the heatmap
pdf(file=paste0('bayes/', prefix, "-heatmap.pdf"), width=9, height=7)
ggplot(data = melted_mtx, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text.align = 1) +
    scale_fill_viridis(name = "K", na.value = 'gainsboro', option='viridis')
                       # rescale the color palette so the zero threshold is obvious
                       # values=rescale(c(-1, 0-.Machine$double.eps, 0, 0+.Machine$double.eps,1)))
dev.off()
