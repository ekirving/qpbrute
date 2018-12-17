#!/usr/bin/env Rscript
library(admixturegraph)
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

# load all the thinned chains
mcmc.regex <- paste0(prefix, "-(.+)-thinned.csv")
mcmc.files <- list.files(path='bayes', pattern=mcmc.regex, full.names=TRUE)
graphs <- str_match(mcmc.files, mcmc.regex)[,2]
names(mcmc.files) <- graphs
chains <- lapply(mcmc.files, function(x) data.matrix(read.csv(x)))

# compute the likelihoods for all models
ll <- data.frame()
for (graph in graphs) {
    ll <- rbind(ll, model_likelihood_n(chains[[graph]][,'likelihood'], 100))
}
rownames(ll) <- graphs

# reverse sort by likelihood
ll <- ll[order(-ll$mean),]
chains <- chains[row.names(ll)]

# save the model likelihoods
write.csv(ll, file=paste0('bayes/', prefix, "-likelihoods.csv"))

# make a matrix to hold all the pairwise combinations
num_chains <- length(chains)
perms <- permutations(num_chains, 2, names(chains))
mtx <- matrix(nrow=num_chains, ncol=num_chains)
rownames(mtx) <- names(chains)
colnames(mtx) <- names(chains)

# compare the Bayes factors for all model pairs
for(i in 1:nrow(perms)) {
    x <- perms[i,1]
    y <- perms[i,2]
    K <- model_bayes_factor_n(chains[[x]][, "likelihood"],
                              chains[[y]][, "likelihood"], 100)
    mtx[x,y] <- K[,'mean']
}

# convert the matrix into the 3-column format required by geom_tile
melted_mtx <- melt(mtx)

# we consider K > 150 to be very strong statistical support, so cap at that threshold
# see https://en.wikipedia.org/wiki/Bayes_factor#Interpretation
melted_mtx$value <- clamp(melted_mtx$value, lower=-150, upper=150)

# set the sort oder for the heatmap by reordering the factors
melted_mtx$Var1 <- factor(melted_mtx$Var1, levels=rev(names(chains)))
melted_mtx$Var2 <- factor(melted_mtx$Var2, levels=names(chains))

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
