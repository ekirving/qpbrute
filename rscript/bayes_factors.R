#!/usr/bin/env Rscript
library(admixturegraph, quietly=T)
library(gtools, quietly=T)
library(stringr, quietly=T)
library(ggplot2, quietly=T)
library(reshape2, quietly=T)
library(viridis, quietly=T)
library(scales, quietly=T)
library(raster, quietly=T)
library(data.table, quietly=T)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
num_burn <- strtoi(args[2])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/qpbrute')
# prefix <- 'pygmyhog'
# num_burn <- 1e5

# load all the thinned chains
mcmc.regex <- paste0(prefix, "-(.+)-chain-(\\d).csv")
mcmc.files <- list.files(path='bayes', pattern=mcmc.regex, full.names=TRUE)
names(mcmc.files) <- str_match(mcmc.files, mcmc.regex)[,2]
graphs <- unique(names(mcmc.files))

# load all the chains, and burn them in
chains.all <- lapply(mcmc.files, function(x) {
    cat("Loading chain: ", x, "\n")
    burn_in(fread(x, header = T, sep = ','), k=num_burn)
})

chains <- list()
cat("Merging replicate chains.", "\n")
for (graph in graphs) {
    chains[[graph]] <- rbindlist(chains.all[names(chains.all) == graph])
}
remove(chains.all)

cat("Computing the likelihoods for all models.", "\n")
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

cat("Comparing the Bayes factors for all model pairs.", "\n")
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
