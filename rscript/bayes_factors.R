#!/usr/bin/env Rscript
require(admixturegraph)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_prefix <- args[2]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/qpbrute')
output_prefix <- 'pygmyhog'

# load all the thinned matrices
matrices <- list.files(pattern = "\\.dbf$")




# model_likelihood_n(thinned[, "likelihood"], 100)
# model_bayes_factor_n(thinned[, "likelihood"], thinned[, "likelihood"], 100)
