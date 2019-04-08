#!/usr/bin/env bash

# fit qpGraph models
python qpbrute.py \
    --par test/sim1.par \
    --prefix sim1 \
    --pops A B C X \
    --out Out

# calculate Bayes Factors
python qpbayes.py \
    --geno test/sim1.geno \
    --ind test/sim1.ind \
    --snp test/sim1.snp \
    --prefix sim1 \
    --pops A B C X \
    --out Out