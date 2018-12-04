#!/usr/bin/env bash

python qpbayes.py \
    --geno pygmyhog/autosomes.filtered.geno \
    --ind pygmyhog/autosomes.filtered.ind \
    --snp pygmyhog/autosomes.filtered.snp \
    --prefix pygmyhog \
    --pops LIB EUWB NCWB SCWB ISEA_SVSV \
    --out AF
    