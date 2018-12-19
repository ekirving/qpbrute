#!/usr/bin/env bash

export OMP_NUM_THREADS=1

nohup time Rscript rscript/model_likelihood.R pygmyhog test1 dstats/pygmyhog.csv 2  5  100000 &> nohup-test1.out &
nohup time Rscript rscript/model_likelihood.R pygmyhog test2 dstats/pygmyhog.csv 2 10  100000 &> nohup-test2.out &
nohup time Rscript rscript/model_likelihood.R pygmyhog test3 dstats/pygmyhog.csv 2  5  500000 &> nohup-test3.out &
nohup time Rscript rscript/model_likelihood.R pygmyhog test4 dstats/pygmyhog.csv 2 10  500000 &> nohup-test4.out &
nohup time Rscript rscript/model_likelihood.R pygmyhog test5 dstats/pygmyhog.csv 2  5 1000000 &> nohup-test5.out &
