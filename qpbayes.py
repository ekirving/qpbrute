#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Compare all fitted models to each other using Bayes factors from admixture_graph

import sys
import os
import argparse
import glob
import re
import itertools

import pandas as pd

from time import time
from datetime import timedelta
from collections import defaultdict

# use the Pathos library for improved multi-processing
import pathos.multiprocessing as mp

# import the custom modules
from utils import run_cmd

# import all the constants
from consts import *


class QPBayes:

    def __init__(self, geno, snp, ind, prefix, nodes, outgroup, verbose, nthreads):
        """
        Initialise the object attributes
        """
        self.geno_file = geno
        self.snp_file = snp
        self.ind_file = ind
        self.prefix = prefix
        self.nodes = nodes
        self.outgroup = outgroup
        self.verbose = verbose
        self.nthreads = nthreads

        # sanity check the outgroup is not in the node list
        if outgroup in nodes:
            self.nodes.remove(outgroup)

        # find all the PDFs, and extract their graph names
        self.graphs = [re.search(r'a[0-9]-(.+).pdf', pdf).group(1) for pdf in glob.glob('pdf/{}*'.format(prefix))]

        # self.dot_path = 'graphs/{}'.format(prefix)
        self.par_file = "dstats/{}.par".format(prefix)
        self.csv_file = "dstats/{}.csv".format(prefix)
        self.dstat_file = "dstats/{}.log".format(prefix)
        self.tests_file = "dstats/{}.tests".format(prefix)
        self.log_file = '{}.bayes.log'.format(prefix)

        # clean up the log file
        if os.path.exists(self.log_file):
            os.remove(self.log_file)

        # open the log file for writing
        self.log_handle = open(self.log_file, 'a')

    def log(self, message):
        """
        Handle message logging to file/stdout.
        """
        # send message to the log file
        print >> self.log_handle, message
        self.log_handle.flush()

        if self.verbose:
            # echo to stdout
            print message
            sys.stdout.flush()

    def calculate_dstats(self):
        """
        Use `qpDstat` to calculate D-stats for all possible permutations of populations.

        See https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics
        """

        # get all the samples, grouped by population
        samples = defaultdict(list)
        with open(self.ind_file, 'r') as fin:
            for line in fin.readlines():
                sample, gender, population = line.split()
                samples[population].append(sample)

        # compose the list of all 3-way tests (we're doing outgroup D-stats)
        tests = set()
        for x, y, z in itertools.permutations(self.nodes, 3):
            tests.add((self.outgroup, x, y, z))

        # write all the tests to disk
        with open(self.tests_file, 'w') as fout:
            fout.writelines(' '.join(test) + '\n' for test in tests)

        # compose the config settings
        config = [
            "genotypename: {}".format(self.geno_file),
            "snpname:      {}".format(self.snp_file),
            "indivname:    {}".format(self.ind_file),
            "popfilename:  {}".format(self.tests_file),  # Program will run the method for all listed 4-way tests
            "blgsize:      0.005"                        # TODO parameterise
            "f4mode:       YES"                          # TODO f4 statistics not D-stats are computed
        ]

        # the params to be defined in a .par file
        with open(self.par_file, 'w') as fout:
            fout.write("\n".join(config))

        self.log("INFO: There are {:,} D-stat tests to compute for {} populations.".format(len(tests), len(self.nodes)))

        # run qpDstat
        log = run_cmd(["qpDstat", "-p", self.par_file])

        # save the log file
        with open(self.dstat_file, 'w') as fout:
            fout.write(log)

        results = list()
        columns = ['W', 'X', 'Y', 'Z', 'D', 'Z.value']

        # parse the results from the log file
        for line in log.splitlines():
            if 'result:' in line:
                results.append(dict(zip(columns, line.split()[1:7])))

        # convert to DataFrame and save to disk
        pd.DataFrame(results, columns=columns).to_csv(self.csv_file, index=False)

    def calculate_bayes_factors(self):
        """
        Use `admixturegraph` to calculate Bayes factors for all fitted graphs.

        See https://github.com/mailund/admixture_graph
        """
        self.log("INFO: There are {:,} graphs to compute Bayes factors for.".format(len(self.graphs)))

        if self.nthreads > 1:
            # compute the model likelihoods
            pool = mp.ProcessingPool(self.nthreads)
            pool.map(self.model_likelihood, self.graphs)
        else:
            # compute likelihoods without multi-threading
            for graph in self.graphs:
                self.model_likelihood(graph)

    def model_likelihood(self, graph):
        """
        Run the MCMC to calculate the model likelihoods
        """
        run_cmd(["Rscript",
                 "rscript/model_likelihood.R",
                 self.prefix,
                 graph,
                 self.csv_file])

    def find_best_model(self):
        """
        Compare Bayes factors to find the best fitting model.
        """
        run_cmd(["Rscript",
                 "rscript/bayes_factors.R",
                 self.prefix])


def calculate_bayes_factors(geno, snp, ind, prefix, nodes, outgroup, verbose=True, nthreads=CPU_CORES_MAX):
    """
    Find the best fitting graph by calculating the Bayes factors for each model found by QPBrute.
    """
    start = time()

    # instantiate the class
    qpb = QPBayes(geno, snp, ind, prefix, nodes, outgroup, verbose, nthreads)

    # calculate outgroup D-stats for all 3-way permutations of the populations
    qpb.calculate_dstats()

    # calculate Bayes factors for all the fitted graphs, using the pre-computed D-stats
    qpb.calculate_bayes_factors()

    # compare Bayes factors to find the best fitting model
    qpb.find_best_model()

    qpb.log("INFO: Calculating Bayes factors took: {}".format(timedelta(seconds=time() - start)))


if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description='Compare all fitted models to each other using Bayes factors.')
    parser.add_argument("--geno", help="Input genotype file (eigenstrat format)", metavar='example.geno', required=True)
    parser.add_argument("--snp", help="Input snp file (eigenstrat format)", metavar='example.snp', required=True)
    parser.add_argument("--ind", help="Input indiv file (eigenstrat format)", metavar='example.ind', required=True)
    parser.add_argument("--prefix", help="Output prefix", metavar='example', required=True)
    parser.add_argument("--pops", nargs='+', help='List of populations', metavar=('A', 'B'), required=True)
    parser.add_argument("--out", help="Outgroup population", metavar='OUT', required=True)

    argv = parser.parse_args()

    calculate_bayes_factors(argv.geno, argv.snp, argv.ind, argv.prefix, argv.pops, argv.out)
