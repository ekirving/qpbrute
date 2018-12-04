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

# import the custom modules
from utils import run_cmd

# import all the constants
from consts import *


class QPBayes:

    def __init__(self, graph_names, log_file, dot_path, nodes, outgroup, verbose, nthreads):
        """
        Initialise the object attributes
        """
        self.graph_names = graph_names
        self.log_file = log_file
        self.dot_path = dot_path
        self.verbose = verbose
        self.nthreads = nthreads

        # open the file for writing
        self.log_handle = open(log_file, 'a')

        if outgroup in nodes:
            nodes.remove(outgroup)

        self.nodes = nodes
        self.outgroup = outgroup

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
        Use AdmixTools to calculate D-stats for all possible permutations of populations.

        See https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics
        """

        geno_file = "pygmyhog/autosomes.filtered.geno"
        snp_file = "pygmyhog/autosomes.filtered.snp"
        ind_file = "pygmyhog/autosomes.filtered.ind"
        par_file = "pygmyhog/autosomes.filtered.dstat.par"

        tests_file = 'pygmyhog/autosomes.filtered.tests'
        dstat_file = 'pygmyhog.dstat'
        csv_file = 'pygmyhog.csv'

        # get all the samples, grouped by population
        samples = defaultdict(list)
        with open(ind_file, 'r') as fin:
            for line in fin.readlines():
                sample, gender, population = line.split()
                samples[population].append(sample)

        tests = set()

        # compose the list of all 4-way tests
        for x, y, z in itertools.permutations(self.nodes, 3):
            tests.add((self.outgroup, x, y, z))

        # write all the tests to disk
        with open(tests_file, 'w') as fout:
            fout.writelines(' '.join(test) + '\n' for test in tests)

        # compose the config settings
        config = [
            "genotypename: {}".format(geno_file),
            "snpname:      {}".format(snp_file),
            "indivname:    {}".format(ind_file),
            "popfilename:  {}".format(tests_file),  # Program will run the method for all listed 4-way tests
            "blgsize:      0.005"                   # TODO parameterise
        ]

        # the params to be defined in a .par file
        with open(par_file, 'w') as fout:
            fout.write("\n".join(config))

        self.log("INFO: There are {:,} D-stat tests to compute".format(len(tests)))

        # run qpDstat
        log = run_cmd(["qpDstat", "-p", par_file])

        # save the log file
        with open(dstat_file, 'w') as fout:
            fout.write(log)

        results = list()
        columns = ['W', 'X', 'Y', 'Z', 'D', 'Z.value']

        # parse the results from the log file
        for line in log.splitlines():
            if 'result:' in line:
                results.append(dict(zip(columns, line.split()[1:7])))

        # convert to DataFrame and save to disk
        pd.DataFrame(results, columns=columns).to_csv(csv_file, index=False)


def calculate_bayes_factors(graph_names, prefix, nodes, outgroup, verbose=True, nthreads=CPU_CORES_MAX):
    """
    Find the best fitting graph by calculating the Bayes factors for each model found by QPBrute.
    """
    start = time()

    log_file = '{}.bayes.log'.format(prefix)
    dot_path = 'graphs/{}'.format(prefix)

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    qpb = QPBayes(graph_names, log_file, dot_path, nodes, outgroup, verbose, nthreads)

    qpb.log("INFO: There are {:,} graphs to compute Bayes factors for".format(len(set(graph_names))))

    # get all the D-stat permutations
    qpb.calculate_dstats()

    qpb.log("INFO: Calculating Bayes factors took: %s" % timedelta(seconds=time() - start))


if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description='Compare all fitted models to each other using Bayes factors.')
    parser.add_argument("--prefix", help="Output prefix", metavar='example', required=True)
    parser.add_argument("--pops", nargs='+', help='List of populations', metavar=('A', 'B'), required=True)
    parser.add_argument("--out", help="Outgroup population", metavar='OUT', required=True)

    argv = parser.parse_args()

    # find all the PDFs, and extract the graph names
    files = glob.glob('pdf/{}*'.format(argv.prefix))
    graphs = [re.search(r'a[0-9]-(.+).pdf', filename).group(1) for filename in files]

    calculate_bayes_factors(graphs, argv.prefix, argv.pops, argv.out)
