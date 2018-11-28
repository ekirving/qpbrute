#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Compare all fitted models to each other using Bayes factors from admixture_graph

import sys
import os
import argparse
import glob
import re

from time import time
from datetime import timedelta

# import the custom modules
from utils import run_cmd

# import all the constants
from consts import *


class QPBayes:

    def __init__(self, graph_names, log_file, dot_path, verbose, nthreads):
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


def calculate_bayes_factors(graph_names, prefix, verbose=True, nthreads=CPU_CORES_MAX):
    """
    Find the best fitting graph by calculating the Bayes factors for each model found by QPBrute.
    """
    log_file = '{}.bayes.log'.format(prefix)
    dot_path = 'graphs/{}'.format(prefix)

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    qpb = QPBayes(graph_names, log_file, dot_path, verbose, nthreads)

    qpb.log("INFO: There are {:,} graphs to compute Bayes factors for".format(len(set(graph_names))))

    pass


if __name__ == "__main__":
    start = time()

    # parse the command line arguments
    parser = argparse.ArgumentParser(description='Compare all fitted models to each other using Bayes factors.')
    parser.add_argument("--prefix", help="Output prefix", metavar='example', required=True)

    argv = parser.parse_args()

    # find all the PDFs, and extract the graph names
    files = glob.glob('pdf/{}*'.format(argv.prefix))
    graphs = [re.search(r'a[0-9]-(.+).pdf', filename).group(1) for filename in files]

    calculate_bayes_factors(graphs, argv.prefix)

    print "INFO: Calculating Bayes factors took: %s" % timedelta(seconds=time() - start)
