#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare all fitted models to each other using Bayes factors from admixture_graph
"""
__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2018"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import argparse
import glob
import itertools
import os
import re
import sys
from collections import defaultdict
from datetime import timedelta
from time import time

import pandas as pd
import pathos.multiprocessing as mp

from qpbrute.consts import (
    MCMC_NUM_BURN,
    CPU_CORES_HIGH,
    MCMC_NUM_CHAINS,
    MCMC_NUM_TEMPS,
    MCMC_NUM_ITERS,
    CPU_CORES_MAX,
)
from qpbrute.utils import run_cmd


class QPBayes:
    def __init__(
        self,
        geno,
        snp,
        ind,
        prefix,
        nodes,
        outgroup,
        chains,
        heated,
        iterations,
        burnin,
        verbose,
        threads,
    ):
        """
        Initialise the object attributes
        """
        self.geno_file = geno
        self.snp_file = snp
        self.ind_file = ind
        self.prefix = prefix
        self.nodes = nodes
        self.outgroup = outgroup
        self.mcmc_chains = int(chains)
        self.mcmc_heated = int(heated)
        self.mcmc_iters = int(float(iterations))
        self.mcmc_burn = int(float(burnin))
        self.verbose = verbose
        self.threads = int(threads)

        # sanity check the outgroup is not in the node list
        if outgroup in nodes:
            self.nodes.remove(outgroup)

        # find all the PDFs, and extract their graph names
        self.graphs = [
            re.search(r"a[0-9]-(.+).pdf", pdf).group(1)
            for pdf in glob.glob("pdf/{}*".format(prefix))
        ]

        # self.dot_path = 'graphs/{}'.format(prefix)
        self.dstat_par = "dstats/{}.par".format(prefix)
        self.dstat_csv = "dstats/{}.csv".format(prefix)
        self.dstat_log = "dstats/{}.log".format(prefix)
        self.dstat_tests = "dstats/{}.tests".format(prefix)
        self.bayes_log = "{}.bayes.log".format(prefix)

        # clean up the log file
        if os.path.exists(self.bayes_log):
            os.remove(self.bayes_log)

        # open the log file for writing
        self.log_handle = open(self.bayes_log, "a")

    def log(self, message):
        """
        Handle message logging to file/stdout.
        """
        # send message to the log file
        print(message, file=self.log_handle)
        self.log_handle.flush()

        if self.verbose:
            # echo to stdout
            print(message)
            sys.stdout.flush()

    def calculate_dstats(self):
        """
        Use `qpDstat` to calculate D-stats for all possible permutations of populations.

        See https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics
        """

        if os.path.isfile(self.dstat_csv):
            # only run once
            return

        # get all the samples, grouped by population
        samples = defaultdict(list)
        with open(self.ind_file, "r") as fin:
            for line in fin.readlines():
                sample, gender, population = line.split()
                samples[population].append(sample)

        # compose the list of all 3-way tests (we're doing outgroup D-stats)
        tests = set()
        for x, y, z in itertools.permutations(self.nodes, 3):
            tests.add((self.outgroup, x, y, z))

        # write all the tests to disk
        with open(self.dstat_tests, "w") as fout:
            fout.writelines(" ".join(test) + "\n" for test in tests)

        # compose the config settings
        config = [
            "genotypename: {}".format(self.geno_file),
            "snpname:      {}".format(self.snp_file),
            "indivname:    {}".format(self.ind_file),
            "popfilename:  {}".format(
                self.dstat_tests
            ),  # Program will run the method for all listed 4-way tests
            "blgsize:      0.005"  # TODO parameterize
            "f4mode:       YES",  # TODO f4 statistics not D-stats are computed
        ]

        # the params to be defined in a .par file
        with open(self.dstat_par, "w") as fout:
            fout.write("\n".join(config))

        self.log(
            "INFO: There are {:,} D-stat tests to compute for {} populations.".format(
                len(tests), len(self.nodes)
            )
        )

        # run qpDstat
        log = run_cmd(["qpDstat", "-p", self.dstat_par])

        # save the log file
        with open(self.dstat_log, "w") as fout:
            fout.write(log)

        results = list()
        columns = ["W", "X", "Y", "Z", "D", "Z.value"]

        # parse the results from the log file
        for line in log.splitlines():
            if "result:" in line:
                results.append(dict(zip(columns, line.split()[1:7])))

        # convert to DataFrame and save to disk
        pd.DataFrame(results, columns=columns).to_csv(self.dstat_csv, index=False)

    def calculate_bayes_factors(self):
        """
        Use `admixturegraph` to calculate Bayes factors for all fitted graphs.

        See https://github.com/mailund/admixture_graph
        """
        self.log(
            "INFO: There are {:,} graphs to compute Bayes factors for.".format(
                len(self.graphs)
            )
        )

        if self.threads > 1:
            # compute the model likelihoods
            pool = mp.ProcessingPool(self.threads)
            pool.map(self.model_likelihood, self.graphs)
        else:
            # compute likelihoods without multi-threading
            for graph in self.graphs:
                self.model_likelihood(graph)

    def model_likelihood(self, graph):
        """
        Run the MCMC to calculate the model likelihoods
        """
        log_file = "bayes/{}-{}-likelihood.log".format(self.prefix, graph)

        if not os.path.isfile("bayes/{}-{}-burn-gelman.pdf".format(self.prefix, graph)):
            # only run once
            run_cmd(
                [
                    "Rscript",
                    "rscript/model_likelihood.R",
                    self.prefix,
                    graph,
                    self.dstat_csv,
                    self.mcmc_chains,
                    self.mcmc_heated,
                    self.mcmc_iters,
                    self.mcmc_burn,
                ],
                env={"OMP_NUM_THREADS": "1"},
                stdout=open(log_file, "w"),
            )

        self.log("INFO: Bayes factor done for graph {}".format(graph))

    def find_best_model(self):
        """
        Compare Bayes factors to find the best fitting model.
        """
        log_file = "bayes/{}-bayes.log".format(self.prefix)

        run_cmd(
            ["Rscript", "rscript/bayes_factors.R", self.prefix, MCMC_NUM_BURN],
            stdout=open(log_file, "w"),
        )


def calculate_bayes_factors(
    geno,
    snp,
    ind,
    prefix,
    nodes,
    outgroup,
    chains,
    heated,
    iterations,
    burnin,
    verbose=True,
    threads=CPU_CORES_HIGH,
):
    """
    Find the best fitting graph by calculating the Bayes factors for each model found by QPBrute.
    """
    start = time()

    # instantiate the class
    qpb = QPBayes(
        geno,
        snp,
        ind,
        prefix,
        nodes,
        outgroup,
        chains,
        heated,
        iterations,
        burnin,
        verbose,
        threads,
    )

    if burnin >= iterations:
        raise RuntimeError(
            "ERROR: MCMC burn in must be less than the number of iterations"
        )

    # calculate outgroup D-stats for all 3-way permutations of the populations
    qpb.calculate_dstats()

    # calculate Bayes factors for all the fitted graphs, using the pre-computed D-stats
    qpb.calculate_bayes_factors()

    # compare Bayes factors to find the best fitting model
    qpb.find_best_model()

    qpb.log(
        "INFO: Calculating Bayes factors took: {}".format(
            timedelta(seconds=time() - start)
        )
    )


def qpbayes():
    # parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Compare all fitted models to each other using Bayes factors."
    )
    parser.add_argument(
        "--geno",
        help="Input genotype file (eigenstrat format)",
        metavar="example.geno",
        required=True,
    )
    parser.add_argument(
        "--snp",
        help="Input snp file (eigenstrat format)",
        metavar="example.snp",
        required=True,
    )
    parser.add_argument(
        "--ind",
        help="Input indiv file (eigenstrat format)",
        metavar="example.ind",
        required=True,
    )
    parser.add_argument(
        "--prefix", help="Output prefix", metavar="example", required=True
    )
    parser.add_argument(
        "--pops",
        nargs="+",
        help="List of populations",
        metavar=("A", "B"),
        required=True,
    )
    parser.add_argument(
        "--out", help="Outgroup population", metavar="OUT", required=True
    )
    parser.add_argument(
        "--chains",
        help="Number of replicate MCMC chains to run (default: %s)" % MCMC_NUM_CHAINS,
        metavar="NUM",
        default=MCMC_NUM_CHAINS,
    )
    parser.add_argument(
        "--heated",
        help="Number of heated  chains per replicate (default: %s)" % MCMC_NUM_TEMPS,
        metavar="NUM",
        default=MCMC_NUM_TEMPS,
    )
    parser.add_argument(
        "--iterations",
        help="Number of MCMC iterations per chain (default: %.1e)" % MCMC_NUM_ITERS,
        metavar="NUM",
        default=MCMC_NUM_ITERS,
    )
    parser.add_argument(
        "--burnin",
        help="Number of MCMC iterations to burnin (default: %.1e)" % MCMC_NUM_BURN,
        metavar="NUM",
        default=MCMC_NUM_BURN,
    )
    parser.add_argument(
        "--threads",
        help="Number of threads to use (default: %s)" % CPU_CORES_MAX,
        metavar="NUM",
        default=CPU_CORES_MAX,
    )

    argv = parser.parse_args()

    calculate_bayes_factors(
        argv.geno,
        argv.snp,
        argv.ind,
        argv.prefix,
        argv.pops,
        argv.out,
        argv.chains,
        argv.heated,
        argv.iterations,
        argv.burnin,
        threads=argv.threads,
    )


if __name__ == "__main__":
    qpbayes()
