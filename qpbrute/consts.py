#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Static constants for qpBrute
"""
__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2018"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import multiprocessing

# the maximum number of populations that are computationally tractatable
MAX_POPULATIONS = 7

# how many CPU cores does this machine have
TOTAL_CORES = multiprocessing.cpu_count()

# no single worker should use more than 30% of the available cores
CPU_CORES_MED = int(TOTAL_CORES * 0.3)  # 30%
CPU_CORES_HIGH = int(TOTAL_CORES * 0.5)  # 50%
CPU_CORES_MAX = TOTAL_CORES

# default colours for print the graphs
COLOURS = [
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
    "#b15928",
    "#a6cee3",
    "#b2df8a",
    "#fb9a99",
    "#fdbf6f",
    "#cab2d6",
    "#ffff99",
]

# config settings for the admixturegraph MCMC
MCMC_NUM_CHAINS = 2
MCMC_NUM_TEMPS = 5
MCMC_NUM_ITERS = int(2e6)
MCMC_NUM_BURN = int(1.1e6)

ROOT_NODE = "R"

FOLDERS = ["bayes", "dstats", "graphs", "pdf"]  # "clusters",
