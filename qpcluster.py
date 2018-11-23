#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import csv
import os
import numpy as np

from itertools import izip
from cStringIO import StringIO

# use the Pathos library for improved multi-processing
import pathos.multiprocessing as mp

# TODO this throws a warning...
# import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

from graph_tool import *
from graph_tool.topology import *


class ClusterQpgraph():

    def __init__(self, graph_names, log_file, dot_path, csv_file, mtx_file, verbose, nthreads):
        """
        Initialise the object attributes
        """
        self.graph_names = graph_names
        self.log_file = log_file
        self.dot_path = dot_path
        self.csv_file = csv_file
        self.mtx_file = mtx_file
        self.verbose = verbose
        self.nthreads = nthreads

        self.graphs = []

        self.verbose = True

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

    @staticmethod
    def parse_dot_file(path):
        """
        The graph-tool library doesn't like the header attributes used by qpGraph, so we need to filter them out
        """
        with open(path, 'r') as fin:
            rows = fin.readlines()

        # exclude lines 2-4, which contain the problematic metadata
        text = "".join(rows[:1] + rows[5:])

        return StringIO(text)

    def calculate_distance(self, args):
        """
        Calculate the similarity distance for two graphs.

        See https://graph-tool.skewed.de/static/doc/topology.html#graph_tool.topology.similarity
        """

        # extract the tuple of arguments
        i, j = args

        # calculate the distance scores between graph pairs (scores are not symmetric; i.e. A->B != B->A)
        d1 = similarity(self.graphs[i], self.graphs[j], distance=True)
        d2 = similarity(self.graphs[j], self.graphs[i], distance=True)

        # enforce symmetry in the matrix by taking the max distance
        dist = max(d1, d2)

        return i, j, dist

    def build_matrix(self):
        """
        Build a symmetrical distance matrix for all graphs.
        """

        # instantiate all the graph objects
        for graph_name in self.graph_names:
            dot_file = self.dot_path + '-{name}.dot'.format(name=graph_name)
            graph = load_graph(self.parse_dot_file(dot_file), fmt='dot')
            self.graphs.append(graph)

        # how many graphs are we comparing
        size = len(self.graph_names)

        # initialise the distance matrix
        dist_matrix = np.zeros([size, size])

        # get all the i,j pairs for one diagonal half
        idxs = [(i, j) for i in range(1, size) for j in range(i)]

        self.log("INFO: Calculating distance matrix for {:,} graph pairs".format(len(idxs)))

        if self.nthreads > 1:
            # we need to buffer the results to use multi-threading
            pool = mp.ProcessingPool(self.nthreads)
            results = pool.map(self.calculate_distance, idxs)
        else:
            # compute distances without multi-threading
            results = []
            for i, j in idxs:
                result = self.calculate_distance((i, j))
                results.append(result)

        # populate the distance matrix
        for i, j, dist in results:
            dist_matrix[i, j] = dist_matrix[j, i] = dist

        # save the matrix
        np.save(self.mtx_file, dist_matrix)

        return dist_matrix

    def get_matrix(self):
        """
        Load the distance matix from file, or build it if necessary.
        """
        try:
            # load the distance matrix from file
            dist_matrix = np.load(self.mtx_file)

            self.log("INFO: Loaded distance matrix from file %s" % self.mtx_file)

        except IOError:
            # file doesn't exist, so build it
            dist_matrix = self.build_matrix()

        return dist_matrix


def cluster_qpgraph(graph_names, dot_path, log_file, pdf_file, csv_file, mtx_file, verbose=False, nthreads=1):
    """
    Compare all fitting graphs and compute the number of clusters.
    """

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    cq = ClusterQpgraph(graph_names, log_file, dot_path, csv_file, mtx_file, verbose, nthreads)

    cq.log("INFO: There are {:,} graphs to compare".format(len(set(graph_names))))

    # get the distance matrix
    dist_matrix = cq.get_matrix()

    cq.log("INFO: Calculating the hierarchical clusters (linkage matrix)")

    # calculate the hierarchical clusters, using Ward's minimum variance method
    # https://en.wikipedia.org/wiki/Ward%27s_method
    Z = linkage(dist_matrix, method='ward')

    # print a dendrogram of the clusters
    pprint_dendrogram(Z, truncate_mode='lastp', p=10, leaf_rotation=90.,
                      leaf_font_size=12., show_contracted=True, pdf=pdf_file)

    cq.log("INFO: Printed hierarchical clustering dendrogram %s" % pdf_file)

    # automatically assign graphs to clusters
    # https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/#Inconsistency-Method
    clusters = fcluster(Z, t=10, criterion='inconsistent', depth=10)

    cq.log("INFO: Found %s clusters using inconsistency criterion (t=%s)" % (len(set(clusters)), 10))

    with open(csv_file, 'wb') as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(['Graph', 'Cluster'])
        for graph, cluster in izip(graph_names, clusters):
            csv_writer.writerow([graph, cluster])

    cq.log("INFO: Saved clusters to file %s" % csv_file)




def pprint_dendrogram(*args, **kwargs):
    """
    Helper function for plotting a dendrogram.

    See https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    """
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)
    pdf = kwargs.pop('pdf', 'clustering.pdf')

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        # plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')

    plt.savefig(pdf, format='pdf')

    return ddata
