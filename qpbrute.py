#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Build all possible graphs using a stepwise addition order algorithm, bounded by |Z| < 3

import re
import sys
import copy
import itertools
import hashlib
import os
import random
import argparse

from time import time
from datetime import timedelta
from collections import OrderedDict, defaultdict
from cStringIO import StringIO
from Bio import Phylo

import xml.etree.ElementTree as ElemTree

# use the Pathos library for improved multi-processing
import pathos.multiprocessing as mp

# import the custom modules
from utils import run_cmd, pprint_qpgraph

# import all the constants
from consts import *


class QPBrute:
    def __init__(
        self,
        par_file,
        log_file,
        dot_path,
        pdf_path,
        nodes,
        outgroup,
        exhaustive,
        verbose,
        threads,
        skeleton,
        qpgraph,
        no_admix,
        max_outlier,
        print_offset,
    ):
        """
        Initialise the object attributes 
        """
        self.par_file = par_file
        self.dot_path = dot_path
        self.pdf_path = pdf_path
        self.verbose = verbose
        self.threads = int(threads)
        self.skeleton = skeleton
        self.qpgraph = qpgraph
        self.no_admix = no_admix
        self.max_outlier = int(max_outlier)
        self.print_offset = int(print_offset)

        # should we try all possible graphs, or should we stop when we find something reasonable
        self.exhaustive_search = exhaustive

        # open the file for writing
        self.log_handle = open(log_file, "a")

        if outgroup in nodes:
            nodes.remove(outgroup)

        self.nodes = nodes
        self.outgroup = outgroup

        self.root_node = ROOT_NODE
        self.problem_nodes = []
        self.tested_graphs = set()
        self.solutions = set()

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

    def recurse_tree(self, root_tree, new_tag, remaining, depth=0):
        """
        Permute all possible new trees, by adding the new node to all branches.
    
        If no resulting tree passes the outlier threshold then try adding the node to all possible pairs of branches.
        """
        new_trees = []

        # get all the nodes in the tree (skip the outgroup)
        target_nodes = [
            node for node in root_tree.findall(".//*") if node.tag != self.outgroup
        ]

        # add the new node to every branch in the tree
        for target_node in target_nodes:

            # clone the current tree and add the new node
            new_tree = copy.deepcopy(root_tree)
            self.insert_node(new_tree, target_node, new_tag)
            new_trees.append(new_tree)

        # test all the trees
        results = self.test_trees(new_trees, depth)

        # process the results
        node_placed = self.check_results(results, remaining, depth)

        # test all the admixture possibilities
        if not node_placed and not self.no_admix:

            admix_trees = []

            # permute all the two parent admixture possibilities
            pairs = list(itertools.combinations(target_nodes, 2))

            for target1, target2 in pairs:
                # skip duplicate targets (this happens when there is already an admixture node in the tree)
                if target1.tag == target2.tag:
                    continue

                # clone the current tree
                new_tree = copy.deepcopy(root_tree)

                # make a new intermediate node
                admix_label = self.new_label(new_tree, admix=True)

                # add two admix nodes as the children of both targets
                admix_nodes = [
                    self.insert_node(
                        new_tree,
                        target1,
                        admix_label,
                        attrs={"internal": "1", "admix": "1", "side": "l"},
                    ),
                    self.insert_node(
                        new_tree,
                        target2,
                        admix_label,
                        attrs={"internal": "1", "admix": "1", "side": "r"},
                    ),
                ]

                # choose the actual parent based on the sort order of the tag name (needed for unique tree hashing)
                admix_node = (
                    admix_nodes[0] if target1.tag < target2.tag else admix_nodes[1]
                )

                # add the new node as the child of the preferred admix node
                self.insert_node(new_tree, admix_node, new_tag, append=True)

                admix_trees.append(new_tree)

            # test all the admixture trees
            results = self.test_trees(admix_trees, depth)

            # process the results
            node_placed = self.check_results(results, remaining, depth)

        if not node_placed:

            # we could not place the node via either method :(
            if (
                new_tag not in self.problem_nodes
                and remaining
                and not self.exhaustive_search
            ):
                self.log("WARNING: Unable to place node '%s' at this time." % new_tag)

                self.problem_nodes.append(new_tag)

                # add the problem node to end of the list, as we may be able to add it later on
                remaining.append(new_tag)

                # try and add the other nodes
                self.recurse_tree(root_tree, remaining[0], remaining[1:], depth)

            else:
                raise NodeUnplaceable(
                    "ERROR: Cannot place node '%s' in the graph." % new_tag
                )

    def test_trees(self, new_trees, depth):
        """
        Run qpGraph on a list of trees
        """
        if self.threads > 1:
            # we need to buffer the results to use multi-threading
            pool = mp.ProcessingPool(self.threads)
            results = pool.map(
                self.run_qpgraph, itertools.izip(new_trees, itertools.repeat(depth))
            )
        else:
            # test the trees without multi-threading
            results = []
            for new_tree in new_trees:
                result = self.run_qpgraph((new_tree, depth))
                results.append(result)

        return results

    def check_results(self, results, remaining, depth):
        """
        Check the results from qpGraph
        """
        # were we able to place the new node
        placed_node = False

        for new_tree, outliers, graph_name in results:

            # add this graph to the list of those we've tested
            self.tested_graphs.add(graph_name)

            # did our new trees pass the threshold
            if len(outliers) <= self.max_outlier:

                # recursively add any remaining nodes
                if remaining:
                    self.recurse_tree(new_tree, remaining[0], remaining[1:], depth + 1)
                else:
                    if len(outliers) == 0:
                        self.log(
                            "SUCCESS: Placed all nodes on a graph without outliers!"
                        )
                    else:
                        self.log(
                            "INFO: Placed all nodes on a graph with {} outlier(s).".format(
                                len(outliers)
                            )
                        )

                    # add this graph to the list of solutions
                    self.solutions.add(graph_name)

                # we successfully placed the new node!
                placed_node = True

        return placed_node

    def insert_node(self, new_tree, target_node, new_tag, attrs=None, append=False):
        """
        Helper function to add a new node on the branch leading to the target node.
        """
        # get the target node in the new tree
        target_xpath = ".//" + target_node.tag

        if target_node.get("side"):
            target_xpath += '[@side="%s"]' % target_node.get("side")

        target_node = new_tree.find(target_xpath)

        if append:
            # append the new node directly to the target
            parent_node = target_node

        else:
            # get the parent of the target
            parent_node = new_tree.find(target_xpath + "/..")

            # does the target node have a sibling
            if len(parent_node) > 1:
                label = self.new_label(new_tree)
                parent_node.remove(target_node)

                # add an intermediate node, to act as the parent for the new node
                parent_node = ElemTree.SubElement(parent_node, label)
                parent_node.set("internal", "1")

                # re add the target node
                # ElemTree.SubElement(parent_node, target_node.tag)
                parent_node.append(target_node)

        # add the new node as a sibling to the target
        new_node = ElemTree.SubElement(parent_node, new_tag)

        if attrs:
            # add any node attributes
            for key, value in attrs.iteritems():
                new_node.set(key, value)

        return new_node

    def run_qpgraph(self, params):
        """
        Run qpGraph on the given tree 
        """
        # extract the tuple of arguments
        new_tree, depth = params

        # convert the tree to newick format
        newick = self.print_newick_tree(new_tree)

        # get unique names for the output files
        graph_name = self.hash_text(newick)
        grp_file = self.dot_path + "-{name}.graph".format(name=graph_name)
        dot_file = self.dot_path + "-{name}.dot".format(name=graph_name)
        log_file = self.dot_path + "-{name}.log".format(name=graph_name)
        xml_file = self.dot_path + "-{name}.xml".format(name=graph_name)

        try:
            # if the log file exists then we've run the analysis already
            with open(log_file, "r") as fin:
                log = fin.read()

        except IOError:
            # save the xml file
            new_tree.write(
                xml_file
            )  # TODO find a better internal representation for this

            # convert the tree to qpGraph format
            graph = self.export_qpgraph(new_tree)

            # save the graph file
            with open(grp_file, "w") as fout:
                fout.write(graph)

            # run qpGraph
            log = run_cmd(
                ["qpGraph", "-p", self.par_file, "-g", grp_file, "-d", dot_file],
                env={"OMP_NUM_THREADS": "1"},
            )

            # save the log file
            with open(log_file, "w") as fout:
                fout.write(log)

        # parse the log and extract the outliers
        outliers, worst_fstat = self.extract_outliers(log.splitlines())

        # count the leaf nodes
        all_nodes = new_tree.findall(".//*")
        num_nodes = len([node for node in all_nodes if node.get("internal") != "1"])
        num_admix = len([node for node in all_nodes if node.get("admix") == "1"]) / 2
        num_outliers = len(outliers)

        # only print PDFs for graphs that pass the threshold
        if num_outliers <= self.max_outlier and num_nodes > (
            len(self.nodes) - self.print_offset
        ):

            # embed some useful metadata info in the PDF name
            pdf_file = self.pdf_path + "-n{nodes}-o{out}-a{admix}-{name}.pdf".format(
                nodes=num_nodes, out=num_outliers, admix=num_admix, name=graph_name
            )

            # pretty print the qpGraph dot file
            pprint_qpgraph(dot_file, pdf_file)

        # output some summary stats
        self.log(
            "{padding}{tree} \tnodes={nodes}\t admix={admix}\t outliers={out}\t worst={worst}\t {name}".format(
                padding="  " * depth,
                name=graph_name,
                tree=newick.ljust(80),
                nodes=num_nodes,
                admix=num_admix,
                out=len(outliers),
                worst=worst_fstat[-1],
            )
        )

        return new_tree, outliers, graph_name

    @staticmethod
    def extract_outliers(log):
        """
        Parse the log file and extract the outliers
        """
        outliers = []
        read_log = False
        worst_fstat = []

        header = ["Fit", "Obs", "Diff", "Std.", "error", "Z"]

        for line in log:
            if "outliers" in line:
                read_log = True
                continue
            elif "worst f-stat" in line:
                worst_fstat = line.split()
                read_log = False
                continue

            if read_log and len(line.strip()) > 0:
                outlier = line.split()
                if outlier != header:
                    # save all the outliers
                    outliers.append(outlier)

        return outliers, worst_fstat

    def export_qpgraph(self, root_tree):
        """
        Convert the ElementTree into qpGraph format
        """
        # clone the tree because this process is destructive
        local_tree = copy.deepcopy(root_tree)

        graph = "root\t{root}\n".format(root=self.root_node)

        # get all the nodes
        nodes = local_tree.findall(".//*")

        for node in nodes:
            # skip non-leaf nodes
            if len(node) == 0 and node.get("admix") != "1":
                graph += "label\t{node}\t{node}\n".format(node=node.tag)

        # build the list of edges
        graph += self.export_qpgraph_node(local_tree)

        return graph

    def export_qpgraph_node(self, root_tree, parent_node=None):
        """
        Recursively export all the edges in the graph
        """
        graph = ""

        if parent_node is None:
            parent_node = root_tree.getroot()

        for child_node in list(parent_node):

            if child_node.get("printed") != "1":

                # is this an admixture node or a normal node
                matches = root_tree.findall(".//" + child_node.tag + "/..")

                if len(matches) > 1:
                    # admixture branch
                    parent1, parent2 = matches

                    middle1 = child_node.tag + "a"
                    middle2 = child_node.tag + "b"
                    code1 = self.hash_text(middle1)
                    code2 = self.hash_text(middle2)

                    # don't admix from a bifurcating node; intermediate nodes to accommodate drift
                    graph += "edge\t{code}\t{parent}\t{middle}\n".format(
                        code=code1, parent=parent1.tag, middle=middle1
                    )
                    graph += "edge\t{code}\t{parent}\t{middle}\n".format(
                        code=code2, parent=parent2.tag, middle=middle2
                    )

                    # now admix from the two middle nodes
                    graph += "admix\t{child}\t{parent1}\t{parent2}\t50\t50\n".format(
                        parent1=middle1, parent2=middle2, child=child_node.tag
                    )

                    # flag both nodes so we don't export them twice
                    parent1.find(child_node.tag).set("printed", "1")
                    parent2.find(child_node.tag).set("printed", "1")

                else:
                    # regular branch
                    code = self.hash_text(child_node.tag)
                    graph += "edge\t{code}\t{parent}\t{child}\n".format(
                        code=code, parent=parent_node.tag, child=child_node.tag
                    )

            # leaf nodes
            if len(child_node) > 0:
                # now convert the children
                graph += self.export_qpgraph_node(root_tree, child_node)

        return graph

    def import_qpgraph(self, graph_file):
        """
        Convert a qpGraph model into ElementTree format
        """
        root, labels, edges, admixes = self.parse_qpgraph(graph_file)

        added_a = []
        added_i = []

        # construct the ElementTree
        root_node = ElemTree.Element(self.root_node)
        root_tree = ElemTree.ElementTree(root_node)

        # recursively insert the nodes
        self.insert_qpgraph_nodes(
            root, root_node, labels, edges, admixes, added_a, added_i
        )

        return root_tree

    def insert_qpgraph_nodes(
        self, parent, parent_elem, labels, edges, admixes, added_a, added_i
    ):
        """
        Recursively insert the qpgraph nodes
        """
        for child in edges[parent]:
            # add all the child nodes
            if child in labels:
                if child in edges:
                    raise RuntimeError(
                        "ERROR: Sampled populations cannot have direct descendents ({} -> {})".format(
                            child, edges[child]
                        )
                    )
                # leaf node
                ElemTree.SubElement(parent_elem, labels[child])
            else:
                if child in admixes:
                    # admixture node
                    if child not in added_a:
                        added_a.append(child)
                        side = "l"
                    else:
                        side = "r"

                    name = "a{}".format(added_a.index(child) + 1)
                    attrs = [("internal", "1"), ("admix", "1"), ("side", side)]
                else:
                    # internal node
                    if len(edges[child]) == 1:
                        # skip drift branches in the qpGraph model, as we add these automatically on export
                        self.insert_qpgraph_nodes(
                            child, parent_elem, labels, edges, admixes, added_a, added_i
                        )
                        continue

                    added_i.append(child)
                    name = "n{}".format(added_i.index(child) + 1)
                    attrs = [("internal", "1")]

                child_elem = ElemTree.SubElement(parent_elem, name)
                for key, val in attrs:
                    child_elem.set(key, val)

                # don't recurse on duplicated admixture nodes
                if child_elem.get("side") != "r":
                    self.insert_qpgraph_nodes(
                        child, child_elem, labels, edges, admixes, added_a, added_i
                    )

    @staticmethod
    def parse_qpgraph(graph_file):
        """
        Parse the elements of a qpGraph model
        """
        labels = {}
        edges = defaultdict(list)
        admixes = []

        # parse the contents of the graph file
        with open(graph_file, "r") as fin:
            for line in fin:
                # root node
                m = re.match(r"^root\s+(\w+)$", line)
                if m is not None:
                    root = m.group(1)
                    continue

                # tip labels
                m = re.match(r"^label\s+(\w+)\s+(\w+)$", line)
                if m is not None:
                    label, alias = m.groups()
                    labels[alias] = label
                    continue

                # normal edges
                m = re.match(r"^edge\s+(\w+)\s+(\w+)\s+(\w+)$", line)
                if m is not None:
                    _, parent, child = m.groups()
                    edges[parent].append(child)
                    continue

                # admixture edges
                m = re.match(
                    r"^admix\s+(\w+)\s+(\w+)\s+(\w+)(\s+(\d+)\s+(\d+))?$", line
                )
                if m is not None:
                    child, parent1, parent2 = m.groups()[0:3]
                    edges[parent1].append(child)
                    edges[parent2].append(child)
                    admixes.append(child)
                    continue

                if line != "":
                    raise RuntimeError(
                        "Invalid line found in {}\n{}".format(graph_file, line)
                    )

        return root, labels, edges, admixes

    # TODO use the birthday problem formula to calculate the chance of collisions
    #  https://en.wikipedia.org/wiki/Birthday_problem
    @staticmethod
    def hash_text(text, length=12):
        """
        Generate a unique key by hashing a string
        """
        return hashlib.sha1(text).hexdigest()[0:length]

    @staticmethod
    def new_label(root_tree, admix=False):
        """
        Return a new label for a node
        """
        all_nodes = root_tree.findall(".//*")

        if admix:
            num = len([node for node in all_nodes if node.get("admix") == "1"]) / 2
        else:
            num = len(
                [
                    node
                    for node in all_nodes
                    if node.get("internal") == "1" and node.get("admix") != "1"
                ]
            )

        return "{pre}{num}".format(pre="a" if admix else "n", num=num + 1)

    def print_newick_tree(self, root_tee):
        """
        Convert an ElementTree into a ladderized Newick tree.
        """
        # TODO output as a series of newick trees, not as a hybrid (and use weighted pairwise distance to do clustering)
        # TODO use the hash of all the newick trees, to fix the duplicates problem
        newick = self.export_newick_tree(root_tee.getroot())

        # load into Phylo so we can sort the tree (i.e. ladderize)
        tree = Phylo.read(StringIO(newick), "newick")
        tree.ladderize()

        # export the tree back to a string
        fout = StringIO()
        Phylo.write(tree, fout, "newick")
        newick = fout.getvalue()

        # remove the branch lengths
        newick = newick.replace(":0.00000", "").strip()

        # get the order of admix nodes in the tree
        order = list(OrderedDict.fromkeys(re.findall("a[0-9]+", newick)))

        # normalise the node numbering
        for i, old in enumerate(order):
            newick = newick.replace(old, "n%s" % (i + 1))

        # replace n0 with a0 (to preserve the existing cache)
        newick = re.sub(r"n([0-9]+)", r"a\1", newick)

        return newick

    @staticmethod
    def export_newick_tree(parent_node):
        """
        Convert an ElementTree tree into Newick format
        """
        if len(parent_node) == 0:
            return parent_node.tag
        else:
            children = [
                (child_node.tag, QPBrute.export_newick_tree(child_node))
                for child_node in parent_node
            ]
            children.sort()
            tag_name = (
                ""
                if re.match("n[0-9]+|" + ROOT_NODE, parent_node.tag)
                else parent_node.tag
            )
            return "(" + ",".join(node for tag, node in children) + ")%s" % tag_name

    def find_graph(self):
        """
        Build and test all possible trees and graphs
        """
        self.log("INFO: Starting list %s" % self.nodes)

        if self.skeleton or self.qpgraph:
            # load the skeleton tree / graph
            root_tree = (
                ElemTree.parse(self.skeleton)
                if self.skeleton
                else self.import_qpgraph(self.qpgraph)
            )

            # recursively add all the other nodes
            self.recurse_tree(
                root_tree, self.nodes[0], self.nodes[1:] if len(self.nodes) > 1 else []
            )
        else:
            # setup a simple 2-node tree
            root_node = ElemTree.Element(self.root_node)
            root_tree = ElemTree.ElementTree(root_node)

            ElemTree.SubElement(root_node, self.outgroup)
            ElemTree.SubElement(root_node, self.nodes[0])

            # recursively add all the other nodes
            self.recurse_tree(root_tree, self.nodes[1], self.nodes[2:])


class NodeUnplaceable(Exception):
    """
    Node cannot be placed in the graph without exceeding outlier threshold
    """

    pass


def permute_qpgraph(
    par_file,
    prefix,
    nodes,
    outgroup,
    exhaustive=True,
    verbose=True,
    threads=CPU_CORES_MAX,
    skeleton=None,
    qpgraph=None,
    no_admix=False,
    max_outlier=0,
    print_offset=0,
):
    """
    Find the best fitting graph for a given set of nodes, by permuting all possible graphs.
    """

    log_file = prefix + ".log"
    dot_path = "graphs/" + prefix
    pdf_path = "pdf/" + prefix

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    pq = QPBrute(
        par_file,
        log_file,
        dot_path,
        pdf_path,
        nodes,
        outgroup,
        exhaustive,
        verbose,
        threads,
        skeleton,
        qpgraph,
        no_admix,
        max_outlier,
        print_offset,
    )

    if len(nodes) > 7:
        raise RuntimeError(
            "ERROR: Maximum number of populations is 7 (n={})".format(len(nodes))
        )

    if ROOT_NODE in nodes:
        raise RuntimeError(
            "ERROR: '{}' is a reserved name for the root of the graph".format(ROOT_NODE)
        )

    # get all the permutations of possible node orders
    all_nodes_perms = list(itertools.permutations(nodes, len(nodes)))

    # randomise the list of starting orders
    random.shuffle(all_nodes_perms)

    pq.log(
        "INFO: There are {:,} possible starting orders for the given nodes.".format(
            len(all_nodes_perms)
        )
    )
    pq.log(
        "INFO: Performing %s search."
        % ("an exhaustive" if pq.exhaustive_search else "a heuristic")
    )

    # remove the starting set (so we don't do it twice)
    all_nodes_perms.remove(tuple(nodes))

    # keep looping until we find a solution, or until we've exhausted all possible starting orders
    while not pq.solutions or pq.exhaustive_search:

        try:
            # find the best fitting graph for this starting order
            pq.find_graph()

        except NodeUnplaceable as error:
            # log the error
            pq.log(error)

        try:
            # try starting with a different node order
            pq.nodes = list(all_nodes_perms.pop())

        except IndexError:
            # we've run out of node orders to try
            if not pq.solutions:
                pq.log(
                    "ERROR: Cannot resolve the graph from any permutation of the given nodes."
                )

            break

    pq.log(
        "FINISHED: Found {:,} unique solution(s) from a total of {:,} unique graphs!".format(
            len(pq.solutions), len(pq.tested_graphs)
        )
    )

    return pq.solutions


if __name__ == "__main__":

    start = time()

    # parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Test all possible qpgraph models for a given set of populations."
    )
    parser.add_argument(
        "--par", help="Param file for qpGraph", metavar="example.par", required=True
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
        "--skeleton", help="A skeleton model to add all other nodes to", metavar="FILE"
    )
    parser.add_argument(
        "--qpgraph", help="A qpgraph model to add all other nodes to", metavar="FILE"
    )
    parser.add_argument(
        "--no_admix",
        help="Do not insert new nodes with admixture branches",
        action="store_true",
    )
    parser.add_argument(
        "--max_outlier",
        help="How many outliers are permitted before pruning a graph",
        metavar="NUM",
        default=0,
    )
    parser.add_argument(
        "--print_offset",
        help="Print PDFs for graphs with (N - offset) nodes",
        metavar="NUM",
        default=0,
    )
    parser.add_argument(
        "--threads",
        help="Number of threads to use (default: %s)" % CPU_CORES_MAX,
        metavar="NUM",
        default=CPU_CORES_MAX,
    )

    argv = parser.parse_args()

    # test all the models
    permute_qpgraph(
        argv.par,
        argv.prefix,
        argv.pops,
        argv.out,
        threads=argv.threads,
        skeleton=argv.skeleton,
        qpgraph=argv.qpgraph,
        no_admix=argv.no_admix,
        max_outlier=argv.max_outlier,
        print_offset=argv.print_offset,
    )

    print "INFO: Permute execution took: %s" % timedelta(seconds=time() - start)
