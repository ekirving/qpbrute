#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shared utilities for qpBrute
"""
__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2018"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import csv
import os
import re
import subprocess

import graphviz

from qpbrute.consts import COLOURS


def run_cmd(cmd, shell=False, stdout=None, stderr=None, env=None, verbose=False):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param stdout: File handle to redirect stdout
    :param stderr: File handle to redirect stderr
    :param env: Optional dictionary of local environment settings
    :param verbose: Print the command being run
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # print the command
    if verbose:
        print(" ".join(cmd))

    # handle custom environment variables
    local_env = os.environ
    if env:
        local_env.update(env)

    stdout = subprocess.PIPE if not stdout else stdout
    stderr = subprocess.PIPE if not stderr else stderr

    # run the command
    proc = subprocess.Popen(
        cmd, shell=shell, stdout=stdout, stderr=stderr, env=local_env
    )

    # fetch any output and error
    (out, err) = proc.communicate()

    # bail if something went wrong
    if proc.returncode != 0:

        # decode return codes
        if proc.returncode == 139:
            err = f"Segmentation fault (core dumped) {err}"

        raise RuntimeError(err)

    return out


def trim_ext(full_path, n=1):
    return ".".join(full_path.split(".")[:-n])


def pprint_qpgraph(dot_file, pdf_file):
    """
    Pretty print a qpGraph dot file using the graphviz library.
    """

    # extract the body contents from the dot file
    with open(dot_file, "rU") as fin:
        body = re.search("{(.+)}", fin.read(), re.DOTALL).group(1).strip().split("\n")

    # make a new direct graph
    dot = graphviz.Digraph(filename=trim_ext(pdf_file), body=body, format="pdf")

    # remove the messy graph label
    dot.attr("graph", label="")

    # set Node defaults
    dot.node_attr["shape"] = "point"
    dot.node_attr["fontname"] = "arial"
    dot.node_attr["fontsize"] = "11"

    # set Edge defaults
    dot.edge_attr["arrowhead"] = "vee"
    dot.edge_attr["fontcolor"] = "#838b8b"  # grey
    dot.edge_attr["fontname"] = "arial"
    dot.edge_attr["fontsize"] = "11"

    nodes = []

    # extract the leaf nodes from the body of the graph
    for line in dot:
        match = re.search(r"^ *([a-z_0-9.]+) +\[", line, re.IGNORECASE)
        if match:
            nodes.append(match.group(1))

    # sort the nodes
    nodes.sort()

    try:
        # TODO make node colours a CLI param
        # load the optional colour file
        colours = dict(csv.reader(open("nodes.list", "r"), delimiter="\t"))
    except IOError:
        colours = {}

    # set leaf node attributes
    for idx, node in enumerate(nodes):
        colour = colours.get(node, COLOURS[idx])
        dot.node(node, shape="ellipse", color=colour, fontcolor=colour)

    # render the graph (basename.pdf)
    dot.render(cleanup=True)
