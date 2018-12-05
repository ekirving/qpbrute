#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

# import all the constants
from consts import *


def run_cmd(cmd, returnout=True, shell=False, verbose=True):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param returnout: Return stdout as a string
    :param shell: Use the native shell
    :param verbose: Print the command before execution
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # print the command
    if verbose:
        print(u' '.join(cmd))

    # run the command
    proc = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

    if returnout:
        return stdout


def trim_ext(fullpath, n=1):
    return ('.').join(fullpath.split('.')[:-n])


def pprint_qpgraph(dot_file, pdf_file):
    """
    Pretty print a qpGraph dot file using the graphviz library.
    """
    import graphviz
    import re

    # extract the body contents from the dot file
    with open(dot_file, 'rU') as fin:
        body = re.search("{(.+)}", fin.read(), re.DOTALL).group(1).strip().split("\n")

    # make a new direct graph
    dot = graphviz.Digraph(filename=trim_ext(pdf_file), body=body, format='pdf')

    # remove the messy graph label
    dot.attr('graph', label='')

    # set Node defaults
    dot.node_attr['shape'] = 'point'
    dot.node_attr['fontname'] = 'arial'
    dot.node_attr['fontsize'] = '11'

    # set Edge defaults
    dot.edge_attr['arrowhead'] = 'vee'
    dot.edge_attr['fontcolor'] = '#838b8b'  # grey
    dot.edge_attr['fontname'] = 'arial'
    dot.edge_attr['fontsize'] = '11'

    nodes = []

    # extract the leaf nodes from the body of the graph
    for line in dot:
        match = re.search("^ +([a-z_]+) +\[", line, re.IGNORECASE)
        if match:
            nodes.append(match.group(1))

    # sort the nod
    nodes.sort()

    # set leaf node attributes
    for idx, node in enumerate(nodes):
        dot.node(node, shape='ellipse', color=COLOURS[idx], fontcolor=COLOURS[idx])

    # render the graph (basename.pdf)
    dot.render(cleanup=True)
