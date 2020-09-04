#!/usr/bin/env python
# -*- coding: utf-8 -*-

### Import modules ###
from Bio import Phylo
import argparse
import sys
import os
######################

def get_node2Ntips(tree):
    # Caution: all internal nodes and tips must be named!!!!
    node2Ntips = {}
    for tip in tree.get_terminals():
        node2Ntips[tip.name] = 1
    
    idx = 0
    for internal in tree.get_nonterminals():
        if internal.name is None:
            internal.name = "internal_{}".format(int(idx))
            idx += 1
    for internal in reversed(list(tree.get_nonterminals())):
        node2Ntips[internal.name] = sum([node2Ntips[child.name] for child in internal.clades])
        
    return node2Ntips

# func
def decompose(
    treefile,
    number,
    output
    ):
    # All nodes must be named!!!

    tree = Phylo.read(treefile, 'newick')

    # count downstream tips

    node2Ntips = get_node2Ntips(tree)

    # print(node2Ntips)

    # decompose tree

    stack = [tree.clade]
    counter = 1
    while counter != 0:
        node = stack.pop()
        if ( node2Ntips[node.name] > number ):
            stack.extend(node.clades)
            counter += 1
        else:
            newtree    = Phylo.BaseTree.Tree(node)
            Phylo.write(newtree, output + "/Down_"+node.name+".nwk", 'newick')
            node.clades = []
            counter -= 1

    Phylo.write(tree, output + "/deUp.nwk", 'newick')
        
    
# interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=SOFTWARE_NAME+' v1.0.0', add_help = True)
    
    parser.add_argument(
        "-v", "--version",
        help    = "Print "+SOFTWARE_NAME+" version",
        action  = 'store_true',
        default = False
        )
    
    parser.add_argument(
        "-t", "--tree",
        help    = "Input tree file path",
        )
    
    parser.add_argument(
        "-n", "--number",
        help = "Number of tips",
        type = int,
        )
    
    parser.add_argument(
        "-o", "--output",
        help = "Output folder",
        type = str,
        )

    # parse arguments
    args = parser.parse_args()

    # execute main_function()
    decompose(
        treefile = args.tree,
        number   = args.number,
        output   = args.output
        )