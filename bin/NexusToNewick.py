#!/usr/bin/env python

import dendropy
import argparse

parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-t', "--treeFile", type=str, help='tree')
parser.add_argument('-o', "--outputtree", type=str, help='tree')

args = parser.parse_args()
tree_path = args.treeFile
outputtree = args.outputtree

tree_path = '/home/nehleh/PhiloBacteria/Results/num_4/num_4_beasttree.nexus'
outputtree = '/home/nehleh/PhiloBacteria/Results/num_4/num_4_beasttree.newick'

tree = dendropy.Tree.get_from_path(tree_path, 'nexus')
tree.write(path = outputtree, schema="newick")