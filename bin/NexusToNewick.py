#!/usr/bin/env python

import dendropy
import argparse

parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-t', "--treeFile", type=str, required=True, help='tree')
parser.add_argument('-o', "--outputtree", type=str, required=True, help='tree')

args = parser.parse_args()
tree_path = args.treeFile
outputtree = args.outputtree

tree = dendropy.Tree.get_from_path(tree_path, 'nexus')
tree.write(path = outputtree, schema="newick")