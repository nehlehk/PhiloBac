#!/usr/bin/env python

import argparse


parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-t', "--treeFile", type=str, required=True, help='tree')
parser.add_argument('-o', "--outputtree", type=str, required=True, help='tree')

args = parser.parse_args()
tree_path = args.treeFile
outputtree = args.outputtree

# tree_path = "/home/nehleh/physher/release/physher_PB2_json_zero.txt"
# outputtree = "/home/nehleh/physher/release/physher_PB2_zero.newick"


with open(tree_path) as f:
    for line in f:
        if line[0] == '(':
            myfile = open(outputtree, 'w')
            myfile.write(line)
            myfile.close()