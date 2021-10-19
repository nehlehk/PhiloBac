#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
from utility import *




def give_descendents_CFML(tree,node_label,result):
    if "NODE" in str(node_label):
        internal_recom_node = tree.find_node_with_label(node_label)
        children = internal_recom_node.child_nodes()
        for n in range(len(children)):
          r_node= children[n].label
          if "NODE" in str(r_node):
            give_descendents_CFML(tree,r_node,result)
          else:
            result.add(r_node)
    return result
# **********************************************************************************************************************
def CFML_resultFig(tree,CFMLData):
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    taxa = CFMLData.shape[1]
    for i in range(taxa):
        ax = fig.add_subplot(taxa, 1, i + 1)
        if i >= tips_num:
            node_label = str('NODE '+ str(i+1))
            desc = set()
            d = give_descendents_CFML(tree, node_label, desc)
            ax.plot(CFMLData[:, i], label=str(i+1) + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(CFMLData[:, i], label=i, color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        # ax.plot(CFMLData[:, i],label=i, color=color[i % 5])
        # ax.legend(bbox_to_anchor=(0.04, 1.33) ,prop={'size':10} )
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    # plt.show()
    plt.savefig("CFML_Recombination.jpeg")
# **********************************************************************************************************************
def CFML_recombination(CFML_recomLog,cfml_tree,tips_num):
    CFMLData = np.zeros((alignment_len, nodes_number))
    df = pd.read_csv(CFML_recomLog, sep='\t', engine='python')
    # print(df)
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        if "NODE_" in str(node):
            node = node[5:]
            mynode = int(give_taxon_index(cfml_tree, node,tips_num))
        else:
            mynode = int(node)
        CFMLData[s:e, mynode] = 1
        # CFMLData[s:e,int(node)] = 1

    return CFMLData
# **********************************************************************************************************************


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str, required=True, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    parser.add_argument('-cl', "--cfmllogFile", type=str, help='cfmlFile')
    parser.add_argument('-ct', "--cfmltreefile", type=str, help='cfmltreefile')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    args = parser.parse_args()

    clonal_path = args.clonaltreeFile
    cfml_log = args.cfmllogFile
    cfml_tree = args.cfmltreefile
    genomefile = args.alignmentFile
    baciSimLog = args.recomlogFile

    clonal_tree = Tree.get_from_path(clonal_path, 'newick')
    cfml_tree = Tree.get_from_path(cfml_tree, 'newick')

    # CFML_tree has label
    # set_label(cfml_tree)


    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    nodes_number = len(clonal_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    set_index(clonal_tree)

    realData = real_recombination(baciSimLog,clonal_tree,nodes_number,alignment_len,tips_num)
    CFMLData = CFML_recombination(cfml_log,cfml_tree,tips_num)
    CFML_resultFig(cfml_tree, CFMLData)
    rmse_real_CFML = mean_squared_error(realData, CFMLData, squared=False)
    write_rmse(rmse_real_CFML,'RMSE_CFML.csv')
