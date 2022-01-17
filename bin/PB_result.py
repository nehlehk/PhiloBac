#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
from utility import *
from BCBio import GFF


def set_index(tree):
    for node in tree.postorder_node_iter():
      node.index = -1
      node.annotations.add_bound_attribute("index")

    s = len(tree.leaf_nodes())
    for node in tree.postorder_node_iter():
      if not node.is_leaf():
          node.index = s
          node.label = str(node.index)
          s += 1
      else:
          node.index = int(node.taxon.label)
          node.label = str(node.index)
# **********************************************************************************************************************
def Gubbins_recombination(gubbins_log,gubbins_tree,nodes_number,alignment_len):
    starts = []
    ends = []
    desc = []
    mrca = []
    gubb_handle = open(gubbins_log)
    for rec in GFF.parse(gubb_handle):
        for feature in rec.features:
            starts.append(feature.location.start)
            ends.append(feature.location.end)
            d = str(feature.qualifiers['taxa'][0])
            kids = [int(s) for s in d.split() if s.isdigit()]
            desc.append(kids)
            mrca.append(my_mrca(gubbins_tree, kids))
    gubb_handle.close()

    all_data = {'nodes': desc, 'start': starts, 'end': ends , 'mrca' :mrca }
    df = pd.DataFrame(all_data)
    # print(df)

    GubbData = np.zeros((alignment_len, nodes_number))
    for i in range(len(df)):
        s = int(df['start'][i])
        e = int(df['end'][i])
        node = int(df['mrca'][i])
        GubbData[s:e, node] = 1

    return GubbData,df
# **********************************************************************************************************************
def Gubbins_resultFig(gubbins_tree,GubbData,tips_num,nodes_number,df):
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    for i in range(nodes_number):
        ax = fig.add_subplot(nodes_number, 1, i + 1)
        if i >= tips_num:
            d = (df.loc[df['mrca'] == i]['nodes'])
            d = d.drop_duplicates()
            if len(d.values) > 0 :
                txt = str(d.values[0])
            else:
                desc = set()
                txt = str(give_descendents(gubbins_tree, i, desc, tips_num))
                # print(txt)
            ax.plot(GubbData[:, i], label=str(i) + ' is mrca:' + txt, color=color[i % 5])
        else:
            ax.plot(GubbData[:, i], label=i, color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    # plt.show()
    plt.savefig("Gubbins_Recombination.jpeg")
# **********************************************************************************************************************

if __name__ == "__main__":

    # clonal_path = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/num_1_Clonaltree.tree'
    # genomefile = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/num_1_recom_1_Wholegenome_1_1.fasta'
    # baciSimLog = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/num_1_recom_1_BaciSim_Log.txt'
    # gubbins_log = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/gubbins.recombination_predictions.gff'
    # gubb_tree = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/gubbins.node_labelled.final_tree.tre'
    # gubb_csv = '/home/nehleh/PhyloCode/RecomPhyloHMM/Results/num_1/gubbins.per_branch_statistics.csv'

    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str,  help='clona tree from BaciSim')
    parser.add_argument('-a', "--alignmentFile", type=str,  help='fasta file')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-gl', "--PBlogFile", type=str, help='PB Log File')
    parser.add_argument('-gt', "--PBtreefile", type=str, help='PB tree File')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    parser.add_argument('-st', "--status", type=str, default='2',help='2 for the two states hmm and 8 for eight states of hmm , 2,8 for both ')
    parser.add_argument('-p', "--threshold", type=float, default=0.9, help='threshold')

    args = parser.parse_args()


    genomefile = args.alignmentFile
    pb_log = args.PBlogFile
    pb_tree = args.PBtreefile
    simulation = args.simulation
    initialstat = args.status
    threshold = args.threshold

    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    PB_tree = Tree.get_from_path(pb_tree, 'newick')
    set_index(PB_tree)
    nodes_num_g = len(PB_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size


    if initialstat.find('2') != -1:
        status = 2
        phyloHMMData2 = recom_resultFig_dm(recom_prob, tips_num, threshold, status, 'PB_Recom_two.jpeg')



    if initialstat.find('8') != -1:
        status = 8
        phyloHMMData8 = recom_resultFig_dm(recom_prob,tips_num,threshold,status,'PB_Recom_eight.jpeg')


    if simulation == 1 :
        clonal_path = args.clonaltreeFile
        baciSimLog = args.recomlogFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick')
        nodes_number_c = len(clonal_tree.nodes())
        set_index(clonal_tree,alignment)
        realData = real_recombination(baciSimLog, clonal_tree, nodes_number_c, alignment_len, tips_num)
        print(realData.shape)
        # print(nodes_number_c)
        print(PB_tree.as_ascii_plot(show_internal_node_labels=True))
        if initialstat.find('2') != -1:
            rmse_real_philo2 = mean_squared_error(realData,phyloHMMData2,squared=False)
            write_rmse(rmse_real_philo2, 'RMSE_PB_two.csv')
        if initialstat.find('8') != -1:
            print(phyloHMMData8.shape)
            rmse_real_philo8 = mean_squared_error(realData,phyloHMMData8,squared=False)
            write_rmse(rmse_real_philo8, 'RMSE_PB_eight.csv')









