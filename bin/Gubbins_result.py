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
def rescale_gubbtree(gubbins_tree,gubb_csv,alignment_len):
    gubbins_tree.scale_edges(1/alignment_len)
    g_tree = gubbins_tree.as_string(schema="newick")
    g = g_tree.replace('\n',"")
    myfile = open('./Gubbinstree_rescale.tree', 'w')
    myfile.write(g)
    myfile.close()
    #
    # df = pd.read_csv(gubb_csv, sep='\t', engine='python')
    # print(df)
    # for node in gubbins_tree.levelorder_node_iter():
    #     print(node.label)
    #     if int(node.label) < tips_num:
    #         temp = (df.loc[df['Node'] == str(node.label)]['Bases in Clonal Frame']).values[0]
    #         node.edge_length = node.edge_length /temp
    #         print(temp)
    #     else:
    #         node.edge_length = node.edge_length / alignment_len
    #     myfile = open('./Gubbinstree_rescale.tree', 'w')
    #     g_tree = gubbins_tree.as_string(schema="newick")
    #     g = g_tree.replace('\n',"")
    #     myfile.write(g_tree)
    #     myfile.close()
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
    parser.add_argument('-gl', "--gubblogFile", type=str, help='Gubbins Log File')
    parser.add_argument('-gt', "--gubbtreefile", type=str, help='Gubbins tree File')
    parser.add_argument('-gs', "--gubbcsvfile", type=str, help='Gubbins per_branch_statistics csv File')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')

    args = parser.parse_args()


    genomefile = args.alignmentFile
    gubbins_log = args.gubblogFile
    gubb_tree = args.gubbtreefile
    gubb_csv = args.gubbcsvfile
    simulation = args.simulation

    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    gubbins_tree = Tree.get_from_path(gubb_tree, 'newick')
    set_index(gubbins_tree,alignment)
    nodes_num_g = len(gubbins_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size

    GubbData, df = Gubbins_recombination(gubbins_log, gubbins_tree, nodes_num_g, alignment_len)
    Gubbins_resultFig(gubbins_tree, GubbData, tips_num, nodes_num_g, df)
    rescale_gubbtree(gubbins_tree, gubb_csv, alignment_len)


    if simulation == 1:
        clonal_path = args.clonaltreeFile
        baciSimLog = args.recomlogFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick')
        nodes_num_c = len(clonal_tree.nodes())
        set_index(clonal_tree,alignment)
        realData = real_recombination(baciSimLog, clonal_tree, nodes_num_c, alignment_len, tips_num)
        rmse_real_CFML = mean_squared_error(realData, GubbData, squared=False)
        write_rmse(rmse_real_CFML, 'RMSE_Gubbins.csv')









