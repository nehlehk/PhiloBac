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
import operator
import itertools

def recom_output(recom_prob,tips_num,threshold,status,nodes_number):
    output = np.zeros((alignment_len,nodes_number))
    if status == 2:
        for i in range(len(recom_prob)):
            if (recom_prob['recom_nodes'][i] < tips_num):
                for j in range(alignment_len):
                    if (float(recom_prob['posterior1'][i][j]) >= threshold):
                        output[j, recom_prob['recom_nodes'][i]] = 1
            else:
                # for j in range(alignment_len):
                #     if (recom_prob['posterior'][i][j][1] >= mixtureProb):
                #         output[j, recom_prob['target_node'][i]] = 1
                for j in range(i + 1, len(recom_prob)):
                    if (recom_prob['recom_nodes'][i] == recom_prob['target_node'][j]) and (recom_prob['recom_nodes'][j] == recom_prob['target_node'][i]):
                        for k in range(alignment_len):
                            if ((float(recom_prob['posterior1'][i][k]) >= threshold) and (float(recom_prob['posterior1'][j][k]) >= threshold)):
                                output[k, recom_prob['target_node'][i]] = 1
                            # if (recom_prob['posterior'][i][k] < recom_prob['posterior'][j][k]):
                            #   recom_prob['posterior'][i][k] = recom_prob['posterior'][j][k]
                            # if (recom_prob['posterior'][i][k] >= mixtureProb):
                            #     output[k, recom_prob['target_node'][i]] = 1

    return output
# **********************************************************************************************************************
def recom_resultFig_dm(pb_tree,recom_prob,tips_num,mixtureProb,status,outputname,nodes_number):
    output = recom_output(recom_prob,tips_num,mixtureProb,status,nodes_number)
    # output = np.zeros((alignment_len, nodes_number))
    # for i in range(len(recom_prob)):
    #     r = recom_prob['recom_nodes'][i]
    #     if isinstance(r, int):
    #       if (r < tips_num):
    #           for j in range(alignment_len):
    #               if (recom_prob['posterior'][i][j] >= mixtureProb):
    #                 output[j, r] = 1
    #       else:
    #           for j in range(i+1,len(recom_prob)):
    #             r2 = recom_prob['recom_nodes'][j]
    #             if isinstance(r2, int):
    #               if (r == recom_prob['target_node'][j]) and (r2 == recom_prob['target_node'][i]) :
    #                 for k in range(alignment_len):
    #                   if ((recom_prob['posterior'][i][k] >= mixtureProb) and (recom_prob['posterior'][j][k] >= mixtureProb)):
    #                       output[k, recom_prob['target_node'][i]] = 1
    #     else:
    #         for w in range(len(r)):
    #             m_node = int(r[w])
    #             if (m_node < tips_num):
    #                 for k in range(alignment_len):
    #                     if (recom_prob['posterior'][i][k] >= mixtureProb):
    #                         output[k, m_node] = 1
    #             # else:
    #             #   for j in range(i+1,len(recom_prob)):
    #             #     r2 = recom_prob['recom_nodes'][j]
    #             #     if isinstance(r2, int) and (r2 > tips_num) and (r2 != m_node) :
    #             #       if (m_node == recom_prob['target_node'][j]) and (r2 == recom_prob['target_node'][i]) :
    #             #         print("m_node:",m_node)
    #             #         print("r2:",r2)
    #             #         for k in range(alignment_len):
    #             #           # if (recom_prob['posterior'][i][k] < recom_prob['posterior'][j][k]):
    #             #           #   recom_prob['posterior'][i][k] = recom_prob['posterior'][j][k]
    #             #           if (recom_prob['posterior'][i][k] >= mixtureProb):
    #             #             output[k, recom_prob['target_node'][i]] = 1


    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    clonaltree = Tree.get_from_path(pb_tree, 'newick')
    set_index(clonaltree,alignment)
    for i in range(nodes_number):
        ax = fig.add_subplot(nodes_number, 1, i + 1)
        if i >= tips_num:
            desc = set()
            d = give_descendents(clonaltree, i, desc,tips_num)
            ax.plot(output[:, i], label=str(i) + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(output[:, i], label=give_taxon(clonaltree, i), color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        ax.set_frame_on(False)
        ax.axis('off')

    ax.axis('on')
    ax.set_yticklabels([])
    plt.savefig(outputname)

    return output
# **********************************************************************************************************************
def phyloHMM_Log(c_tree,output,outputname):
    nodes = []
    starts = []
    ends = []
    recomlens = []

    for i in range(output.shape[1]):
        non_zeros = [[i for i, value in it] for key, it in itertools.groupby(enumerate(output[:, i]), key=operator.itemgetter(1)) if key != 0]
        for j in range(len(non_zeros)):
            if i < tips_num:
                n = give_taxon(c_tree, i)
            else:
                n = i
            nodes.append(n)
            starts.append(non_zeros[j][0])
            ends.append(non_zeros[j][len(non_zeros[j]) - 1])
            recomlens.append(non_zeros[j][len(non_zeros[j]) - 1] - non_zeros[j][0])

    all_data = {'nodes': nodes, 'start': starts, 'end': ends, 'len': recomlens}
    df = pd.DataFrame(all_data)
    df = df.sort_values(by=['nodes'], ascending=[True])
    df.to_csv(outputname, sep='\t', header=True , index = False)

    return df
# **********************************************************************************************************************


if __name__ == "__main__":

    # path = '/home/nehleh/PhiloBacteria/Results/num_3'
    # clonal_path = path+'/num_3_Clonaltree.tree'
    # genomefile = path+'/num_3_recom_1_Wholegenome_3_1.fasta'
    # baciSimLog = path+'/num_3_recom_1_BaciSim_Log.txt'
    # pb_tree = path+'/num_3_recom_1_physherTree_PB.newick'
    # recomProb = path+'/num_3_recom_1_Recom_prob_two.csv'


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str,  help='clona tree from BaciSim')
    parser.add_argument('-a', "--alignmentFile", type=str,  help='fasta file')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-pb', "--PBtreefile", type=str, help='PBtree')
    parser.add_argument('-rp', "--recomProb", type=str, help='recomProb')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    parser.add_argument('-st', "--status", type=str, default='2',help='2 for the two states hmm and 8 for eight states of hmm , 2,8 for both ')
    parser.add_argument('-p', "--threshold", type=float, default=0.5, help='threshold')

    args = parser.parse_args()


    genomefile = args.alignmentFile
    pb_tree = args.PBtreefile
    recomProb = args.recomProb
    simulation = args.simulation
    initialstat = args.status
    threshold = args.threshold


    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    PB_tree = Tree.get_from_path(pb_tree, 'newick')
    set_index(PB_tree,alignment)
    nodes_num_pb = len(PB_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size


    # print(PB_tree.as_ascii_plot(show_internal_node_labels=True))


    # f = open('/home/nehleh/Desktop/sisters/mutiple_sisters/Recom_prob_two.csv', "r")

    f = open(recomProb,"r")
    recom_prob = pd.read_csv(f, sep=',')
    recom_prob['posterior1'] = recom_prob['posterior1'].str.replace("[","")
    recom_prob['posterior1'] = recom_prob['posterior1'].str.replace("]","")
    recom_prob['posterior1'] = recom_prob['posterior1'].str.split()
    # print(type(recom_prob['posterior1'][0][0]))



    if initialstat.find('2') != -1:
        status = int(2)
        phyloHMMData2 = recom_resultFig_dm(pb_tree, recom_prob, tips_num, threshold, status, 'PB_Recom_two.jpeg', nodes_num_pb)
        # print(phyloHMMData2.shape)
        phyloHMM_log = phyloHMM_Log(PB_tree, phyloHMMData2, 'PB_Log_two.txt')
        # print("phyloHMMData2")
        # print(phyloHMMData2[780:1000,1])


    if initialstat.find('8') != -1:
        status = 8
        phyloHMMData8 = recom_resultFig_dm(pb_tree, recom_prob, tips_num, threshold, status, 'PB_Recom_eight.jpeg', nodes_num_pb)
        phyloHMM_log = phyloHMM_Log(PB_tree, phyloHMMData8, 'PB_Log_eight.txt')


    if simulation == 1 :
        clonal_path = args.clonaltreeFile
        baciSimLog = args.recomlogFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick')
        nodes_number_c = len(clonal_tree.nodes())
        set_index(clonal_tree,alignment)
        realData = real_recombination(baciSimLog, clonal_tree, nodes_number_c, alignment_len, tips_num)
        # print("realData")
        # print(realData.shape)
        # print(realData[780:1000,1])
        # print(clonal_tree.as_ascii_plot(show_internal_node_labels=True))
        if initialstat.find('2') != -1:
            rmse_real_philo2 = mean_squared_error(realData,phyloHMMData2,squared=True)
            write_rmse(rmse_real_philo2, 'RMSE_PB_two.csv')
        if initialstat.find('8') != -1:
            # print(phyloHMMData8.shape)
            rmse_real_philo8 = mean_squared_error(realData,phyloHMMData8,squared=False)
            write_rmse(rmse_real_philo8, 'RMSE_PB_eight.csv')









