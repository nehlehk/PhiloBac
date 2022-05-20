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

def recom_output(pb_tree,recom_prob,tips_num,threshold,status,nodes_number):
    output = np.zeros((alignment_len,nodes_number))
    rmsedata = np.zeros((alignment_len,tips_num))
    if status == 2:
        for i in range(len(recom_prob)):
            if (int(recom_prob['recom_nodes'][i]) < tips_num):
                for j in range(alignment_len):
                    if (float(recom_prob['posterior_rec'][i][j]) >= threshold):
                        output[j, recom_prob['recom_nodes'][i]] = 1
                        rmsedata[j, recom_prob['recom_nodes'][i]] = 1
            else:
                for k in range(alignment_len):
                    if (float(recom_prob['posterior_rec'][i][k]) >= threshold):
                        output[k, recom_prob['recom_nodes'][i]] = 1
                        desc = set()
                        d = give_descendents(pb_tree, int(recom_prob['recom_nodes'][i]), desc,tips_num)
                        # print(d)
                        for elm in d:
                            rmsedata[k, int(elm)] = 1
    return output,rmsedata
# **********************************************************************************************************************
def recom_resultFig_dm(pb_tree,recom_prob,tips_num,mixtureProb,status,outputname,nodes_number):
    output,rmsedata = recom_output(PB_tree,recom_prob,tips_num,mixtureProb,status,nodes_number)
    # print(output)
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    clonaltree = Tree.get_from_path(pb_tree, 'newick',rooting='force-rooted')
    # clonaltree.reroot_at_midpoint(update_bipartitions=True)
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
    plt.text(0.045,2.7, r'test test test test test test test test test', fontsize=14,
             horizontalalignment='left', verticalalignment='top',
             bbox=dict(boxstyle="sawtooth", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8),)
             )
    plt.savefig(outputname)
    # plt.show()
    return output,rmsedata
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

    # path = '/home/nehleh/Desktop/sisters/mutiple_sisters/'
    # clonal_path = path+'/clonaltree.tree'
    # genomefile = path+'/num_1_wholegenome_1.fasta'
    # baciSimLog = path+'/BaciSim_Log.txt'
    # pb_tree = path+'/PB_two.newick'
    # recomProb = '/home/nehleh/PhiloBacteria/bin/Recom_prob_two.h5'



    path = '/home/nehleh/PhiloBacteria/Results/num_1'
    clonal_path = path+'/num_1_Clonaltree.tree'
    genomefile = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_Wholegenome.fasta'
    baciSimLog = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_BaciSim_Log.txt'
    # pb_tree = path+'/num_1_nu_0.05_Rlen_1000_Rrate_0.01_physherTree_PB.newick'
    pb_tree = '/home/nehleh/PhiloBacteria/bin/pb_tree.newick'
    recomProb = '/home/nehleh/PhiloBacteria/bin/Recom_prob_two.h5'
    # recomProb = path+'/num_1_nu_0.05_Rlen_1000_Rrate_0.01_Recom_prob_two.h5'
    baciSimStat = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_Recom_stat.csv'
    json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str,  help='clona tree from BaciSim')
    parser.add_argument('-a', "--alignmentFile", type=str,  help='fasta file')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-pb', "--PBtreefile", type=str, help='PBtree')
    parser.add_argument('-rp', "--recomProb", type=str, help='recomProb')
    parser.add_argument('-rs', "--recomstat", type=str, help='recomstat')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    parser.add_argument('-st', "--status", type=str, default='2',help='2 for the two states hmm and 8 for eight states of hmm , 2,8 for both ')
    parser.add_argument('-p', "--threshold", type=float, default=0.8, help='threshold')

    args = parser.parse_args()


    genomefile = args.alignmentFile
    pb_tree = args.PBtreefile
    recomProb = args.recomProb
    baciSimStat = args.recomstat
    simulation = args.simulation
    initialstat = args.status
    threshold = args.threshold


    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    PB_tree = Tree.get_from_path(pb_tree, 'newick' ,  rooting='force-rooted')
    set_index(PB_tree,alignment)
    print(PB_tree.as_ascii_plot(show_internal_node_labels=True))
    nodes_num_pb = len(PB_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    recom_prob = pd.read_hdf(recomProb, key='recom_prob')


    if initialstat.find('2') != -1:
        status = int(2)
        phyloHMMData2,rmse_PB2 = recom_resultFig_dm(pb_tree, recom_prob, tips_num, threshold, status, 'PB_Recom_two.jpeg', nodes_num_pb)
        phyloHMM_log = phyloHMM_Log(PB_tree, phyloHMMData2, 'PB_Log_two.txt')
        # print(phyloHMM_log)
        recom_count_pb2= len(phyloHMM_log)
        write_value(len(phyloHMM_log),'PB_rcount_two.csv')
        phyloHMM_log[['len']].to_csv('PB_delta_two.csv', index=False)
        mean_dalta = phyloHMM_log[['len']].mean()



    if initialstat.find('8') != -1:
        status = 8
        phyloHMMData8 = recom_resultFig_dm(pb_tree, recom_prob, tips_num, threshold, status, 'PB_Recom_eight.jpeg', nodes_num_pb)
        phyloHMM_log = phyloHMM_Log(PB_tree, phyloHMMData8, 'PB_Log_eight.txt')



    if simulation == 1 :
        clonal_path = args.clonaltreeFile
        baciSimLog = args.recomlogFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick',rooting='force-rooted' )
        clonal_tree.resolve_polytomies(update_bipartitions=True)
        nodes_number_c = len(clonal_tree.nodes())
        set_index(clonal_tree,alignment)
        # print(clonal_tree.as_ascii_plot(show_internal_node_labels=True))
        realData,rmse_real = real_recombination(baciSimLog, clonal_tree, nodes_number_c, alignment_len, tips_num)


        recom_stat = pd.read_csv(open(baciSimStat, "r"), sep=',')
        # print(recom_stat)
        num_recom_real = len(recom_stat)
        write_value(len(recom_stat), 'baci_rcount.csv')
        recom_stat[['len']].to_csv('baci_delta.csv' , index=False)



        if initialstat.find('2') != -1:
            rmse_real_philo2 = mean_squared_error(rmse_real,rmse_PB2,squared=False)
            write_value(rmse_real_philo2, 'RMSE_PB_two.csv')
            # print(rmse_real_philo2)
        if initialstat.find('8') != -1:
            # print(phyloHMMData8.shape)
            rmse_real_philo8 = mean_squared_error(rmse_real,phyloHMMData8,squared=False)
            write_value(rmse_real_philo8, 'RMSE_PB_eight.csv')