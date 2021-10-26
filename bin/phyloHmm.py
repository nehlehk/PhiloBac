#!/usr/bin/env python

import numpy as np
from utility import *
import pandas as pd
import matplotlib.pyplot as plt
import hmmlearn.base
import hmmlearn._utils
from dendropy import Tree, DnaCharacterMatrix
import dendropy
import xml.etree.ElementTree as ET
import argparse
import csv
import operator
import itertools
import json
import ast



class phyloLL_HMM(hmmlearn.base._BaseHMM):
    """
    params : string, optional
        Controls which parameters are updated in the training
        process.  Can contain any combination of 's' for startprob,
        't' for transmat, 'e' for emissionprob.
        Defaults to all parameters.
    init_params : string, optional
        Controls which parameters are initialized prior to
        training.  Can contain any combination of 's' for
        startprob, 't' for transmat, 'e' for emissionprob.
        Defaults to all parameters.

    """
    def __init__(self, n_components, trees, model ,child_order,X_child_order,n_iter=10, tol=1e-2, params="st", init_params="st"):
        super().__init__(n_components)
        self.trees = trees
        self.model = model
        self.child_order = child_order
        self.X_child_order = X_child_order

    def _init(self, X, lengths):
        init = 1. / self.n_components
        # if 's' in self.init_params or not hasattr(self, "startprob_"):
        #     self.startprob_ = np.full(self.n_components, init)
        # if 't' in self.init_params or not hasattr(self, "transmat_"):
        #     self.transmat_ = np.full((self.n_components, self.n_components), init)
        # n_fit_scalars_per_param = self._get_n_fit_scalars_per_param()
        # n_fit_scalars = sum(n_fit_scalars_per_param[p] for p in self.params)
        # if X.size < n_fit_scalars:
        #     _log.warning("Fitting a model with {} free scalar parameters with "
        #                  "only {} data points will result in a degenerate "
        #                  "solution.".format(n_fit_scalars, X.size))


    #     ==========================================================================
    def _check(self):
        self.startprob_ = np.asarray(self.startprob_)
        if len(self.startprob_) != self.n_components:
            raise ValueError("startprob_ must have length n_components")
        if not np.allclose(self.startprob_.sum(), 1.0):
            raise ValueError("startprob_ must sum to 1.0 (got {:.4f})"
                             .format(self.startprob_.sum()))

        self.transmat_ = np.asarray(self.transmat_)
        if self.transmat_.shape != (self.n_components, self.n_components):
            raise ValueError(
                "transmat_ must have shape (n_components, n_components)")
        if not np.allclose(self.transmat_.sum(axis=1), 1.0):
            raise ValueError("rows of transmat_ must sum to 1.0 (got {})"
                             .format(self.transmat_.sum(axis=1)))

    #     ==========================================================================
    def _compute_log_likelihood(self, X):
        return compute_logprob_phylo(X, self.trees, self.model, self.child_order, self.X_child_order,self.n_components)

    #     ==========================================================================

    def _initialize_sufficient_statistics(self):
        stats = super()._initialize_sufficient_statistics()
        return stats

    #     ==========================================================================

    def _accumulate_sufficient_statistics(self, stats, X, framelogprob,
                                          posteriors, fwdlattice, bwdlattice):

        super()._accumulate_sufficient_statistics(stats, X, framelogprob, posteriors, fwdlattice, bwdlattice)


    # ==========================================================================
    def _do_mstep(self, stats):
        super()._do_mstep(stats)

# **********************************************************************************************************************
def compute_logprob_phylo(X,recom_trees,model,child_order,X_child_order,status):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        # print(children)
        for site_id, partial in enumerate(X):
            # order = child_order.index(X_child_order[tree_id * len(children)])
            if status == 2:
                order = X_child_order.index(child_order[0])
            if status == 8:
                if tree_id == 0:
                    order = child_order.index(X_child_order[0])
                else:
                    order = X_child_order.index(child_order[0])
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[order * 4:(order + 1) * 4])
            for i in range(1, len(children)):
                # order = child_order.index(X_child_order[(tree_id* len(children)) + i])
                if status == 2:
                    order = X_child_order.index(child_order[i])
                if status == 8:
                    if tree_id == 0:
                        order = child_order.index(X_child_order[0])
                    else:
                        order = X_child_order.index(child_order[0])
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[order * 4:(order + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    return result
# **********************************************************************************************************************
def make_beast_xml_partial(tipdata,tree,xml_path,outputname):
    my_tipdata = tipdata.transpose(1, 0, 2)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(my_tipdata.shape[0]):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon" , str(give_taxon(tree,i)))
        c.set("uncertain" , "true")
        for j in range(my_tipdata.shape[1]):
          x = x + str(repr(my_tipdata[i,j,:]))[7:-2] + ';'
        c.text = '\n' + x +'\n'
        data.insert(i,c)
        c.tail = "\n"

    my_xml.write(outputname ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_beast_xml_gap(tipdata,tree,xml_path,threshold,outputname):
    my_tipdata = tipdata.transpose(1, 0, 2)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(my_tipdata.shape[0]):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon" , str(give_taxon(tree,i)))
        for j in range(my_tipdata.shape[1]):
            if (my_tipdata[i,j,0] >= threshold) and (my_tipdata[i,j,1] >= threshold) and (my_tipdata[i,j,2] >= threshold) and (my_tipdata[i,j,3] >= threshold):
                x = x + '-'
            else:
                x = x + str(alignment[i][j])
        c.text = '\n' + x +'\n'
        data.insert(i,c)
        c.tail = "\n"

    my_xml.write(outputname,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_beast_xml_delCol(recom_prob, threshold,outputname):
    r_output= recom_output(recom_prob, threshold)
    del_col = []
    for i in range(r_output.shape[0]):
        for j in range(r_output.shape[1]):
            if r_output[i][j] == 1 :
              del_col.append(i)

    my_tipdata = tipdata.transpose(1, 0, 2)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(my_tipdata.shape[0]):
            x = ''
            c = ET.Element("sequence")
            c.set("taxon" , str(give_taxon(tree,i)))
            for j in range(my_tipdata.shape[1]):
              if not (j in del_col):
                # print(j)
                x = x + str(alignment[i][j])
            c.text = '\n' + x +'\n'
            data.insert(i,c)
            c.tail = "\n"

    my_xml.write(outputname ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def set_tips_partial(column,tips_num):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        partial[site][tip][i] = 1
    return partial
# **********************************************************************************************************************
def computelikelihood_mixture(tree,alignment_len,column,tip_partial,model,tips_num):
    partial = np.zeros(((alignment_len,(2 * tips_num) -1, 4)))
    partial[:,0:tips_num,:] = tip_partial
    persite_ll = np.zeros(alignment_len)
    uniqueCol = list(set(column))
    for u in range(len(uniqueCol)):
      indexes = []
      indexes = [id for id, x in enumerate(column) if x == uniqueCol[u]]
      site = indexes[0]
      for node in tree.postorder_node_iter():
          if not node.is_leaf():
              children = node.child_nodes()
              partial[site,node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[site,int(children[0].index)])
              for i in range(1, len(children)):
                  partial[site, node.index] *= np.dot(model.p_matrix(children[i].edge_length),partial[site,int(children[i].index)])
      partial[indexes,:,:] = partial[site,:,:]
      p = np.dot(partial[site,tree.seed_node.index] , model.get_pi())
      persite_ll[indexes] = np.log(p)

    return persite_ll, partial
# **********************************************************************************************************************
def make_hmm_input_mixture(tree,alignment_len,column,tip_partial,model,tips_num):
    sitell, partial = computelikelihood_mixture(tree,alignment_len,column,tip_partial,model,tips_num)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment_len, children_count * 4))
    for id, child in enumerate(children):
        # print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, int(child.index), :]
    return x
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
    # co_recom = nu/2
    co_recom = nu
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")
    return recombination_tree
# **********************************************************************************************************************
def recom_maker(r_tree,index,nu):
    filter_fn = lambda n: hasattr(n, 'index') and n.index == index
    recombined_node = r_tree.find_node(filter_fn=filter_fn)
    return tree_evolver_rerooted(r_tree, recombined_node, nu)
# **********************************************************************************************************************
def give_rho(node,recom_prob,site,status,node_order):
  if status == 2:
    rho = recom_prob[site][1]
  if status == 8:
    rho = recom_prob[site][node_order]
  return rho
# **********************************************************************************************************************
def update_mixture_partial(column,node,tipdata,posterior,status,node_order):
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    if status == 2:
        rho = give_rho(node,posterior,site,status,node_order)
    if status == 8:
        rho = give_rho(node,posterior,site,status,node_order)
    for i in range(4):
        if i == my_number:
            tipdata[site,node.index,i] = 1
        else:
            tipdata[site,node.index,i] = rho
  return tipdata
# **********************************************************************************************************************
def recom_output(recom_prob, threshold,status):
    output = np.zeros((alignment_len, nodes_number))
    if status == 2:
        for i in range(len(recom_prob)):
            if (recom_prob['recom_nodes'][i] < tips_num):
                for j in range(alignment_len):
                    if (recom_prob['posterior'][i][j][1] >= threshold):
                        output[j, recom_prob['recom_nodes'][i]] = 1
            else:
                # for j in range(alignment_len):
                #     if (recom_prob['posterior'][i][j][1] >= mixtureProb):
                #         output[j, recom_prob['target_node'][i]] = 1
                for j in range(i + 1, len(recom_prob)):
                    if (recom_prob['recom_nodes'][i] == recom_prob['target_node'][j]) and (
                            recom_prob['recom_nodes'][j] == recom_prob['target_node'][i]):
                        for k in range(alignment_len):
                            if ((recom_prob['posterior'][i][k][1] >= threshold) and (
                                    recom_prob['posterior'][j][k][1] >= threshold)):
                                output[k, recom_prob['target_node'][i]] = 1
                            # if (recom_prob['posterior'][i][k] < recom_prob['posterior'][j][k]):
                            #   recom_prob['posterior'][i][k] = recom_prob['posterior'][j][k]
                            # if (recom_prob['posterior'][i][k] >= mixtureProb):
                            #     output[k, recom_prob['target_node'][i]] = 1
    if status == 8:
        for i in range(len(recom_prob)):
            r = recom_prob['recom_nodes'][i]
            if isinstance(r, int):
                if (r < tips_num):
                    for j in range(alignment_len):
                        if (recom_prob['posterior'][i][j] >= threshold):
                            output[j, r] = 1
                else:
                    for j in range(i + 1, len(recom_prob)):
                        r2 = recom_prob['recom_nodes'][j]
                        if isinstance(r2, int):
                            if (r == recom_prob['target_node'][j]) and (r2 == recom_prob['target_node'][i]):
                                for k in range(alignment_len):
                                    if ((recom_prob['posterior'][i][k] >= threshold) and (
                                            recom_prob['posterior'][j][k] >= threshold)):
                                        output[k, recom_prob['target_node'][i]] = 1
            else:
                for w in range(len(r)):
                    m_node = int(r[w])
                    if (m_node < tips_num):
                        for k in range(alignment_len):
                            if (recom_prob['posterior'][i][k] >= threshold):
                                output[k, m_node] = 1
                    # else:
                    #   for j in range(i+1,len(recom_prob)):
                    #     r2 = recom_prob['recom_nodes'][j]
                    #     if isinstance(r2, int) and (r2 > tips_num) and (r2 != m_node) :
                    #       if (m_node == recom_prob['target_node'][j]) and (r2 == recom_prob['target_node'][i]) :
                    #         print("m_node:",m_node)
                    #         print("r2:",r2)
                    #         for k in range(alignment_len):
                    #           # if (recom_prob['posterior'][i][k] < recom_prob['posterior'][j][k]):
                    #           #   recom_prob['posterior'][i][k] = recom_prob['posterior'][j][k]
                    #           if (recom_prob['posterior'][i][k] >= mixtureProb):
                    #             output[k, recom_prob['target_node'][i]] = 1
    return output
# **********************************************************************************************************************
def internal_plot(c_tree,posterior,hiddenStates,score,r_node,t_node,status):
    if status == 2 :
        for i in range(len(posterior)):
            poster = posterior[i]
            # hidden = hiddenStates[i]
            sc = score[i]
            fig = plt.figure(figsize=(10, 3))
            # ax = fig.add_subplot(2, 1, 1)
            # ax.plot(hidden)
            # ax.set_title("Node"+str(t_node[i]))
            # ax.set_ylabel("Clonal - Recombination State")
            ax2 = fig.add_subplot(1, 1, 1)
            ax2.set_title("Node" + str(t_node[i]))
            ax2.plot(poster[:, 0], label="ClonalFrame")
            if r_node[i] < tips_num:
                label = give_taxon(c_tree, r_node[i])
            else:
                label = r_node[i]
            ax2.plot(poster[:, 1], label="Node"+str(label))
            ax2.set_ylabel("posterior probability for each state")
            ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))
            plt.savefig("posterior_two_node"+str(t_node[i])+"_node"+str(label)+".jpeg")
            # plt.show()
    if status == 8 :
        for i in range(len(posterior)):
            poster = posterior[i]
            hidden = hiddenStates[i]
            sc = score[i]
            fig = plt.figure(figsize=(10, 3))
            # ax = fig.add_subplot(2, 1, 1)
            # ax.plot(hidden)
            # ax.set_title("Node" + str(tips_num+i))
            # ax.set_ylabel("Clonal - Recombination State")
            ax2 = fig.add_subplot(1, 1, 1)
            ax2.set_title("Node" + str(tips_num + i))
            ax2.plot(poster[:, 0], label="ClonalFrame")
            if (isinstance(r_node[i * 7], int)) and (r_node[i * 7]) < tips_num:
                label1 = give_taxon(c_tree, r_node[i * 7])
            else:
                label1 = r_node[i * 7]
            ax2.plot(poster[:, 1], label=str(label1))
            if (isinstance(r_node[i * 7 + 1], int)) and r_node[i * 7 + 1] < tips_num:
                label2 = give_taxon(c_tree, r_node[i * 7 + 1])
            else:
                label2 = r_node[i * 7 + 1]
            ax2.plot(poster[:, 2], label=str(label2))
            if (isinstance(r_node[i * 7 + 2], int)) and r_node[i * 7 + 2] < tips_num:
                label3 = give_taxon(c_tree, r_node[i * 7 + 2])
            else:
                label3 = r_node[i * 7 + 2]
            ax2.plot(poster[:, 3], label=str(label3))
            if (isinstance(r_node[i * 7 + 3], int)) and r_node[i * 7 + 3] < tips_num:
                label4 = give_taxon(c_tree, r_node[i * 7 + 3])
            else:
                label4 = r_node[i * 7 + 3]
            ax2.plot(poster[:, 4], label=str(label4))
            if (isinstance(r_node[i * 7 + 4], int)) and r_node[i * 7 + 4] < tips_num:
                label5 = give_taxon(c_tree, r_node[i * 7 + 4])
            else:
                label5 = r_node[i * 7 + 4]
            ax2.plot(poster[:, 5], label=str(r_node[i * 7 + 4]))
            if (isinstance(r_node[i * 7 + 5], int)) and r_node[i * 7 + 5] < tips_num:
                label6 = give_taxon(c_tree, r_node[i * 7 + 5])
            else:
                label6 = r_node[i * 7 + 5]
            ax2.plot(poster[:, 6], label=str(r_node[i * 7 + 5]))
            if (isinstance(r_node[i * 7 + 6], int)) and r_node[i * 7 + 6] < tips_num:
                label7 = give_taxon(c_tree, r_node[i * 7 + 6])
            else:
                label7 = r_node[i * 7 + 6]
            ax2.plot(poster[:, 7], label=str(r_node[i * 7 + 6]))
            ax2.set_ylabel("posterior probability for each state")
            ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))

            plt.savefig("posterior_eight_node" + str(tips_num + i) + ".jpeg")
# **********************************************************************************************************************
def recom_resultFig_dm(recom_prob,mixtureProb,outputname):
    output = recom_output(recom_prob,mixtureProb,status)
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
    clonaltree = Tree.get_from_path(tree_path, 'newick')
    set_index(clonaltree)
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
# *********************************************************************************************************************
def write_best_nu(best_nu,outputname):
    with open(outputname, mode='w') as bestnu_file:
        nu_writer = csv.writer(bestnu_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        nu_writer.writerow(best_nu)
# **********************************************************************************************************************
def phylohmm(tree,alignment_len,column,nu,p_start,p_trans,tips_num,status):
    mytree = []
    posterior = []
    hiddenStates = []
    score = []
    tipdata = set_tips_partial(column,tips_num)
    r_node = []
    t_node = []
    single_posterior = []

    print(tree.as_ascii_plot(show_internal_node_labels=True))

    best_nu = []

    for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
      # if target_node.index < 15:
        print(target_node.index)
        recombination_trees = []
        recombination_nodes = []
        child_order = []
        mytree.append(Tree.get_from_path(tree_path, 'newick'))
        set_index(mytree[id_tree])

        # ----------- Step 1 : Make input for hmm ------------------------------------------------------
        # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

        mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
        recombination_trees.append(mytree[id_tree].as_string(schema="newick"))

        # --------------  Step 1.2: Calculate X based on this re-rooted tree

        X = make_hmm_input_mixture(mytree[id_tree],alignment_len,column,tipdata,GTR_sample,tips_num)
        # to keep the order of clonal tree children
        X_child_order = []
        for id, child in enumerate(target_node.child_node_iter()):
            X_child_order.append(child.index)


        if status == 2 :
            # -------------- find best nu ----------------------
            # my_nu = np.arange(0.001, 0.0715, 0.01)
            #
            # nu = give_best_nu(X, my_nu, tree_path, recombination_trees[0], target_node, X_child_order)
            # best_nu.append([target_node.index, nu])

            # # nu = best[id_tree]
            # print(nu)
            # ----------- Step 2: make recombination trees -----------------------------------------------

            temptree = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree)

            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
            kids = temptree.seed_node.child_nodes()

            for k1, kid1 in enumerate(kids):
                child_order.append(kid1.index)  # keep the order of children after reroot
                # if (kid1.index < tips_num) or (target_node.index < kid1.index):  # to avoid having repeated internal nodes
                temptree = Tree.get_from_path(tree_path, 'newick')
                set_index(temptree)
                filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
                target_node_temp = temptree.find_node(filter_fn=filter_fn)
                temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)

                # t_max = temptree.max_distance_from_root()
                # print("t_max:",t_max,"kid1.edge_length:",kid1.edge_length)
                # my_nu = ((t_max  - kid1.edge_length) / t_max) * nu
                # print("nu:   ",my_nu)
                # recombination_trees.append(recom_maker(temptree, kid1.index, my_nu))

                recombination_trees.append(recom_maker(temptree, kid1.index, nu))
                # recombination_trees.append(recom_maker(temptree, kid1.index, nu[k1]))
                recombination_nodes.append(kid1)

            # print(recombination_trees)
            # ----------- Step 3: Call phyloHMM -----------------------------------------------------
            for h in range(1, len(recombination_trees)):
                model = phyloLL_HMM(n_components=status, trees=[recombination_trees[0], recombination_trees[h]],model=GTR_sample, child_order=child_order, X_child_order=X_child_order)
                model.startprob_ = p_start

                p_trans_nu0 = np.array([[1, 0], [1, 0]])

                # if nu[h-1] <= my_nu[0]:
                #     model.transmat_ = p_trans_nu0
                # else:
                #     model.transmat_ = p_trans

                model.transmat_ = p_trans

                p = model.predict_proba(X)
                hidden = model.predict(X)

                posterior.append(p)
                hiddenStates.append(hidden)
                score.append(model.score(X))

                r_node.append(recombination_nodes[h - 1].index)
                t_node.append(target_node.index)

            tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
            set_index(tree_updatePartial)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == X_child_order[h-1]
            update_child = tree_updatePartial.find_node(filter_fn=filter_fn)
            if update_child.is_leaf():
                update_mixture_partial(column, update_child, tipdata, p , status , 1)

        if status == 8 :
            # -------------- find best nu ----------------------
            # my_nu = np.arange(0.001, 0.0715, 0.01)
            #
            # nu = give_best_nu(X, my_nu, tree_path, recombination_trees[0], target_node, X_child_order)
            # best_nu.append([target_node.index,nu])
            # nu = best[id_tree]
            # print(nu)

            # ----------- Step 2: make recombination trees -----------------------------------------------
            # recom_child_order = []
            temptree = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree)

            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
            kids = temptree.seed_node.child_nodes()

            for k1, kid1 in enumerate(kids):
                child_order.append(kid1.index)

            # print(kids)
            def myFunc(e):
                return e.index

            kids.sort(key=myFunc)
            for k1, kid1 in enumerate(kids):
                temptree = Tree.get_from_path(tree_path, 'newick')
                set_index(temptree)
                filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
                target_node_temp = temptree.find_node(filter_fn=filter_fn)
                temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
                recombination_trees.append(recom_maker(temptree, kid1.index, nu))
                # recom_child_order.append(kid1.index)
                for k2, kid2 in enumerate(kids):
                    if kid1.index < kid2.index:
                        recombination_trees.append(recom_maker(temptree, kid2.index, nu))
                        # recom_child_order.append(kid2.index)
                    if k1 == k2 == 2:
                        recombination_trees.append(recom_maker(temptree, kids[0].index, nu))
                        # recom_child_order.append(kids[0].index)
                # recombination_trees.append(recom_maker(temptree, kid1.index, kid1.edge_length * nu))
                # recom_child_order.append(kid1.index)
                # for k2, kid2 in enumerate(kids):
                #     if kid1.index < kid2.index:
                #         recombination_trees.append(recom_maker(temptree, kid2.index, kid2.edge_length * nu))
                #         recom_child_order.append(kid2.index)
                #     if k1 == k2 == 2:
                #         recombination_trees.append(recom_maker(temptree, kids[0].index, kids[0].edge_length * nu))
                #         recom_child_order.append(kids[0].index)

            # ----------- Step 3: Call phyloHMM -----------------------------------------------------
            # print("recombination_trees after best nu:")
            # print(recombination_trees)
            model = phyloLL_HMM(n_components=status, trees=recombination_trees, model=GTR_sample, child_order=child_order,recom_child_order=X_child_order)
            model.startprob_ = p_start

            p_trans_nu0 = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0], ])

            # if nu <= my_nu[0]:
            #     model.transmat_ = p_trans_nu0
            # else:
            #     model.transmat_ = p_trans

            model.transmat_ = p_trans

            p = model.predict_proba(X)
            hidden = model.predict(X)

            posterior.append(p)
            hiddenStates.append(hidden)

            ss = model.score(X)
            score.append(ss)
            # print("final score is: ",ss)

            r_node.append(kids[0].index)
            t_node.append(target_node.index)
            single_posterior.append(p[:, 1])

            r_node.append([kids[0].index, kids[1].index])
            t_node.append(target_node.index)
            single_posterior.append(p[:, 2])

            r_node.append([kids[0].index, kids[1].index, kids[2].index])
            t_node.append(target_node.index)
            single_posterior.append(p[:, 3])

            r_node.append(kids[1].index)
            t_node.append(target_node.index)
            single_posterior.append(p[:, 4])

            r_node.append([kids[1].index, kids[2].index])
            t_node.append(target_node.index)
            single_posterior.append(p[:, 5])

            r_node.append(kids[2].index)
            t_node.append(target_node.index)
            single_posterior.append(p[:, 6])

            r_node.append([kids[0].index, kids[2].index])
            t_node.append(target_node.index)
            single_posterior.append(p[:, 7])

            tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
            set_index(tree_updatePartial)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_partial = tree_updatePartial.find_node(filter_fn=filter_fn)
            for id, child in enumerate(target_node_partial.child_node_iter()):
                # print(child.index)
                if child.is_leaf():
                    order = X_child_order.index(child.index)
                    # print("my beloved child:", child.index , child.taxon , "order:" , order+1)
                    update_mixture_partial(column,child, tipdata, p,status, order + 1)


    recom_prob = pd.DataFrame({'recom_nodes':r_node, 'target_node':t_node,  'posterior':posterior})
    return tipdata,posterior,hiddenStates,score,recom_prob,r_node,t_node,best_nu
# **********************************************************************************************************************
if __name__ == "__main__":

    tree_path = '/home/nehleh/Desktop/sisters/1/num_1_RAxML_bestTree.tree'
    genomefile = '/home/nehleh/Desktop/sisters/1/num_1_wholegenome_1.fasta'
    baciSimLog = '/home/nehleh/Desktop/sisters/1/BaciSim_Log.txt'
    clonal_path = '/home/nehleh/Desktop/sisters/1/clonaltree.tree'



    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    # parser.add_argument('-t', "--raxmltree", type=str, required= True, help='tree')
    # parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    # parser.add_argument('-cl', "--clonaltreeFile", type=str, help='clonalclonaltreeFile tree from BaciSim')
    # parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-nu', "--nuHmm", type=float,default=0.03,help='nuHmm')
    parser.add_argument('-p', "--threshold", type=float, default=0.9, help='threshold')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-s', "--startProb", type=list, default= [0.99, 0.01],help='frequencies')
    parser.add_argument('-m', "--transmat", type=list, default= [[0.999, 0.001],  [0.001, 0.999]], help='rates')
    parser.add_argument('-st', "--status", type=str,  default='2', help='2 for the two states hmm and 8 for eight states of hmm , [2,8] for both ')
    parser.add_argument('-x', "--xmlFile", type=str, help='xmlFile')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    args = parser.parse_args()

    # tree_path = args.raxmltree
    # genomefile = args.alignmentFile
    pi = args.frequencies
    rates = args.rates
    nu = args.nuHmm
    p_start = args.startProb
    p_trans = args.transmat
    threshold = args.threshold
    xml_path = args.xmlFile
    initialstat = args.status
    simulation = args.simulation


    # ============================================ setting parameters ==================================================
    tree = Tree.get_from_path(tree_path, 'newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    nodes_number = len(tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    set_index(tree)
    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)



    if initialstat.find('2') != -1:
        status = 2

        p_start = np.array([0.99, 0.01])
        p_trans = np.array([[0.999, 0.001],
                            [0.001, 0.999]])
        print('status = 2')
        tipdata, posterior, hiddenStates, score, recom_prob, r_node, t_node, best_nu = phylohmm(tree,alignment_len,column,nu, p_start, p_trans,tips_num,status)
        c_tree = Tree.get_from_path(tree_path, 'newick')
        set_index(c_tree)
        phyloHMMData2 = recom_resultFig_dm(recom_prob, threshold,'PhyloHMM_Recombination_two.jpeg')
        phyloHMM_log = phyloHMM_Log(c_tree, phyloHMMData2,'Recom_phyloHMM_Log_two.txt')
        write_best_nu(best_nu,'best_nu_two.txt')
        # # # ======================================= providing xml files for beast ============================================
        make_beast_xml_partial(tipdata, c_tree, xml_path,'RecomPartial_two.xml')
        make_beast_xml_gap(tipdata, tree, xml_path, 0.5,'Recom_Gap_two.xml')
        make_beast_xml_delCol(recom_prob, 0.5,'Recom_DelCol_two.xml')

    #----------------------------------------------------------------------------------------------------------------------------------------------------
    if initialstat.find('8') != -1:
        status = 8
        p_start = np.array([0.93, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
        p_trans = np.array([[0.9993, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001],
                            [0.0001, 0.9993, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001],
                            [0.0001, 0.0001, 0.9993, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001],
                            [0.0001, 0.0001, 0.0001, 0.9993, 0.0001, 0.0001, 0.0001, 0.0001],
                            [0.0001, 0.0001, 0.0001, 0.0001, 0.9993, 0.0001, 0.0001, 0.0001],
                            [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.9993, 0.0001, 0.0001],
                            [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.9993, 0.0001],
                            [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.9993], ])

        print('status = 8')

        tipdata, posterior, hiddenStates, score, recom_prob, r_node, t_node, best_nu = phylohmm(tree,alignment_len,nu, p_start, p_trans,tips_num,status)
        c_tree = Tree.get_from_path(tree_path, 'newick')
        set_index(c_tree)
        phyloHMMData8 = recom_resultFig_dm(recom_prob,threshold,'PhyloHMM_Recombination_eight.jpeg')
        phyloHMM_log = phyloHMM_Log(c_tree, phyloHMMData8,'Recom_phyloHMM_Log_eight.txt')
        write_best_nu(best_nu,'best_nu_eight.txt')
        # # ======================================= providing xml files for beast ============================================
        make_beast_xml_partial(tipdata, c_tree, xml_path,'RecomPartial_eight.xml')
        make_beast_xml_gap(tipdata, tree, xml_path, 0.5,'Recom_Gap_eight.xml')
        make_beast_xml_delCol(recom_prob, 0.5,'Recom_DelCol_eight.xml')




    if simulation == 1 :
        # clonal_path = args.clonaltreeFile
        # baciSimLog = args.recomlogFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick')
        nodes_number = len(clonal_tree.nodes())
        set_index(clonal_tree)
        realData = real_recombination(baciSimLog, clonal_tree, nodes_number, alignment_len, tips_num)
        if status == 2:
            rmse_real_philo= mean_squared_error(realData,phyloHMMData2,squared=False)
            write_rmse(rmse_real_philo, 'RMSE_PhiloBac_two.csv')
        if status == 8:
            rmse_real_philo= mean_squared_error(realData,phyloHMMData8,squared=False)
            write_rmse(rmse_real_philo, 'RMSE_PhiloBac_eight.csv')