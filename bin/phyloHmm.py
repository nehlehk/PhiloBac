#!/usr/bin/env python
from builtins import print

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
import scipy.optimize as spo
from scipy.optimize import Bounds
import os
import sys
import time


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
    XX = np.transpose(X.reshape((X.shape[0], 3, 4)), (0, 2, 1))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        orders = [X_child_order.index(child_order[i]) for i in range(len(children))]
        branch_lengths = np.array( [children[i].edge_length for i in range(len(children))])
        matrices = model.p_t(branch_lengths)
        result[:, tree_id] = np.log(model.get_pi() @ (matrices @ XX[..., orders].T).prod(0))

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
def make_beast_xml_delCol(recom_prob,tips_num,threshold,outputname):
    r_output= recom_output(recom_prob,tips_num,threshold,status)
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
    partial_new =  np.rollaxis(partial, 2, 0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[..., node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[..., children[0].index])
            for i in range(1, len(children)):
                partial_new[..., node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[..., children[i].index])

    persite_ll = np.log(model.get_pi() @ partial_new[..., tree.seed_node.index])

    return persite_ll, partial
# **********************************************************************************************************************
def make_hmm_input_mixture(tree,alignment_len,column,tip_partial,model,tips_num):
    sitell, partial = computelikelihood_mixture(tree, alignment_len, column, tip_partial, model, tips_num)

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
def update_mixture_partial_PSE(column,node,tipdata,posterior,status,node_order,nu):
  # print("update_mixture_partial_PSE:",node.index)
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    if status == 2:
        rho = give_rho(node,posterior,site,status,node_order)
    if status == 8:
        rho = give_rho(node,posterior,site,status,node_order)
    for i in range(4):
        if i == my_number:
            tipdata[site,node.index,i] = 1 - (nu * rho)
        else:
            tipdata[site,node.index,i] = (nu * rho)/3
  return tipdata
# **********************************************************************************************************************
def update_mixture_partial_PSE_AllNodes(column,node,tipdata,posterior,status,node_order,nu):
  for site in range(alignment_len):
    if node.index < tips_num:
        dna = column[site]
        my_number = give_index(dna[node.index])
        if status == 2:
            rho = give_rho(node,posterior,site,status,node_order)
        if status == 8:
            rho = give_rho(node,posterior,site,status,node_order)
        for i in range(4):
            if i == my_number:
                tipdata[site,node.index,i] = 1 - (nu * rho)
            else:
                tipdata[site,node.index,i] = (nu * rho)/3
    else:
        tipdata[site, node.index] = (nu * rho) / 3
  return tipdata
# **********************************************************************************************************************
def recom_output(recom_prob,tips_num,threshold,status):
    output = np.zeros((alignment_len,nodes_number))
    if status == 2:
        for i in range(len(recom_prob)):
            if (recom_prob['recom_nodes'][i] < tips_num):
                for j in range(alignment_len):
                    if (recom_prob['posterior'][i][j][1] >= threshold):
                        output[j, recom_prob['recom_nodes'][i]] = 1
            else:
                for j in range(i + 1, len(recom_prob)):
                    if (recom_prob['recom_nodes'][i] == recom_prob['target_node'][j]) and (recom_prob['recom_nodes'][j] == recom_prob['target_node'][i]):
                        for k in range(alignment_len):
                            if ((recom_prob['posterior'][i][k][1] >= threshold) and (recom_prob['posterior'][j][k][1] >= threshold)):
                                output[k, recom_prob['target_node'][i]] = 1
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
                                    if ((recom_prob['posterior'][i][k] >= threshold) and (recom_prob['posterior'][j][k] >= threshold)):
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
def internal_plot(c_tree,posterior0,posterior1,hiddenStates,score,r_node,t_node,status):
    if status == 2 :
        for i in range(len(posterior0)):
            # poster = posterior[i]
            hidden = hiddenStates[i]
            sc = score[i]
            fig = plt.figure(figsize=(10, 3))
            ax = fig.add_subplot(2, 1, 1)
            ax.plot(hidden)
            ax.set_title("Node"+str(t_node[i]))
            ax.set_ylabel("Clonal - Recombination State")
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.set_title("Node" + str(t_node[i]))
            ax2.plot(posterior0[i], label="ClonalFrame")
            if r_node[i] < tips_num:
                label = give_taxon(c_tree, r_node[i])
            else:
                label = r_node[i]
            ax2.plot(posterior1[i], label="Node"+str(label))
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
            ax = fig.add_subplot(2, 1, 1)
            ax.plot(hidden)
            ax.set_title("Node" + str(tips_num+i))
            ax.set_ylabel("Clonal - Recombination State")
            ax2 = fig.add_subplot(2, 1, 2)
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
def score_plot_one(nu_score, my_nu,i,status,r_node):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(my_nu, nu_score)
    ax.axvline(my_nu[nu_score.index(max(nu_score))], color='r', ls='-.')
    ax.annotate(str(my_nu[nu_score.index(max(nu_score))]), xy=(my_nu[nu_score.index(max(nu_score))], max(nu_score)),
                xytext=(my_nu[nu_score.index(max(nu_score))] + 0.015, max(nu_score) - 1),
                arrowprops=dict(facecolor='black', shrink=0.05))
    ax.set_ylabel("score")
    ax.set_xlabel("nu")
    ax.xaxis.set_label_coords(0.98, -0.025)
    if status == 2:
        ax.set_title("Node" + str(i)+"_child"+str(r_node.index))
        plt.savefig("two_nu_score_node" + str(i) + "_child" + str(r_node.index) + ".jpeg")
        ax.set_title("Node" + str(i))
        # plt.savefig("two_nu_score_node" + str(i) + ".jpeg")
    if status == 8:
        ax.set_title("Node" + str(i))
        plt.savefig("eight_nu_score_node" + str(i) + ".jpeg")
    # plt.show()
# **********************************************************************************************************************
def nu_trees(nu,target_node,tree_path,status):
    recom_child_order = []
    recombination_trees = []
    child_order = []
    recombination_nodes = []
    temptree = Tree.get_from_path(tree_path, 'newick')
    # temptree.reroot_at_midpoint(update_bipartitions=False)
    set_index(temptree, alignment)
    filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
    target_node_temp = temptree.find_node(filter_fn=filter_fn)
    temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
    kids = temptree.seed_node.child_nodes()

    if status == 2:
        for k1, kid1 in enumerate(kids):
            child_order.append(kid1.index)  # keep the order of children after reroot
            temptree = Tree.get_from_path(tree_path, 'newick')
            # temptree.reroot_at_midpoint(update_bipartitions=False)
            set_index(temptree, alignment)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
            recombination_trees.append(recom_maker(temptree, kid1.index, nu))
            recombination_nodes.append(kid1)
    if status == 8:
        for k1, kid1 in enumerate(kids):
            child_order.append(kid1.index)

        # print(kids)
        def myFunc(e):
            return e.index

        kids.sort(key=myFunc)
        for k1, kid1 in enumerate(kids):
            temptree = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree, alignment)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
            recombination_trees.append(recom_maker(temptree, kid1.index, nu))
            recom_child_order.append(kid1.index)
            for k2, kid2 in enumerate(kids):
                if kid1.index < kid2.index:
                    recombination_trees.append(recom_maker(temptree, kid2.index, nu))
                    recom_child_order.append(kid2.index)
                if k1 == k2 == 2:
                    recombination_trees.append(recom_maker(temptree, kids[0].index, nu))
                    recom_child_order.append(kids[0].index)

    return recombination_trees,child_order,recombination_nodes
# **********************************************************************************************************************
def give_best_nu(X,tree_path,clonal,target_node,X_child_order,status,my_nu):
    if status == 2:
        score = np.zeros((3, len(my_nu)))
        for id, nu in enumerate(my_nu):
            # print(nu)
            r_trees, child_order, r_node = nu_trees(nu, target_node, tree_path, status)
            # print(r_node)
            # print(X_child_order)
            # print("child_order:",child_order)
            for h in range(0, len(r_trees)):
                model = phyloLL_HMM(n_components=2, trees=[clonal, r_trees[h]], model=GTR_sample,child_order=child_order, X_child_order=X_child_order)
                model.startprob_ = p_start
                model.transmat_ = p_trans
                s = model.score(X)
                score[h][id] = s

        best_nu = []
        for h in range(0, len(r_trees)):
            kid_score = list(score[h])
            # print("score:", kid_score)
            score_plot_one(kid_score, my_nu, target_node.index,status,r_node[h])
            best_nu.append(my_nu[kid_score.index(max(kid_score))])

        print("best_nu:" ,best_nu)


        # best_nu = []
        # def fn(h,nu):
        #     r_trees, child_order , r_node = nu_trees(nu, target_node, tree_path, status)
        #     model = phyloLL_HMM(n_components=status, trees= [clonal, r_trees[h]], model=GTR_sample, child_order=child_order, X_child_order=X_child_order)
        #     # model = phyloLL_HMM(n_components=status, trees=[clonal, r_trees[child_order.index(X_child_order[h])]], model=GTR_sample,  child_order=child_order, X_child_order=X_child_order)
        #     model.startprob_ = p_start
        #     model.transmat_ = p_trans
        #     s = model.score(X)
        #     # print(h, nu, s)
        #     return -s
        #
        # def fn0(nu):
        #     return fn(0,nu)
        # def fn1(nu):
        #     return fn(1,nu)
        # def fn2(nu):
        #     return fn(2,nu)
        #
        # result = spo.minimize_scalar(fn0, method="bounded", bounds=(0.001, 0.1) , options={'disp': 1})
        # best_nu.append(result.x)
        #
        # result1 = spo.minimize_scalar(fn1, method="bounded", bounds=(0.001, 0.1) , options={'disp': 1})
        # best_nu.append(result1.x)
        #
        # result2 = spo.minimize_scalar(fn2, method="bounded", bounds=(0.001, 0.1) , options={'disp': 1}) #bounded
        # best_nu.append(result2.x)


    if status == 8:
        # score = []
        def fn(nu):
            r_trees, child_order = nu_trees(nu, target_node, tree_path, status)
            recombination_trees = [clonal, r_trees[0], r_trees[1], r_trees[2], r_trees[3], r_trees[4], r_trees[5], r_trees[6]]
            model = phyloLL_HMM(n_components=status, trees=recombination_trees, model=GTR_sample, child_order=child_order, X_child_order=X_child_order)
            model.startprob_ = p_start
            model.transmat_ = p_trans
            score = model.score(X)
            print(nu,score)
            return -(score)


        # x = np.linspace(0, 1, 100)
        result = spo.minimize_scalar(fn, method="bounded", bounds=(0, 0.09), options={'disp': 1})
        # result = spo.minimize(fn, 0.001, method='Nelder-Mead' , bounds=my_nu) #, bounds=my_nu
        # result = spo.minimize(fn, 0.0001, method='trust-constr', options={'disp': True}) #, bounds=my_nu
        print("nu = {} , Score = {}".format(result.x, result.fun))
        best_nu = result.x

    print(best_nu)
    return best_nu
# *********************************************************************************************************************
def write_best_nu(best_nu,outputname):
    with open(outputname, mode='w') as bestnu_file:
        nu_writer = csv.writer(bestnu_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        nu_writer.writerow(best_nu)
# **********************************************************************************************************************
def phylohmm(tree,alignment_len,column,nu,p_start,p_trans,tips_num,status):
    mytree = []
    posterior0 = []
    posterior1 = []
    hiddenStates = []
    score = []
    tipdata = set_tips_partial(column,tips_num)
    r_node = []
    t_node = []
    single_posterior = []

    # print(tree.as_ascii_plot(show_internal_node_labels=True))

    best_nu = []

    for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
      # if target_node.index < 15:
        print(target_node.index)
        recombination_trees = []
        recombination_nodes = []
        child_order = []
        mytree.append(Tree.get_from_path(tree_path, 'newick'))
        # mytree[id_tree].reroot_at_midpoint(update_bipartitions=False)
        set_index(mytree[id_tree],alignment)

        # ----------- Step 1 : Make input for hmm ------------------------------------------------------
        # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

        mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
        recombination_trees.append(mytree[id_tree].as_string(schema="newick"))

        # --------------  Step 1.2: Calculate X based on this re-rooted tree
        X = make_hmm_input_mixture(mytree[id_tree],alignment_len,column,tipdata,GTR_sample,tips_num)
        # print(X[57500])


        # to keep the order of clonal tree children
        X_child_order = []
        for id, child in enumerate(target_node.child_node_iter()):
            X_child_order.append(child.index)

        my_nu = np.arange(0.01, 0.1, 0.01)
        if status == 2 :
            # -------------- find best nu ----------------------
            nu = give_best_nu(X, tree_path, recombination_trees[0], target_node, X_child_order,status,my_nu)
            best_nu.append([target_node.index, nu])
            # ----------- Step 2: make recombination trees -----------------------------------------------

            temptree = Tree.get_from_path(tree_path, 'newick')
            # temptree.reroot_at_midpoint(update_bipartitions=False)
            set_index(temptree,alignment)


            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)
            kids = temptree.seed_node.child_nodes()
            # print(kids)

            for k1, kid1 in enumerate(kids):
                child_order.append(kid1.index)  # keep the order of children after reroot
                temptree = Tree.get_from_path(tree_path, 'newick')
                # temptree.reroot_at_midpoint(update_bipartitions=False)
                set_index(temptree,alignment)
                filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
                target_node_temp = temptree.find_node(filter_fn=filter_fn)
                temptree.reroot_at_node(target_node_temp, update_bipartitions=False, suppress_unifurcations=True)

                # recombination_trees.append(recom_maker(temptree, kid1.index, nu))
                recombination_trees.append(recom_maker(temptree, kid1.index, nu[k1]))
                recombination_nodes.append(kid1)

            print("main function:",recombination_nodes)
            # ----------- Step 3: Call phyloHMM -----------------------------------------------------
            for h in range(1, len(recombination_trees)):
                model = phyloLL_HMM(n_components=status, trees=[recombination_trees[0], recombination_trees[h]],model=GTR_sample, child_order=child_order, X_child_order=X_child_order)
                model.startprob_ = p_start

                p_trans_nu0 = np.array([[1, 0],
                                        [1, 0]])

                if nu[h-1] <= my_nu[0]:
                    model.transmat_ = p_trans_nu0
                else:
                    model.transmat_ = p_trans

                # model.transmat_ = p_trans

                p = model.predict_proba(X)
                hidden = model.predict(X)
                posterior0.append(p[:,0])
                posterior1.append(p[:,1])
                hiddenStates.append(hidden)
                score.append(model.score(X))

                r_node.append(recombination_nodes[h - 1].index)
                t_node.append(target_node.index)

                tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
                # tree_updatePartial.reroot_at_midpoint(update_bipartitions=False)
                set_index(tree_updatePartial,alignment)
                # filter_fn = lambda n: hasattr(n, 'index') and n.index == X_child_order[h-1]
                filter_fn = lambda n: hasattr(n, 'index') and n.index == recombination_nodes[h - 1].index
                update_child = tree_updatePartial.find_node(filter_fn=filter_fn)
                if update_child.is_leaf():
                    update_mixture_partial_PSE(column, update_child, tipdata, p, status, 1,nu[h-1])
                    # update_mixture_partial(column,update_child, tipdata, p,status, 1)
                    # update_mixture_partial_PSE(column, update_child, tipdata, p, status, 1, nu)


        if status == 8 :
            # -------------- find best nu ----------------------
            nu = give_best_nu(X, my_nu, tree_path, recombination_trees[0], target_node, X_child_order,status)
            best_nu.append([target_node.index,nu])

            # ----------- Step 2: make recombination trees -----------------------------------------------
            temptree = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree,alignment)

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
                set_index(temptree,alignment)
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


            # ----------- Step 3: Call phyloHMM -----------------------------------------------------
            model = phyloLL_HMM(n_components=status, trees=recombination_trees, model=GTR_sample, child_order=child_order,X_child_order=X_child_order)
            model.startprob_ = p_start

            p_trans_nu0 = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0],
                                    [1, 0, 0, 0, 0, 0, 0, 0], ])


            if nu <= my_nu[0]:
                model.transmat_ = p_trans_nu0
            else:
                model.transmat_ = p_trans

            print(model.transmat_)

            # model.transmat_ = p_trans

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
            set_index(tree_updatePartial,alignment)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_partial = tree_updatePartial.find_node(filter_fn=filter_fn)
            for id, child in enumerate(target_node_partial.child_node_iter()):
                # print(child.index)
                if child.is_leaf():
                    order = X_child_order.index(child.index)
                    # print("my beloved child:", child.index , child.taxon , "order:" , order+1)
                    update_mixture_partial(column,child, tipdata, p,status, order + 1)


    if status == 2:
        np.set_printoptions(threshold=np.inf)
        myposterior0 = np.array(posterior0, dtype='double')
        myposterior1 = np.array(posterior1, dtype='double')
        recom_prob = pd.DataFrame({'recom_nodes': r_node, 'target_node': t_node, 'posterior0': pd.Series(list(myposterior0)) , 'posterior1': pd.Series(list(myposterior1))})


    if status == 8:
        recom_prob = pd.DataFrame({'recom_nodes': r_node, 'target_node': t_node, 'posterior': single_posterior})

    return tipdata,posterior0,posterior1,hiddenStates,score,recom_prob,r_node,t_node,best_nu
# **********************************************************************************************************************
def make_CATG_file(tips_num,alignment,alignment_len,tipdata,column,tree,outputname,error_flag):
    taxon = tree.taxon_namespace
    myfile = open(outputname, 'w')
    myfile.write(str(tips_num)+'\t'+str(alignment_len)+'\n')
    for idx, name in enumerate(alignment):
        myfile.write(str(name).replace("'", ""))
        myfile.write('\t')
    myfile.write('\n')
    for i in range(alignment_len):
        myfile.write(str(column[i]))
        myfile.write('\t')
        for j in range(tips_num):
            if error_flag == 1:
                if (tipdata[i,j,0] < 0.5):
                    tipdata[i, j, 0] = 0
                if (tipdata[i,j,1] < 0.5):
                    tipdata[i, j, 1] = 0
                if (tipdata[i,j,2] < 0.5):
                    tipdata[i, j, 2] = 0
                if (tipdata[i,j,3] < 0.5):
                    tipdata[i, j, 3] = 0
            temp = str(','.join(map(str, tipdata[i, j])))
            myfile.write(temp.strip())
            myfile.write('\t')
        myfile.write('\n')
    myfile.close()
# **********************************************************************************************************************
def make_physher_json_partial(tipdata,tree,json_path,outputname):
    my_tipdata = tipdata.transpose(1, 0, 2)
    with open(json_path) as json_file:
        data = json.load(json_file)
        taxon = []
        seq = []
        for i in range(my_tipdata.shape[0]):
            taxon.append(str(give_taxon(tree, i)))
            # seq.append(np.array2string(my_tipdata[i]))
            # seq.append(str(my_tipdata[i]).replace("[", "").replace("]", "").replace("\n", "").replace(" ", ","))
            seq.append(np.array2string(my_tipdata[i], separator=',').replace("[", "").replace("]", "").replace("\n", "").replace(" ", ""))

        partial = dict(zip(taxon, seq))
        data['model']['sitepattern']['partials'] = partial
        data['model']['tree']['newick'] = str(tree)

    jsonFile = open(outputname, "w")
    json.dump(data, jsonFile, indent=4)
    jsonFile.close()
# **********************************************************************************************************************


if __name__ == "__main__":

    # path = os.path.dirname(os.path.abspath(__file__))


    # path = '/home/nehleh/Desktop/sisters/mutiple_sisters/'
    # tree_path = path+'/num_1_RAxML_bestTree.tree'
    # genomefile = path+'/num_1_wholegenome_1.fasta'
    # baciSimLog = path+'/BaciSim_Log.txt'
    # clonal_path = path+'/clonaltree.tree'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'

    # path = '/home/nehleh/PhiloBacteria/Results_slides/num_4'
    # tree_path = path+'/num_4_recom_1_RAxML_bestTree.tree'
    # # tree_path = path+'/num_1_beasttree.newick'
    # clonal_path = path+'/num_4_Clonaltree.tree'
    # genomefile = path+'/num_4_recom_1_Wholegenome_4_1.fasta'
    # baciSimLog = path+'/num_4_recom_1_BaciSim_Log.txt'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--raxmltree", type=str, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, help='fasta file')
    parser.add_argument('-cl', "--clonaltreeFile", type=str, help='clonalclonaltreeFile tree from BaciSim')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-nu', "--nuHmm", type=float,default=0.033,help='nuHmm')
    parser.add_argument('-p', "--threshold", type=float, default=0.3, help='threshold')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-s', "--startProb", type=list, default= [0.99, 0.01],help='frequencies')
    parser.add_argument('-m', "--transmat", type=list, default= [[0.999, 0.001],  [0.001, 0.999]], help='rates')
    parser.add_argument('-st', "--status", type=str,  default='2', help='2 for the two states hmm and 8 for eight states of hmm , 2,8 for both ')
    parser.add_argument('-js', "--jsonFile", type=str, default='/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json', help='jsonFile')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    args = parser.parse_args()

    # parser.add_argument('-xml', "--xmlFile", type=str, default='/home/nehleh/PhiloBacteria/bin/template/GTR_template.xml' ,help='xmlFile')

    tree_path = args.raxmltree
    genomefile = args.alignmentFile
    json_path = args.jsonFile
    pi = args.frequencies
    rates = args.rates
    nu = args.nuHmm
    p_start = args.startProb
    p_trans = args.transmat
    threshold = args.threshold
    # xml_path = args.xmlFile
    initialstat = args.status
    simulation = args.simulation

    # ============================================ setting parameters ==================================================
    tree = Tree.get_from_path(tree_path, 'newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")

    # tree.reroot_at_midpoint(update_bipartitions=False)

    set_index(tree, alignment)

    nodes_number = len(tree.nodes())
    # print("nodes_number:",nodes_number)
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    set_index(tree,alignment)
    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)

    # print(column)
    # print(nodes_number)
    print(tree.as_ascii_plot(show_internal_node_labels=True))

    if initialstat.find('2') != -1:
        status = 2
        p_start = np.array([0.9, 0.1])
        p_trans = np.array([[1-(1/alignment_len), 1/alignment_len],
                            [1/alignment_len, 1-(1/alignment_len)]])
        # p_trans = np.array([[0.9999,0.0001],
        #                    [0.0001,0.9999]])

        start = time.time()
        tipdata, posterior0, posterior1, hiddenStates, score, recom_prob, r_node, t_node, best_nu = phylohmm(tree,alignment_len,column,nu, p_start, p_trans,tips_num,status)
        end = time.time()
        print("time phylohmm",end - start)

        # print(recom_prob)
        # print(tipdata[40000])

        # pd.options.display.max_colwidth = None
        start = time.time()
        recom_prob.to_hdf('Recom_prob_two.h5', key='prob' ,mode='w')
        end = time.time()
        print("time recom_prob.to_hdf",end - start)


        c_tree = Tree.get_from_path(tree_path, 'newick')
        set_index(c_tree,alignment)
        # c_tree.reroot_at_midpoint(update_bipartitions=False)
        internal_plot(c_tree,posterior0,posterior1,hiddenStates,score,r_node,t_node,status)
        write_best_nu(best_nu,'PB_nu_two.txt')
        # make_CATG_file(tips_num, alignment, alignment_len, tipdata, column, tree, 'PB_Two.catg', 0)
        start= time.time()
        make_physher_json_partial(tipdata, tree, json_path, 'PB_two.json')
        end = time.time()
        print("time make_physher_json_partial",end - start)


        # # # ======================================= providing xml files for beast ============================================
        # make_beast_xml_partial(tipdata, c_tree, xml_path,'PB_Partial_two.xml')
        # make_beast_xml_gap(tipdata, tree, xml_path, 0.5,'PB_Gap_two.xml')
        # make_beast_xml_delCol(recom_prob,tips_num,0.5,'PB_Del_two.xml')

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

        # print('status = 8')

        tipdata, posterior, hiddenStates, score, recom_prob, r_node, t_node, best_nu = phylohmm(tree,alignment_len,column,nu,p_start,p_trans,tips_num,status)



        c_tree = Tree.get_from_path(tree_path, 'newick')
        set_index(c_tree,alignment)
        # internal_plot(c_tree, posterior, hiddenStates, score, r_node, t_node, status)
        # phyloHMMData8 = recom_resultFig_dm(recom_prob,tips_num,threshold,status,'PB_Recom_eight.jpeg')
        # phyloHMM_log = phyloHMM_Log(c_tree, phyloHMMData8,'PB_Log_eight.txt')
        write_best_nu(best_nu,'PB_nu_eight.txt')
        # # ======================================= providing xml files for beast ============================================
        # make_beast_xml_partial(tipdata, c_tree, xml_path,'PB_Partial_eight.xml')
        # make_beast_xml_gap(tipdata, tree, xml_path, 0.5,'PB_Gap_eight.xml')
        # make_beast_xml_delCol(recom_prob,tips_num,0.5,'PB_Del_eight.xml')