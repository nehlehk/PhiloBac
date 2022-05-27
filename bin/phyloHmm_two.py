#!/usr/bin/env python
import concurrent.futures
import numpy as np
from utility import *
import pandas as pd
import dendropy
import argparse
import csv
import json
import scipy.optimize as spo
import time
from scipy.optimize import Bounds






def compute_logprob_phylo_bw(X,recom_trees,model,tip_partial,alignment_len):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        set_index(state_tree,alignment)
        persite_ll, partial = computelikelihood_mixture_bw(state_tree, alignment_len, tip_partial, model, tips_num)
        result[:, tree_id] = persite_ll

    return result
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
def computelikelihood_mixture(tree,alignment_len,tip_partial,model,tips_num):
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
def computelikelihood_mixture_bw(tree,alignment_len,tip_partial,model,tips_num):
    partial = np.zeros(((alignment_len,(2 * tips_num) - 1,4)))
    partial[:,0:tips_num,:] = tip_partial
    persite_ll = np.zeros(alignment_len)
    partial_new =  np.rollaxis(partial,2,0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[..., node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[..., children[0].index])
            for i in range(1, len(children)):
                partial_new[..., node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[..., children[i].index])

    persite_ll = (model.get_pi() @ partial_new[..., tree.seed_node.index])

    return persite_ll, partial
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
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
def update_mixture_partial_PSE(column,node,tipdata,posterior_recom,nu):
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    epsilon = nu * posterior_recom[site]
    # print(epsilon)
    for i in range(4):
        if i == my_number:
            tipdata[site,node.index,i] = 1 - (epsilon)
        else:
            tipdata[site,node.index,i] = (epsilon)/3
  return tipdata
# **********************************************************************************************************************
def my_best_nu(X,tree_path,clonal,target_node,tipdata,p_trans,p_start):
    def fn(nu):
        temptree = Tree.get_from_path(tree_path, 'newick', rooting='force-rooted')
        set_index(temptree, alignment)
        r_trees = recom_maker(temptree, target_node.index, nu)
        emission = compute_logprob_phylo_bw(X, [clonal, r_trees], GTR_sample, tipdata, alignment_len)
        alpha , c = forward(X, p_trans, emission, p_start)
        s = np.sum(np.log(c))
        # print(nu,s)
        return s

    result = spo.minimize_scalar(fn, method="bounded", bounds=(0.0, 0.5) , options={'disp': 1})
    return result.x
# *********************************************************************************************************************
def write_best_nu(best_nu,outputname):
    with open(outputname, mode='w') as bestnu_file:
        nu_writer = csv.writer(bestnu_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        nu_writer.writerow(best_nu)

# **********************************************************************************************************************
def baumwelch_parallel(target_node):
    print(target_node)
    recombination_trees = []
    clonaltree = Tree.get_from_path(tree_path,'newick',rooting='force-rooted')
    set_index(clonaltree, alignment)
    my_nu = np.arange(0.0001, 0.1, 0.01)
    recombination_trees.append(clonaltree.as_string(schema="newick"))
    X = partial[:, target_node.index, :]
    nu = my_best_nu(X,tree_path,recombination_trees[0],target_node,tipdata,p_trans,p_start)
    # print("nu:",nu)
    recombination_trees.append(recom_maker(clonaltree, target_node.index, nu))
    emission = compute_logprob_phylo_bw(X, recombination_trees, GTR_sample, tipdata, alignment_len)
    best_trans , p = baum_welch(X,p_trans,emission,p_start,n_iter=1)
    p = p.T
    # print("best_trans:",best_trans)
    # hidden = viterbi(X,best_trans,emission,p_start)
    p_trans_nu0 = np.array([[1, 0],[1, 0]])
    if nu <= my_nu[0]: # if the best nu is smaller than a threshod, we consider there is no recombiantion on that branch
        trans = p_trans_nu0
    else:
        trans = best_trans
    R_over_theta = -(np.log(trans[0][0]) - target_node.edge_length)

    if trans[1][1] == 0:
        delta = 0
    else:
        delta = -1/np.log(trans[1][1])
    # print("delta:", delta)
    return p,target_node,nu
# **********************************************************************************************************************
def phylohmm_baumwelch(tree,alignment_len,column,nu,p_start,p_trans,tips_num):
    mytree = []
    tipdata = set_tips_partial(column,tips_num)
    posterior0 = []
    posterior1 = []
    t_node = []
    best_nu = []*18

    # make per_site likelihood and per_site partials for all nodes including tips and internal nodes
    persite_ll, partial = computelikelihood_mixture(tree, alignment_len, tipdata, GTR_sample, tips_num)

    # each node play the role of target nodes
    for id_tree, target_node in enumerate(tree.postorder_node_iter()):
        if target_node != tree.seed_node:
        #     print(target_node)
            recombination_trees = []
            mytree.append(Tree.get_from_path(tree_path, 'newick', rooting='force-rooted'))
            set_index(mytree[id_tree],alignment)
            #     step 1 --- make hmm input
            #  take the partials of the target node as input of hmm
            X = partial[:,target_node.index,:]

            #     step 2 --- make recombination tree
            my_nu = np.arange(0.0001, 0.1, 0.01)
            recombination_trees.append(mytree[id_tree].as_string(schema="newick"))
            # find the best nu for target branch based on the maximizing score value of hmm
            nu = my_best_nu(X,tree_path,recombination_trees[0],target_node,tipdata,p_trans,p_start)
            best_nu[target_node.index] = (target_node.index, nu)
            # print("nu_first:",nu)

            # make recombination tree for target node using best nu
            recombination_trees.append(recom_maker(mytree[id_tree],target_node.index,nu))
            emission = compute_logprob_phylo_bw(X, recombination_trees, GTR_sample, tipdata, alignment_len)
            best_trans , p = baum_welch(X, p_trans, emission, p_start, n_iter=1)
            p = p.T
            # print("first:",best_trans)

            hidden = viterbi(X,best_trans,emission,p_start)
            p_trans_nu0 = np.array([[1, 0],
                                    [1, 0]])
            if nu <= my_nu[0]: # if the best nu is smaller than a threshold, we consider there is no recombination on that branch
                trans = p_trans_nu0
            else:
                trans = best_trans

            R_over_theta = -(np.log(trans[0][0]) - target_node.edge_length)
            # print("R_over_theta:",R_over_theta)
            if trans[1][1] == 0:
                delta = 0
            else:
                delta = -1/np.log(trans[1][1])
            # print("delta:",delta)

            posterior0.append(p[:, 0])  # posterior probability for clonal tree
            posterior1.append(p[:, 1])  # posterior probability for recombination tree
            t_node.append(target_node.index)
            # Update tip partials based on the posterior probability
            if target_node.is_leaf():
                update_mixture_partial_PSE(column, target_node,tipdata,p[:, 1],nu)



    np.set_printoptions(threshold=np.inf)
    myposterior0 = np.array(posterior0, dtype='double')
    myposterior1 = np.array(posterior1, dtype='double')
    # print(myposterior1)
    # recom_prob = pd.DataFrame( {'recom_nodes': t_node, 'posterior0': pd.Series(list(myposterior0)),'posterior1': pd.Series(list(myposterior1))})
    recom_prob = pd.DataFrame({'recom_nodes': t_node, 'posterior_rec': pd.Series(list(myposterior1))})
    write_best_nu(best_nu,'PB_nu_two.txt')

    return tipdata,recom_prob,myposterior1,best_nu
# **********************************************************************************************************************
def make_physher_json_partial(tipdata,tree,json_path,outputname):
    my_tipdata = tipdata.transpose(1, 0, 2)
    with open(json_path) as json_file:
        data = json.load(json_file)
        taxon = []
        seq = []
        for i in range(my_tipdata.shape[0]):
            taxon.append(str(give_taxon(tree, i)))
            seq.append(np.array2string(my_tipdata[i], separator=',').replace("[", "").replace("]", "").replace("\n", "").replace(" ", ""))

        partial = dict(zip(taxon, seq))
        data['model']['sitepattern']['partials'] = partial
        data['model']['tree']['newick'] = str(tree)

    jsonFile = open(outputname, "w")
    json.dump(data, jsonFile, indent=4)
    jsonFile.close()
# **********************************************************************************************************************
def make_physher_json_partial_1(partial_dict,tree,json_path,outputname):
    with open(json_path) as json_file:
        data = json.load(json_file)
        data['model']['sitepattern']['partials'] = partial_dict
        data['model']['tree']['newick'] = str(tree)

    jsonFile = open(outputname, "w")
    json.dump(data, jsonFile, indent=4)
    jsonFile.close()
# **********************************************************************************************************************
def forward(X, trans, emission, initial_distribution):
    alpha = np.zeros((X.shape[0], trans.shape[0]))
    c = np.zeros(X.shape[0])
    alpha[0, :] = initial_distribution * emission[0]
    c[0] = 1/alpha[0].sum(axis=0)
    alpha[0] *= c[0]
    for t in range(1, X.shape[0]):
        for j in range(trans.shape[0]):
            alpha[t, j] = alpha[t - 1].dot(trans[:, j]) * emission[t, j]
        c[t] = 1/alpha[t].sum(axis=0)
        alpha[t] *= c[t]

    return alpha , c
# **********************************************************************************************************************
def backward(X, trans, emission,c):
    beta = np.zeros((X.shape[0], trans.shape[0]))
    # setting beta(T) = 1
    beta[X.shape[0] - 1] = np.ones((trans.shape[0]))

    # Loop in backward way from T-1 to
    for t in range(X.shape[0] - 2, -1, -1):
        for j in range(trans.shape[0]):
            beta[t, j] = (beta[t + 1] * emission[t + 1, :]).dot(trans[j, :])
        beta[t] *= c[t]

    return beta
# **********************************************************************************************************************
def baum_welch(X, trans, emission, initial_distribution, n_iter=1):
    M = trans.shape[0]
    T = len(X)

    for n in range(n_iter):
        recombination_trees = []
        alpha,c  = forward(X, trans, emission, initial_distribution)
        beta = backward(X, trans, emission,c)

        gamma = np.zeros((M, T - 1))
        xi = np.zeros((M, M, T - 1))
        for t in range(T - 1):
            gamma[:, t] = (alpha[t, :] * beta[t, :]) / np.dot(alpha[t, :], beta[t, :])
            denominator = np.dot(np.dot(alpha[t, :].T, trans) * emission[t + 1].T, beta[t + 1, :])
            for i in range(M):
                numerator = alpha[t, i] * trans[i, :] * emission[t + 1].T * beta[t + 1, :].T
                xi[i, :, t] = numerator / denominator

        trans = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))

        # Add additional T'th element in gamma
        gamma = np.hstack((gamma , (alpha[T-1, :] * beta[T-1, :] / np.dot(alpha[T-1, :] , beta[T-1, :])).reshape((-1, 1)) ))



    return trans,gamma
# **********************************************************************************************************************
def viterbi(X, trans, emission, initial_distribution):
    T = X.shape[0]
    M = trans.shape[0]

    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * emission[0])

    prev = np.zeros((T - 1, M))

    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(trans[:, j]) + np.log(emission[t,j])

            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)

            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)

    # Path Array
    S = np.zeros(T)

    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])

    S[0] = last_state

    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1

    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)

    return S
# **********************************************************************************************************************
def computelikelihood(brs, nu, posterior):
    # brs: [B] B: number of branches
    # nu: [B] B: number of branches
    # posterior: [B,N] N: number of sites
    # P [B,4,4]

    partial = np.empty((4,nodes_number, alignment_len))
    partial[:,0:tips_num,:] = tipdata.T
    P = GTR_sample.p_t(brs)
    P_nu = GTR_sample.p_t(brs + nu)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            if children[0].is_leaf():
                partial[:, node.index, :] = np.dot(P[children[0].index], partial[:,children[0].index, :])
            else:
                post = np.expand_dims(np.expand_dims(posterior[children[0].index], -1), -1)
                P_mixture = ((1 - post) * P[children[0].index] + post * P_nu[children[0].index])
                partial[:, node.index, :] = np.squeeze(P_mixture @ np.expand_dims(partial[:, children[0].index, :].T, -1), -1).T
            for i in range(1, len(children)):
                if children[i].is_leaf():
                    partial[:, node.index, :] *= np.dot(P[children[i].index], partial[:,children[i].index, :])
                else:
                    post = np.expand_dims(np.expand_dims(posterior[children[i].index], -1), -1)
                    P_mixture = ((1 - post) * P[children[i].index] + post * P_nu[children[i].index])
                    partial[:, node.index, :] *= np.squeeze(P_mixture @ np.expand_dims(partial[:, children[i].index, :].T, -1), -1).T

    persite_ll = np.log(GTR_sample.get_pi() @ partial[:, tree.seed_node.index, :])
    return persite_ll.sum()
# **********************************************************************************************************************
if __name__ == "__main__":

    # path = '/home/nehleh/Desktop/sisters/mutiple_sisters/'
    # tree_path = path+'/num_1_RAxML_bestTree.tree'
    # genomefile = path+'/num_1_wholegenome_1.fasta'
    # baciSimLog = path+'/BaciSim_Log.txt'
    # clonal_path = path+'/clonaltree.tree'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'

    path = '/Users/mathieu/Downloads'
    tree_path = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_RAxML_bestTree.tree'
    # tree_path = path+'/num_1_beasttree.newick'
    clonal_path = path+'/num_1_Clonaltree.tree'
    genomefile = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_Wholegenome.fasta'
    # baciSimLog = path+'/num_1_nu_0.07_Rlen_500_Rrate_0.02_BaciSim_Log.txt'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--raxmltree", type=str, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, help='fasta file')
    parser.add_argument('-cl', "--clonaltreeFile", type=str, help='clonalclonaltreeFile tree from BaciSim')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-nu', "--nuHmm", type=float,default=0.033,help='nuHmm')
    parser.add_argument('-p', "--threshold", type=float, default=0.5, help='threshold')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-s', "--startProb", type=list, default= [0.99, 0.01],help='frequencies')
    parser.add_argument('-m', "--transmat", type=list, default= [[0.999, 0.001],  [0.001, 0.999]], help='rates')
    parser.add_argument('-st', "--status", type=str,  default='2', help='2 for the two states hmm and 8 for eight states of hmm , 2,8 for both ')
    parser.add_argument('-js', "--jsonFile", type=str, default='/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json', help='jsonFile')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    args = parser.parse_args()

    # tree_path = args.raxmltree
    # genomefile = args.alignmentFile
    # json_path = args.jsonFile
    pi = args.frequencies
    rates = args.rates
    nu = args.nuHmm
    p_start = args.startProb
    p_trans = args.transmat
    threshold = args.threshold
    simulation = args.simulation


    # ============================================ setting parameters ==================================================
    tree = Tree.get_from_path(tree_path,schema='newick',rooting='force-rooted')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile),schema="fasta")
    set_index(tree, alignment)
    print(tree.as_ascii_plot(show_internal_node_labels=True))

    nodes_number = len(tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)


    p_start = np.array([0.9, 0.1])
    p_trans = np.array([[1-(1/alignment_len), 1/alignment_len],
                        [1/alignment_len, 1-(1/alignment_len)]])


    # doing hmm
    start = time.time()
    tipdata, recom_prob, posterior, best_nu = phylohmm_baumwelch(tree, alignment_len, column, nu, p_start, p_trans, tips_num)

    brs = [None]*18
    for _, node in enumerate(tree.postorder_node_iter()):
        brs[node.index] = node.edge_length
    brs.pop()

    # def max_LL():
    #     initial_guess = [0.01]
    #     bounds = Bounds([0], [0.1])
    #     result = spo.minimize(computelikelihood, initial_guess, method='trust-constr', bounds=bounds, options={'verbose': 1})
    #     print("edge_length = {} " "likelihood = {}".format(result.x[0], -result.fun))
    #     return result
    #
    #
    # max_LL()
        
    initial_guess = np.asarray(brs)
    # posterior = np.zeros((17, 5000))
    # best_nu = np.zeros(17)
    posterior = np.asarray(posterior)
    best_nu = np.asarray(best_nu)[:,1]
    for node in tree.postorder_node_iter():
        if node.index < len(best_nu):
            print(f"{node.index} {best_nu[node.index]}")
    # tipdata = set_tips_partial(column, tips_num)
    print(brs)
    print(best_nu)
    print(initial_guess.shape)
    bounds = [[1.e-20, 0.1]]*len(brs)
    print(best_nu.shape, posterior.shape)
    fn = lambda x: -computelikelihood(np.exp(x), best_nu, posterior)
    print('LnL', -fn(np.log(initial_guess)))
    # print('LnL', -fn(np.log(initial_guess+10.)))
    # print(GTR_sample.p_t(np.array([1.0])))
    # print(GTR_sample.p_t(np.array([.1])))
    # print(GTR_sample.p_t(np.array([.000000001])))
    # print(GTR_sample.p_t(np.array([1., 0.1, .000000001])))
    # exit(2)
    # result = spo.minimize(fn, np.log(initial_guess), method='CG', options={'ftol': 0.0000001, 'maxiter': 1000})
    print(bounds)
    result = spo.minimize(lambda x: -computelikelihood(x, best_nu, posterior),
        initial_guess, method='TNC', bounds=bounds, options={'verbose': 1})

    print("edge_length = {} " "likelihood = {}".format(result.x[0], -result.fun))
    # print(np.exp(result.x))
    print(result.x)
    for i, j in zip(initial_guess, result.x):
        print(f"{i} {j} {i-j} {i/j}")


    # computelikelihood(0.1)


    # posterior_rec = []
    # nodes = []
    # tnodes= []
    # best_nu = []
    # np.set_printoptions(threshold=np.inf)
    # tipdata = set_tips_partial(column, tips_num)
    # persite_ll, partial = computelikelihood_mixture(tree, alignment_len, tipdata, GTR_sample, tips_num)
    #
    # for node in tree.postorder_node_iter():
    #     if node != tree.seed_node:
    #         nodes.append(node)
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = [executor.submit(baumwelch_parallel, target_node) for target_node in nodes]
    #     for res in concurrent.futures.as_completed(results):
    #         p = res.result()[0]
    #         target_node = res.result()[1]
    #         nu = res.result()[2]
    #         posterior_rec.append(np.array(p[:, 1], dtype='double')) # posterior probality for recombination tree
    #         tnodes.append(target_node.index)
    #         best_nu.append([target_node.index, nu])
    #         if target_node.is_leaf():
    #             update_mixture_partial_PSE(column,target_node,tipdata,p[:, 1],nu)
    #
    #
    # # print(tipdata[36700])
    #
    # recom_prob = pd.DataFrame( {'recom_nodes': tnodes,'posterior_rec': pd.Series(list(posterior_rec))})
    # write_best_nu(best_nu, 'PB_nu_two.txt')
    # # end = time.time()
    # # # print(recom_prob)
    # # print("time phylohmm", end - start)
    #
    #
    # start = time.time()
    # my_tipdata = tipdata.transpose(1, 0, 2)
    # taxon = []
    # seq = []
    # def partial_dic(i):
    #     taxon = str(give_taxon(tree, i))
    #     seq = np.array2string(my_tipdata[i], separator=',').replace("[", "").replace("]", "").replace("\n","").replace(" ","")
    #     return taxon,seq
    #
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     f1 = [executor.submit(partial_dic, i) for i in range(my_tipdata.shape[0])]
    #     for res in concurrent.futures.as_completed(f1):
    #         taxon.append(res.result()[0])
    #         seq.append(res.result()[1])
    # partial_dict = dict(zip(taxon,seq))
    # make_physher_json_partial_1(partial_dict, tree, json_path, 'PB_two.json')
    # end = time.time()
    # print("time_make_physher_json_partial:", end - start)
    #
    #
    # start = time.time()
    # #  write updated tip partail to make recombination result based on that
    # recom_prob.to_hdf('Recom_prob_two.h5',key='recom_prob',mode='w')
    # end = time.time()
    # print("time recom_prob.to_hdf",end - start)