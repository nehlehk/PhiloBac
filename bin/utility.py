#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
from dendropy.simulate import treesim
import csv


def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
# **********************************************************************************************************************
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
          node.index = node.taxon.label
          node.label = str(node.index)
# **********************************************************************************************************************
def set_label(tree):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
# **********************************************************************************************************************
def give_taxon(tree,index):
    for node in tree.postorder_node_iter():
        if int(node.index) == index:
            return int(node.taxon.label)
# **********************************************************************************************************************
def give_taxon_index(tree,taxa,tips_num):
    node_mapping = np.zeros((tips_num,2))
    i = 0
    for node in tree.postorder_node_iter():
      if node.is_leaf():
        node_mapping[i][0] = int(str(node.taxon.label))
      # else:
      #   node_mapping[i][0] = int(str(node.index))
        node_mapping[i][1] = int(str(node.index))
        i = i+1

    # print(node_mapping)
    for i in range(len(node_mapping)):
      if int(node_mapping[i][0]) == int(taxa):
        return node_mapping[i][1]
# **********************************************************************************************************************
def give_descendents(tree,node_index,result,tips_num):
  if node_index >= tips_num:
    internal_recom_node = tree.find_node_with_label(str(node_index))
    if not (internal_recom_node is None):
        children = internal_recom_node.child_nodes()
        for n in range(len(children)):
          r_node= int(children[n].index)
          if r_node >= tips_num:
            give_descendents(tree,r_node,result,tips_num)
          else:
            result.add(give_taxon(tree,r_node))
  return result
# **********************************************************************************************************************
def nodes_separation(nodes):
  nodes = str(nodes)
  nodes = nodes[1:-1]
  mynodes = nodes.split(",")
  return mynodes
# **********************************************************************************************************************
def give_equivalent_node(recomtree):
  e = []
  for edge in recomtree.postorder_edge_iter():
    if edge.length is None:
      edge.length = 0
    e.append(edge.length)

  m = np.max(e)
  for node in recomtree.postorder_node_iter():
    if node.edge_length == m and node.is_leaf():
      return node.label,m
    elif node.edge_length == m and node.is_internal():
      return node.label,node.distance_from_tip() + m
# **********************************************************************************************************************
def remove_internal_labels(strtree,tips_num):
  tree = Tree.get_from_string(strtree,schema='newick')
  for node in tree.postorder_node_iter():
    if not node.label is None:
      if (int(node.label)>= tips_num):
        node.label = None

  return tree.as_string(schema="newick")
# **********************************************************************************************************************
def make_taxa(tips_num):
    taxon_list= []
    for i in range(tips_num):
      taxon_list.append(str(i))
    taxa = dendropy.TaxonNamespace(taxon_list)
    return taxa
# **********************************************************************************************************************
def make_clonaltree(tips_num,max_tMRCA):
    taxon_list= []
    for i in range(tips_num):
      taxon_list.append(str(i))
    taxa = dendropy.TaxonNamespace(taxon_list)
    tree = treesim.pure_kingman_tree(taxon_namespace=taxa,pop_size=3)
    t = tree.max_distance_from_root()
    normal_co = t/ max_tMRCA
    tree.scale_edges(1/normal_co)
    clonal_tree = tree.as_string(schema="newick")
    clonal_tree = clonal_tree.replace('\n',"")
    myfile = open('./Clonaltree.tree', 'w')
    myfile.write(clonal_tree)
    myfile.close()

    return clonal_tree
# **********************************************************************************************************************
def real_recombination(recomLog,clonaltree,nodes_number,alignment_len,tips_num):
    realData = np.zeros((alignment_len,nodes_number))
    df = pd.read_csv(recomLog,sep='\t', engine='python')
    recom = df.loc[df['status'] != 'clonal']
    recom = recom.reset_index(drop=True)
    for i in range(len(recom)):
        s = recom['start'][i]
        t = recom['end'][i]
        nodes = nodes_separation(recom['nodes'][i])
        for i in range(len(nodes)):
            mynode = int(nodes[i])
            if mynode < tips_num:
                mynode = int(give_taxon_index(clonaltree, mynode,tips_num))
            realData[s:t,mynode] = 1

    return realData
# **********************************************************************************************************************
def write_rmse(rmse_value,filename):
    with open(filename, mode='w') as rmse_file:
        rmse_writer = csv.writer(rmse_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        rmse_writer.writerow([rmse_value])
# **********************************************************************************************************************
def my_mrca(tree,tips):
  pdm = tree.phylogenetic_distance_matrix()
  taxon = tree.taxon_namespace
  node0 = [i for i,x in enumerate(taxon) if x.label==str(tips[0])]
  node1 = [i for i,x in enumerate(taxon) if x.label==str(tips[len(tips)-1])]
  myMrca = pdm.mrca(taxon[node0[0]], taxon[node1[0]])

  return myMrca.index
# **********************************************************************************************************************