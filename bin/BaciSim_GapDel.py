#!/usr/bin/env python

import numpy as np
import pandas as pd
from utility import *
import argparse
from dendropy import Tree
import dendropy
import xml.etree.ElementTree as ET
import json
import os.path


def make_fasta_del(realData,clonaltree,alignment):
    del_col = []
    for i in range(realData.shape[0]):
        for j in range(realData.shape[1]):
            if realData[i][j] == 1 :
              del_col.append(j)

    taxon = []
    seq = []
    for i in range(realData.shape[0]):
        x = ''
        # taxon.append(str(i))
        if i < tips_num:
            taxon.append(str(give_taxon(clonaltree, i)))
            for j in range(realData.shape[1]):
              if not (j in del_col):
                x = x + str(alignment[i][j])
            seq.append(x)

    with open("Del_alignment.fasta", "w") as fasta_file:
        for i in range(len(taxon)):
            fasta_file.write(">" + taxon[i] + "\n" + seq[i] + "\n")
    fasta_file.close()
# **********************************************************************************************************************
def make_fasta_gap(realData,clonaltree,alignment):
    taxon = []
    seq = []
    for i in range(realData.shape[0]):
        if i < tips_num:
            x = ''
            taxon.append(str(give_taxon(clonaltree, i)))
            for j in range(realData.shape[1]):
                if realData[i][j]==1:
                    x = x + '-'
                else:
                    x = x + str(alignment[i][j])
            seq.append(x)

    for i in range(tips_num,realData.shape[0]):
        desc = set()
        desc_nodes = list(give_descendents(clonaltree, i, desc,tips_num))
        # print(i, ":", desc_nodes)
        for k in range(len(desc_nodes)):
            d = int(desc_nodes[k])
            # print(d)
            kid = int(give_taxon_index(clonaltree, d,tips_num))
            # print(kid)
            x = list(seq[kid])
            for j in range(realData.shape[1]):
                if realData[i][j]==1:
                    x[j] = '-'
            seq[kid] = (''.join(x))


    with open("Gap_alignment.fasta", "w") as fasta_file:
        for i in range(len(taxon)):
            fasta_file.write(">" + taxon[i] + "\n" + seq[i] + "\n")
    fasta_file.close()
# **********************************************************************************************************************
def make_xml_seq(xml_path,clonaltree):
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(tips_num):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon", str(give_taxon(clonaltree, i)))
        c.text = '\n' + str(alignment[i]) + '\n'
        data.insert(i, c)
        c.tail = "\n"

    my_xml.write('OriginalSeq.xml', encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_xml_gap(realData,xml_path,clonaltree):
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(realData.shape[0]):
        x = ''
        c = ET.Element("sequence")
        # c.set("taxon", str(i))
        if i < tips_num:
            c.set("taxon", str(give_taxon(clonaltree, i)))
            for j in range(realData.shape[1]):
                if realData[i][j] == 1:
                    x = x + '-'
                else:
                    x = x + str(alignment[i][j])
            c.text = '\n' + x +'\n'
            data.insert(i,c)
            c.tail = "\n"

    my_xml.write('BaciSim_Gap.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_xml_del(realData,xml_path,clonaltree):
    del_col = []
    for i in range(realData.shape[0]):
        for j in range(realData.shape[1]):
            if realData[i][j] == 1 :
              del_col.append(j)

    # print(del_col)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(realData.shape[0]):
        x = ''
        c = ET.Element("sequence")
        if i < tips_num:
            c.set("taxon" , str(give_taxon(clonaltree,i)))
            # c.set("taxon", str(i))
            for j in range(realData.shape[1]):
              if not (j in del_col):
                x = x + str(alignment[i][j])
            c.text = '\n' + x +'\n'
            data.insert(i,c)
            c.tail = "\n"

    my_xml.write('BaciSim_Del.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_xml_partial(realData,xml_path,clonaltree,gap_prob,normal_prob):
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for i in range(realData.shape[0]):
        x = ''
        c = ET.Element("sequence")
        if i < tips_num:
            c.set("taxon" , str(give_taxon(clonaltree,i)))
            c.set("uncertain" , "true")
            for j in range(realData.shape[1]):
                if realData[i][j] == 1:
                    partial[j][i][0:4] = gap_prob
                else:
                    partial[j][i] = normal_prob
                k = give_index(str(alignment[i][j]))
                partial[j][i][k] = 1
                x = x + str(repr(partial[j][i]))[7:-2] + ';'
            c.text = '\n' + x +'\n'
            data.insert(i,c)
            c.tail = "\n"

    my_xml.write('BaciSim_partial.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_xml_partial_certain(realData,xml_path,clonaltree):
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for i in range(realData.shape[0]):
        x = ''
        c = ET.Element("sequence")
        if i < tips_num:
            c.set("taxon" , str(give_taxon(clonaltree,i)))
            c.set("uncertain" , "true")
            for j in range(realData.shape[1]):
                if realData[i][j] == 1:
                    partial[j][i][0:4] = [1,1,1,1]
                else:
                    partial[j][i] = [0,0,0,0]
                k = give_index(str(alignment[i][j]))
                partial[j][i][k] = 1
                x = x + str(repr(partial[j][i]))[7:-2] + ';'
            c.text = '\n' + x +'\n'
            data.insert(i,c)
            c.tail = "\n"

    my_xml.write('BaciSim_partial_certian.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_physher_json_partial(tipdata,tree,json_path,outputname,error_flag):
    my_tipdata = tipdata.transpose(1, 0, 2)
    # print(my_tipdata.shape)
    with open(json_path) as json_file:
        data = json.load(json_file)
        taxon = []
        seq = []
        for i in range(my_tipdata.shape[0]):
            x = ''
            taxon.append(str(give_taxon(tree, i)))
            for j in range(my_tipdata.shape[1]):
                temp = ",".join(str(x) for x in (my_tipdata[i, j, :]))
                if j == my_tipdata.shape[1]:
                    x = x + temp
                else:
                    x = x + temp + ','
            seq.append(x)

        partial = dict(zip(taxon, seq))
        data['model']['sitepattern']['partials'] = partial
        data['model']['tree']['newick'] = str(tree)

    jsonString = json.dumps(data, indent=4)
    jsonFile = open(outputname, "w")
    jsonFile.write(jsonString)
    jsonFile.close()
# **********************************************************************************************************************
def make_json_partial(realData,json_path,raxml_tree,gap_prob,normal_prob,output):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    with open(json_path) as json_file:
        data = json.load(json_file)
        taxon = []
        seq = []
        for i in range(realData.shape[0]):
            x = ''
            if i < tips_num:
                taxon.append(str(give_taxon(raxml_tree, i)))
                for j in range(realData.shape[1]):
                    if realData[i][j] == 1:
                        partial[j][i][0:4] = gap_prob
                    else:
                        partial[j][i] = normal_prob
                    k = give_index(str(alignment[i][j]))
                    partial[j][i][k] = 1
                    x = x + str(repr(partial[j][i]))[7:-2] + ','
                seq.append(x)

        partial = dict(zip(taxon, seq))
        data['model']['sitepattern']['partials'] = partial
        data['model']['tree']['newick'] = str(raxml_tree)

        # align = dict(zip(taxon, seq))
        # data['model']['sitepattern']['alignment']['sequences'] =  align
        # data['model']['tree']['newick'] = str(raxml_tree)

    jsonString = json.dumps(data, indent=4)
    jsonFile = open(output, "w")
    jsonFile.write(jsonString)
    jsonFile.close()
# **********************************************************************************************************************
def make_json_partial_PSE(realData,json_path,raxml_tree,gap_prob,normal_prob,nu,output):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    with open(json_path) as json_file:
        data = json.load(json_file)
        taxon = []
        seq = []
        for i in range(realData.shape[0]):
            x = ''
            if i < tips_num:
                taxon.append(str(give_taxon(raxml_tree, i)))
                for j in range(realData.shape[1]):
                    if realData[i][j] == 1:
                        partial[j][i][0:4] = (gap_prob[0:4] * nu)/3
                        k = give_index(str(alignment[i][j]))
                        partial[j][i][k] = 1 - (gap_prob[0] * nu)
                    else:
                        partial[j][i] = (normal_prob[0:4] * nu)/3
                        k = give_index(str(alignment[i][j]))
                        partial[j][i][k] = 1 - (normal_prob[0] * nu)
                    x = x + str(repr(partial[j][i]))[7:-2] + ','
                seq.append(x)

        partial = dict(zip(taxon, seq))
        data['model']['sitepattern']['partials'] = partial
        data['model']['tree']['newick'] = str(raxml_tree)
        # align = dict(zip(taxon, seq))
        # data['model']['sitepattern']['alignment']['sequences'] =  align
        # data['model']['tree']['newick'] = str(raxml_tree)

    jsonString = json.dumps(data, indent=4)
    jsonFile = open(output, "w")
    jsonFile.write(jsonString)
    jsonFile.close()
# **********************************************************************************************************************
def make_CATG_file_certain(realData,tips_num,alignment_len,column,tree,outputname):
    taxon = tree.taxon_namespace
    myfile = open(outputname, 'w')
    myfile.write(str(tips_num)+'\t'+str(alignment_len)+'\n')
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for i in range(len(taxon)):
        myfile.write(str(taxon[i].label))
        myfile.write('\t')
    myfile.write('\n')
    for i in range(alignment_len):
        # print(column[i])
        myfile.write(str(column[i]))
        myfile.write('\t')
        for j in range(tips_num):
            if realData[j][i] == 1:
                temp = str(','.join(map(str, [1,1,1,1])))
            else:
                k = give_index(str(column[i][j]))
                partial[i][j][k] = 1
                # print(partial[i][j])
                temp = str(','.join(map(str, partial[i][j])))
            myfile.write(temp.strip())
            myfile.write('\t')
        myfile.write('\n')
    myfile.close()
# **********************************************************************************************************************
def make_CATG_file_partial(realData,tips_num,alignment_len,column,tree,outputname,gap_prob,normal_prob):
    taxon = tree.taxon_namespace
    myfile = open(outputname, 'w')
    myfile.write(str(tips_num)+'\t'+str(alignment_len)+'\n')
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for i in range(len(taxon)):
        myfile.write(str(taxon[i].label))
        myfile.write('\t')
    myfile.write('\n')
    for i in range(alignment_len):
        # print(column[i])
        myfile.write(str(column[i]))
        myfile.write('\t')
        for j in range(tips_num):
            if realData[j][i] == 1:
                # temp = str(','.join(map(str, gap_prob)))
                partial[i][j] = gap_prob
            else:
                partial[i][j] = normal_prob
                k = give_index(str(column[i][j]))
                partial[i][j][k] = 1
            # print(partial[i][j])
            temp = str(','.join(map(str, partial[i][j])))
            # print(temp)
            myfile.write(temp.strip())
            myfile.write('\t')
        myfile.write('\n')
    myfile.close()
# **********************************************************************************************************************




if __name__ == "__main__":


    # tree_path = '/home/nehleh/Desktop/examples/num_4/num_4_Clonaltree.tree'
    # genomefile = '/home/nehleh/Desktop/examples/num_4/num_4_recom_1_Wholegenome_4_1.fasta'
    # xml_path = '/home/nehleh/PhyloCode/RecomPhyloHMM/bin/GTR_template.xml'
    # recomLog = '/home/nehleh/Desktop/examples/num_4/num_4_recom_1_BaciSim_Log.txt'
    # json_path = '/home/nehleh/PhyloCode/RecomPhyloHMM/bin/GTR_template.json'
    # raxml_path = '/home/nehleh/Desktop/examples/num_4/num_4_recom_1_RAxML_bestTree.tree'
    #
    path = os.getcwd()

    parser=argparse.ArgumentParser(description='''You did not specify any parameters. ''',epilog="""All's well that ends well.""")



    parser.add_argument('-t', "--clonaltree", type=str, default=path+'/Clonaltree.tree' , help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, default= path+'/', help='fasta file')
    parser.add_argument('-r', "--raxmltree", type=str, default= path+'/', help='raxmltree')
    parser.add_argument('-l', "--recomlogFile", type=str, default=path+'/BaciSim_Log.txt'  ,help='recombination log file')
    parser.add_argument('-x', "--xmlFile", default= path+'/template/GTR_template.xml' , type=str, help='xmlFile')
    parser.add_argument('-j', "--jsonFile", default= path+'/template/GTR_temp_partial.json', type=str, help='jsonFile')

    args = parser.parse_args()

    tree_path = args.clonaltree
    recomLog = args.recomlogFile
    genomefile = args.alignmentFile
    raxml_path = args.raxmltree
    xml_path = args.xmlFile
    json_path = args.jsonFile



    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    clonal_tree = Tree.get_from_path(tree_path, 'newick')
    raxml_tree = Tree.get_from_path(raxml_path, 'newick')
    tips_num = len(alignment)
    nodes_number = len(clonal_tree.nodes())
    column = get_DNA_fromAlignment(alignment)


    alignment_len = alignment.sequence_size
    set_index(clonal_tree,alignment)
    set_index(raxml_tree,alignment)

    # print(clonal_tree.as_ascii_plot(show_internal_node_labels=True))

    realData = real_recombination(recomLog,clonal_tree,nodes_number,alignment_len,tips_num)
    realData = realData.transpose()


    make_fasta_gap(realData, clonal_tree, alignment)
    make_fasta_del(realData, clonal_tree, alignment)

    # make_xml_seq(xml_path,clonal_tree)
    # make_xml_gap(realData, xml_path, clonal_tree)
    # make_xml_del(realData, xml_path, clonal_tree)
    # make_CATG_file_certain(realData, tips_num, alignment_len, column, clonal_tree, 'BaciSim_Best_Two.catg')


    gap_prob = np.array([0.999,0.999,0.999,0.999])
    normal_prob = np.array([0.0001,0.0001,0.0001,0.0001])
    nu = 0.08

    # make_xml_partial(realData, xml_path, clonal_tree,gap_prob,normal_prob)
    # make_CATG_file_partial(realData, tips_num, alignment_len, column, clonal_tree, 'BaciSim_partial_Two_00001_99_.catg', gap_prob, normal_prob)
    # make_xml_partial_certain(realData, xml_path, clonal_tree)

    # make_json_partial(realData,json_path,raxml_tree,gap_prob,normal_prob,"BaciSim_partial.json")

    make_json_partial_PSE(realData,json_path,raxml_tree,np.array([0.999,0.999,0.999,0.999]),np.array([0.0001,0.0001,0.0001,0.0001]),nu,"BaciSim_PSE_999.json")


