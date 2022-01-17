#!/usr/bin/env python

import pandas as pd
import argparse
from dendropy import Tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltree", type=str, help='clonaltree')
    parser.add_argument('-r', "--raxml", type=str, help='raxml')
    # parser.add_argument('-rg', "--raxGap", type=str, help='raxGap')
    # parser.add_argument('-rd', "--raxDel", type=str, help='raxDel')
    parser.add_argument('-cfml', "--cfml", type=str, help='cfml')
    parser.add_argument('-g', "--gubbins", type=str, help='gubbins')
    # parser.add_argument('-b', "--beast", type=str, help='gubbins')
    # parser.add_argument('-bp', "--beastpartial", type=str, help='beastpartial')
    # parser.add_argument('-n', "--recomnum", type=str, help='recomnum')
    parser.add_argument('-o', "--output", type=str, help='output')
    parser.add_argument('-o2', "--output2", type=str, help='output_root')

    args = parser.parse_args()

    cl_tree = args.clonaltree
    r_tree = args.raxml
    # rgap_tree = args.raxGap
    # rdel_tree = args.raxDel
    cf_tree = args.cfml
    g_tree = args.gubbins
    # b_tree = args.beast
    # bp_tree = args.beastpartial
    # recom_count = args.recomnum
    outputfile = args.output
    outputfile2 = args.output2



    clonal_tree = Tree.get_from_path(cl_tree, 'newick')
    raxml_tree = Tree.get_from_path(r_tree, 'newick')
    # raxgap_tree = Tree.get_from_path(rgap_tree, 'newick')
    # raxdel_tree = Tree.get_from_path(rdel_tree, 'newick')
    cfml_tree = Tree.get_from_path(cf_tree, 'newick')
    gubbins_tree = Tree.get_from_path(g_tree, 'newick')
    # beast_tree = Tree.get_from_path(b_tree, 'newick')
    # beast_par_tree = Tree.get_from_path(bp_tree, 'newick')

    cl_nodes = []
    cl_edge_length = []
    r_edge_length = []
    r_edge_length_diff = []
    r_edge_length_diff_r = []
    rg_edge_length = []
    rg_edge_length_diff = []
    rg_edge_length_diff_r = []
    rd_edge_length = []
    rd_edge_length_diff = []
    rd_edge_length_diff_r = []
    cfml_edge_length = []
    cfml_edge_length_diff = []
    cfml_edge_length_diff_r = []
    g_edge_length = []
    g_edge_length_diff = []
    g_edge_length_diff_r = []
    b_edge_length = []
    b_edge_length_diff = []
    b_edge_length_diff_r = []
    bp_edge_length = []
    bp_edge_length_diff = []
    bp_edge_length_diff_r = []
    count = []

    # df_count = pd.read_csv(recom_count,sep='\t', engine='python')
    # print(df_count)



    for i in range(len(clonal_tree)):
        cl_node = clonal_tree.find_node_with_taxon_label(str(i))

        cl_nodes.append(cl_node.taxon.label)
        cl_edge_length.append(cl_node.edge_length)

        r_node = raxml_tree.find_node_with_taxon_label(str(i))
        r_edge_length.append(r_node.edge_length)
        r_edge_length_diff.append(abs(r_node.edge_length - cl_node.edge_length))
        r_edge_length_diff_r.append(abs(r_node.distance_from_root() - cl_node.distance_from_root()))

        # rg_node = raxgap_tree.find_node_with_taxon_label(str(i))
        # rg_edge_length.append(rg_node.edge_length)
        # rg_edge_length_diff.append(abs(rg_node.edge_length - cl_node.edge_length))
        # rg_edge_length_diff_r.append(abs(rg_node.distance_from_root() - cl_node.distance_from_root()))
        #
        # rd_node = raxdel_tree.find_node_with_taxon_label(str(i))
        # rd_edge_length.append(rd_node.edge_length)
        # rd_edge_length_diff.append(abs(rd_node.edge_length - cl_node.edge_length))
        # rd_edge_length_diff_r.append(abs(rd_node.distance_from_root() - cl_node.distance_from_root()))

        cfml_node = cfml_tree.find_node_with_taxon_label(str(i))
        cfml_edge_length.append(cfml_node.edge_length)
        cfml_edge_length_diff.append(abs(cfml_node.edge_length - cl_node.edge_length))
        cfml_edge_length_diff_r.append(abs(cfml_node.distance_from_root() - cl_node.distance_from_root()))

        g_node = gubbins_tree.find_node_with_taxon_label(str(i))
        g_edge_length.append(g_node.edge_length)
        g_edge_length_diff.append(abs(g_node.edge_length - cl_node.edge_length))
        g_edge_length_diff_r.append(abs(g_node.distance_from_root() - cl_node.distance_from_root()))

        # b_node = beast_tree.find_node_with_taxon_label(str(i))
        # b_edge_length.append(b_node.edge_length)
        # b_edge_length_diff.append(abs(b_node.edge_length - cl_node.edge_length))
        # b_edge_length_diff_r.append(abs(b_node.distance_from_root() - cl_node.distance_from_root()))
        #
        # bp_node = beast_par_tree.find_node_with_taxon_label(str(i))
        # bp_edge_length.append(bp_node.edge_length)
        # bp_edge_length_diff.append(abs(bp_node.edge_length - cl_node.edge_length))
        # bp_edge_length_diff_r.append(abs(bp_node.distance_from_root() - cl_node.distance_from_root()))

        # r_num  = (df_count.loc[df_count['nodes'] == int(cl_node.taxon.label)]['count'].values)
        # if len(r_num) == 0:
        #     r_num = 0
        # else:
        #     r_num = int(r_num)
        # count.append(r_num)


    all_data = {'node': cl_nodes,'gap_r':rg_edge_length_diff, 'del_r':rd_edge_length_diff,'cfml':cfml_edge_length_diff ,'gubbins': g_edge_length_diff, 'count': count}
    df = pd.DataFrame(all_data)
    # print(df)

    all_data_root = {'node': cl_nodes,'gap_r':rg_edge_length_diff_r, 'del_r':rd_edge_length_diff_r,'cfml':cfml_edge_length_diff_r ,'gubbins': bp_edge_length_diff_r,  'count': count}
    df_root = pd.DataFrame(all_data_root)


    df.to_csv(outputfile, sep='\t', index = False , header= False)
    df_root.to_csv(outputfile2, sep='\t', index=False, header=False)


