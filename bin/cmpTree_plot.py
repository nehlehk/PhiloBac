#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np


parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-c', "--cmpTreesResult", type=str, help='cmpTreesResult')
parser.add_argument('-p', "--PB_dist",type=str, help='PB_dist')
parser.add_argument('-g', "--Gubb_dist",type=str, help='Gubb_dist')
parser.add_argument('-m', "--CFML_dist",type=str, help='CFML_dist')
args = parser.parse_args()

cmpTrees = args.cmpTreesResult


# cmpTrees = '/home/nehleh/PhiloBacteria/hpc/nu_09/Summary_Results/all_cmpTrees.result'


f = open(cmpTrees, "r")
df = pd.read_csv(f,sep='\t', names=['No', 'clonalTree','otherTrees', 'Quartet' , 'PathDiffernce' , 'RF' , 'MatchingSplit' , 'UMAST' , 'RFWeighted' , 'GeoUnrooted' ] ,header=None)



new_df= df.loc[df['otherTrees'].between('0' , '4' , inclusive=False )].apply(pd.to_numeric)
new_df['otherTrees'].replace(to_replace=[1,2,3],value=['CFML','Gubbins','PhiloBacter'] ,inplace=True)
# print(new_df)
new_df = new_df.loc[(new_df['otherTrees'] != 'Gubbins') ]





fig = plt.figure(figsize=(10,5))


ax6 = fig.add_subplot(1, 2, 1)
ax6 = sns.boxplot(x = 'otherTrees', y="RFWeighted" ,data=new_df )
ax6 = sns.stripplot(x = 'otherTrees', y="RFWeighted" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax6.set_title('RFWeighted' , fontsize=9)

ax7 = fig.add_subplot(1, 2, 2)
ax7 = sns.boxplot(x = 'otherTrees', y="GeoUnrooted" ,data=new_df )
ax7 = sns.stripplot(x = 'otherTrees', y="GeoUnrooted" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax7.set_title('GeoUnrooted' , fontsize=9)


plt.ylabel('Value')
# plt.xlabel('')
# plt.xticks([])
plt.savefig("TreeCmp_summary.jpeg")
# plt.show()


dist_PBtree = args.PB_dist
dist_Gubb = args.Gubb_dist
dist_CFML = args.CFML_dist

f = open(dist_PBtree, "r")
df = pd.read_csv(f, sep=';', names=['PB'], header=None)

f1 = open(dist_Gubb, "r")
df1 = pd.read_csv(f1, sep=';', names=['Gubbins'], header=None)

f2 = open(dist_CFML, "r")
df2 = pd.read_csv(f2, sep=';', names=['CFML'], header=None)


new = pd.concat([df2,df1,df], axis=1)

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(1, 1, 1)
ax1 = sns.boxplot( data=new )
ax1 = sns.stripplot(data=new,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('tree comparison_ euclidean_distance' , fontsize=9)
plt.ylabel('euclidean_distance')


# plt.show()
plt.savefig("Dist_summary.jpeg")