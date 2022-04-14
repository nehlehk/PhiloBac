#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np


parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-c', "--cmpTreesResult", type=str, help='cmpTreesResult')
args = parser.parse_args()

cmpTrees = args.cmpTreesResult

# cmpTrees = '/home/nehleh/PhiloBacteria/hpc/nu_09/Summary_Results/all_cmpTrees.result'


f = open(cmpTrees, "r")
df = pd.read_csv(f,sep='\t', names=['No', 'clonalTree','otherTrees', 'Quartet' , 'PathDiffernce' , 'RF' , 'MatchingSplit' , 'UMAST' , 'RFWeighted' , 'GeoUnrooted' ] ,header=None)

# print(df)

new_df= df.loc[df['otherTrees'].between('0' , '4' , inclusive=False )].apply(pd.to_numeric)
new_df['otherTrees'].replace(to_replace=[1,2,3],value=['CFML','Gubbins','PB2'] ,inplace=True)
# print(new_df)
# g_df = new_df.loc[(new_df['otherTrees'] != 'Gubbins') ]





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



