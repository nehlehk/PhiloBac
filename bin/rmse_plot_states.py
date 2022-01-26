#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns



parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-t', "--PB_two", type=str, help='rmse_philobacter_two')
parser.add_argument('-g', "--Gubbins", type=str, help='rmse_gubbins')
parser.add_argument('-c', "--CFML", type=str, help='rmse_CFML')
args = parser.parse_args()


rmse_PB_two = args.PB_two
rmse_Gubbins = args.Gubbins
rmse_CFML = args.CFML


# rmse_phylohmm_two = '/home/nehleh/Desktop/temp_reslut/rmse_phylohmm_two.csv'
# rmse_phylohmm_eight = '/home/nehleh/Desktop/temp_reslut/rmse_phylohmm_eight.csv'
# rmse_CFML = '/home/nehleh/Desktop/temp_reslut/rmse_CFML.csv'




f = open(rmse_PB_two, "r")
df = pd.read_csv(f, sep=';', names=['PB_two'], header=None)

f1 = open(rmse_Gubbins, "r")
df1 = pd.read_csv(f1, sep=';', names=['Gubbins'], header=None)

f2 = open(rmse_CFML, "r")
df2 = pd.read_csv(f2, sep=';', names=['CFML'], header=None)




new = pd.concat([df,df1,df2], axis=1)

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(1, 1, 1)
ax1 = sns.boxplot( data=new )
ax1 = sns.stripplot(data=new,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('RMSE values' , fontsize=9)
plt.ylabel('RMSE')


# plt.show()
plt.savefig("RMSE_summary.jpeg")



# Pr = np.unique(df['prob'])
# Nu = np.unique(df['nu'])
# pr_num = len(Pr)
# nu_num = len(Nu)
# i = 1
# fig = plt.figure(figsize=(15,15))
# for n in Nu:
#     for p in Pr:
#         s_df = df.loc[df['nu'] == n]
#         temp = s_df.loc[s_df['prob'] == p]
#         new = pd.concat([temp['Phylohmm_two'],df2],axis=1)
#
#         s_df1 = df1.loc[df1['nu'] == n]
#         temp1 = s_df1.loc[s_df1['prob'] == p]
#         new2 = pd.concat([temp1['Phylohmm_eight'], new], axis=1)
#
#         ax1 = fig.add_subplot(pr_num, nu_num, i)
#         ax1 = sns.boxplot( data=new2  )
#         ax1 = sns.stripplot(data=new2,  jitter=True, dodge=True, marker='o', color=".1")
#         ax1.set_title('RMSE values when nu =' + str(n) + ' and prob =' + str(p) , fontsize=9)
#         plt.ylabel('RMSE')
#         i = i+1
#
#
#
# # plt.show()
# plt.savefig("RMSE_comparison_states.jpeg")