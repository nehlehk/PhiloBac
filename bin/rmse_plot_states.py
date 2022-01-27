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




new = pd.concat([df2,df1,df], axis=1)

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(1, 1, 1)
ax1 = sns.boxplot( data=new )
ax1 = sns.stripplot(data=new,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('RMSE values' , fontsize=9)
plt.ylabel('RMSE')


# plt.show()
plt.savefig("RMSE_summary.jpeg")