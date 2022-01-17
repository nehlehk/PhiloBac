#!/usr/bin/env python

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":

    # parser=argparse.ArgumentParser(description='''Enter you input ''',  epilog="""All's well that ends well.""")
    # parser.add_argument('-f', "--file", type=str  , help='input file')
    #
    # args = parser.parse_args()
    # input = args.file

    input = '/home/nehleh/PhiloBacteria/hpc_beagle/100_100/All_result.csv'

    df = pd.read_csv(input, sep=',', engine='python')
    print(df)
    df['log_numpy'] = np.log2(df['time_lp_numpy'])
    df['log_cpu'] = np.log2(df['time_lp_beagle'])
    df['log_gpu'] = np.log2(df['time_lp_beagle_GPU'])
    print(df)
    # df.plot.bar(x="seq_len", y=["time_lp_numpy","time_lp_beagle","time_lp_beagle_GPU"], figsize=(8, 8), fontsize=10, rot=20, title="time to calculate Tree likelihood with Numpy and BeaglePy, tips number =100")

    # df.plot(x="seq_len", y=["time_lp_numpy","time_lp_beagle","time_lp_beagle_GPU"], kind='line' , figsize=(8, 8), fontsize=10 )
    df.plot(x="seq_len", y=["log_numpy", "log_cpu", "log_gpu"], kind='line', figsize=(8, 8), fontsize=10)
    plt.xlabel("Sequence length")
    plt.ylabel("log (time - second) ")
    plt.legend(["Numpy", "BeaglePy_CPU" , "BeaglePy_GPU"])
    plt.savefig("All_result.jpeg")
    plt.show()


