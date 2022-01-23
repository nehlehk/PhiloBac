#!/usr/bin/env python

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='''Enter you input ''',  epilog="""All's well that ends well.""")
    parser.add_argument('-f', "--file", type=str  , help='input file')

    args = parser.parse_args()
    input = args.file

    input = '/home/nehleh/PhiloBacteria/hpc_beagle/100_100/All_result.csv'

    df = pd.read_csv(input, sep=',', engine='python')


    df.plot(x="seq_len", y=["time_lp_numpy","time_lp_beagle","time_lp_beagle_GPU"], kind='line' , figsize=(8, 8), fontsize=10 )
    plt.xlabel("Sequence length")
    plt.ylabel("time (second) -- log scale ")
    plt.legend(["Numpy", "BeaglePy_CPU" , "BeaglePy_GPU"])
    plt.yscale('log')
    plt.savefig("All_result.jpeg")
    plt.show()


