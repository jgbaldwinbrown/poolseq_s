#!/usr/bin/env python3

import sys
import pandas as pd
from scipy import stats

from statsmodels.stats.multicomp import pairwise_tukeyhsd

def get_data(inconn):
    data = pd.read_csv(inconn, sep="\t")
    # data.drop(['Unnamed: 0'], axis =1, inplace = True)
    return(data)

def do_tukey(data):
    # perform multiple pairwise comparison (Tukey HSD)
    m_comp = pairwise_tukeyhsd(endog=data['p.value'], groups=data['full_treatment'], alpha=0.05)
    df = pd.DataFrame(data=m_comp._results_table.data[1:], columns=m_comp._results_table.data[0])
    return(df)

def do_tukeys(data):
    unique_chrpositions = list(data["chrpos"].unique())
    out = []
    for chrpos in unique_chrpositions:
        chrpos_data = data[data["chrpos"] == chrpos]
        out.append(do_tukey(chrpos_data))
    return(unique_chrpositions, out)

def write_tukey(tukey_results, outconn):
    outconn.write(str(tukey_results) + "\n")

def write_tukey_pair(pair, pair_data, out_prefix):
    outpath = out_prefix + "-".join(pair) + "_tukey.txt"
    pair_data.to_csv(outpath, sep="\t")

def write_tukeys(chrpos, tukey_results, out_prefix):
    unique_group_pairs = zip(list(tukey_results[0]["group1"]), list(tukey_results[0]["group2"]))
    for pair in unique_group_pairs:
        print(pair)
        pair_data = pd.concat([x[(x["group1"] == pair[0]) & (x["group2"] == pair[1])] for x in tukey_results])
        pair_data["chrpos"] = chrpos
        write_tukey_pair(pair, pair_data, out_prefix)

def main():
    data = get_data(sys.stdin)
    out_prefix = sys.argv[1]
    # print(data)
    chrpos, tukey_results = do_tukeys(data)
    write_tukeys(chrpos, tukey_results, out_prefix)

if __name__ == "__main__":
    main()

