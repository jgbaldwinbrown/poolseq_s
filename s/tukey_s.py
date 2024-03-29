#!/usr/bin/env python3

import sys
import pandas as pd
from scipy import stats

from statsmodels.stats.multicomp import pairwise_tukeyhsd
# from statsmodels.stats.multicomp import TukeyHSDResults

from typing import TextIO, Tuple, List, Dict, Union, Optional, Any, NewType

def get_data(inconn: TextIO) -> pd.DataFrame:
    data = pd.read_csv(inconn, sep="\t")
    # data.drop(['Unnamed: 0'], axis =1, inplace = True)
    return(data)

def do_tukey(data):
    # perform multiple pairwise comparison (Tukey HSD)
    m_comp = pairwise_tukeyhsd(endog=data['s'], groups=data['full_treatment'], alpha=0.05)
    df = pd.DataFrame(data=m_comp._results_table.data[1:], columns=m_comp._results_table.data[0])
    return(df)

def do_tukeys(data: pd.DataFrame):
    unique_chrpositions: List[int] = list(data["chrpos"].unique())
    out: Any = []
    used_chrpositions: List[int] = []
    chrpos: int = 0
    for chrpos in unique_chrpositions:
        chrpos_data = data[data["chrpos"] == chrpos]
        try:
            tukey: Any = do_tukey(chrpos_data)
            out.append(tukey)
            used_chrpositions.append(chrpos)
        except:
            pass
    return(unique_chrpositions, out, used_chrpositions)

def write_tukey(tukey_results: Any, outconn: TextIO):
    outconn.write(str(tukey_results) + "\n")

def write_tukey_pair(pair: Tuple[str, str], pair_data: pd.DataFrame, out_prefix: str):
    outpath: str = out_prefix + "-".join(pair) + "_tukey.txt"
    pair_data.to_csv(outpath, sep="\t")

def write_tukeys(chrpos: List[int], tukey_results: Any, out_prefix: str):
    unique_group_pairs = zip(list(tukey_results[0]["group1"]), list(tukey_results[0]["group2"]))
    a_chrpos: int = 0
    tukey_result: Any = None
    for a_chrpos, tukey_result in zip(chrpos, tukey_results):
        tukey_result["chrpos"] = a_chrpos
    pair: Tuple[str, str]
    for pair in unique_group_pairs:
        pair_data: pd.DataFrame = pd.concat([x[(x["group1"] == pair[0]) & (x["group2"] == pair[1])] for x in tukey_results])
        # print(pair)
        # print(pair_data)
        # print(chrpos)
        # pair_data["chrpos"] = chrpos
        write_tukey_pair(pair, pair_data, out_prefix)

def main() -> None:
    data: pd.DataFrame = get_data(sys.stdin)
    out_prefix: str = sys.argv[1]
    chrpos, tukey_results, used_chrpos = do_tukeys(data)
    write_tukeys(used_chrpos, tukey_results, out_prefix)

if __name__ == "__main__":
    main()

