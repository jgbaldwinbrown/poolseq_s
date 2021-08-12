#!/usr/bin/env python3

import sys
import manhatify as mh
import pandas as pd
import math
import argparse
from typing import Optional, Union, List, Dict, Tuple

def plot_p(genedata: pd.DataFrame, chrom_offsets: pd.DataFrame, outpath: str, ymin: Optional[float], ymax: Optional[float]) -> None:
    if (not ymin) or (not ymax):
        ymin = max(min(genedata["p"]), 0)
        ymax = max(genedata["p"])
    mh.plot_manhat(genedata,
        outpath,
        chrom_offsets,
        "p",
        title="Tukey test probability of\ndifference in selection coefficient",
        yname = "-log10(p)",
        dims = (20, 6),
        scale = 1.5,
        log=True,
        named_xticks = True,
        chrom_col = "Scaffold",
        geom = "point",
        ylim = (ymin,ymax),
        color_col = "reject"
    )

def get_args() -> Tuple[str, str, str, Optional[float], Optional[float], Optional[float], Optional[float]]:
    parser: argparse.ArgumentParser = argparse.ArgumentParser("Generate Manhattan plots of S and P based on the output of poolseq_s.")
    parser.add_argument("data", help="Data to plot (required).")
    parser.add_argument("chrlens", help="Bed files specifying chromosome lengths (required).")
    parser.add_argument("-o", "--out_prefix", help="Output file path (default=\"out\").", default = "out")
    parser.add_argument("-y", "--ymin", help="Y minimum for plotting selection coefficient (default = minimum of data).", type=float, default=None)
    parser.add_argument("-Y", "--ymax", help="Y maximum for plotting selection coefficient (default = minimum of data).", type=float, default=None)
    parser.add_argument("-p", "--pmin", help="Y minimum for plotting p value (default = minimum of data).", type=float, default=None)
    parser.add_argument("-P", "--pmax", help="Y maximum for plotting p value (default = minimum of data).", type=float, default=None)

    args: argparse.Namespace = parser.parse_args()

    return (args.data, args.chrlens, args.out_prefix, args.ymin, args.ymax, args.pmin, args.pmax)

def get_chr(s: str) -> str:
    return s.split(":")[0]

def get_pos(s: str) -> int:
    return int(s.split(":")[1])

def nlog10(f: float) -> float:
	return -(math.log10(f))

def main() -> None:
    data_inpath: str
    chrlen_bed_inpath: str
    out_prefix: str
    ymin: Optional[float]
    ymax: Optional[float]
    pmin: Optional[float]
    pmax: Optional[float]
    data_inpath, chrlen_bed_inpath, out_prefix, ymin, ymax, pmin, pmax = get_args()
    # data_inpath = sys.argv[1]
    # chrlen_bed_inpath = sys.argv[2]
    # out_prefix = sys.argv[3]
    with open(chrlen_bed_inpath, "r") as inconn:
        chrlens = mh.get_chrom_lens_from_bed(inconn)
    gdata: pd.DataFrame
    gdata = mh.get_data_from_table(data_inpath, header_row = 0)
    gdata["Scaffold"] = gdata["chrpos"].apply(get_chr)
    gdata["Position"] = gdata["chrpos"].apply(get_pos)
    gdata["p"] = gdata["p-adj"].apply(nlog10)
    # try:
    #     gdata["nlogp"] = gdata["p"].apply(lambda x: -math.log10(x))
    # except:
    #     pass
    # print(gdata)
    genedata: pd.DataFrame
    chrom_offsets: pd.DataFrame
    # print(gdata.dtypes)
    genedata, chrom_offsets, extra_data = mh.manhatify(gdata, chrlens, chrom_col = "Scaffold", bp_col = "Position", val_col = "p", offset = 10000)
    # print(genedata)
    # print(chrom_offsets)
    # try:
    plot_p(genedata, chrom_offsets, out_prefix + "_p.pdf", pmin, pmax)
    # except:
    #     pass

if __name__ == "__main__":
    main()
