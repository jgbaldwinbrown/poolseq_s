#!/usr/bin/env python3

import sys
import manhatify as mh
import pandas as pd
import math
import argparse

def plot_s(genedata, chrom_offsets, outpath, ymin, ymax):
    if (not ymin) or (not ymax):
        ymin = max(min(genedata["s"]), -1)
        ymax = min(max(genedata["s"]), 1)
    mh.plot_manhat(genedata,
        outpath,
        chrom_offsets,
        "s",
        title="Selection coefficient",
        yname = "s",
        dims = (20, 6),
        scale = 1.5,
        log=True,
        named_xticks = True,
        chrom_col = "Scaffold",
        geom = "line",
        ylim = (ymin,ymax)
    )

def plot_p(genedata, chrom_offsets, outpath, ymin, ymax):
    if (not ymin) or (not ymax):
        ymin = max(min(genedata["p"]), 0)
        ymax = max(genedata["p"])
    mh.plot_manhat(genedata,
        outpath,
        chrom_offsets,
        "p",
        title="Probability of selection coefficient occurring randomly",
        yname = "-log10(p)",
        dims = (20, 6),
        scale = 1.5,
        log=True,
        named_xticks = True,
        chrom_col = "Scaffold",
        geom = "line",
        ylim = (ymin,ymax)
    )

def get_args():
    parser = argparse.ArgumentParser("Generate Manhattan plots of S and P based on the output of poolseq_s.")
    parser.add_argument("data", help="Data to plot (required).")
    parser.add_argument("chrlens", help="Bed files specifying chromosome lengths (required).")
    parser.add_argument("-o", "--out_prefix", help="Bed files specifying chromosome length (default=\"out\").", default = "out")
    parser.add_argument("-y", "--ymin", help="Y minimum for plotting selection coefficient (default = minimum of data).", type=float, default=None)
    parser.add_argument("-Y", "--ymax", help="Y maximum for plotting selection coefficient (default = minimum of data).", type=float, default=None)
    parser.add_argument("-p", "--pmin", help="Y minimum for plotting p value (default = minimum of data).", type=float, default=None)
    parser.add_argument("-P", "--pmax", help="Y maximum for plotting p value (default = minimum of data).", type=float, default=None)

    args = parser.parse_args()

    return (args.data, args.chrlens, args.out_prefix, args.ymin, args.ymax, args.pmin, args.pmax)

def main():
    data_inpath, chrlen_bed_inpath, out_prefix, ymin, ymax, pmin, pmax = get_args()
    # data_inpath = sys.argv[1]
    # chrlen_bed_inpath = sys.argv[2]
    # out_prefix = sys.argv[3]
    with open(chrlen_bed_inpath, "r") as inconn:
        chrlens = mh.get_chrom_lens_from_bed(inconn)
    gdata = mh.get_data_from_table(data_inpath, column_names = ["Scaffold", "Position", "s", "p"], header_row = 0)
    gdata["s"] = gdata["s"].apply(pd.to_numeric, errors = "coerce")
    gdata["p"] = gdata["p"].apply(pd.to_numeric, errors = "coerce")
    # try:
    #     gdata["nlogp"] = gdata["p"].apply(lambda x: -math.log10(x))
    # except:
    #     pass
    # print(gdata)
    genedata, chrom_offsets = mh.manhatify(gdata, chrlens, chrom_col = "Scaffold", bp_col = "Position", val_col = "s", offset = 10000)
    # print(genedata)
    # print(chrom_offsets)
    try:
        plot_s(genedata, chrom_offsets, out_prefix + "_s.pdf", ymin, ymax)
    except:
        pass
    try:
        plot_p(genedata, chrom_offsets, out_prefix + "_p.pdf", pmin, pmax)
    except:
        pass

if __name__ == "__main__":
    main()
