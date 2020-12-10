#!/usr/bin/env python3

import vcf
import numpy as np
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append("/home/jgbaldwinbrown/Documents/git_repositories/geneview")
#/home/jgbaldwinbrown/Documents/git_repositories/geneview/geneview
import geneview as gv

def parse_my_args():
    parser = argparse.ArgumentParser("Plot the chosen field from a vcf file.")
    parser.add_argument("vcf", nargs="?", help="Input VCF file; default stdin")
    parser.add_argument("-f", "--field",  help="VCF INFO field to plot; required.", required=True)
    parser.add_argument("-o", "--outfile", help="file to plot in; default = show")

    args = parser.parse_args()
    return(args)

def get_arg_vars(args):
    if args.vcf:
        inconn = open(args.vcf, "r")
    else:
        inconn = sys.stdin
    field = args.field
    if args.outfile:
        outfile = args.outfile
        pdf = True
    else:
        outfile = ""
        pdf = False
    return((inconn, field, outfile, pdf))

def vcf2pd(vcfin, field):
    data = vcf2list(vcfin, field)
    out = pd.DataFrame.from_records(data)
    out.columns = ["CHROM", "POS", field]
    return(out)

def vcf2list(vcfin, field):
    out = []
    for record in vcfin:
        record_data = []
        record_data.append(record.CHROM)
        record_data.append(record.POS)
        record_data.append(float(record.INFO[field][0]))
        out.append(record_data)
    return(out)

def plot_data(data, outfile, pdf, field):

    #xtick = ['chr'+c for c in map(str, range(1, 15) + ['16', '18', '20', '22'])]
    gv.gwas.manhattanplot(data,  
                             xlabel="Chromosome", 
                             ylabel=field, 
                             xticklabel_kws={'rotation': 'vertical'}
                             )
    #g = sns.FacetGrid(data)
    #g.map(sns.relplot, x='POS', y=field, data=data, col='CHROM')
    #g.set(xlim = 
    #g.set(ylim=(-1, 11), yticks=[0, 5, 10]);
    #myplot=sns.relplot(x='POS', y=field, data=data, col = 'CHROM')
    if pdf:
        plt.savefig(outfile)
    else:
        plt.show()

def main():
    args = parse_my_args()

    inconn, field , outfile, pdf = get_arg_vars(args)

    vcfin = vcf.Reader(inconn)
    data = vcf2pd(vcfin, field)
    inconn.close()

    plot_data(data, outfile, pdf, field)
    
    print(data)

if __name__ == "__main__":
    main()

#>>> for i in b:
#...     for j in i:
#...         try:print(j.data.AD)
