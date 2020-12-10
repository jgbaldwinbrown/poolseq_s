#!/usr/bin/env python3

import math
import vcf
import numpy as np
import sys
import argparse
import copy

def parse_my_args():
    parser = argparse.ArgumentParser("Compute the Cochran-Mantel-Haenszel test on a VCF file")
    parser.add_argument("vcf", nargs="?", help="Input VCF file; default stdin")

    args = parser.parse_args()
    return(args)

def get_arg_vars(args):
    if args.vcf:
        inconn = open(args.vcf, "r")
    else:
        inconn = sys.stdin

    return(inconn)

def write_afs(vcfin, outconn):
    out = []
    for index, record in enumerate(vcfin):
        calls = [call for call in record]
        if index == 0:
            for call in calls:
                out.append(call.sample + "_af")
                out.append(call.sample + "_count")
            outconn.write("\t".join(map(str, out)) + "\n")
        out = []
        for call in calls:
            count = call.data.AD[0] + call.data.AD[1]
            af = float(call.data.AD[0]) / float(count)
            out.append(af)
            out.append(count)
        outconn.write("\t".join(map(str, out)) + "\n")

def main():
    args = parse_my_args()

    inconn = get_arg_vars(args)
    outconn = sys.stdout

    vcfin = vcf.Reader(inconn)
    write_afs(vcfin, outconn)

    inconn.close()

if __name__ == "__main__":
    main()
