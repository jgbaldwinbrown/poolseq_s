#!/usr/bin/env python3

import math
import vcf
import numpy as np
import statsmodels.api as sm
import sys
import argparse
import copy

def parse_my_args():
    parser = argparse.ArgumentParser("Compute the Cochran-Mantel-Haenszel test on a VCF file")
    parser.add_argument("vcf", nargs="?", help="Input VCF file; default stdin")
    parser.add_argument("-c", "--control",  help="comma separated, 0-indexed VCF columns to use as controls; required.", required=True)
    parser.add_argument("-t", "--test",  help="comma separated, 0-indexed VCF columns to use as test data; required.", required=True)

    args = parser.parse_args()
    return(args)

def get_arg_vars(args):
    if args.vcf:
        inconn = open(args.vcf, "r")
    else:
        inconn = sys.stdin

    control = [int(x) for x in args.control.split(",")]
    test = [int(x) for x in args.test.split(",")]
    if not len(control) == len(test):
        exit("unequal number of control and test pops!")
    return((inconn, control, test))

def cmh_vcf(vcfin, control, test, outwriter):
    tester = np.ndarray(shape=(2,2,len(control)))
    for record in vcfin:
        calls = [call for call in record]
        for i in range(len(control)):
            tester[i,0,0] = calls[control[i]].data.AD[0]
            tester[i,0,1] = calls[control[i]].data.AD[1]
            tester[i,1,0] = calls[test[i]].data.AD[0]
            tester[i,1,1] = calls[test[i]].data.AD[1]
        cmh = sm.stats.StratifiedTable(tester)
        writeout(cmh, record, outwriter)
        #print(cmh.summary())

        #control_ad0, control_ad1 = [calls.AD[:2] for i in control]
        #test_ad0, test_ad1 = [calls.AD[:2] for i in test]

def writeout(cmh, record, writer):
    #newrecord = copy.deepcopy(record)
    record.INFO["CMH"] = str(cmh.test_null_odds().pvalue)
    record.INFO["NLOG10CMH"] = str(-math.log10(cmh.test_null_odds().pvalue))
    writer.write_record(record)

def main():
    args = parse_my_args()

    inconn, control, test = get_arg_vars(args)

    vcfin = vcf.Reader(inconn)
    outwriter = vcf.Writer(sys.stdout, vcfin)
    cmh_vcf(vcfin, control, test, outwriter)

    inconn.close()

if __name__ == "__main__":
    main()

#>>> for i in b:
#...     for j in i:
#...         try:print(j.data.AD)
