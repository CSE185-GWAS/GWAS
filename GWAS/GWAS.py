#!/usr/bin/env python
import argparse
from GWAS import __version__
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import sys
from . import p_val 
from . import readvcf
# parse argument from command line, and call function to execute the analysis 
def main():
    parser = argparse.ArgumentParser(
        prog='GWAS-py',
        description='Python script to perform GWAS analysis that by default outputs a csv file for further analysis with optional arguments to save graphs'
    )

    parser.add_argument('--vcf', help='input vcf file that contains information about genotype and phenotype for analysis',
                        metavar='FILE', type=str, required=True)
    parser.add_argument('-p', '--phen', help='corresponding phenotype file', 
                        metavar='FILE', type=str, required=True)
    parser.add_argument('--maf', help='optional threshold for filtering out SNPs with MAF below this threshold; default is 0.05',
                        type=float, required=False)
    parser.add_argument('-o', '--out', help='optional file to output GWAS results',
                        metavar='FILE', type=str, required=False)
    parser.add_argument('-q', '--qq', help='optional qqman plot saved as "qqplot.png"', required=False, action='store_true')
    parser.add_argument('-m', '--manhattan', help='optional manhattan plot saved as "manhattanplot.png"', action='store_true', required=False)
    
    args = parser.parse_args()
    print(args)
    if not args.maf:
        args.maf = 0.05
    performAnalysis(args.vcf, args.phen, args.maf, args.out, args.qq, args.manhattan) 

# perform GWAS based on input and print error message if input file is invalid 
def performAnalysis(vcf, phen, maf=0.05, outPath='out.csv', willOutputQQ=False, willOutputManhattan=False):
    if not os.path.isfile(vcf) or not os.path.isfile(phen):
        sys.stderr.write('vcf and phenotype files missing\n')
        sys.exit(1)
    
    # generate the output df based on vcf file
    geno_df = readvcf.genoDf(vcf, phen, outPath, maf)

    # convert result file into dataframe for plottings 
    df = pd.read_csv(outPath)
    if willOutputQQ:
        p_val.QQPlot(df)
    if willOutputManhattan:
        p_val.manhattanPlot(df)

# return messages based on command line input 
if __name__=='__main__':
    main()
