#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import sys
import readvcf
import p_val

def performAnalysis(vcf, phen, outPath='out.csv', graphType='Both'):
    if not os.path.isfile(vcf) or not os.path.isfile(phen):
        sys.stderr.write('vcf and phenotype files missing\n')
        sys.exit(1)
    
    # generate the output df based on vcf file
    readvcf.genoDf(vcf, phen, outPath)
    geno_df = pd.read_csv(outPath)
    
    # TODO: create plots for specified types
    if graphType == 'qq':
        p_val.QQPlot(geno_df['P'].tolist())

    elif graphType == 'manhattan':
        p_val.manhattanPlot(geno_df)

    # if no input, print both graph by default
    else:
        p_val.QQPlot(geno_df['P'].tolist())
        p_val.manhattanPlot(geno_df)


# return messages based on command line input 
if __name__=='__main__':
    parser = argparse.ArgumentParser(
        prog='GWAS-py',
        description='Python script to perform GWAS analysis that by default outputs a csv file for further analysis with optional arguments to save graphs'
    )

    parser.add_argument('--vcf', help='input vcf file that contains information about genotype and phenotype for analysis',
                        metavar='FILE', type=str, required=True)
    parser.add_argument('-p', '--phen', help='corresponding phenotype file', 
                        metavar='FILE', type=str, required=True)
    parser.add_argument('--maf', help='optional threshold for filtering out SNPs with MAF below this threshold',
                        type=float, required=False)
    parser.add_argument('-o', '--out', help='optional file to output GWAS results',
                        metavar='FILE', type=str, required=False)
    parser.add_argument('-g', '--graph_type', help='optional type of graph to save; graph will not be saved otherwise',
                        metavar='FILE', type=str, required=False)
    
    args = parser.parse_args()
    performAnalysis(args.vcf, args.phen, args.out, args.graph_type) 