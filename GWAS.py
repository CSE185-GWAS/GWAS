#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import sys
import readvcf
import p_val

def performAnalysis(vcf, phen, outPath='out.csv', willOutputQQ=False, willOutputManhattan=False):
    if not os.path.isfile(vcf) or not os.path.isfile(phen):
        sys.stderr.write('vcf and phenotype files missing\n')
        sys.exit(1)
    
    # generate the genotype df based on vcf file if the dataframe is not processed yet 
    geno_df = readvcf.genoDf(vcf, phen, outPath)
    #geno_df = pd.read_csv(vcf)
    #p_values, beta_values = p_val.calculatePVal(phen, geno_df)
    #print(geno_df)
    #print(p_values)
    #print(beta_values)
    #out_df = pd.DataFrame({'CHR': geno_df['CHROM'], 'SNP': geno_df['ID'], 'BP': geno_df['POS'], 
    #                       'BETA': beta_values, 'P': p_values})
    #out_df.to_csv('out.txt', sep='\t', index=False)
    # print(f'{p_value[:10]}')
    # print(f'{beta_values[:10]}')
    
    df = pd.read_csv(outPath)
    if willOutputQQ:
        p_val.QQPlot(df)
    if willOutputManhattan:
        p_val.manhattanPlot(df)

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
    parser.add_argument('-q', '--qq', help='optional qqman plot saved as "qqplot.png"', required=False, action='store_true')
    parser.add_argument('-m', '--manhattan', help='optional manhattan plot saved as "manhattanplot.png"', action='store_true', required=False)
    
    args = parser.parse_args()
    performAnalysis(args.vcf, args.phen, args.out, args.qq, args.manhattan) 