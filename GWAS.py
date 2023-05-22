#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import sys

def performAnalysis(vcf, phen, graphType=None, graphPath=None):
    if not os.path.isfile(vcf) or not os.path.isfile(phen):
        sys.stderr.write('vcf and phenotype files missing')
        sys.exit(1)

    # TODO: create plots for specified types
    if graphType == 'qq':
        print()
    elif graphType == 'manhattan':
        print()

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
    parser.add_argument('-o', '--out', help='optional file to output GWAS graph results; default file name is graph.png if not specified',
                        metavar='FILE', type=str, required=False)
    parser.add_argument('-g', '--graph_type', help='optional type of graph to save; graph will not be saved otherwise',
                        metavar='FILE', type=str, required=False)
    
    args = parser.parse_args()
    performAnalysis(args.vcf, args.phen, args.out, args.graph_type) 