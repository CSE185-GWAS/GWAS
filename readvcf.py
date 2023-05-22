# this file can read in vcf or vcf.gzip file, arrange it into a csv file 
# reference: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

import gzip 
import io
import os
import pandas as pd
# used to store content in vcf file
def read_vcf(path, file_format):
    lines = []
    if (file_format == 'gzip'):  
        with gzip.open(path, 'rt') as f:
            for l in f:
                if (not l.startswith('##')):
                    lines.append(l)
        # read the info from vcf file into a pandas dataframe
        # by using read_csv
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
    else:
        with open(path, 'r') as f:
            for l in f:
                if (not l.startswith('##')):
                    lines.append(l)
        # read the info from vcf file into a pandas dataframe
        # by using 
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

# this is used to test vcf file by lab 3 file
def testReadVCF():
    path = "lab3_gwas.vcf.gz"
    vcf_df = read_vcf(path, 'gzip')
    vcf_df.to_csv('out_vcf.csv')

# alt can be multiple separated by ','
def assignGenoType(ref, alt, query):
    query = query.split('|')
    alt = alt.split(',')
    genotype = ''
    dictionary = {
        "0" : ref
    }
    index = 1
    for alternative in alt:
        dictionary[str(index)] = alternative
        index+=1
    
    for num in query:
        if num in dictionary.keys():
            genotype += dictionary[num]
    
    return genotype


def testGenotype():
    vcf_df = pd.read_csv('out_vcf.csv')
    ref_allele = vcf_df['REF']
    alt_allele = vcf_df['ALT']
    
    size = vcf_df.shape[0]
    # sample names
    column_name = vcf_df.columns
    # base on vcf format, extract 12th column till the end for sample names
    sample_name = column_name[10:]

    for sample in sample_name: 
        genotype_list = []
        sample_allele = vcf_df[sample]
        for i in range(size):
            genotype_list.append(assignGenoType(ref_allele[i], alt_allele[i], sample_allele[i]))
        vcf_df[sample] = genotype_list
    print(vcf_df)
    vcf_df.to_csv('final_vcf.csv')

testGenotype()