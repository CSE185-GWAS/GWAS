# this file can read in vcf or vcf.gzip file, arrange it into a csv file 
# reference: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

import gzip 
import io
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import p_val as pval
# this will convert list of queries to a genotype list 
# it will only look at one alt allele
def convert_to_genoTypes(ref, alt, queries):
    genotypes = []
    dictionary = {
        "0" : ref,
        "1" : alt
    }

    for query in queries: 
        genotype = ''
        num1  = query[0]
        num2 = query[2]
        if num1 in dictionary.keys():
            genotype += dictionary[num1]
        if num2 in dictionary.keys():
            genotype += dictionary[num2]

        genotypes.append(genotype)
    
    return genotypes

# used to store content & analyze genotype data in vcf file
# this will return a df with genotype info grabbed from vcf file
def read_vcf(path, file_format):
    lines = []
    if (file_format == 'gzip'):  
        with gzip.open(path, 'rt') as f:
            for l in f:
                if (l.startswith('#CHROM')):
                    lines.append(l)
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split('\t')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]
                    # as it is rare to have multiple allele, 
                    # only keep first alternative allele for analysis
                    alt_allele = alt_allele.split(',')[0]
                    genotypes = convert_to_genoTypes(ref_allele, alt_allele, queires)
                    info = line[:9] + genotypes
                    mod_l = '\t'.join(info)
                    last_char = mod_l[-1]
                    # replace \t at end with \n
                    mod_l = mod_l[:-1] + last_char + '\n'
                    lines.append(mod_l)
    else: 
        with open(path, 'r') as f:
            for l in f:
                if (l.startswith('#CHROM')):
                    lines.append(l)
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split('\t')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]
                    # as it is rare to have multiple allele, 
                    # only keep first alternative allele for analysis
                    alt_allele = alt_allele.split(',')[0]
                    genotypes = convert_to_genoTypes(ref_allele, alt_allele, queires)
                    info = line[:9] + genotypes
                    mod_l = '\t'.join(info)
                    last_char = mod_l[-1]
                    # replace \t at end with \n
                    mod_l = mod_l[:-1] + last_char + '\n'
                    lines.append(mod_l)
    # read the info from vcf file into a pandas dataframe
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

# this function will read in path of vcf file and convert it 
# to a df with genotype info, stored as geno.csv
def genoDf(path):
    if ('gz' in path):
        vcf_df = read_vcf(path, 'gzip')
    else:
        vcf_df = read_vcf(path, 'vcf')
    vcf_df = vcf_df.drop(columns=['QUAL','FILTER','INFO','FORMAT'])
    vcf_df.to_csv('geno_vcf.csv', index=False)
    
    print('Geno Dafarame is created')
    return vcf_df
    
# uncomment below to test this out using lab3 data
#genoDf('lab3_gwas.vcf.gz')









#********************** Unit Testing for function *****************************#
# this is used to read vcf file by lab 3 file
def testReadVCF():
    path = "lab3_gwas.vcf.gz"
    vcf_df = read_vcf(path, 'gzip')
    vcf_df.to_csv('out_vcf.csv')

# this is used to test if geno df function 
def testGenoDF():
    genoDf('lab3_gwas.vcf.gz')

# test code for assigning genotype by reading from csv file extracted from vcf file
def analyze_labb3_genotype():
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
            genotype_list.append(assign_multi_genoType(ref_allele[i], alt_allele[i], sample_allele[i]))
        vcf_df[sample] = genotype_list
    print(vcf_df)
    vcf_df.to_csv('final_vcf.csv')


#**Function that will output all genotype, even by >1 alternative allele*****#
# alt can be multiple separated by ','
# this will assign genotype to sample based on query 
def assign_multi_genoType(ref, alt, query):
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

