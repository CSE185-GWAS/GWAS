# this file can read in vcf or vcf.gzip file, arrange it into a csv file 
# reference: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

import gzip 
import io
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import p_val as pval

# weird alt allele example: <CN0> CAGC GT very long sequences
#**Function that will output all genotype, even by >1 alternative allele*****#
# alt can be multiple separated by ','
# this will assign genotype to sample based on query 
def assign_multi_genoType(ref, alt, queries):
    genotypes = []
    # create dictionary based on input alt
    alt = alt.split(',')
    genotype = ''
    dictionary = {
        "0" : ref
    }
    
    alt_allele_counts = []
    # if there is multi alt allele, map them to dictionary one by one 
    index = 1
    for alternative in alt:
        dictionary[str(index)] = alternative
        index+=1
        alt_allele_counts.append(0)
    # fill in genotypes list based on query 
    for query in queries: 
        genotype = ''
        query = query.split('|')
        for num in query:
            if num in dictionary.keys():
                # if it is a alternaitve allele
                if (int(num) > 0):
                    # update its count in gt_allele_counts
                    alt_allele_counts[int(num)-1] += 1
                genotype += dictionary[num]
        genotypes.append(genotype)
    # if more than one alt allele exists
    sig_alt_allele = ''
    if (len(alt_allele_counts) > 1):
        # find the significant allele
        sig_alt_allele = alt_allele_counts.index(max(alt_allele_counts)) + 1
        sig_alt_allele = dictionary[str(sig_alt_allele)]
        #print('most significant alternative allele is: {}'.format(sig_alt_allele))
        # reassign dictionary to only contains most significant alt allele 
        dictionary = [ref, sig_alt_allele]
        for i in range(len(genotypes)):
            geno = genotypes[i]
            for allele in geno:
                # if it has minor alleles 
                if (allele not in dictionary):
                    # mark it as empty because we will ignore it 
                    genotypes[i] = 'N' 

    return genotypes, sig_alt_allele

# used to store content & analyze genotype data in vcf file
# this will return a df with genotype info grabbed from vcf file
def read_vcf(path, file_format):
    lines = []
    if (file_format == 'gzip'):  
        with gzip.open(path, 'rt') as f:
            for l in f:
                if (l.startswith('#CHROM')):
                    # remove useless info from column names (QUAL,FILTER,INFO,FORMAT)
                    l = l.split('\t')
                    l_part1 = l[:5]
                    l_part2 = l[9:]
                    l_list = l_part1 + l_part2
                    l = '\t'.join(l_list)
                    lines.append(l)
              
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split('\t')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]
                    # ignore invalid SNP with >1 length ref_allele
                    if (len(ref_allele) > 1):
                        print('Invalid SNP, ignore this line')
                        continue 

                    # assign genotypes based on input
                    genotypes, sig_alt_allele = assign_multi_genoType(ref_allele, alt_allele, queires)
                    # if there are multi alt allele exist
                    # reformat alt alelle to most significant alt allele
                    if sig_alt_allele != '':
                        line[4] = sig_alt_allele
                    # ignore invalid SNP with >1 length alternative allele 
                    if (len(line[4]) > 1):
                        print('Invalid SNP, ignore this line')
                        continue 

                    info = line[:9] + genotypes
                    # remove useless info (QUAL, INFO, FILTER, FORMAT)from line 
                    info_pt1 = info[:5]
                    info_pt2 = info[9:]
                    info = info_pt1 + info_pt2
                    mod_l = '\t'.join(info)
                    last_char = mod_l[-1]
                    # replace \t at end with \n
                    mod_l = mod_l[:-1] + last_char + '\n'
                    lines.append(mod_l)
    else: 
        with open(path, 'r') as f:
            for l in f:
                if (l.startswith('#CHROM')):
                    # remove useless info from column names (QUAL,FILTER,INFO,FORMAT)
                    l = l.split(' ')
                    l_part1 = l[:5]
                    l_part2 = l[9:]
                    l_list = l_part1 + l_part2
                    l = '\t'.join(l_list)
                    lines.append(l)
              
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split(' ')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]
                    # ignore invalid SNP with >1 length ref_allele
                    if (len(ref_allele) > 1):
                        print('Invalid SNP, ignore this line')
                        continue 

                    # assign genotypes based on input
                    genotypes, sig_alt_allele = assign_multi_genoType(ref_allele, alt_allele, queires)
                    # if there are multi alt allele exist
                    # reformat alt alelle to most significant alt allele
                    if sig_alt_allele != '':
                        line[4] = sig_alt_allele

                    if (len(line[4]) > 1):
                        print('Invalid SNP, ignore this line')
                        continue 
                    
                    info = line[:9] + genotypes
                    # remove useless info (QUAL, INFO, FILTER, FORMAT)from line 
                    info_pt1 = info[:5]
                    info_pt2 = info[9:]
                    info = info_pt1 + info_pt2
                    mod_l = '\t'.join(info)
                    last_char = mod_l[-1]
                    # replace \t at end with \n
                    mod_l = mod_l[:-1] + last_char + '\n'
                    lines.append(mod_l)
    #print(''.join(lines[:10]))
    # read the info from vcf file into a pandas dataframe
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

# this function will read in path of vcf file and convert it 
# to a df with genotype info, stored as geno.csv
def genoDf(path):
    print('Creating Geno Dafarame...')
    if ('gz' in path):
        vcf_df = read_vcf(path, 'gzip')
    else:
        vcf_df = read_vcf(path, 'vcf')
        
    vcf_df.to_csv('test_multi.csv', index=False)
    
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

def test_multi_genotype():
    print(assign_multi_genoType('G', 'A,T,C', ['0|0', '1|0', '0|2', '2|0', '1|2','0|3']))

# test code for assigning genotype by reading from csv file extracted from vcf file
def analyze_multi_genotype():
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