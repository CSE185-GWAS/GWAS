# this file can read in vcf or vcf.gzip file, read in, analyze, and output p-value and beta to output files
# reference: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

import gzip 
import io
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import p_val as pval
# this function destruct geno and return if it needs to be deleted for containing minor alleles 
def destructGeno(geno, allowed_alle, minor_alleles):
    for i in range(2):
        # find allowed allele in geno
        ref_pos = geno.find(allowed_alle[0])
        alt_pos = geno.find(allowed_alle[1])
    
        # if none of allowed geno found, return deleted = TRUE
        comp_len = 0                     
        if (ref_pos != 0 and alt_pos != 0):
            return True
        # if both refers to first position, use the longer one 
        elif(ref_pos == 0 and alt_pos == 0):
            comp_len = max(len(allowed_alle[0]), len(allowed_alle[1]))
            # but if using larger length make it = 0
            if ((i == 0 and geno[comp_len:] == '')):
                comp_len = min(len(allowed_alle[0]), len(allowed_alle[1]))
        # otherwise, record its length to compare with potential minor allele 
        elif(ref_pos == 0):
            comp_len = len(allowed_alle[0])
        else:
            comp_len = len(allowed_alle[1])
    
        # store minor alleles position
        minor_pos = []
        for minor in minor_alleles:
            minor_pos.append(geno.find(minor))
           
        minor_alt = ''
        for j in range(len(minor_pos)):
            # once we recognize other minor allele existence
            if (minor_pos[j] == 0):
                minor_alt = minor_alleles[j]
   
        # if there is a minor allele start at pos 0, and length higher than ref/alt alllele start at pos 0
        # geno contains minor allele that needs to be removed
        if len(minor_alt) > comp_len:
            # when checking for first allele, make sure using minor alt does not use up all nucleotide
            if ((i == 0 and geno[len(minor_alt):] != '')):
                return True
            
            # when it's second allele, it is definitely the minor alt so return 1
            if (i == 1 and geno[len(minor_alt):] == ''):
                return True 
        # update geno
        geno = geno[comp_len:]
    return False  
    
#**Function that will output all genotype, even by >1 alternative allele*****#
# alt can be multiple separated by ','
# this will assign genotype to sample based on query '0|1', '1|2'
# then it computes most significant alt allele and changed all other genotype
# with other alt allele to 'N'
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
        if '|' in query:
            query = query.split('|')
        elif '/' in query:
            query = query.split('/')
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
        minor_allele = list(dictionary.values())
        minor_allele.remove(ref)
        minor_allele.remove(sig_alt_allele)

        dictionary = [ref, sig_alt_allele]
        # check genotype with other insignificant alt allele by
        # checking if it has invalid length
        maxLen = max(len(ref), len(sig_alt_allele))
        maxLen *= 2
        minLen = min(len(ref), len(sig_alt_allele))
        minLen *= 2
        # or checking it has allele nonexistent in updated dictionary 
        for i in range(len(genotypes)):
            geno = genotypes[i]
            # check for invalid length genos
            if len(geno) > maxLen or len(geno) < minLen:
                genotypes[i] = 'N'
                continue
            # check for invalid allele 
            if destructGeno(geno, dictionary, minor_allele):
                genotypes[i] = 'N' 
                continue 
          
    return genotypes, sig_alt_allele

# used to read in each line of vcf, find genotypes and phenotypes,
# filter out invalid snp with low maf, and conduct linear regression
# to find the beta and pvalue, returning output gwas df 
def read_vcf(path, file_format, phen, maf_threhold=0.05):
    # record output file lines
    lines = []
    # used to record removed snps
    removed_snps = []
    # map sample id to phenotype value
    pheno_dict = pval.createDictFromPhenotype(phen)
    # store sample id that helps to map to phenotype  
    geno_cols = []
    result_file = open('out.txt', 'w')
    # if input fie is zipped, unzipped and analyze it 
    if (file_format == 'gzip'):  
        with gzip.open(path, 'rt') as f:
            for l in f:
                # when column headers found
                if (l.startswith('#CHROM')):
                    # find sample ids 
                    l = l.split('\t')
                    # keep the geno_cols 
                    geno_cols = l[9:]
                    # remove \n for last sample id
                    geno_cols[-1] = geno_cols[-1].strip()
                    # rename column names for output files 
                    l_list = ['CHR', 'BP','SNP','REF','A1','BETA', 'P\n']
                    l = '\t'.join(l_list)
                    result_file.write(l)
                    lines.append(l)
                # when snp rows are found
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split('\t')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]

                    # assign genotypes based on input
                    genotypes, sig_alt_allele = assign_multi_genoType(ref_allele, alt_allele, queires)
                    # if there are multi alt allele exist
                    # reformat alt alelle to most significant alt allele
                    if sig_alt_allele != '':
                        line[4] = sig_alt_allele
                        alt_allele = sig_alt_allele
                    
                    # calculate maf
                    maf = pval.calculate_maf(ref_allele, alt_allele, genotypes)
                    # remove snp if maf is too low 
                    if (maf < maf_threhold):
                        removed_snps.append(line[2])
                        print('maf is too low, removed from list')
                        continue
       
                    print(f'snp {line[2]}')
                    # print genotype mapping for linear regression
                    genotype_mapping = {ref_allele+ref_allele: 0, ref_allele+alt_allele: 1, alt_allele+ref_allele: 1, alt_allele+alt_allele: 2}
                    # find beta and pvalue by linear regression 
                    obs_beta, p_value = pval.getSingleP(genotypes, pheno_dict, geno_cols, genotype_mapping)
                    # remove useless info (QUAL, INFO, FILTER, FORMAT)from line, append beta and pval
                    info = line[:5] + [str(obs_beta), str(p_value) + '\n']
                    mod_l = '\t'.join(info)
                    # write to output files
                    result_file.write(mod_l)
                    lines.append(mod_l)
    # same as above, but directly use file, no step to unzip 
    else: 
        with open(path, 'r') as f:
            for l in f:
                if (l.startswith('#CHROM')):
                    # changed column names for output files
                    l = l.split('\t')
                    geno_cols = l[9:]
                    geno_cols[-1] = geno_cols[-1].strip()
                    l_list = ['CHR', 'BP','SNP','REF','A1','BETA', 'P\n']
                    l = '\t'.join(l_list)
                    print(l)
                    result_file.write(l)
                    lines.append(l)
              
                elif (not l.startswith('##') and not l.startswith('#CHROM')):
                    line = l.split('\t')
                    # remove \n character for last element 
                    line[-1] = line[-1].strip()
                    # grab queries, ref allele and alt allele
                    queires =  line[9:]
                    ref_allele  =  line[3]
                    alt_allele  =  line[4]

                    # assign genotypes based on input
                    genotypes, sig_alt_allele = assign_multi_genoType(ref_allele, alt_allele, queires)
                    # if there are multi alt allele exist
                    # reformat alt alelle to most significant alt allele
                    if sig_alt_allele != '':
                        line[4] = sig_alt_allele
                        alt_allele = sig_alt_allele
                    
                    # calculate maf
                    maf = pval.calculate_maf(ref_allele, alt_allele, genotypes)
                    if (maf < maf_threhold):
                        removed_snps.append(line[2])
                        print('maf is too low, removed from list')
                        continue

                    print(f'snp {line[2]}')
                    genotype_mapping = {ref_allele+ref_allele: 0, ref_allele+alt_allele: 1, alt_allele+ref_allele: 1, alt_allele+alt_allele: 2}
                    obs_beta, p_value = pval.getSingleP(genotypes, pheno_dict, geno_cols, genotype_mapping)
                    
                    info = line[:5] + [str(obs_beta), str(p_value) + '\n']
                    # remove useless info (QUAL, INFO, FILTER, FORMAT)from line 
                    mod_l = '\t'.join(info)
                    result_file.write(mod_l)
                    lines.append(mod_l)
    # return result as a pandas dataframe 
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'CHR': str, 'BP': int, 'SNP': str, 'REF': str, 'A1': str, 'BETA':float, 'P':float},
        sep='\t'
    )

# this function will read in path of vcf file and convert it 
# to a df with gwas linear output 
def genoDf(path, phen, outPath, maf=0.05):
    print('Creating Geno Dafarame...')
    if ('gz' in path):
        vcf_df = read_vcf(path, 'gzip', phen, maf_threhold=maf)
    else:
        vcf_df = read_vcf(path, 'vcf', phen, maf_threhold=maf)
        
    vcf_df.to_csv(outPath, index=False)
    
    print('Analysis complete, please check output file for details')
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
    # print(vcf_df)
    vcf_df.to_csv('final_vcf.csv')



