# this file include statistical script for compute maf, assign genotype, and conduct linear regressions 
# reference: lab 3  
import math 
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd
from qqman import qqman


# this function will compute maf based on ref allele, alt allele, and genotypes
# genotypes = list of genotypes ['AA', 'TT']
def calculate_maf(ref_allele, alt_allele, genotypes):
    # calculate maf
    maf = -1
    ref_allele_freq = 0
    alt_allele_freq = 0
    total = 0
    ref_l = len(ref_allele)
    alt_l = len(alt_allele)
    for genotype in genotypes:
        # skip empty genotypes 
        if (genotype == 'N'):
            continue
        # update total # of alleles by 2
        total += 2
        for i in range(2): 
            # find index of allele position, check first allele position 
            ref_index = genotype.find(ref_allele)
            alt_index = genotype.find(alt_allele)
            if (ref_index == 0 and alt_index == 0):
                if (ref_l > alt_l):
                    ref_allele_freq += 1
                    genotype = genotype[ref_l:]
                else:
                    alt_allele_freq+=1
                    genotype = genotype[alt_l:]
            elif(ref_index == 0):
            # if first allele is ref allele, split it using ref_l 
                ref_allele_freq += 1
                genotype = genotype[ref_l:]
            # otherwise, split it using alt_l
            else:
                alt_allele_freq+=1
                genotype = genotype[alt_l:]

    # if genotype all N, return 0
    if (total == 0):
        return 0
    # find minor allele via proportion 
    maf = min(ref_allele_freq/total, alt_allele_freq/total)
    return maf

# this function will return a dictionary where the keys are the sample IDs 
# and the values are their corresponding phenotype values
# phenofile = path to csv / tsv file use sampleID as row and phenotype as columns
# pheno_col = phenotype we want to test with 
def createDictFromPhenotype(pheno_file, pheno_col=''):
    pheno_df = pd.read_csv(pheno_file)
    # omit other columns
    cols_to_drop = []
    for col in pheno_df.columns:
        if col != pheno_col:
            cols_to_drop.append(col)
    pheno_df.drop(cols_to_drop, axis=1)
    dictionary = dict(pheno_df.values)

    if set(map(type, dictionary)) != {str}:
        dictionary = {str(int(key)): val for key, val in dictionary.items()}
    return dictionary

# conduct linear regression and get corresponding beta and pval 
def LinReg(gts, pts):
    X = sm.add_constant(gts)
    model = sm.OLS(pts, X)
    results = model.fit()
    beta, pvalue = results.params[1], results.pvalues[1]
    return beta, pvalue

# compute current snp's p value 
def getSingleP(genotypes, pheno_dict, geno_cols, genotype_mapping):
    # record gts nums and pheno nums for linear regression 
    pheno_vars = []
    geno_vars = []
    for i in range(len(genotypes)):
        genotype = genotypes[i]
        # if geno empty, skip 
        if (genotype == 'N'):
            continue 
        # find genotype's sample id to find corresponding phenotype
        # also convert genotype to number by genotype mapping
        else:
            if geno_cols[i] in pheno_dict:
                pheno_vars.append(pheno_dict[geno_cols[i]])
                geno_vars.append(genotype_mapping[genotype])
    # compute beta and pvalues        
    obs_beta, obs_pval = LinReg(geno_vars, pheno_vars)
    print("observed beta is: {}".format(obs_beta))
    print("observed pval is: {}".format(obs_pval))

    return obs_beta, obs_pval
# function to plot qq plot based on result df 
def QQPlot(df):
    pvals = df['P'].tolist()
    pvals.sort()
    unif = list(np.random.uniform(0, 1, size=len(pvals)))
    unif.sort()

    # Make a QQ plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(-1*np.log10(unif), -1*np.log10(pvals), s=5, color="black")
    largestX = math.ceil(max(-1*np.log10([item for item in unif if item > 0])))
    largestY = math.ceil(max(-1*np.log10([item for item in pvals if item > 0])))
    smallerDimension = min(largestX, largestY)
    ax.plot([0, smallerDimension], [0, smallerDimension])
    ax.set_xlim(left=0, right=largestX)
    ax.set_ylim(bottom=0, top=largestY)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    plt.savefig('qqplot.png')
# function to plot manhattan plot based on result df 
def manhattanPlot(df):
    #df needs to have columns 'CHR'(chromosome), 'BP'(basepair), 'P'(p-value)
    qqman.manhattan(df, out='manhattanplot.png')

