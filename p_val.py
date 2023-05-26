# this file 
# reference: lab 3  
import random
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as stats
import pandas as pd

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
    return dict(pheno_df.values)

# thsi function will return two lists where one will be the genotype values 
# and the other will be their corresponding phenotype values; used to association testing 
def generateGenoTypeAndPhenotype(pheno_dict, geno_df):
    pheno_val = []
    geno_cols = list(geno_df.columns)

    for geneIds in geno_cols:
        if geneIds == 'REF' or geneIds == 'ALT':
            continue
        else:
            pheno_val.append(pheno_dict[geneIds])
    # for now, only look at first SNP
    snps = geno_df.iloc[0]

    ref_allele = geno_df['REF'].iloc[0]
    alt_allele = geno_df['ALT'].iloc[0]
    genotype_mapping = {ref_allele+ref_allele: 0, ref_allele+alt_allele: 1, alt_allele+ref_allele: 1, alt_allele+alt_allele: 2}

    geno_val = []
    gt_counts = {ref_allele+ref_allele: 0, ref_allele+alt_allele: 0, alt_allele+ref_allele: 0, alt_allele+alt_allele: 0}
    for genotype in snps:
        # omit invalid genotype
        if (len(genotype) > 1):
            geno_val.append(genotype_mapping[genotype])
            gt_counts[genotype] += 1
            
    return geno_val, pheno_val, gt_counts


def generateMAF(gt_counts):
    print('Finding minor allele frequencies of current SNP...')
    gts = gt_counts.keys()
    alleles = set("".join(gts))
    # Compute minor allele frequency
    maf = -1
    # calculate two allele freuqnecy separately
    allele_list = list(alleles)
    allele_1 = allele_list[0]
    allele1_freq = 0
    allele2_freq = 0
    total = 0
    for gt in gts:
        freq = gt_counts[gt]
        total = total + freq
        for allele in gt:
            if (allele == allele_1):
                allele1_freq = allele1_freq + freq
            else:
                allele2_freq = allele2_freq + freq
    allele1_freq = allele1_freq / (total * 2)
    allele2_freq = allele2_freq / (total * 2)
    
    # intialize maf
    maf = min(allele1_freq,allele2_freq)
    print('maf is {}'.format(maf))
    return maf

def generateNullDistribution(maf, N):
    # generate null genotypes
    null_geno = []
    for i in range(N):
        gt = int(random.random() < maf)+int(random.random() < maf)
        null_geno.append(gt)
    # Normalize geno
    null_geno = np.array(null_geno)
    null_geno = (null_geno-np.mean(null_geno))/np.sqrt(np.var(null_geno))
    
    # generate null phenotypes (normalized)
    null_pheno = np.random.normal(0, np.sqrt(1), size=len(null_geno))

    return null_geno, null_pheno


#input should be a list of beta values, calculate p value
def caculateSingleP(gts, pts, gt_counts):
    size = sum(gt_counts.values())
    maf = generateMAF(gt_counts)
    
    # Get null list of beta values
    num_tests = 10000
    beta_null_dist = []
    for i in range(num_tests):    
        null_geno,null_pheno= generateNullDistribution(maf, size)
        gwas_beta = LinReg(null_geno, null_pheno)
        beta_null_dist.append(abs(gwas_beta))
    beta_null_dist = np.array(beta_null_dist)
    
    # Observed beta value
    gwas_beta = LinReg(gts, pts)
    print('Observed beta value is: {}'.format(gwas_beta))
    print("Computing pval by simulated null distribution...")
    print(beta_null_dist)
    pval = findPval(gwas_beta, beta_null_dist)
    return pval

def findPval(observed, null_values):
    # Calculate the t-statistic and p-value
    pval = None
    # your code here
    total = len(null_values)
    count = 0
    # calculate pval by |beta|
    observed = abs(observed)
    for beta in null_values:
        if (beta >= observed):
            count = count + 1

    pval = count / total
    return pval

def LinReg(gts, pts):
    X = sm.add_constant(gts)
    model = sm.OLS(pts, X)
    results = model.fit()
    beta = results.params[1]
    return beta


def calculatePVal(pheno, geno_df):
    pheno_dict = createDictFromPhenotype(pheno)
    # only keep essential info from geno_vcf for t-test 
    geno_df = geno_df.iloc[:, 3: ].drop(columns=['QUAL','FILTER','INFO','FORMAT'])
    geno_val, pheno_val, gt_counts = generateGenoTypeAndPhenotype(pheno_dict, geno_df)
    p_value = caculateSingleP(geno_val, pheno_val, gt_counts)
    print('p_value for first snp is :')
    print(p_value)
    return p_value


# Not called yet
def QQPlot(pvals):
    pvals.sort()
    unif = list(np.random.uniform(0, 1, size=len(pvals)))
    unif.sort()

    # Make a QQ plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(-1*np.log10(unif), -1*np.log10(pvals), s=5, color="black")
    ax.plot([0, 3], [0,3])
    ax.set_xlim(left=0, right=3)
    ax.set_ylim(bottom=0, top=max(-1*np.log10([item for item in pvals if item >0])))
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")












#**********Optional TODO: Account for multiple alter allele*******************#
def multi_maf(gt_counts):
    gts = gt_counts.keys()
    print("Genotype frequencies:")
    for gt in gts:
        print("%s: %.2f"%(gt, gt_counts[gt]/sum(gt_counts.values())))

    alleles = set("".join(gts))
    print("Possible alleles: %s"%alleles)
    # Compute minor allele frequency
    a_list = list(alleles)
    # count the first possible allele
    freq = [0, 0, 0, 0]
        
    for gt in gts:
        base = list(gt)
        for i in range(len(a_list)):
            if base[0] == a_list[i]:
                freq[i] += gt_counts[gt]
            if base[1] == a_list[0]:
                fre_1 += gt_counts[gt]
        
        max_index = freq.index(max(freq))
        freq[max_index] = float('-inf')
        second_max_index = freq.index(max(freq))
        allele_list = list(alleles)
        top_two = [allele_list[max_index], allele_list[second_max_index]]
        filtered_gt_counts = {key: value for key, value in gt_counts.items() if all(char in top_two for char in key)}

        # look for maf in the top two bases
        total_a = 2 * sum(filtered_gt_counts.values())
        fre_1 = 0
        fre_2 = 0
        filtered_gts = filtered_gt_counts.keys()
        for filtered_gt in filtered_gts:
            ind = list(filtered_gts)
            if ind[0] == top_two[0]:
                fre_1 += filtered_gt_counts[filtered_gt]/total_a
            elif ind[0] == top_two[1]:
                fre_2 += filtered_gt_counts[filtered_gt]/total_a
            if ind[1] == top_two[0]:
                fre_1 += filtered_gt_counts[filtered_gt]/total_a
            elif ind[1] == top_two[1]:
                fre_2 += filtered_gt_counts[filtered_gt]/total_a
                
        maf = min(fre_1, fre_2)
        return maf

