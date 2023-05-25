import random
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as stats

#input should be a list of beta values, calculate p value
def ComputeP(gts, pts, gt_counts):
    sizeN = sum(gt_counts.values())
    maf = ComputeMinorAlleleFrequency(gt_counts)
    
    # Get null list of beta values
    gts_null = SimulateGenotypes(maf, sizeN)
    num_tests = 10000
    beta_null_dist = []
    for i in range(num_tests):    
        pts_null = SimulatePhenotype(gts_null, 0)
        gwas_beta = LinReg(gts_null, pts_null)
        beta_null_dist.append(abs(gwas_beta))
    
    # Observed beta value
    gwas_beta = LinReg(gts, pts)
    pval = one_sample_ttest(gwas_beta, beta_null_dist)
    return pval

def one_sample_ttest(observed, null_values):
    # Calculate the t-statistic and p-value
    p_value = stats.ttest_1samp(observed, null_values.mean())
    return p_value

def ComputeMinorAlleleFrequency(gt_counts):
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

# simulates an array of randomly generated genotypes for N people using a specified minor allele frequency maf
def SimulateGenotypes(maf, N):
    gts = []
    for i in range(N):
        gt = sum(random.random() < maf)+sum(random.random() < maf)
        gts.append(gt)
    # Technical note: scale to have mean 0 var 1
    gts = np.array(gts)
    gts = (gts-np.mean(gts))/np.sqrt(np.var(gts))
    return gts

# Beta should be 0 for null distribution
def SimulatePhenotype(gts, Beta):
    if Beta<-1 or Beta>1:
        print("Error: Beta should be between -1 and 1")
        return [None]*len(pts)
    pts = Beta*gts + np.random.normal(0, np.sqrt(1-Beta**2), size=len(gts))
    return pts

def LinReg(gts, pts):
    X = sm.add_constant(gts)
    model = sm.OLS(pts, X)
    results = model.fit()
    beta = results.params[1]
    return beta

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