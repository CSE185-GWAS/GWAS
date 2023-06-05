import pandas as pd

if __name__=='__main__':
    #renamed the lab 3 phenotype file from lab3_gwas.phen to lab3_gwas.txt
    df = pd.read_csv('lab3_gwas.txt', sep='\t', index_col=0, header=None)
    df.columns = ['id', 'LDL']
    df.to_csv('../benchmark/lab3_phen.csv', index=False)