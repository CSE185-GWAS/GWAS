import pandas as pd

if __name__=='__main__':
    df = pd.read_csv('Palmer_Lab_Pheno.csv')
    df = df[df['Vendor'] == 'Charles River']
    print(len(df))
    print(df.isna())
    df = df.dropna(subset='Weight')
    print(len(df))
    cols_to_drop = list(df.columns)
    cols_to_drop.remove('Sample')
    cols_to_drop.remove('Weight')
    df = df.drop(cols_to_drop, axis=1)
    df.to_csv('Palmer_dataset_Pheno_Test.tsv', index=False)
