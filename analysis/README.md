This directory contains files for preprocessing:
* lab3_gwas.txt: the original simulated LDL cholesterol file from Lab 3
* Palmer_Lab_Pheno.csv: the .csv file from the sheet "Filtered Samps - Phenos & Covar" within “SD_PavCA_SampleInfo_PhenotypeInfo.xlsx” from the Palmer Lab dataset
* preprocess_lab3_pheno.py: script used to convert the Lab 3 phenotype file into the correct formatting for use with our tool
* preprocess_palmer_pheno.py: script used to extract weight and convert the Palmer Lab phenotype file into the correct formatting for use with our tool


This directory also contains files for comparing output from GWAS and plink
* compare_result.ipynb

Please change input files paths at first two cells to compare results as you own need. The script will compare two output files' p value column, beta value column, and their snp column from GWAS linear analysis. 
