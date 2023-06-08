### Datasets
#### From Lab 3
Run our tool and PLINK using lab 3 1000 Genome Project for LDL Cholesterol level data as input

`lab3_gwas.vcf.gz` is a VCF file containing the LD-pruned SNPs for a subset of the 1000 Genomes samples

`lab3_phen.csv` is the normalized LDL values for each sample

Run this command to generate results. Here, we set maf to 0.05, which is the same as we did in lab 3. use `-m` and `-q` flags to generate the Manhattan and QQ plot:
```
GWAS-py --vcf benchmark/lab3_gwas.vcf.gz --phen benchmark/lab3_phen.csv -o lab3_results.csv -m -q
```
Generate output from PLINK:
```
PLINK —-linear --vcf lab3_gwas.vcf.gz --phen lab3_gwas.phen –maf 0.05 -o lab3_gwas.assoc.linear
```
The plots and output files are generated in this directory.

#### From Palmer Lab
Run our tool on Palmer Lab data. From the Palmer Lab website, we got SNP data and phenotype information for a Pavlovian conditioned (PavCA) approach was applied to a sample of 4,608 Sprague Dawley Rats. 

`allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf.gz` is the genotype data consists of SNPs from rats produced by the vendor Charles River

`palmer_dataset_phen.csv` is the weight values of these rats extracted from the “Filtered Samps - Phenos & Covar” sheet within the file “SD_PavCA_SampleInfo_PhenotypeInfo.xlsx”

Again, run this command to generate results:
```
GWAS-py --vcf benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf --phen benchmark/palmer_dataset_phen.csv -o out2.txt -m -q
```
Running with PLINK:
```
plink --vcf benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf --phen pheno/Palmer_Lab_Pheno_Test.tsv —maf 0.05 -o out2.txt
```
Same as above, output files are generated in this directory.
