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
plink --linear --vcf lab3_gwas.vcf.gz --pheno lab3_gwas.phen --maf 0.05  --out lab3_gwas --allow-no-sex 
```
The plots and output files are generated in this directory.

#### From Palmer Lab
Run our tool on Palmer Lab data. From the Palmer Lab website, we got SNP data and phenotype information for a Pavlovian conditioned (PavCA) approach was applied to a sample of 4,608 Sprague Dawley Rats. 

`allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf.gz` is the genotype data consists of SNPs from rats produced by the vendor Charles River

`palmer_dataset_phen.csv` is the weight values of these rats extracted from the “Filtered Samps - Phenos & Covar” sheet within the file “SD_PavCA_SampleInfo_PhenotypeInfo.xlsx”

Again, run this command to generate results:
```
GWAS-py --vcf benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf.gz --phen benchmark/palmer_dataset_phen.csv -o out2.csv -m -q
```
Running with PLINK:
```
plink --linear --vcf benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf.gz --pheno pheno/Palmer_Lab.phen  --maf 0.05 --out Palmer_Lab --allow-no-sex
```
Same as above, output files are generated in this directory.
