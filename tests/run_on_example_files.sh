#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
echo "$parent_path"
cd "$parent_path"
python "../GWAS/GWAS.py" --vcf "../example-files/chr10.recode.vcf.gz" --phen "../benchmark/lab3_phen.csv" --out "example_file_results.csv"