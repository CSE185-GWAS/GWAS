# The GWAS Graph Tool: Performing and Graphing GWAS
### Tool Description:
Our GWAS Tool can read in a VCF file, a phenotype file (in csv or tsv), performing genome-wide
association studies by linear regression and outputting a csv file of GWAS results. Optionally,
users can pass in arguments to choose the type of plots (e.g. qq plot) would like to generate as
a .png file. We would compare it with the plink GWAS tool that is used in Lab 3. The tool will be
implemented in Python and hosted on Github such that users can clone the repository and run
the associated script through the command-line, like:

```
git clone https://github.com/CSE185-GWAS/GWAS.git # clone script to your local computer  
python GWAS.py --vcf <input_file> --phen <phenotype_file> –o <output_file_name> -g <graph_name> # run the command when you're in the directory of gwas folder you clone 
```

```--vcf <input_file>``` : used to mark **path to the input vcf file** that contains information about SNP and genotype information of samples

```--phen <phenotype_file>``` or ```-p <phenotype_file>```: used to mark the **path to the phenotype tsv or csv file**

```--out <output_file_name>``` or ```-o <output_file_name>```: used to mark **name of GWAS output file**

```--graph_type <graph_name>``` or```-g <graph_name>```: used to mark the name of plot user want. It can be either **'qq'** or **'manhattan'**. This is optional tag that user can add to see output graphs. 


### Benchmarking:
We plan to benchmark the tool by applying our tool to the subset of 1000 Genome Project data
for LDL cholesterol level provided in Lab 3 and comparing our GWAS result with the output in
Lab 3 using plink. In addition, we will use a public dataset, which contains information about the
phenotype of hippocampus, prefrontal cortex, and striatum of outbred CFW mice, and check if
our output graph is similar to the GWAS graphs published in the relevant paper that uses the
same dataset.

We would also compare the memory consumption of the plink tool used in Lab 3 versus our tool
using the terminal command ‘ps aux | grep <program name>’. Since our group members are
also using different operating systems (MacOS and Windows), we also want to examine the
memory performance differences between our machines.
  
### Datasets:
We plan to run our tools on the Lab 3 dataset consisting of genotype data from the 1000
Genome Project and the simulated phenotype data for LDL cholesterol level. We also plan to
use this public dataset, which was used to analyze the phenotype of hippocampus, prefrontal
cortex, and striatum of outbred CFW mice:
- Dataset: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2rs41
- Associated paper: https://www.nature.com/articles/ng.3609
