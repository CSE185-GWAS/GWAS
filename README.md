# The GWAS Graph Tool: Performing and Graphing GWAS
### Tool Description & Set Up:
Our GWAS Tool can read in a VCF file, a phenotype file (in csv or tsv), performing genome-wide
association studies by linear regression and outputting a csv file of GWAS results. Optionally,
users can pass in arguments to choose the type of plots (e.g. qq plot) would like to generate as
a .png file. The tool is be implemented in Python and hosted on Github that you can clone the repository and run
the associated script through the command-line argument, specifically:

```
git clone https://github.com/CSE185-GWAS/GWAS.git # clone script to your local computer  
python GWAS.py --vcf <input_file> --phen <phenotype_file> –o <output_file_name> -g <graph_name> # run the command when you're in the directory of gwas folder you clone 
```

```--vcf <input_file>``` : used to mark **path to the input vcf file** that contains information about SNP and genotype information of samples

```--phen <phenotype_file>``` or ```-p <phenotype_file>```: used to mark the **path to the phenotype tsv or csv file**

```--out <output_file_name>``` or ```-o <output_file_name>```: used to mark **name of GWAS output file**

```--graph_type <graph_name>``` or```-g <graph_name>```: used to mark the name of plot user want. It can be either **'qq'** or **'manhattan'**. This is optional tag that user can add to see output graphs. 

If you are using VScode for running our tool: 
As our tool uses packages like pandas, below are command-line argument you can refer to for set up a python environment:

```
python3 -m venv studysession
source /Users/apple/Desktop/GWAS #this will vary depend on directory that you download our tool 
pip install pandas
pip install matplotlib
pip install statsmodel
```

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
  
For current progress, after you download the dataset from lab3 and store them in same directory of out script, you can use:
  
```
python GWAS.py --vcf lab3_gwas.vcf.gz --phen lab3_gwas.phen --o out.txt
```
  
to see the genotype file from vcf and p-value for first SNP.  
  
### Datasets:
We plan to run our tools on the Lab 3 dataset consisting of genotype data from the 1000
Genome Project and the simulated phenotype data for LDL cholesterol level. We also plan to
use this public dataset, which was used to analyze the phenotype of hippocampus, prefrontal
cortex, and striatum with genotype of outbred CFW mice:
- Dataset: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2rs41
- Associated paper: https://www.nature.com/articles/ng.3609
