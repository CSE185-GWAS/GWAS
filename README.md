# The GWAS Graph Tool: Performing and Graphing GWAS
### Tool Description & Set Up:
Our GWAS Tool can read in a VCF file, a phenotype file (in csv format), performing genome-wide
association studies by linear regression and outputting a csv file of GWAS results. Optionally,
users can pass in arguments to choose the type of plots (e.g. qq plot) would like to generate as
a .png file. The tool is be implemented in Python and hosted on Github that you can clone the repository and run
the associated script through the command-line argument (when you're in our GWAS tool folder you cloned), specifically:

```
git clone https://github.com/CSE185-GWAS/GWAS.git # clone script to your local computer  
python GWAS.py --vcf <input_file> --phen <phenotype_file> –o <output_file_name> -g <graph_name> # run the command when you're in the directory of gwas folder you clone 
```

```--vcf <input_file>``` : used to mark **path to the input vcf file** that contains information about SNP and genotype information of samples

```--phen <phenotype_file>``` or ```-p <phenotype_file>```: used to mark the **path to the phenotype csv file** The phenotype file should have rows as sample ID and columns as phenotype name (ex. LDL)

```--out <output_file_name>``` or ```-o <output_file_name>```: used to mark **name of GWAS output file**


### On progress:
```--maf <maf_threhold>```: used to mark the threhold for minor allele frequency of SNP. By default (if you don't specify this tag), this threhold value is set to 0.05 to filter out SNP with low maf. For now, users' input won't affect maf we use yet. 

```--graph_type <graph_name>``` or```-g <graph_name>```: used to mark the name of plot user want. It can be either **'qq'** or **'manhattan'**. This is optional tag that user can add to see output graphs. We have functions available for plotting, but have not connected with our main script. 


### If you are using VScode for running our tool: 
As our tool uses packages like pandas, below are command-line argument you can refer to for set up a python environment in the GWAS folder you cloned from this repository:

```
python3 -m venv studysession
source /Users/apple/Desktop/GWAS/studysession/bin/activate #this will vary depend on directory that you download our tool 
pip install pandas
pip install matplotlib
pip install statsmodels
pip install qqman
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
  
For current progress, after you download the dataset from lab3 (both lab3_gwas.vcf.gz and lab3_gwas.vcf.gz.tbi) and store them in same directory of out script (GWAS folder you cloned), you can use:
  
```
python GWAS.py --vcf lab3_gwas.vcf.gz --phen testcases/lab3_phen.csv -o lab3_results.csv
```

Alternatively, you can use path to the lab3 vcf file you downloaded to replace ```lab3_gwas.vcf.gz``` if you don't have it on the folder of GWAS. Please make sure you have ```lab3_gwas.vcf.gz``` and ```lab3_gwas.vcf.gz.tbi``` in the same directory though.
  
In the process of analyzing snps, our tool will also generate a text output file called ```out.txt```. As analyzing output of lab3 would take for a while (~1h), to see its progress, you can refer to the ```out.txt```. 

Also, if you don't like seeing messages we print on terminal, feel free to use  
  
```
nohup python GWAS.py --vcf lab3_gwas.vcf.gz --phen testcases/lab3_phen.csv -o lab3_results.csv &
```
instead.
  
We also put the lab3 output csv file generated by following command in the ```lab3_results.csv.gz``` . Feel free to take a look over it by using:
  ```
  gzip -d lab3_results.csv.gz
  ```
when you're in our GWAS folder. 

### On Progress:  
We have not resolved some issue when directly passing in vcf file, so for now please make sure the your input file is in gz format along with its index file. *The input vcf file format should have .gz extension and along with its index gz.tbi*. 

### Datasets:
We plan to run our tools on the Lab 3 dataset consisting of genotype data from the 1000
Genome Project and the simulated phenotype data for LDL cholesterol level. We also plan to
use this public dataset, which was used to analyze the phenotype of hippocampus, prefrontal
cortex, and striatum with genotype of outbred CFW mice:
- Dataset: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2rs41
- Associated paper: https://www.nature.com/articles/ng.3609
