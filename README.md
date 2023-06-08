# The GWAS Graph Tool: Performing and Graphing GWAS
### Tool Description:
Our GWAS Tool can read in a VCF file, a phenotype file (in csv format), performing genome-wide
association studies by linear regression and outputting a csv file of GWAS results. Optionally,
users can pass in arguments to choose the type of plots (e.g. qq plot) would like to generate as
a .png file. 

It is a python package hosted on Github that you can clone the repository and run
the associated script through the command-line argument (when you're in our GWAS tool folder you cloned). 

### Download & Set up Environment 
```
git clone https://github.com/CSE185-GWAS/GWAS.git # clone our tool to your local computer 
cd GWAS # go to our tools folder 
pip install pandas # set up package environment, you can also refers to requirements.txt
pip install matplotlib
pip install statsmodels
pip install qqman
python setup.py install --prefix=$HOME # set up our package 
pwd # this will return the absolute path to our tool on local computer
# For Linux / mac users, replace /path/to/your/project/ with what pwd command returns above 
export PYTHONPATH="${PYTHONPATH}:/path/to/your/project/"

# For Windows non-Linux users, replace /path/to/your/project/ with what pwd command returns above 
set PYTHONPATH=%PYTHONPATH%;C:\path\to\your\project\
```

### If you are using VSCode, and see error messages that *pip* command is not found?
To set up a python environment on your local computer, you can use command like
below:
```
python3 -m venv studysession
source /Users/apple/Desktop/GWAS/studysession/bin/activate #this will vary depend on directory that you download our tool 
# download all required python modules  
pip install pandas
pip install matplotlib
pip install statsmodels
pip install qqman
```

### If you are download the package by instruction, but it is showing *GWAS-py command not found*?
Our apology that we have not figure out why it happens for some of our users. We are still working on it and would update instruction once we understand the issue better. 
Alternatively, please directly use our main script as below:

```python GWAS/GWAS.py --vcf <input_file> --phen <phenotype_file> –o <output_file_name> -m -q  ```

Essentially, just replace our package name with main script. A simple use case would be: 

```python GWAS/GWAS.py --vcf testcases/test.vcf --phen benchmark/lab3_phen.csv -o testcases/test.csv -q -m ```


### Tool Usage: 

```GWAS-py --vcf <input_file> --phen <phenotype_file> –o <output_file_name> -m -q  ```

### Required input tags 

```--vcf <input_file>``` : used to mark **path to the input vcf file** that contains information about SNP and genotype information of samples. (ex. input.vcf.gz / input.vcf)

```--phen <phenotype_file>``` or ```-p <phenotype_file>```: used to mark the **path to the phenotype csv file** (ex. lab3_phen.csv)
*The phenotype file should have rows as sample ID and columns as phenotype name (ex. LDL)*

***need help with preprocess phenotype file?***

As suggested, our tool require phenotype file to have sample ID as rows and phenotype values as columns. 
If you are kind of confused on how to pre-process your data, feel free to look through the script in the ```analysis``` directory and look at how we preprocess different phenotype text files into a csv pheno file our tool would take. 

### Optional output and filtering tags 

```--out <output_file_name>``` or ```-o <output_file_name>```: used to mark **directory and the name of GWAS output file** By default, we set it to *out.csv*. (ex. result.csv)


```--maf <maf_threhold>```: used to mark the threhold for minor allele frequency of SNP. By default, this threhold value is set to 0.05 to filter out SNP with low maf. (ex. 0.03)

### Graphing utility tags

```-q```: used to mark that user want a  **'qq'** plot of gwas result.

```-m```: used to mark that user want a **'manhattan'** plot of gwas result 

You can use both tag at same command to mark that you want both graphics. 


### Simple Use Case:
A simple use case of our tool with small dataset you can test would be:
```
GWAS-py --vcf testcases/test.vcf --phen benchmark/lab3_phen.csv -o testcases/test.csv -q -m
```
We store the example output files in the README in testcases. Feel free to take a look over it.

### Test scripts: 
You can also use 
 ```
  chmod 777 ./tests/run_on_example_files.sh
  ./tests/run_on_example_files.sh
```

to test a medium-size dataset. We put the output csv file in the tests/ directory so you can take a look over it before proceeding. 

### Benchmark:
We benchmark our tool by applying our tool to the subset of 1000 Genome Project data
for LDL cholesterol level provided in Lab 3 and comparing our GWAS result with the output in Lab 3 using plink --linear with maf of 0.05. 

To benchmark our tool, type 
  
```
GWAS-py --vcf benchmark/lab3_gwas.vcf.gz --phen benchmark/lab3_phen.csv -o lab3_results.csv -m -q 
```
  
In the process of analyzing snps, our tool will also generate a text output file called ```out.txt```. As analyzing output of lab3 would take for a while (~30 minutes), to see its progress, you can refer to the ```out.txt```. 

Also, if you don't like seeing messages we print on terminal, feel free to use  
  
```
nohup GWAS-py --vcf benchmark/lab3_gwas.vcf.gz --phen benchmark/lab3_phen.csv -o lab3_results.csv -m -q  &
```
instead.
  
We also put the lab3 output csv file generated by following command in the ```benchmark/lab3_results.csv.gz``` . Feel free to take a look over it by using:
  ```
  gzip -d benchmark/lab3_results.csv.gz
  ```
when you're in our GWAS folder. You can also find results for the qqplot and marathanplot, as well as result for another public dataset from Palmar Lab that we tested with in the benchmark directory.


### other Dataset:
To test our tool with other dataset, type 
  
```
GWAS-py --vcf benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf --phen benchmark/palmer_dataset_phen.csv -o palmer_daraset.csv -m -q
```
We also put the output in the ```palmer_dataset.csv```. If you want to run plink on this dataset too, please use ```benchmark/allChr.allSamps.90DR2.maf01.hweE7.noIBD.CharlesRiverOnly.vcf``` as well as ```pheno/```
