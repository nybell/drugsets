# DRUGSETS
DRUg Gene SET AnalysiS (DRUGSETS) is a command line interface (CLI) tool implemented in python to perform genetically informed drug repositioning using drug gene set analysis using MAGMA.

Nathaniel Bell (n.y.bell@vu.nl)

## Prerequisites    
DRUGSETS was developed using [Python 3.8.5](https://www.python.org/) and [R 4.1.0](https://www.r-project.org/) with the following packages. Older versions of these packages will not work:    

Python  
[tqdm](https://tqdm.github.io) == 4.62.3   

R  
[tidyr](https://tidyr.tidyverse.org/) == 1.1.3  
[dyplr](https://dplyr.tidyverse.org/) == 1.0.7  
   
Additionally, the [MAGMA](https://ctg.cncr.nl/software/magma)[^1] software tool for gene and gene-set analysis is needed to run DRUGSETS. Full documentation of MAGMA can be found [here](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf). MAGMA needs to be able to be executed globally, which can be setup by adding the following line to the end of the user's .bashrc or .zshrc file.    
   
`export PATH=$PATH:/INSERT/PATH/TO/FILE/magma_v1.09b_mac/`   
   
## Setup   
   
If [Python 3.8.5](https://www.python.org/) is installed, all prerequisite Python packages can be installed using the following code:    
   
```
pip install --upgrade pip  
pip install tqdm==4.62.3  
```
   
Users can also create a virtual environment in which to run DRUGSEA using [Conda](https://www.anaconda.com/products/individual) or [pyenv](https://github.com/pyenv/pyenv). A guide to installing pyenv can be found [here](https://github.com/pyenv/pyenv), and a guide to installing Conda can be found [here](https://docs.anaconda.com/anaconda/navigator/tutorials/index.html).    
   
A virtual environment with all prerequisites can be created using pyenv via:    
    
```
pyenv virtualenv 3.8.5 venv_drugsets  
pip install --upgrade pip  
pip install tqdm==4.62.3  
pyenv activate venv_drugsets   
```   
A virtual environment with all prerequisites can be created using Conda via:    
   
```   
conda create --name venv_drugsea python=3.8.5 tqdm=4.62.3 numpy=1.21.1 pandas=1.3.1 scipy=1.7.1 scikit-learn=1.0 matplotlib=3.4.2
conda activate venv_drugsea
```   

The necessary R packages can be installed using the following code:

```
require(remotes)
install_version("tidyr", version = "1.1.3", repos = "http://cran.us.r-project.org")
install_version("dplyr", version = "1.0.7", repos = "http://cran.us.r-project.org")
```

## Setup & Input data     
   
1. Clone repository    
`$ git clone https://github.com/nybell/drugsea`   
   
2. Input data  
The only input data needed by the user is the [GENE_RESULTS].genes.raw file from MAGMA gene analysis. All drug data used in the drug gene set analysis is included in the download, and is stored in the directory `/DATA/GENESETS/`. ***The data needs to be stored in this directory in order to work.***    
   
## Running drug gene set analysis   
    
Drug gene set analysis is done by executing the script `drugsets.py`. The only flags required for running drug gene set analysis are `--geneassoc`, `--drugsets`, and `--out`. Enrichment analyses can be added using the optional flags `--enrich` and `--nsize`.     

* `--geneassoc` or `-g`: specifies the filepath and filename of the user's [GENE_RESULTS].genes.raw file (i.e., `/INSERT/PATH/TO/[GENE_RESULTS].genes.raw`)
* `--drugsets` or `-d`: specifies which type of drug genesets to use. There are four options:
    * `solo`: drug gene sets for each individual drug in our data. **NOTE:** if you want to test for enrichment you much use `solo`
    * `atc`: drug gene sets for each ATC III code category in our data 
    * `moa`: drug gene sets for each mechanism of action category in our data 
    * `ind`: drug gene sets for each clinical indication category in our data 
    * `all`: run drug gene-set analysis for everytype of drug gene set  
* `--enrich` or `-e`: tests all groups in the given category to see if drugs in each group have higher Z-statistics from the drug gene-set analysis than drugs not in the group (e.g., "do the drugs with ATC code N03A have higher Z-statistics from the drug gene-set analysis than drugs without ATC CODE N03A?"). Tests all drugs in the category type that is input (e.g., for input 'atc', all ATC code groups will be tested) and then is Bonferroni corrected for the number of groups tested. There are three input options:
    * `atc`: test drugs grouped by ATC code 
    * `moa`: test drugs grouped by mechanism of action 
    * `ind`: test drugs grouped by clinical indication 
    * `all`: test all types of drug groups
* `--conditional` or `-c`: Specifies whether or not to run competitive gene-set analysis while conditioning on a gene set of all druggable genes. Input options are `yes` and `no` (default = `yes`).
* `--id` or `-i`: Indicate which gene naming convention is used for your .genes.raw file. Options are "entrez" and "ensembl v105", and "ensembl v92". If you ran MAGMA using FUMA, then use "ensembl92".
* `--setsize` or `-s`: specfify minimum N for drug gene sets. 
* `--nsize` or `-n`: minimum sample size of drug categories to use when testing for enrichment (e.g., when set to 5, ATC code categories with less than 5 drugs will not be tested for enrichment). 
* `--correct` or `-p`: select Bonferroni or FDR correction for the drug group associations tests. Input options are `bonf` and `fdr`.     
* `--out` or `-o`: specify prefix for output files  
    
### Example usage    
    
The following code tests individual drug gene set analysis for associated with schizophrenia. The SCZ_SAMPLE.genes.raw was created using MAGMA and the 2021 schizophrenia GWAS summary statistics, downloaded [here](https://www.med.unc.edu/pgc/download-results/). Additionally, the code tests for ATC III codes for enrichment of drugs that are highly associated with schizophrenia, with a minimum sample size of 5 drugs for each ATC code. Lastly, it specifies to not show output from MAGMA (to see this go to the .log file in the directory `/drugsets-main/OUTPUT/`.       
  
1. Go to directory      
`$ cd /PATH/TO/drugsets-main/`     
  
2. Execute script    
`python ./drugsets.py --geneassoc /[PATH]/[TO]/[FILE]/SCZ_SAMPLE.genes.raw --drugsets solo --conditional yes --setsize 5 --enrich atc --nsize 5 --out /[PATH]/[TO]/[FILE]/SCZ_output`  

or the following if the working directory is not the drugsets-main directory:  

`python /[PATH]/[TO]/[DRUGSETS-MAIN]/drugsets.py --geneassoc /[PATH]/[TO]/[FILE]/SCZ_SAMPLE.genes.raw --drugsets solo --conditional yes --setsize 5 --enrich atc --nsize 5 --out /[PATH]/[TO]/[FILE]/SCZ_output`   
    
## Output    
   
The formats of the output files from the drug gene set analysis in MAGMA (.gsa.out, .gsa.genes.out .gsa.sets.genes.out) is described [here](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf) on page 23.    
    
#### Drug gene set analysis
* `.gsa.genes.out` contains information for each gene used in the MAGMA analysis   
* `.gsa.out` contains information for each drug gene set used in the MAGMA analysis    
* `.gsa.sets.genes.out` contains information on significant drug gene sets after Bonferroni correction, and the genes in each gene set  
    
#### Enrichment analysis 
* `lnreg.[GROUP].[OUT].csv` contains the results from the enrichment analysis for every drug category tested
* `lnreg.BONF.[GROUP].[OUT].csv` contains the results from the enrichment analysis for drug significant after Bonferroni correction   
  
Both enrichment output files contain the same columns:
* `GROUP` the drug category tested 
* `BETA` Beta value from the dependent linear regression model of the effect of group (e.g., which ATC III code a drug has) on the Z-statistic of the drug from the drug gene-set analysis run in MAGMA
* `SIGMA` Sigma value from the dependent linear regression model
* `T` Test statistic from the hypothesis H0: Bgroup = 0 (i.e., testing the null hypothesis that the Beta for a drug group is equal to 0)
* `P` Two-tailed p-value from the T statistic


# Support   
   
Any questions, issues, or feedback can be posted to the DRUGSETS repository's issue tracker.   


## References 
 [^1]: de Leeuw C, Mooij J, Heskes T, Posthuma D (2015): MAGMA: Generalized gene-set analysis of GWAS data. PLoS Comput Biol 11(4): e1004219. doi:10.1371/journal.pcbi.1004219
