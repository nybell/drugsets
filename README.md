# DRUGSEA
DRUg Gene SEt Analysis (DRUGSEA) is a command line interface (CLI) tool implemented in python to perform genetically informed drug repositioning using drug gene set analysis using MAGMA.    

## Prerequisites    
DRUGSEA was developed using [Python 3.8.5](https://www.python.org/) following packages. Older versions of these packages may work but have not been tested:    
   
[tqdm](https://tqdm.github.io) == 4.62.3   
[numpy](https://www.numpy.org) == 1.21.1    
[pandas](https://pandas.pydata.org) == 1.3.1   
[scipy](https://www.scipy.org) == 1.7.1   
[scikit-learn](http://scikit-learn.org) == 1.0   
[matplotlib](https://matplotlib.org) == 3.4.2   
   
Additionally, the [MAGMA](https://ctg.cncr.nl/software/magma)[^1] software tool for gene and gene-set analysis is needed to run DRUGSEA. Full documentation of MAGMA can be found [here](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf). MAGMA needs to be able to be executed globally, which can be setup by adding the following line to the end of the user's .bashrc or .zshrc file.    
   
`export PATH=$PATH:/INSERT/PATH/TO/FILE/magma_v1.09b_mac/`   
   
## Setup   
   
If [Python 3.8.5](https://www.python.org/) is installed, all prerequisite Python packages can be installed using the following code:    
   
```
pip install --upgrade pip
pip install tqdm == 4.62.3
pip install numpy == 1.21.1
pip install pandas == 1.3.1
pip install scipy == 1.7.1
pip install scikit-learn == 1.0
pip install matplotlib == 3.4.2
```
   
Users can also create a virutal environment in which to run DRUGSEA using [Conda](https://www.anaconda.com/products/individual) or [pyenv](https://github.com/pyenv/pyenv). A guide to installing pyenv can be found [here](https://github.com/pyenv/pyenv), and a guide to installing Conda can be found [here](https://docs.anaconda.com/anaconda/navigator/tutorials/index.html).    
   
A virtual environment with all prerequisites can be created using pyenv via:    
    
```
pyenv virtualenv 3.8.5 venv_drugsea
pip install --upgrade pip
pip install tqdm == 4.62.3
pip install numpy == 1.21.1
pip install pandas == 1.3.1
pip install scipy == 1.7.1
pip install scikit-learn == 1.0
pip install matplotlib == 3.4.2
pyenv activate venv_drugsea
```   
A virtual environment with all prerequisites can be created using Conda via:    
   
```   
conda create --name venv_drugsea python=3.8.5 tqdm=4.62.3 numpy=1.21.1 pandas=1.3.1 scipy=1.7.1 scikit-learn=1.0 matplotlib=3.4.2
conda activate venv_drugsea
```   
   
## Setup & Input data     
   
1. Clone repository    
`$ git clone https://github.com/nybell/drugsea`   
   
2. Add input data to directory `/DATA/MAGMA_ANNOT/`    
The only input data needed by the user is the [GENE_RESULTS].genes.raw file from MAGMA gene analysis. ***In order for DRUGSEA to run, users must put their [GENE_RESULTS].genes.raw into the `/DATA/MAGMA_ANNOT/` directory.*** All drug data used in the drug gene set analysis is included in the download, and is stored in the directory `/DATA/GENESETS/`.    
   
## Running drug gene set analysis   
    
Drug gene set analysis is done by executing the script `drugsea.py`. The only flags required for running drug gene set analysis are `--geneassoc`, `--drugsets`, and `--out`. Enrichment analyses can be added using the optional flags `--enrich` and `--nsize`.     

* `--geneassoc` or `-g` specifies the filepath and filename of the user's [GENE_RESULTS].genes.raw file (i.e., `/INSERT/PATH/TO/[GENE_RESULTS].genes.raw`)
* `--drugsets` or `-d` specifies which type of drug genesets to use. There are four options:
    * `solo` drug gene sets for each individual drug in our data. **NOTE:** if you want to test for enrichment you much use `solo`
    * `atc` drug genesets for each ATC III code category in our data 
    * `moa` drug genesets for each mechanism of action category in our data 
    * `ind` drug gene sets for each clinical indication category in our data 
* `--out` or `-o` specific prefix for output files
* `--enrich` or `-e` test a type of category for enrichment of the drugs with the strongest associations to the phenotype measured using drug gene set analysis. There are three input options:
    * `atc` test for enrichment by ATC code 
    * `moa` test for enrichment by mechanism of action 
    * `ind` test for enrichment by clinical indication 
* `--nsize` or `-n` minimum sample size of drug categories to use when testing for enrichment (e.g., when set to 5, ATC code categories with less than 5 drugs will not be tested for enrichment). 
* `--showlog` or `-l` specifies whether or not print the output from MAGMA gene-set analysis to screen. Options are `no` and `yes` (default = `no`) 
    
### Example usage    
    
The following code tests individual drug gene set analysis for associated with schizophrenia. The SCZ_SAMPLE.genes.raw was created using MAGMA and the 2021 schizophrenia GWAS summary statistics, downloaded [here](https://www.med.unc.edu/pgc/download-results/). Additionally, the code tests for ATC III codes for enrichment of drugs that are highly associated with schizophrenia, with a minimum sample size of 5 drugs for each ATC code. Lastly, it specifies to not show output from MAGMA (to see this go to the .log file in the directory `/drugsea/OUTPUT/`.       
  
1. Go to directory      
`$ cd drugsea/`     
  
2. Execute script    
`python ./drugsea.py --geneassoc SCZ_SAMPLE.genes.raw --drugsets solo --out SCZ --enrich atc --nsize 5 --showlog no`    
    
## Output    
   
The formats of the output files from the drug gene set analysis in MAGMA (.gsa.out, .gsa.genes.out .gsa.sets.genes.out) is described [here](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.09.pdf) on page 23.    
    
##### Drug gene set analysis
* `.gsa.genes.out` contains information for each gene used in the MAGMA analysis   
* `.gsa.out` contains information for each drug gene set used in the MAGMA analysis    
* `.gsa.sets.genes.out` contains information on significant drug gene sets after Bonferroni correction, and the genes in each gene set  
    
##### Enrichment analysis 
* `enrich.OUT.txt` contains the results from the enrichment analysis for every drug category tested
* `enrich.bonf.OUT.txt` contains the results from the enrichment analysis for drug significant after Bonferroni correction   
  
Both enrichment output files contain the same columns:
* `GROUP` the drug category tested 
* `MWU` Wilcoxon Mann Whitney U statistic from a one-tailed test whether the p-values in the group are lower then the p-values for drugs not in the group
* `P` p-value from the Wilcoxon Mann Whitney U test
* `AUC` Area under the enrichment curve for drugs in the group tested for enrichment, as described [here](https://www.nature.com/articles/s41598-017-12325-3)


# Support   
   
Any questions, issues, or feedback can be posted to the DRUGSEA repository's issue tracker.   


## References 
 [^1]: de Leeuw C, Mooij J, Heskes T, Posthuma D (2015): MAGMA: Generalized gene-set analysis of GWAS data. PLoS Comput Biol 11(4): e1004219. doi:10.1371/journal.pcbi.1004219
