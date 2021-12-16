# drugsea
DRUg Gene SEt Analysis (DRUGSEA) is a command line interface (CLI) tool implemented in python to perform genetically informed drug repositioning using drug gene set analysis using MAGMA.    

## Prerequisites    
DRUGSEA was developed using [Python 3.8.5](https://www.python.org/) following packages. Older versions of these packages may work but have not been tested:    
   
[tqdm](https://tqdm.github.io) == 4.62.3   
[numpy](https://www.numpy.org) == 1.21.1    
[pandas](https://pandas.pydata.org) == 1.3.1   
[scipy](https://www.scipy.org) == 1.7.1   
[scikit-learn](http://scikit-learn.org) == 1.0   
[matplotlib](https://matplotlib.org) == 3.4.2   
   
Additionally, the [MAGMA](https://ctg.cncr.nl/software/magma) software tool for gene and gene-set analysis is needed to run DRUGSEA. MAGMA needs to be able to be executed globally, which can be setup by adding the following line to the user's .bashrc or .zshrc file.    
   
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
   
2. Go to directory    
`$ cd drugsea/`   
   
3. Add input data to directory `/DATA/MAGMA_ANNOT/`    
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

`python ./drugsea.py --geneassoc SCZ_SAMPLE.genes.raw --drugsets solo --out SCZ --enrich atc --nsize 5 --showlog no`    
  














