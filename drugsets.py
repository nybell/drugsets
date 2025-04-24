# DRUGSETS v2.0
# Author: Nathan B.
# Performs drug gene-set analysis using GWAS sumstats

# import modules
import os
import argparse
import pandas as pd
import logging
from scipy import stats
import drugsets_func as df_func

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# drug set analysis class
class DrugSetAnalysis:
    def __init__(self, raw_file, drugsets, magma_path, out, conditional, setsize, gene_id, enrich, correct_cov, nsize, correct, use_pops=None, pops_mean_cov=True):
        self.raw_file = raw_file
        self.drugsets = drugsets
        self.magma_path = magma_path
        self.out = out
        self.out_dir = os.path.dirname(out)
        self.out_name = os.path.basename(out)
        self.conditional = conditional
        self.setsize = setsize
        self.gene_id = gene_id
        self.enrich = enrich
        self.correct_cov = correct_cov
        self.nsize = nsize
        self.correct = correct
        self.use_pops = use_pops
        self.pops_mean_cov = pops_mean_cov
        self.filepaths = self.set_filepaths()

    def set_filepaths(self):
        cond_prefixes = {'yes': '_cond', 'no': ''}
        gene_id_conditions = ['entrez', 'ensembl', 'ensembl92']
        filepaths = {}
        for gene_id in gene_id_conditions:
            filepaths[gene_id] = {
                'solo': df_func.set_paths(f"solo_{gene_id}", f"{cond_prefixes[self.conditional]}", self.setsize, GENESETDIR),
                'atc': df_func.set_paths(f"atc_{gene_id}", f"{cond_prefixes[self.conditional]}", self.setsize, GENESETDIR),
                'moa': df_func.set_paths(f"moa_{gene_id}", f"{cond_prefixes[self.conditional]}", self.setsize, GENESETDIR),
                'ind': df_func.set_paths(f"ind_{gene_id}", f"{cond_prefixes[self.conditional]}", self.setsize, GENESETDIR)
            }
        return filepaths

    # TO DO: build in gene ID mapping to the .preds file (i.e., replace ENSGID w/ entrez depending on the .raw input file)
    def normalize_pops_and_adjust_genes_out(self, magma_genes_out, pops_preds):
        try:
            pops_df = pd.read_csv(pops_preds, sep='\t')
            pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Score'])
            pops_df['pvalue'] = 2 * stats.norm.sf(abs(pops_df['normed_pops']))
            magma_df = pd.read_csv(magma_genes_out, delimiter= "\s+", comment='#')
            merged_df = magma_df.merge(pops_df[['ENSGID', 'normed_pops', 'pvalue']], left_on='GENE', right_on='ENSGID', how='inner')
            merged_df.drop(columns=['ENSGID', 'ZSTAT', 'P'], inplace=True)
            merged_df.rename(columns={'normed_pops': 'ZSTAT', 'pvalue': 'P'}, inplace=True)
            adjusted_genes_out = os.path.join(self.out_dir, "adjusted.genes.out")
            merged_df.to_csv(adjusted_genes_out, sep='\t', index=False)
            return adjusted_genes_out
        except Exception as e:
            logging.error(f"Error processing PoPs data: {e}")
            raise
            
    def run_analysis(self):
        selected_filepath = self.filepaths[self.gene_id][self.drugsets]
        # Set output labels for each chosen analysis
        if self.drugsets == 'all':
            analysis_labels = ['_SOLO', '_ATC', '_MOA', '_IND']
        else:
            analysis_labels = ['_' + self.drugsets.upper()]
        # log labels
        logging.info(f"Defined analysis label(s): {', '.join(analysis_labels)}")
        # run
        try:
            logging.info("Starting drug gene set analysis...")
            for label in analysis_labels:
                if self.use_pops:
                    magma_genes_out, pops_preds = self.use_pops
                    adjusted_genes_out = self.normalize_pops_and_adjust_genes_out(magma_genes_out, pops_preds)
                    if self.pops_mean_cov:
                        cmd = f"{self.magma_path} --gene-results {self.raw_file} overwrite={adjusted_genes_out} --set-annot {selected_filepath} --gene-covar {pops_cov} --model condition-hide=mean_pops_norm --settings gene-info --out {self.out + label}"
                    else:
                        cmd = f"{self.magma_path} --gene-results {self.raw_file} overwrite={adjusted_genes_out} --set-annot {selected_filepath} --settings gene-info --out {self.out + label}"
                else:
                    cmd = f"{self.magma_path} --gene-results {self.raw_file} --set-annot {selected_filepath} --out {self.out + label}"
                logging.info(f"Running MAGMA command with label {label}: {cmd}")
                os.system(cmd)
            logging.info("Drug gene set analysis completed successfully.")
        except Exception as e:
            logging.error(f"Analysis failed: {e}")

    def enrichment(self):
        selected_filepath = self.filepaths[self.gene_id][self.drugsets]
        try:
            if self.enrich:
                if (self.drugsets == 'solo' or self.drugsets == 'all'):
                    corrdata, metaRdata = df_func.compute_cor_matrix(self.raw_file, selected_filepath, self.out, self.out_name, self.out_dir, self.correct_cov)
                    df_func.test_enrich(corrdata, metaRdata, self.out_name, self.out_dir, self.nsize, self.correct, self.enrich)
                else:
                    logging.error('To test for enrichment "--drugsets" must be set to "solo".')
            else: 
                logging.error('--enrich not specified.')
    
        except Exception as e:
            logging.error(f"Drug enrichment failed: {e}")

# check if running directly
if __name__ == "__main__":
    
    # parse input args
    parser = argparse.ArgumentParser(description="Drug Gene Set Analysis using MAGMA.")
    parser.add_argument("--raw_file", required=False, 
                        help="Full path to MAGMA .raw file")
    parser.add_argument("--drugsets", required=False, choices=['solo', 'atc', 'moa', 'ind', 'all'],
                        help="Type of drug gene set to use (individual drug, ATC code, mechanism of action, clinical indication)")
    parser.add_argument("--magma", required=False, default="magma", 
                        help="Path to MAGMA executable")
    parser.add_argument("--out", required=True, 
                        help="Full path to output directory + file name prefix")
    parser.add_argument("--conditional", required=False, default='yes', choices=['yes','no'], 
                        help="Condition on all all druggable genes (recommended)?")
    parser.add_argument("--enrich", required=False, default=None, choices=[None, 'atc', 'moa', 'ind', 'all'], 
                        help="Specify drug group type to test for enrichment")
    parser.add_argument("--correct_cov", required=False, default="TRUE", choices=["TRUE", "FALSE"], 
                        help="Correct for drug gene set covariance during enrichment testing")
    parser.add_argument("--setsize", type=int, default=5, 
                        help="Minimum gene set size")
    parser.add_argument("--gene_id", required=True, choices=['entrez','ensembl','ensembl92'], 
                        help="Specify gene ID type used in MAGMA .raw file (for FUMA use 'ensembl92')")
    parser.add_argument("--nsize", type=int, default=5, 
                        help="Minimum drug group sample size for enrichment testing")
    parser.add_argument("--correct", default='bonf', choices=['bonf', 'fdr'], 
                        help="Multiple testing correction")
    parser.add_argument("--use_pops", nargs=2, metavar=('MAGMA_GENES_OUT', 'POPS_PREDS'), 
                        help="To use PoPs scores, please supply file paths to MAGMA .genes.out and PoPs .preds files")
    parser.add_argument("--pops_mean_cov", type = bool, choices=[True, False], default=True,
                        help="Condition on normalized PoPS scores across 111 UKB traits (default: True)")
    parser.add_argument("--pops_cov_file", type = str, default='ensembl_mean_pops.cov',
                        help="File with with PoPS score to condition on.")

    # parse input arguments 
    args = parser.parse_args()

    # print welcome and perform input checks
    df_func.welcome_check_inputs(args)

    # Setting directories
    DIR = os.path.dirname(os.path.abspath(__file__))
    DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
    GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))
    ANNOTDIR = os.path.normpath(os.path.join(DATADIR, 'MAGMA_ANNOT'))
    COVDIR = os.path.normpath(os.path.join(DATADIR, 'COV'))
    # OUT = os.path.normpath(args.out)

    # Define PoPS covariate file
    pops_cov = os.path.normpath(os.path.join(COVDIR, args.pops_cov_file))

    # set up analysis class
    analysis = DrugSetAnalysis(
        args.raw_file, 
        args.drugsets, 
        args.magma, 
        args.out, 
        args.conditional, 
        args.setsize, 
        args.gene_id, 
        args.enrich, 
        args.correct_cov, 
        args.nsize, 
        args.correct, 
        args.use_pops,
        args.pops_mean_cov
    )

    # run drug gene set analysis
    analysis.run_analysis()

    # run enrichment (if specified)
    analysis.enrichment()

    # print 
    df_func.print_ending_message()

    


