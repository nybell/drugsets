# DRUGSETS v2.1
# Author: Nathan B.
# Performs drug gene-set analysis using GWAS sumstats

# import modules
import os
import subprocess
import argparse
import pandas as pd
import logging
from scipy import stats
import drugsets_func as df_func

# setup logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S", 
    force = True
)

# drug set analysis class
class DrugSetAnalysis:
    def __init__(self, raw_file, drugsets, magma_path, out, conditional, setsize, gene_id,
                 enrich, enrich_file, correct_cov, nsize, correct, magma_cov=None, use_pops=None, adj_pops=True, dgidb_filter=None):
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
        self.enrich_file = enrich_file
        self.correct_cov = correct_cov
        self.nsize = nsize
        self.correct = correct
        self.magma_cov = magma_cov
        self.use_pops = use_pops
        self.adj_pops = adj_pops
        self.dgidb_filter = dgidb_filter
        self.filepaths = self.set_filepaths()

    def set_filepaths(self):
        cond_prefixes = {'yes': '_cond', 'no': ''}
        gene_id_conditions = ['ensembl'] # 'entrez', 'ensembl', 'ensembl92'
        filepaths = {}
        # for dgidb comp of name
        dgidb_str = f"_dgidb_{self.dgidb_filter.lower()}" if self.dgidb_filter else "_dgidb_median"
        # loop through gene ids 
        for gene_id in gene_id_conditions:
            paths = {}
            # generate drug_mp0 to drug_mp4 paths
            for i in range(5):
                key = f'drug_mp{i}'
                prefix = f"drug_mp{i}{dgidb_str}_{gene_id}"
                paths[key] = df_func.set_paths(prefix, cond_prefixes[self.conditional], self.setsize, GENESETDIR)
            # moa sets
            paths['moa_group'] = df_func.set_paths(f"moa_targ_chembl_group_{gene_id}", cond_prefixes[self.conditional], self.setsize, GENESETDIR)
            paths['moa_nested'] = df_func.set_paths(f"moa_targ_chembl_nested_{gene_id}", cond_prefixes[self.conditional], self.setsize, GENESETDIR)
            paths['moa_specific'] = df_func.set_paths(f"moa_targ_chembl_specific_{gene_id}", cond_prefixes[self.conditional], self.setsize, GENESETDIR)
            
            filepaths[gene_id] = paths
            
        return filepaths
    
    def build_magma_cmd(self, gene_sets, label, adjusted_genes_out=None):
        # base command
        parts = [self.magma_path]

        # gene results (+ PoPS-substituted .genes.out when --use_pops is set)
        if adjusted_genes_out:
            parts.append(f"--gene-results {self.raw_file} overwrite={adjusted_genes_out}")
        else:
            parts.append(f"--gene-results {self.raw_file}")

        # drug gene sets
        parts.append(f"--set-annot {gene_sets}")

        # collect everything to condition on (and hide) into one list
        hide_terms = []

        # optional gene-level covariates (e.g. expression)
        if self.magma_cov:
            cov_file, covars = df_func.get_magma_covars(self.magma_cov)
            covar_str = ",".join(covars)
            parts.append(f"--gene-covar {cov_file} filter-read={covar_str}")
            hide_terms.extend(covars)

        # optional conditioning on the druggable background set
        if self.conditional == 'yes':
            hide_terms.append("druggable")

        # single condition-hide modifier
        if hide_terms:
            parts.append(f"--model condition-hide={','.join(hide_terms)}")

        parts.append("--settings gene-info")
        parts.append(f"--out {self.out + label}")

        return " ".join(parts)

    # TO DO: build in gene ID mapping to the .preds file (i.e., replace ENSGID w/ entrez depending on the .raw input file)
    def normalize_pops_and_adjust_genes_out(self, magma_genes_out, pops_preds, adj_pops):
        try:
            pops_df = pd.read_csv(pops_preds, sep='\t')
            if adj_pops == 'yes':
                # logging.info("Using adjusted PoPS scores.")
                logging.info("[+] INFO: Option currently unavailable. Using unadjusted PoPS scores instead.")
                # pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Z'])
                pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Score'])
            else:
                logging.info("[+] INFO: Using unadjusted PoPS scores.")
                pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Score'])
            # compute p-values & merge with MAGMA genes.out
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
            
    def run_analysis(self, verbose=False):
        # map each requested analysis to its gene-set file key
        if self.drugsets == 'all':
            set_keys = ['drug_mp0', 'drug_mp1', 'drug_mp2', 'drug_mp3', 'drug_mp4', 'moa_group', 'moa_nested', 'moa_specific']
        else:
            set_keys = [self.drugsets]

        logging.info(f"[+] INFO: using parameters: {', '.join('.' + k for k in set_keys)}")

        # PoPS substitution is label-independent: build the adjusted .genes.out once
        adjusted_genes_out = None
        if self.use_pops:
            magma_genes_out, pops_preds = self.use_pops
            adjusted_genes_out = self.normalize_pops_and_adjust_genes_out(
                magma_genes_out, pops_preds, self.adj_pops
            )

        try:
            logging.info("[+] INFO: Defining MAGMA command...")
            for key in set_keys:
                label = '.' + key
                gene_sets = self.filepaths[self.gene_id][key]
                cmd = self.build_magma_cmd(gene_sets, label, adjusted_genes_out)

                logging.info("[+] INFO: Running drug gene-set analysis with MAGMA")
                if verbose:
                    logging.info("[+] INFO: MAGMA Command")
                    logging.info(cmd)
                    subprocess.run(cmd, shell=True, check=True)
                else:
                    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

            logging.info("[+] INFO: Drug gene set analysis completed.")
        except Exception as e:
            logging.error(f"[+] ERROR: Analysis failed: {e}. Check MAGMA log file. Exiting...")
            raise

    def enrichment(self):
        # check 1: if enrich is specified
        if not self.enrich:
            logging.info('[+] INFO: --enrich not specified. Skipping enrichment.')
            return
            
        # check 2: is drugsets set to drug, all, or none
        if self.drugsets not in ('drug_mp0', 'drug_mp1', 'drug_mp2', 'drug_mp3', 'drug_mp4', 'moa_group', 'moa_nested', 'moa_specific', 'all', None, 'None'):
            logging.warning(f'[+] WARNING: Enrichment testing requires drug-level gsa results. Skipping enrichment since --drugsets is set to "{self.drugsets}".')
            return
            
        # check 3: if drug gsa not run, check if enrich_file specified
        if self.drugsets in (None, 'None'):
            if self.enrich_file:
                gsa_out = self.enrich_file
                # out_name = os.path.basename(gsa_out).replace('.gsa.out.txt', '').replace('.gsa.out', '').replace('.drug', '')
                # out_dir = os.path.dirname(gsa_out)
                out_name = self.out_name
                out_dir = self.out_dir
            else:
                logging.warning('[+] WARNING: --enrich_file not specified. Skipping enrichment testing.')
                return
        else:
            gsa_out = self.out
            out_name = self.out_name
            out_dir = self.out_dir

        # set path to gene sets
        gene_sets = self.filepaths[self.gene_id]["drug_mp4"]

        try:
            # execution block (redundant logic removed)
            corrdata, metaRdata = df_func.compute_cor_matrix(self.raw_file, gene_sets, gsa_out, out_name, out_dir, self.correct_cov)
            df_func.test_enrich(corrdata, metaRdata, out_name, out_dir, self.nsize, self.correct, self.enrich)

        except Exception as e:
            logging.error(f"[+] ERROR: Drug enrichment failed: {e}")
            raise

# check if running directly
if __name__ == "__main__":
    
    # parse input args
    parser = argparse.ArgumentParser(description="Drug Gene Set Analysis using MAGMA.")
    parser.add_argument("--raw_file", required=True, 
                        help="Full path to MAGMA .raw file")
    parser.add_argument("--drugsets", required=False, default = None, choices=['drug_mp0', 'drug_mp1', 'drug_mp2', 'drug_mp3', 'drug_mp4', 
                                                                               'moa_group', 'moa_nested', 'moa_specific', 'all'], 
                        help="Type of drug gene set to use (individual drug, ATC code, mechanism of action, clinical indication)")
    parser.add_argument("--dgidb_filter", required=False, default="median", choices=["all", "median", "p75", "none"],
                        help="Filter for DGIdb interaction score to include genes in drug gene sets (default: median)")
    parser.add_argument("--magma", required=False, default="magma", 
                        help="Path to MAGMA executable")
    parser.add_argument("--out", required=True, 
                        help="Full path to output directory + file name prefix")
    parser.add_argument("--conditional", required=False, default='yes', choices=['yes','no'], 
                        help="Condition on all all druggable genes (recommended)?")
    parser.add_argument("--enrich", required=False, default=None, choices=[None, 'atc1', 'atc2', 'atc3', 'atc4'], 
                        help="Specify drug group type to test for enrichment")
    parser.add_argument("--enrich_file", required=False, default=None,
                        help="Full path to gsa.out file to use for enrichment analysis (if different from the one generated in this pipeline)")
    parser.add_argument("--correct_cov", required=False, default="TRUE", choices=["TRUE", "FALSE"], 
                        help="Correct for drug gene set covariance during enrichment testing")
    parser.add_argument("--setsize", type=int, default=2, 
                        help="Minimum gene set size")
    parser.add_argument("--gene_id", required=True, choices=['entrez','ensembl','ensembl92'], 
                        help="Specify gene ID type used in MAGMA .raw file (for FUMA use 'ensembl92')")
    parser.add_argument("--nsize", type=int, default=5, 
                        help="Minimum drug group sample size for enrichment testing")
    parser.add_argument("--correct", default='bonf', choices=['bonf', 'fdr'], 
                        help="Multiple testing correction")
    parser.add_argument("--use_pops", nargs=2, metavar=('MAGMA_GENES_OUT', 'POPS_PREDS'), 
                        help="To use PoPs scores, please supply file paths to MAGMA .genes.out and PoPs .preds files")
    parser.add_argument("--adj_pops", required=False, default='yes', choices=['yes','no'], 
                        help="Are you using corrected PoPS scores?")
    parser.add_argument("--magma_cov", nargs="+", metavar=("COV_FILE", "COVARIATES"),
                        help="First provide the covariate file, then one or more covariate names (i.e., column names) for MAGMA.")

    # parse input arguments 
    args = parser.parse_args()

    # print welcome and perform input checks
    df_func.welcome_check_inputs(args)

    # setting directories
    DIR = os.path.dirname(os.path.abspath(__file__))
    DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
    GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))

    # make gene sets here

    try: 
        # set up analysis class
        analysis = DrugSetAnalysis(
            raw_file=args.raw_file, 
            drugsets=args.drugsets,
            dgidb_filter=args.dgidb_filter,
            magma_path=args.magma, 
            out=args.out, 
            conditional=args.conditional, 
            setsize=args.setsize, 
            gene_id=args.gene_id, 
            enrich=args.enrich, 
            enrich_file=args.enrich_file,
            correct_cov=args.correct_cov, 
            nsize=args.nsize, 
            correct=args.correct, 
            use_pops=args.use_pops,
            adj_pops=args.adj_pops,
            magma_cov=args.magma_cov
        )

        # run drug gene set analysis
        if args.drugsets:
            analysis.run_analysis()
            df_func.print_gsa_message()
        
        if args.enrich:
            # run enrichment (if specified)
            analysis.enrichment()
            df_func.print_group_message()

        # print 
        df_func.print_ending_message()
    except Exception as e:
        logging.error("[!] Analysis failed. Exiting...")
        logging.error(f"Error: {e}")
        raise
    


