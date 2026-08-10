# import packages
import subprocess
import pandas as pd
import logging
import os
import scipy.stats as stats

# setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# function for welcome message and printing input arguments
def welcome_check_inputs(args):
    welcome_text = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃      » Welcome to DRUGSETS v2.1!      ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    # Print welcome message
    print(welcome_text)
    print("🔍 Reading input...\n")
    
    # Print input arguments
    print("📌 Input arguments used:")
    for arg, value in vars(args).items():
        print(f"   » {arg} = {value}")
    
    # Input validation checks
    errors = []

    # Check #1: Validate raw_file extension
    if args.raw_file is not None:
        if not args.raw_file.endswith('.genes.raw'):
            errors.append("[!] ERROR: Gene association file does not end in '.genes.raw'. Please check the MAGMA gene association input file and try again.")

    # Check #2: Validate gene_id when use_pops is enabled
    if args.use_pops is not None and args.gene_id not in ["ensembl", "ensembl92"]:
        errors.append("[!] WARNING: To use PoPs results, please use '.raw' and '.genes.out' files with Ensembl IDs ('ensembl' or 'ensembl92').")

    # Check #3: Check if drugsets is None
    if args.drugsets is None and args.enrich_file is None:
        errors.append("[!] WARNING: No drug gene set analysis or enrichment testing specified. Please specify --drugsets and/or --enrich_file to run the analysis.")
    
    # Check #4: Check if both drugsets and enrich_file are specified
    if args.drugsets is not None and args.enrich_file is not None:
        errors.append("[!] WARNING: Both --drugsets and --enrich_file specified. Please specify only one of these options to avoid confusion between drug gene set analysis and enrichment testing.")

    # Print error messages if any checks failed
    if errors:
        print("\nInput Validation Issues:")
        for error in errors:
            print(error)
        print("\nExiting due to input errors.\n")
        quit()

    # Print initialization message if all checks pass
    print("\nInitialization complete. Starting analysis.\n")


# # function for normalizing PoPs scores and merging with MAGMA genes.out
# def normalize_pops_and_adjust_genes_out(magma_genes_out, pops_preds, output_dir):
#     """Normalize PoPs scores, compute p-values, and merge with MAGMA genes.out."""
#     try:
#         pops_df = pd.read_csv(pops_preds, sep='\t')
#         logging.info(f"Loaded PoPs results from {pops_preds}. Shape: {pops_df.shape}")

#         pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Score'])
#         pops_df['pvalue'] = 2 * stats.norm.sf(abs(pops_df['normed_pops']))

#         magma_df = pd.read_csv(magma_genes_out, sep='\t')
#         logging.info(f"Loaded MAGMA genes.out from {magma_genes_out}. Shape: {magma_df.shape}")

#         merged_df = magma_df.merge(pops_df[['ENSGID', 'normed_pops', 'pvalue']], 
#                                    left_on='GENE', right_on='ENSGID', how='inner')

#         merged_df.drop(columns=['ZSTAT', 'P'], inplace=True)
#         merged_df.rename(columns={'normed_pops': 'ZSTAT', 'pvalue': 'P'}, inplace=True)

#         adjusted_genes_out = os.path.join(output_dir, "pops.adjusted.genes.out")
#         merged_df.to_csv(adjusted_genes_out, sep='\t', index=False)
#         logging.info(f"Saved adjusted MAGMA genes.out file: {adjusted_genes_out}")

#         return adjusted_genes_out

#     except Exception as e:
#         logging.error(f"Error processing PoPs data: {e}")
#         raise

def get_magma_covars(magma_cov):
    if len(magma_cov) < 2:
        raise ValueError("You must provide a covariate file and at least one covariate name.")

    cov_file = magma_cov[0]
    covariates = magma_cov[1:]
    return cov_file, covariates

# define function to create new gene set file with custom size
def setsize(path, file, size):
    
    # create name for new gene set file 
    new = file.replace('.txt', '_min'+str(size)+'.txt')
    new = "/tmp"+new

    # add path
    new = path+new

    # create file 
    with open(os.path.normpath(path+file)) as oldfile, open(os.path.normpath(new), 'w') as newfile:
        for line in oldfile:
                if len(line.split('\t')) -3 >= int(size):
                    newfile.write(line)

# Function to set gene set paths
def set_paths(prefix, cond, set_size, GENESETDIR):
    base = os.path.normpath(os.path.join(GENESETDIR, f"{prefix}{cond}_sets.txt"))
    
    # For default set size
    if set_size == 2:
        return base

    # For custom minimum gene set size
    setsize(GENESETDIR, f"/{prefix}{cond}_sets.txt", set_size)
    return os.path.normpath(os.path.join(GENESETDIR, f"tmp/{prefix}{cond}_sets_min{set_size}.txt"))

# function for computing correlation matrix
def compute_cor_matrix(raw_file, gene_sets, output, out_name, out_dir, correct_cov):
    working = os.path.dirname(os.path.abspath(__file__))
    # working = "/Users/nyb/drugsets"

    # check gene sets file path
    if not os.path.isfile(gene_sets):
        raise FileNotFoundError(f"[!] ERROR: Gene set file not found: {gene_sets}")
    
    # resolve gsa file path
    if output.endswith(('.gsa.out', '.gsa.out.txt')):
        if not os.path.isfile(output):
            raise FileNotFoundError(f"[!] ERROR: Custom MAGMA GSA file not found: {output}")
        gsa = output
    else:
        # search for MAGMA output using prefix
        extensions = ['.gsa.out', '.gsa.out.txt']
        possible_paths = [f"{output}.drug_mp0{ext}" for ext in extensions] + \
                         [f"{output}.drug_mp1{ext}" for ext in extensions] + \
                         [f"{output}.drug_mp2{ext}" for ext in extensions] + \
                         [f"{output}.drug_mp3{ext}" for ext in extensions] + \
                         [f"{output}.drug_mp4{ext}" for ext in extensions] + \
                         [f"{output}{ext}" for ext in extensions]
        
        gsa = next((os.path.normpath(p) for p in possible_paths if os.path.isfile(p)), None)

    if not gsa:
        raise FileNotFoundError(f"[!] ERROR: GSA output file not found for prefix: {output}")

    logging.debug(f"[+] DEBUG: Using GSA file: {gsa}")
    logging.info("[+] INFO: Computing correlation matrix...")

    # run R script to compute correlation matrix
    cmd_corr = f"Rscript --vanilla {os.path.join(working, 'compute_corrs.R')} {raw_file} {gene_sets} {gsa} {out_name} {out_dir} {correct_cov}"
    logging.info(cmd_corr)
    try: 
        subprocess.run(cmd_corr, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"[+] ERROR: Command failed with exit code {e.returncode}")
        logging.error("[+] ERROR: Failing command:")
        logging.error(cmd_corr)
        raise

    corrdata = os.path.normpath(os.path.join(out_dir, f"{out_name}_setcorrs.rdata"))
    metaRdata = os.path.normpath(os.path.join(working, 'DATA', 'meta.data.2026.RData'))

    return corrdata, metaRdata

# function for running linear regression enrichment analysis
def test_enrich(corrdata, metaRdata, out_name, out_dir, nsize, correct, enrich_type):
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    
    def run_lnreg(group):
        logging.info(f"[+] INFO: Running regression for {group}...")
        cmd_lnreg = f"Rscript --vanilla {os.path.join(BASE_DIR, 'compute_lnreg.R')} {corrdata} {metaRdata} {group} {nsize} {out_name}.{group.lower()} {out_dir} {correct}"
        
        # subprocess.run with check=True will raise an exception if the R script fails
        subprocess.run(cmd_lnreg, shell=True, check=True)

    try:
        if enrich_type != 'all':
            logging.info(f"[+] INFO: Testing enrichment for {enrich_type}...")
            run_lnreg(enrich_type)
        else:
            for group in ['atc1', 'atc2', 'atc3', 'atc4', 'moa_group', 'moa_nested', 'moa_specific']:
                logging.info(f"[+] INFO: Testing enrichment for {group}...")
                run_lnreg(group)

        logging.info("[+] INFO: Enrichment analysis finished.")
        
    except subprocess.CalledProcessError as e:
        # e.cmd contains the exact string that was passed to subprocess.run
        logging.error(f"[+] ERROR: Command failed with exit code {e.returncode}")
        logging.error(f"[+] ERROR: Failing command: {e.cmd}")
        raise
    except Exception as e:
        # Catches standard Python errors (e.g., typos, missing variables)
        logging.error(f"[+] ERROR: Analysis failed due to an unexpected error: {e}")
        raise

# print gene set end
def print_gsa_message():
    ending_text = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃     » Gene set analysis complete!     ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    print(ending_text)

def print_group_message():
    ending_text = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃    » Drug group analysis complete!    ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    print(ending_text)

def print_ending_message():
   
    print("Thank you for using DRUGSETS v2.1\n")