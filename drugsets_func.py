import pandas as pd
import logging
import os
import scipy.stats as stats

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# function for welcome message and printing input arguments
def welcome_check_inputs(args):
    welcome_text = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃      🚀 Welcome to DRUGSETS v2.0!     ┃
    ┃      Identify GWAS-enriched drugs     ┃
    ┃       with drug gene-set analysis     ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    # Print welcome message
    print(welcome_text)
    print("🔍 Reading input...\n")
    
    # Print input arguments
    print("📌 Input arguments used:")
    for arg, value in vars(args).items():
        print(f"   🔹 {arg} = {value}")
    
    # Input validation checks
    errors = []

    # Check #1: Validate raw_file extension
    if not args.raw_file.endswith('.genes.raw'):
        errors.append("❌ ERROR: Gene association file does not end in '.genes.raw'. Please check the MAGMA gene association input file and try again.")
    
    # Check #2: Validate gene_id when use_pops is enabled
    if args.use_pops is not None and args.gene_id not in ["ensembl", "ensembl92"]:
        errors.append("⚠️ WARNING: To use PoPs results, please use '.raw' and '.genes.out' files with Ensembl IDs ('ensembl' or 'ensembl92').")

    # Print error messages if any checks failed
    if errors:
        print("\n🚨 Input Validation Issues:")
        for error in errors:
            print(error)
        print("\n🔴 Exiting due to input errors.\n")
        quit()

    # Print initialization message if all checks pass
    print("\n✅ Initialization complete. Starting analysis. 💊💡\n")


# function for normalizing PoPs scores and merging with MAGMA genes.out
def normalize_pops_and_adjust_genes_out(magma_genes_out, pops_preds, output_dir):
    """Normalize PoPs scores, compute p-values, and merge with MAGMA genes.out."""
    try:
        pops_df = pd.read_csv(pops_preds, sep='\t')
        logging.info(f"Loaded PoPs results from {pops_preds}. Shape: {pops_df.shape}")

        pops_df['normed_pops'] = stats.zscore(pops_df['PoPS_Score'])
        pops_df['pvalue'] = 2 * stats.norm.sf(abs(pops_df['normed_pops']))

        magma_df = pd.read_csv(magma_genes_out, sep='\t')
        logging.info(f"Loaded MAGMA genes.out from {magma_genes_out}. Shape: {magma_df.shape}")

        merged_df = magma_df.merge(pops_df[['ENSGID', 'normed_pops', 'pvalue']], 
                                   left_on='GENE', right_on='ENSGID', how='inner')

        merged_df.drop(columns=['ZSTAT', 'P'], inplace=True)
        merged_df.rename(columns={'normed_pops': 'ZSTAT', 'pvalue': 'P'}, inplace=True)

        adjusted_genes_out = os.path.join(output_dir, "pops.adjusted.genes.out")
        merged_df.to_csv(adjusted_genes_out, sep='\t', index=False)
        logging.info(f"Saved adjusted MAGMA genes.out file: {adjusted_genes_out}")

        return adjusted_genes_out

    except Exception as e:
        logging.error(f"Error processing PoPs data: {e}")
        raise

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
def compute_cor_matrix(raw_file, solo, output, out_name, out_dir, correct_cov):
    working = os.path.dirname(os.path.abspath(__file__))
    # working = "/Users/nyb/drugsets"

    gsa = None
    for ext in ['.gsa.out', '.gsa.out.txt']:
        path = os.path.normpath(f"{output}_SOLO{ext}")
        if os.path.exists(path):
            gsa = path
            print(gsa)
            break
    if not gsa:
        raise FileNotFoundError("GSA output file not found.")

    logging.info("Computing correlation matrix...")
    cmd_corr = f"Rscript --vanilla {os.path.join(working, 'compute_corrs.R')} {raw_file} {solo} {gsa} {out_name} {out_dir} {correct_cov}"
    logging.info(f"Running command: {cmd_corr}")
    os.system(cmd_corr)

    corrdata = os.path.normpath(os.path.join(out_dir, f"{out_name}_setcorrs.rdata"))
    metaRdata = os.path.normpath(os.path.join(working, 'DATA', 'metadata.rdata'))

    return corrdata, metaRdata

# function for running linear regression enrichment analysis
def test_enrich(corrdata, metaRdata, out_name, out_dir, nsize, correct, enrich_type):
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    def run_lnreg(group):
        logging.info(f"Running linear regression for {group}...")
        cmd_lnreg = f"Rscript --vanilla {os.path.join(BASE_DIR, 'compute_lnreg.R')} {corrdata} {metaRdata} {group} {nsize} {out_name}_{group.upper()} {out_dir} {correct}"
        logging.info(f"Running command: {cmd_lnreg}")
        os.system(cmd_lnreg)

    if enrich_type != 'all':
        run_lnreg(enrich_type)
    else:
        for group in ['atc', 'moa', 'ind']:
            run_lnreg(group)

    # os.remove(corrdata)
    logging.info("Enrichment analysis finished.")

def print_ending_message():
    ending_text = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃   🎉 Analysis Complete! 🚀            ┃
    ┃   Your drug enrichment results are    ┃
    ┃   ready for interpretation. 💊📊      ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    
    print(ending_text)
    print("\n🔥 Thank you for using DRUGSETS v2.0! 🔥\n")