import pytest
import itertools
from your_script import DrugSetAnalysis  # Import your class

# Define all possible parameter values
geneassoc = ["/Users/nyb/drugsets/DATA/MAGMA_ANNOT/SCZ_SAMPLE.genes.raw"]
drugsets = ["solo", "atc", "moa", "ind"]
magma_path = ["magma"]
out = ["/Users/nyb/Desktop/pytest"]
conditional = ["yes", "no"]
setsize = [2, 5]
gene_id = ["entrez", "ensembl", "ensembl92"]
enrich = [None, "atc", "moa", "ind"]
enrich_auc = [None, "atc", "moa", "ind"]
nsize = [5, 10]
correct = [None, "bonf", "fdr"]
use_pops = [None, ("/Users/nyb/drugsets/DATA/MAGMA_ANNOT/SCZ_SAMPLE.hg37.step2.genes.out", "/Users/nyb/drugsets/DATA/MAGMA_ANNOT/SCZ_SAMPLE.hg37.preds")]

# Generate all combinations of parameters
all_combinations = list(itertools.product(
    geneassoc, drugsets, magma_path, out, conditional, setsize,
    gene_id, enrich, enrich_auc, nsize, correct, use_pops
))

@pytest.mark.parametrize("params", all_combinations)
def test_all_combinations(params):
    geneassoc, drugsets, magma_path, out, conditional, setsize, gene_id, enrich, enrich_auc, nsize, correct, use_pops = params
    
    try:
        analysis = DrugSetAnalysis(
            geneassoc, drugsets, magma_path, out, conditional,
            setsize, gene_id, enrich, enrich_auc, nsize, correct, use_pops
        )
        analysis.run_analysis()
    except Exception as e:
        pytest.fail(f"Failed with parameters: {params}, Error: {e}")
