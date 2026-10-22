"""
Full test of DRUGSETS v2.1 with real MAGMA on the SCZ (PGC3, no MHC) data.

Run from the drugsets folder:
    pytest -v drugsets_pytest.py
    pytest -v drugsets_pytest.py -k "not all_set_types"     # skip the slowest test

PoPS tests use a small null (N_NULL) so they finish quickly. Every test writes to its own
folder under OUT_ROOT, which is emptied at the start of that test, so outputs can be inspected.
"""

import os
import shutil

import pandas as pd
import pytest

from drugsets import DrugSetAnalysis, ALL_SET_KEYS

# ---- inputs ----
RAW = "/Users/nyb/phd_code/drugsets-may2026/output/magma-v110-s1s2-pops-no-mhc/scz.no.mhc.pgc3.step2.genes.raw"
GENES_OUT = "/Users/nyb/phd_code/drugsets-may2026/output/magma-v110-s1s2-pops-no-mhc/scz.no.mhc.pgc3.step2.genes.out"
PREDS = "/Users/nyb/phd_code/drugsets-may2026/pops-out-no-mhc/pops-coexp/scz.no.mhc.pgc3.coexp.preds"
COV = ["/Users/nyb/drugsets/DATA/COV/ensembl_ukb_mean_pops.cov", "median_pops_norm"]   # gene covariate file + column

# ---- settings ----
OUT_ROOT = "/Users/nyb/Desktop/drugsets-testing-2026/pytest"
MAGMA = "magma"
N_NULL = 2_000

SET_KEYS = ALL_SET_KEYS + ["moa_combined"]

pytestmark = pytest.mark.skipif(shutil.which(MAGMA) is None, reason=f"MAGMA not found: {MAGMA}")


# ---- helpers ----

def run(name, **kwargs):
    """Run DRUGSETS as drugsets.py does (analysis -> post-processing -> enrichment) in OUT_ROOT/<name>."""
    out = os.path.join(OUT_ROOT, name, "scz")
    shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    args = dict(raw_file=RAW, drugsets="drug_mp0", magma_path=MAGMA, out=out, conditional="yes", setsize=2,
                gene_id="ensembl", enrich=None, enrich_file=None, enrich_sets="drug_mp0", correct_cov="TRUE",
                nsize=5, correct="bonf", magma_cov=None, use_pops=None, adj_pops="yes", dgidb_filter="median",
                n_null=N_NULL)
    args.update(kwargs)
    analysis = DrugSetAnalysis(**args)

    keys = ALL_SET_KEYS if analysis.drugsets == "all" else [analysis.drugsets]
    for key in keys:
        if not os.path.isfile(analysis.filepaths["ensembl"][key]):
            pytest.skip(f"gene set file not found: {analysis.filepaths['ensembl'][key]}")

    analysis.run_analysis()
    analysis.post_process()
    if analysis.enrich:
        analysis.enrichment()
    return analysis, keys


def header(path):
    """'# KEY = value' lines at the top of a .drugs.out / .moa.out file, as a dict."""
    info = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if " = " in line:
                key, value = line[1:].split(" = ", 1)
                info[key.strip()] = value.strip()
    return info


def check(analysis, keys, p_source):
    """Results, MAGMA files and null files are where they should be, and look right."""
    for key in keys:
        result = f"{analysis.out}.{key}." + ("drugs.out" if key.startswith("drug") else "moa.out")
        assert os.path.isfile(result), f"missing result file {result}"
        assert os.path.isfile(f"{analysis.magma_out}.{key}.gsa.out"), "missing MAGMA .gsa.out in magma/"

        info = header(result)
        assert info["P_SOURCE"] == p_source
        assert info["CORRECTION"].startswith(analysis.gsa_correct)
        assert int(info["MIN_NGENES"]) == analysis.setsize
        assert 0 < int(info["UNIQUE_SETS"]) <= int(info["TOTAL_SETS"])

        res = pd.read_csv(result, sep="\t", comment="#")
        assert len(res) == int(info["UNIQUE_SETS"]), "one row per unique set expected"
        assert (res["NGENES"] >= analysis.setsize).all()
        assert res["P"].between(0, 1).all() and res["P_ADJ"].between(0, 1).all()
        assert "druggable" not in set(res["FULL_NAME"])

        null_prefix = f"{analysis.null_out}.{key}.null"
        if p_source == "empirical":
            null = pd.read_parquet(null_prefix + ".dist.parquet")
            assert len(null) >= 0.9 * analysis.n_null, f"only {len(null)} of {analysis.n_null} null sets in the parquet"
            assert res["P_EMP"].min() >= (1 / (len(null) + 1)) * (1 - 1e-5)   # smallest possible p (file keeps 6 digits)
            assert os.path.isfile(null_prefix + ".sets.key.tsv")
            assert not os.path.exists(null_prefix + ".sets.txt"), "null set file should be deleted"
            assert not os.path.exists(null_prefix + ".gsa.out"), "null .gsa.out should be deleted"
        else:
            assert not os.path.exists(null_prefix + ".dist.parquet"), "no null expected"
            assert res["P_EMP"].isna().all()


# ---- tests ----

@pytest.mark.parametrize("set_type", SET_KEYS)
def test_pops_empirical_each_set_type(set_type):
    """Default PoPS run: empirical p-values for every drug / MoA set type."""
    analysis, keys = run(f"pops_{set_type}", drugsets=set_type, use_pops=(GENES_OUT, PREDS))
    check(analysis, keys, "empirical")


def test_pops_empirical_all_set_types():
    """--drugsets all: one result file and one null per set type."""
    analysis, keys = run("pops_all", drugsets="all", use_pops=(GENES_OUT, PREDS))
    check(analysis, keys, "empirical")


@pytest.mark.parametrize("set_type", ["drug_mp0", "moa_combined"])
def test_pops_unconditioned(set_type):
    """--conditional no: no druggable line in the null sets, still empirical p-values."""
    analysis, keys = run(f"pops_uncond_{set_type}", drugsets=set_type, conditional="no", use_pops=(GENES_OUT, PREDS))
    check(analysis, keys, "empirical")


def test_pops_without_empirical_p():
    """--empirical_p no: MAGMA p-values, no null files."""
    analysis, keys = run("pops_no_empirical", use_pops=(GENES_OUT, PREDS), empirical_p="no")
    check(analysis, keys, "magma")


def test_pops_options():
    """Non-default options: z-scored PoPS, minimum set size 5, Bonferroni, other seed and max_frac."""
    analysis, keys = run("pops_options", use_pops=(GENES_OUT, PREDS), adj_pops="no", setsize=5,
                         gsa_correct="bonf", null_seed=1, max_frac=0.1)
    check(analysis, keys, "empirical")


def test_null_seed_reproducible():
    """Same seed -> same null sets; a different seed -> different null sets."""
    keys = {}
    for name, seed in (("seed_a", 20260823), ("seed_b", 20260823), ("seed_c", 7)):
        analysis, _ = run(name, use_pops=(GENES_OUT, PREDS), null_seed=seed)
        keys[name] = pd.read_csv(f"{analysis.null_out}.drug_mp0.null.sets.key.tsv", sep="\t")
    assert keys["seed_a"].equals(keys["seed_b"])
    assert not keys["seed_a"].equals(keys["seed_c"])


@pytest.mark.parametrize("set_type", ["drug_mp0", "moa_group"])
def test_magma_only(set_type):
    """No PoPS: MAGMA p-values, FDR-corrected, no null."""
    analysis, keys = run(f"magma_{set_type}", drugsets=set_type)
    check(analysis, keys, "magma")


@pytest.mark.parametrize("dgidb", ["all", "median", "p75", "none"])
def test_dgidb_filters(dgidb):
    """Every DGIdb filter level has a drug set file that runs."""
    analysis, keys = run(f"dgidb_{dgidb}", dgidb_filter=dgidb)
    check(analysis, keys, "magma")


@pytest.mark.skipif(not os.path.isfile(COV[0]), reason=f"covariate file not found: {COV[0]}")
def test_pops_with_covariate():
    """--magma_cov: covariate is conditioned on in both the real and the null MAGMA runs."""
    analysis, keys = run("pops_covariate", use_pops=(GENES_OUT, PREDS), magma_cov=COV)
    check(analysis, keys, "empirical")
    with open(f"{analysis.null_out}.drug_mp0.null.log") as fh:
        assert COV[1] in fh.read(), "covariate missing from the null MAGMA run"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not found")
def test_enrichment():
    """--enrich atc1 runs after the drug gene-set analysis."""
    analysis, keys = run("enrich_atc1", enrich="atc1")
    check(analysis, keys, "magma")
