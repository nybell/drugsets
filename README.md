# DRUGSETS

**DRUG gene SET analysiS**: genetically informed drug repositioning with drug and mechanism-of-action (MoA) gene-set analysis in [MAGMA](https://cncr.nl/research/magma/).

DRUGSETS tests whether the genes targeted by a drug (or by drugs sharing a mechanism of action) are more strongly associated with a trait than other genes. It can use standard MAGMA gene results or [PoPS](https://github.com/FinucaneLab/pops) gene scores. For PoPS-based runs it computes empirical p-values from size-matched null gene sets. It can also test whether drugs in an ATC group are enriched for association.

Developed by Nathaniel Bell (n.y.bell@vu.nl).

---

## Installation

**Requirements**
- Python 3 with `numpy`, `pandas`, `scipy` and `pyarrow`
- [MAGMA](https://cncr.nl/research/magma/) v1.10, callable as `magma` (add it to your `PATH`, or give its location with `--magma`)
- R with `tidyr` and `dplyr`, only for drug group enrichment (`--enrich`)

**Setup**
```bash
git clone https://github.com/nybell/drugsets.git
cd drugsets
pip install numpy pandas scipy pyarrow
```

**Repository contents**

| File | Purpose |
|---|---|
| `drugsets.py` | Main script: run this |
| `drugsets_func.py`, `drugsets_null.py`, `drugsets_post.py` | Helper modules (input checks, null gene sets, result processing) |
| `compute_corrs.R`, `compute_lnreg.R` | Drug group enrichment |
| `DATA/GENESETS/` | Drug and MoA gene sets |
| `DATA/meta.data.2026.RData` | Drug annotations (ATC codes) for enrichment |

Keep everything in place: DRUGSETS finds the gene sets and helper scripts relative to `drugsets.py`.

---

## Input files

If running MAGMA gene z-score-based analysis, only the first file is needed. If running the PoPS-based analysis, you need a PoPS .pred file and the MAGMA .genes.out. All files need **Ensembl** gene IDs:

| File | Needed for |
|---|---|
| `[NAME].genes.raw` | every run |
| `[NAME].genes.out` | PoPS-based runs |
| `[NAME].preds` (PoPS output) | PoPS-based runs |

---

## Usage

**MAGMA gene results**
```bash
python drugsets.py \
  --raw_file scz.genes.raw \
  --gene_id ensembl \
  --drugsets drug_mp0 \
  --out results/scz
```

**PoPS gene scores** (empirical p-values are computed by default)
```bash
python drugsets.py \
  --raw_file scz.genes.raw \
  --use_pops scz.genes.out scz.preds \
  --gene_id ensembl \
  --drugsets drug_mp0 \
  --out results/scz
```

The empirical null uses 10.1 million null gene sets per set type by default, which takes time. For a quick test run, add `--n_null 10000`.

**Gene sets** (`--drugsets`)

| Value | Gene sets |
|---|---|
| `drug_mp0` | All drugs |
| `drug_mp1` … `drug_mp4` | Drugs reaching clinical phase 1, 2, 3 or approval (`mp4`) |
| `moa_group` | Broad mechanism-of-action groups |
| `moa_nested` | Specific MoAs nested within broad groups |
| `moa_combined` | `moa_group` + `moa_nested` in one analysis |
| `moa_specific` | All specific MoAs |
| `all` | Every set type above except `moa_combined`, each analysed separately |

**Options**

| Option | Description | Default |
|---|---|---|
| `--raw_file` | MAGMA `.genes.raw` file | required |
| `--gene_id` | Gene ID type of the input (`ensembl`) | required |
| `--out` | Output folder + file name prefix | required |
| `--drugsets` | Gene sets to test (see above) | — |
| `--conditional` | Condition on all druggable genes (`yes`/`no`) | `yes` |
| `--dgidb_filter` | DGIdb interaction score filter for drug sets (`all`, `median`, `p75`, `none`) | `median` |
| `--setsize` | Minimum number of genes per set (genes analysed by MAGMA) | `2` |
| `--gsa_correct` | Multiple-testing correction of the results (`fdr`, `bonf`) | `fdr` |
| `--magma_cov` | Gene covariate file followed by one or more column names, conditioned on in MAGMA | — |
| `--magma` | Path to the MAGMA executable | `magma` |
| `--use_pops` | MAGMA `.genes.out` and PoPS `.preds` file | — |
| `--adj_pops` | PoPS runs: rank-based inverse normal transform of the PoPS scores (`yes`/`no`); otherwise z-scores | `no` |
| `--empirical_p` | PoPS runs: empirical p-values (`yes`/`no`) | `yes` |
| `--n_null` | PoPS runs: number of null gene sets per set type | `10100000` |
| `--null_seed` | PoPS runs: random seed for the null gene sets | `20260823` |
| `--max_frac` | PoPS runs: real sets larger than this fraction of the genes are not used to size null sets | `0.25` |

Run `python drugsets.py --help` for all options.

### Drug group enrichment

Tests whether drugs in an ATC group (`--enrich atc1` … `atc4`) are more strongly associated with the trait than other drugs. Needs R.

**In the same run** as the drug gene-set analysis:
```bash
python drugsets.py \
  --raw_file scz.genes.raw \
  --gene_id ensembl \
  --drugsets drug_mp0 \
  --enrich atc3 \
  --out results/scz
```

**On existing results**: point `--enrich_file` to the `.drugs.out` file of an earlier run (DRUGSETS uses the matching MAGMA `.gsa.out` in its `_magma` folder), without `--drugsets`:
```bash
python drugsets.py \
  --raw_file scz.genes.raw \
  --gene_id ensembl \
  --enrich_file results/scz.drug_mp0.drugs.out \
  --enrich_sets drug_mp0 \
  --enrich atc3 \
  --out results/scz.atc3
```

| Option | Description | Default |
|---|---|---|
| `--enrich` | ATC level to test (`atc1`, `atc2`, `atc3`, `atc4`) | — |
| `--enrich_file` | `.drugs.out` file (or MAGMA `.gsa.out` file) of an earlier run | — |
| `--enrich_sets` | Drug set type of the drug-level results (`drug_mp0` … `drug_mp4`); must match `--drugsets` or the `--enrich_file` | `drug_mp0` |
| `--correct_cov` | Account for correlation between drug gene sets (`TRUE`/`FALSE`) | `TRUE` |
| `--nsize` | Minimum number of drugs per group | `5` |
| `--correct` | Multiple-testing correction (`bonf`, `fdr`) | `bonf` |

Use the same `--conditional` setting as the run that produced the drug-level results. Enrichment uses MAGMA's drug-level test statistics (`BETA / SE`), also for PoPS runs.

---

## Output

Results are written to the `--out` location. Supporting files go into two folders named after the `--out` prefix, so several runs can share one output folder:

```
results/
├── scz.drug_mp0.drugs.out              results, drug sets
├── scz.moa_group.moa.out               results, MoA sets
├── scz_magma/                          MAGMA output (.gsa.out, .log, ...); PoPS runs also scz.pops.genes.out
└── scz_null/                           PoPS runs only: empirical null files
```

For the MAGMA files, see the [MAGMA manual](https://cncr.nl/research/magma/).

### `.drugs.out` / `.moa.out`

One file per set type, one row per unique gene set, sorted by `P`. Tab-delimited, with a summary header (lines starting with `#`): input file, set type, minimum set size, number of sets, number of unique sets, p-value source, correction method and number of significant sets.

| Column | Description |
|---|---|
| `FULL_NAME` | Drug or MoA name |
| `VARIABLE`, `TYPE` | MAGMA set name and type |
| `NGENES` | Number of genes analysed by MAGMA |
| `BETA`, `BETA_STD`, `SE` | MAGMA effect size, standardised effect size and standard error |
| `STAT` | Test statistic, `BETA / SE` |
| `P_MAGMA` | MAGMA p-value |
| `P_EMP` | Empirical p-value (PoPS runs); `NA` otherwise |
| `P` | p-value used for testing: `P_EMP` for PoPS runs with an empirical null, otherwise `P_MAGMA` |
| `P_SOURCE` | `empirical` or `magma` |
| `P_ADJ` | `P` corrected for multiple testing |
| `ADJ_METHOD` | Correction method (`fdr` or `bonf`) |
| `N_TESTS` | Number of tests corrected for (unique sets) |
| `SIG` | `True` if `P_ADJ` < 0.05 |
| `N_DUP` | Number of sets with identical results (1 = unique) |
| `DUP_MEMBERS` | Names of all sets with identical results, `;`-separated |

**Duplicate sets.** Drugs that target exactly the same genes give identical results. These count as one test: one row is kept (the alphabetically first name), and all names are listed in `DUP_MEMBERS`.

**Empirical p-values (PoPS runs).** MAGMA p-values are anti-conservative when PoPS scores replace the gene Z-scores. DRUGSETS therefore draws random gene sets, matching the sizes of the real sets, from the genes in the PoPS analysis. It runs them through MAGMA with the same settings, and computes the one-sided
`P_EMP = (M + 1) / (N + 1)`, where `N` is the number of null sets and `M` the number with a test statistic at least as large as the observed one.

### `[OUT]_null/` (PoPS runs)

| File | Description |
|---|---|
| `[OUT].[SET].null.dist.parquet` | Null distribution: `SET` (null set name), `NGENES`, `STAT` |
| `[OUT].[SET].null.sets.key.tsv` | For each null set, the real set whose size it took: `null_set`, `source_set`, `size` |
| `[OUT].[SET].null.log` | MAGMA log of the null run |

### Drug group enrichment

| File | Description |
|---|---|
| `[OUT].[ATC].groups.out` | All ATC groups tested, sorted by `P` |
| `[OUT].[ATC].[bonf/fdr].groups.out` | Significant groups only; written only if there are any |
| `[OUT]_setcorrs.rdata` | Drug-level statistics and gene-set correlations used for the test |

| Column | Description |
|---|---|
| `GROUP` | ATC code |
| `BETA`, `SE` | Difference in drug-level test statistic between drugs in the group and other drugs (adjusted for gene set size), and its standard error |
| `CI_lower`, `CI_upper` | 95% confidence interval of `BETA` |
| `T`, `DF`, `P` | Test statistic, degrees of freedom and one-sided p-value |
| `GROUP_N` | Number of drugs in the group |
| `fdrP` | FDR-corrected p-value (`--correct fdr` only) |

---

## Standalone scripts

- `drugsets_post.py` re-processes an existing MAGMA `.gsa.out` into a `.drugs.out` / `.moa.out` file, e.g. with a different `--min_size` or `--gsa_correct`.
- `drugsets_null.py` writes size-matched null gene sets for a set file without running MAGMA.

Run either with `--help` for its options.

---

## Support

Please use the [issue tracker](https://github.com/nybell/drugsets/issues) for questions, issues or feedback.

## Citation

If you use DRUGSETS, please also cite MAGMA:

> de Leeuw C, Mooij J, Heskes T, Posthuma D (2015). MAGMA: Generalized gene-set analysis of GWAS data. *PLoS Comput Biol* 11(4): e1004219. doi:10.1371/journal.pcbi.1004219

and, for PoPS-based runs, PoPS (Weeks et al., 2023, *Nature Genetics*).

## License

MIT License
