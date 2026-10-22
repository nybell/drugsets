#!/usr/bin/env python3
"""
drugsets_post.py
Post-processing of MAGMA drug gene-set results for DRUGSETS v2.1.

For one .gsa.out file, process_gsa():
  1. reads the results and drops hidden covariates (e.g. 'druggable', gene covariates)
  2. drops user-excluded sets and sets below a minimum size, using NGENES
     (genes actually used by MAGMA, not genes listed in the set file)
  3. collapses duplicated gene sets: rows with identical NGENES, BETA, BETA_STD, SE and P
     become one row, with all member names kept in DUP_MEMBERS
  4. computes empirical p-values for STAT = BETA / SE against a null distribution (PoPS runs)
  5. corrects for multiple testing (BH-FDR or Bonferroni); after collapsing, every row is one test
  6. writes <prefix>.drugs.out (drug sets) or <prefix>.moa.out (MoA sets): tab-delimited,
     with a '#' comment header summarising set counts and the correction used

Can be imported (process_gsa) or run standalone on existing results:
  python drugsets_post.py --gsa_file out.drug_mp0.gsa.out
"""

import argparse
import logging
import os

import numpy as np
import pandas as pd

GSA_REQUIRED = ["VARIABLE", "NGENES", "BETA", "SE", "P"]
DUP_COLS = ["NGENES", "BETA", "BETA_STD", "SE", "P"]


# ---------------------------------------------------------------------------
# reading inputs
# ---------------------------------------------------------------------------

def read_gsa_out(path):
    """Read a MAGMA .gsa.out file (header comment lines start with '#')."""
    df = pd.read_csv(path, sep=r"\s+", comment="#")
    missing = [c for c in GSA_REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"[!] ERROR: {path} is missing columns: {', '.join(missing)}")
    # MAGMA truncates VARIABLE; FULL_NAME holds the untruncated name when present
    if "FULL_NAME" in df.columns:
        df["FULL_NAME"] = df["FULL_NAME"].fillna(df["VARIABLE"])
    else:
        df["FULL_NAME"] = df["VARIABLE"]
    df["FULL_NAME"] = df["FULL_NAME"].astype(str)
    return df


def read_exclude(path):
    """Sets to drop: a CSV/TSV with a 'full_name' column, or one name per line."""
    try:
        tab = pd.read_csv(path, sep=None, engine="python")
        if "full_name" in tab.columns:
            return set(tab["full_name"].astype(str))
    except Exception:
        pass
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_null(path):
    """Load a null distribution of STAT values (.parquet, or whitespace-delimited text)."""
    if path.endswith(".parquet"):
        stat = pd.read_parquet(path, columns=["STAT"])["STAT"]  # needs pyarrow
    else:
        tab = pd.read_csv(path, sep=r"\s+", comment="#")
        if "STAT" in tab.columns:
            stat = tab["STAT"]
        else:  # headerless single column
            stat = pd.read_csv(path, sep=r"\s+", comment="#", header=None).iloc[:, 0]
    null = pd.to_numeric(stat, errors="coerce").to_numpy(float)
    n_na = int(np.isnan(null).sum())
    null = null[~np.isnan(null)]
    if null.size == 0:
        raise ValueError(f"[!] ERROR: null distribution is empty or all NA: {path}")
    logging.info(f"[+] INFO: Loaded null distribution: N = {null.size} ({n_na} NA dropped)")
    return null


# ---------------------------------------------------------------------------
# core steps
# ---------------------------------------------------------------------------

def collapse_duplicates(df, cols=DUP_COLS, keep="alpha"):
    """Collapse rows with identical values in `cols` into one row per unique test.

    The kept row is the alphabetically first FULL_NAME (keep="alpha", case-insensitive)
    or the shortest FULL_NAME (keep="shortest", ties broken alphabetically); its
    FULL_NAME, VARIABLE and all other columns are kept as they are.
    Adds:
      N_DUP       : number of sets sharing these results (1 = not duplicated)
      DUP_MEMBERS : all member FULL_NAMEs, alphabetically, ';'-separated
    """
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"[!] ERROR: cannot check duplicates, missing columns: {', '.join(missing)}")
    if keep not in ("alpha", "shortest"):
        raise ValueError(f"[!] ERROR: keep must be 'alpha' or 'shortest', not {keep!r}")

    df = df.copy()
    # exact match on the values as read from file (NaN compares equal to NaN here)
    key = df[cols].astype(str).agg("|".join, axis=1)
    grp = df.groupby(key, sort=False)["FULL_NAME"]
    df["N_DUP"] = grp.transform("size")
    df["DUP_MEMBERS"] = grp.transform(lambda s: ";".join(sorted(s, key=lambda n: (n.lower(), n))))

    # order rows so the one to keep comes first within its group
    name_lower = df["FULL_NAME"].str.lower()
    sort_by = [name_lower, df["FULL_NAME"]]
    if keep == "shortest":
        sort_by.insert(0, df["FULL_NAME"].str.len())
    order = np.lexsort([s.to_numpy() for s in reversed(sort_by)])

    out = (df.iloc[order]
             .assign(_key=key.iloc[order].to_numpy())
             .drop_duplicates("_key", keep="first")
             .drop(columns="_key")
             .sort_index())  # back to original file order

    n_groups = int((out["N_DUP"] > 1).sum())
    logging.info(f"[+] INFO: {len(df)} sets collapsed to {len(out)} unique tests "
                 f"({n_groups} duplicate groups, {len(df) - len(out)} rows merged).")
    return out.reset_index(drop=True)


def empirical_p(stat, null):
    """One-sided (upper tail) empirical p = (M + 1) / (N + 1), where M is the number of
    null statistics >= the observed statistic. Matches compute_empirical_p() in the R code."""
    srt = np.sort(np.asarray(null, float))
    stat = np.asarray(stat, float)
    N = srt.size
    M = N - np.searchsorted(srt, stat, side="left")
    p = (M + 1) / (N + 1)
    return np.where(np.isnan(stat), np.nan, p)


# ---------------------------------------------------------------------------
# main entry point
# ---------------------------------------------------------------------------

def infer_set_type(gsa_file):
    """'drug' or 'moa' from the DRUGSETS label in the file name (e.g. .drug_mp0, .moa_group)."""
    name = os.path.basename(gsa_file)
    if ".drug_" in name:
        return "drug"
    if ".moa_" in name:
        return "moa"
    raise ValueError(f"[!] ERROR: cannot tell whether {name} holds drug or MoA sets; "
                     "pass set_type='drug' or set_type='moa'.")


def process_gsa(gsa_file, null_file=None, use_pops=False, min_size=2, exclude=None,
                hide=("druggable",), keep="alpha", method="fdr", alpha=0.05,
                set_type=None, out_file=None):
    """Clean, collapse duplicates, (empirically) re-test and correct one MAGMA .gsa.out file.

    gsa_file  : MAGMA .gsa.out
    null_file : null distribution of STAT; if given, P_EMP becomes the primary p
    use_pops  : True for PoPS runs; warns if no null is supplied
    min_size  : minimum NGENES to keep a set
    exclude   : iterable of FULL_NAMEs to drop before testing
    hide      : covariate names to drop (conditioned-on sets, gene covariates)
    keep      : which duplicate's name to keep: 'alpha' (first alphabetically) or 'shortest'
    method    : multiple-testing correction, 'fdr' (Benjamini-Hochberg, default) or 'bonf'
    alpha     : significance threshold on P_ADJ
    set_type  : 'drug' or 'moa'; inferred from the file name (.drug_* / .moa_*) if not given
    out_file  : output path (default: <prefix>.drugs.out or <prefix>.moa.out)
    """
    if method not in ("fdr", "bonf"):
        raise ValueError(f"[!] ERROR: method must be 'fdr' or 'bonf', not {method!r}")
    if set_type is None:
        set_type = infer_set_type(gsa_file)
    if set_type not in ("drug", "moa"):
        raise ValueError(f"[!] ERROR: set_type must be 'drug' or 'moa', not {set_type!r}")

    logging.info(f"[+] INFO: Post-processing {gsa_file}")
    df = read_gsa_out(gsa_file)
    n_read = len(df)

    # 1. drop covariates / conditioned-on sets
    if "TYPE" in df.columns:
        df = df[df["TYPE"].astype(str).str.upper() != "COVAR"]
    df = df[~df["FULL_NAME"].isin(set(hide or ()))]

    # 2. user exclusions and minimum set size (on NGENES actually analysed)
    if exclude:
        df = df[~df["FULL_NAME"].isin(set(exclude))]
    n_before_size = len(df)
    df = df[df["NGENES"] >= min_size].reset_index(drop=True)
    logging.info(f"[+] INFO: {n_read} rows read; {n_before_size - len(df)} sets dropped with "
                 f"NGENES < {min_size}; {len(df)} sets kept.")
    if df.empty:
        raise ValueError("[!] ERROR: no gene sets left after filtering.")

    # 3. collapse duplicated sets into one row per unique test
    n_total = len(df)
    df = collapse_duplicates(df, keep=keep)

    # 4. test statistic and p-values
    df["STAT"] = df["BETA"] / df["SE"]
    df = df.rename(columns={"P": "P_MAGMA"})
    if null_file:
        df["P_EMP"] = empirical_p(df["STAT"], load_null(null_file))
        df["P"] = df["P_EMP"]
        df["P_SOURCE"] = "empirical"
    else:
        if use_pops:
            logging.warning("[!] WARNING: PoPS run without a null distribution. Reporting MAGMA "
                            "p-values, which are anti-conservative for PoPS runs (P_SOURCE = magma).")
        df["P_EMP"] = np.nan
        df["P"] = df["P_MAGMA"]
        df["P_SOURCE"] = "magma"

    # 5. multiple testing over unique tests (NA p-values stay NA and are not counted)
    p = df["P"].to_numpy(float)
    ok = ~np.isnan(p)
    n_tests = int(ok.sum())
    p_adj = np.full(p.size, np.nan)
    if method == "bonf":
        p_adj[ok] = np.minimum(p[ok] * n_tests, 1.0)
    else:  # Benjamini-Hochberg
        p_ok = p[ok]
        order = np.argsort(p_ok)
        ranked = p_ok[order] * n_tests / np.arange(1, n_tests + 1)
        q = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity from the top
        q_ok = np.empty(n_tests)
        q_ok[order] = np.minimum(q, 1.0)
        p_adj[ok] = q_ok
    df["P_ADJ"] = p_adj
    df["N_TESTS"] = n_tests
    df["ADJ_METHOD"] = method
    df["SIG"] = df["P_ADJ"] < alpha  # NA -> False

    # 6. write
    cols = ["FULL_NAME", "VARIABLE", "TYPE", "NGENES", "BETA", "BETA_STD", "SE", "STAT",
            "P_MAGMA", "P_EMP", "P", "P_SOURCE", "P_ADJ", "ADJ_METHOD", "N_TESTS", "SIG",
            "N_DUP", "DUP_MEMBERS"]
    df = df[[c for c in cols if c in df.columns]].sort_values("P", kind="stable").reset_index(drop=True)

    if out_file is None:
        prefix = gsa_file
        for ext in (".gsa.out.txt", ".gsa.out"):
            if prefix.endswith(ext):
                prefix = prefix[: -len(ext)]
                break
        out_file = prefix + (".drugs.out" if set_type == "drug" else ".moa.out")

    n_sig = int(df["SIG"].sum())                         # significant unique tests
    n_sig_all = int(df.loc[df["SIG"], "N_DUP"].sum())   # ... counting every duplicate member
    method_name = "Benjamini-Hochberg FDR" if method == "fdr" else "Bonferroni"
    header = [
        "# DRUGSETS v2.1 post-processed gene-set results",
        f"# INPUT_FILE = {os.path.basename(gsa_file)}",
        f"# SET_TYPE = {set_type}",
        f"# MIN_NGENES = {min_size}",
        f"# TOTAL_SETS = {n_total}",
        f"# UNIQUE_SETS = {n_tests}",
        f"# P_SOURCE = {df['P_SOURCE'].iloc[0]}",
        f"# CORRECTION = {method} ({method_name})",
        f"# ALPHA = {alpha}",
        f"# SIGNIFICANT_UNIQUE_SETS = {n_sig}",
        f"# SIGNIFICANT_SETS_INCL_DUPLICATES = {n_sig_all}",
    ]
    with open(out_file, "w") as fh:
        fh.write("\n".join(header) + "\n")
        df.to_csv(fh, sep="\t", index=False, na_rep="NA", float_format="%.6g")

    logging.info(f"[+] INFO: {n_total} sets, {n_tests} unique; {n_sig} unique sets with "
                 f"{method} P_ADJ < {alpha}.")
    logging.info(f"[+] INFO: Saved {out_file}")
    return df


# ---------------------------------------------------------------------------
# standalone use
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Post-process a DRUGSETS / MAGMA .gsa.out file.")
    ap.add_argument("--gsa_file", required=True, help="MAGMA .gsa.out file")
    ap.add_argument("--null_dist", help="Null distribution of STAT (.parquet or text)")
    ap.add_argument("--use_pops", action="store_true", help="Results come from a PoPS run")
    ap.add_argument("--min_size", type=int, default=2, help="Minimum NGENES (default 2)")
    ap.add_argument("--exclude", help="File of set names to drop (one per line, or CSV with 'full_name')")
    ap.add_argument("--hide", nargs="*", default=["druggable"], help="Covariate names to drop")
    ap.add_argument("--keep", default="alpha", choices=["alpha", "shortest"],
                    help="Which duplicate's name to keep (default: alphabetically first)")
    ap.add_argument("--gsa_correct", default="fdr", choices=["fdr", "bonf"],
                    help="Multiple-testing correction (default: fdr)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--set_type", choices=["drug", "moa"],
                    help="Drug or MoA sets (default: inferred from the file name)")
    ap.add_argument("--out_file", help="Output path (default: <prefix>.drugs.out / <prefix>.moa.out)")
    a = ap.parse_args()

    process_gsa(a.gsa_file, null_file=a.null_dist, use_pops=a.use_pops, min_size=a.min_size,
                exclude=read_exclude(a.exclude) if a.exclude else None, hide=a.hide,
                keep=a.keep, method=a.gsa_correct, alpha=a.alpha, set_type=a.set_type,
                out_file=a.out_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    main()
