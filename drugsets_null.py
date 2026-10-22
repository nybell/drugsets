#!/usr/bin/env python3
"""
drugsets_null.py
Size-matched null gene sets for DRUGSETS empirical p-values.
Python port of generate_null_sets() in 02.6a_null_gene_sets.R.

Workflow for one null distribution (run by DrugSetAnalysis.run_null() in drugsets.py):
  1. generate_null_sets()   -> MAGMA set file of random, size-matched gene sets
  2. MAGMA on that set file, with exactly the settings of the real run
     (DrugSetAnalysis.build_magma_cmd)
  3. null_distribution()    -> STAT = BETA / SE of every null set, saved as .parquet,
                               which drugsets_post.process_gsa(null_file=...) reads

Genes are drawn from the genes MAGMA actually analysed (the GENE column of a
.genes.out file; for PoPS runs use the PoPS-adjusted .genes.out that DRUGSETS writes).
"""

import argparse                                  # command-line arguments for standalone use
import logging                                   # progress and warning messages
import os                                        # file paths and folders

import numpy as np                               # arrays and random number generation
import pandas as pd                              # reading and writing tables
from scipy.special import gammaln                # log-gamma, used to compute log(n choose k)


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------

def read_gene_pool(genes_out):
    """Unique, non-empty gene IDs from the GENE column of a MAGMA .genes.out file (file order)."""
    go = pd.read_csv(genes_out, sep=r"\s+", comment="#", usecols=["GENE"], dtype=str)  # read only the GENE column, as text
    genes = go["GENE"].dropna()                  # drop missing gene IDs
    genes = genes[genes != ""]                   # drop empty gene IDs
    return np.asarray(pd.unique(genes), dtype=object)  # unique genes, first-seen order, as a numpy array


def read_set_files(set_files):
    """(name, genes) for every line of one or more tab-delimited MAGMA set files."""
    if isinstance(set_files, (str, os.PathLike)):  # a single file was given
        set_files = [set_files]                  # wrap it in a list
    sets = []                                    # collects (set name, gene list) pairs
    for path in set_files:                       # each set file in turn
        with open(path) as fh:                   # open it for reading
            for line in fh:                      # one gene set per line
                parts = line.rstrip("\r\n").split("\t")  # strip the line ending and split on tabs
                if parts[0]:                     # skip lines with no set name (e.g. blank lines)
                    sets.append((parts[0], [g for g in parts[1:] if g]))  # name, then its non-empty genes
    return sets                                  # list of (name, genes) in file order


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------

def _draw_rows(rng, n_pool, s, k, max_cells=20_000_000):
    """k random s-subsets of range(n_pool) as a (k, s) array. Genes are unique within a row.

    Draws m > s indices per row with replacement and keeps the first s distinct ones in
    draw order, which is exactly sampling without replacement. Fast because s is small
    relative to n_pool (generate_null_sets drops sets above max_frac of the pool). The
    rare row with fewer than s distinct draws is redrawn with rng.choice(replace=False).
    """
    out = np.empty((k, s), dtype=np.int64)       # result: one row per null set, one column per gene
    target = min(s + 6 * np.sqrt(s) + 10, 0.9 * n_pool)  # distinct draws to aim for: s plus a wide safety margin
    m = max(s, int(np.ceil(np.log1p(-target / n_pool) / np.log1p(-1.0 / n_pool))))  # draws per row expected to give `target` distinct values
    rows_per_chunk = max(1, max_cells // m)      # rows per batch, to cap memory use

    filled = 0                                   # rows of `out` filled so far
    while filled < k:                            # until every set is drawn
        c = min(rows_per_chunk, k - filled)      # rows in this batch
        A = rng.integers(0, n_pool, size=(c, m))  # m random pool indices per row, with replacement
        order = np.argsort(A, axis=1, kind="stable")  # per row, positions that sort the values (ties keep draw order)
        sA = np.take_along_axis(A, order, axis=1)  # each row sorted
        first_sorted = np.ones(sA.shape, dtype=bool)  # mark first occurrences in the sorted rows...
        first_sorted[:, 1:] = sA[:, 1:] != sA[:, :-1]  # ...a value is new if it differs from the one before it
        first = np.empty_like(first_sorted)      # same marks, but in original draw order
        np.put_along_axis(first, order, first_sorted, axis=1)  # map the marks back to draw positions
        csum = np.cumsum(first, axis=1)          # running count of distinct values along each row
        ok = csum[:, -1] >= s                    # rows that got at least s distinct values
        keep = first & (csum <= s)               # the first s distinct values in each row

        block = out[filled:filled + c]           # the part of `out` this batch fills (a view, so writes go to `out`)
        block[ok] = A[ok][keep[ok]].reshape(-1, s)  # copy the kept values of the good rows
        for i in np.flatnonzero(~ok):            # rare rows with too few distinct values
            block[i] = rng.choice(n_pool, s, replace=False)  # draw them directly instead
        filled += c                              # move on to the next batch
    return out                                   # (k, s) array of pool indices


def _dup_rows(idx):
    """True for every row that repeats an earlier row as a set (gene order ignored)."""
    _, first = np.unique(np.sort(idx, axis=1), axis=0, return_index=True)  # sort genes in each row, find first copy of each distinct row
    dup = np.ones(len(idx), dtype=bool)          # start by marking every row as a duplicate
    dup[first] = False                           # unmark the first copy of each distinct set
    return dup                                   # True = repeats an earlier set


# ---------------------------------------------------------------------------
# main function
# ---------------------------------------------------------------------------

def _make_unique(rng, idx, n_pool, s, max_tries):
    """Redraw rows of idx (null sets of size s) that repeat an earlier row.
    Returns the updated idx and the number of duplicates left."""
    k = len(idx)                                                            # number of sets of this size
    dup = _dup_rows(idx)                                                    # True for each set that repeats an earlier one
    log_choose = gammaln(n_pool + 1) - gammaln(s + 1) - gammaln(n_pool - s + 1)  # log of the number of distinct sets of size s
    if log_choose < np.log(k):                                              # fewer distinct sets exist than were requested
        logging.warning(f"[!] WARNING: uniqueness impossible for size {s} "
                        f"({k} sets requested); {int(dup.sum())} duplicates kept.")  # warn and keep the duplicates
        return idx, int(dup.sum())                                          # nothing can be redrawn
    n_first, tries = int(dup.sum()), 0                                      # duplicates before redrawing; rounds used
    while dup.any() and tries < max_tries:                                  # redraw until unique or out of rounds
        idx[dup] = _draw_rows(rng, n_pool, s, int(dup.sum()))               # new random sets for duplicated rows only
        dup, tries = _dup_rows(idx), tries + 1                              # re-check and count the round
    if n_first:                                                             # report only if there was something to fix
        logging.info(f"  size {s}: {k:,} sets, {n_first:,} duplicate(s) fixed in {tries} round(s)"
                     + (f", {int(dup.sum())} STILL REMAIN" if dup.any() else ""))  # include any that remain
    return idx, int(dup.sum())                                              # updated sets and duplicates left


def generate_null_sets(set_files, genes_out, out_file, n_sets, cond_name="druggable",
                       write_cond=True, min_size=2, max_frac=0.25, unique_sets=False,
                       max_tries=50, seed=20260823):
    """Write a MAGMA set file of size-matched null gene sets, plus a key file.

    set_files   : one or more MAGMA set files. Used ONLY for the size distribution
                  (number of unique genes per set); no gene from these files is sampled.
    genes_out   : MAGMA .genes.out; its GENE column is the sampling pool.
    out_file    : output set file, e.g. 'pheno.drug.null.sets.txt'. The key is written
                  next to it as '<out_file minus .txt>.key.tsv'.
    n_sets      : number of null sets to make.
    cond_name   : name of the background (conditioning) row in set_files; not a real set.
    write_cond  : if True, also write a '<cond_name>' line (its genes restricted to the
                  pool) so MAGMA can condition on it exactly as in the real run.
    min_size    : real sets smaller than this are left out of the size table.
    max_frac    : real sets larger than this fraction of the gene pool are left out of the
                  size table, with a warning (default 0.25, i.e. 25% of the pool).
    unique_sets : if True, redraw null sets that duplicate another set of the same size.
    max_tries   : maximum redraw rounds per size (unique_sets only).
    seed        : RNG seed (numpy's generator, so draws differ from R's for the same seed).

    Returns the key as a DataFrame (null_set, source_set, size).
    """
    rng = np.random.default_rng(seed)                                       # seeded random number generator
    os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)  # make the output folder if missing

    # ---- sampling pool ----
    pool = read_gene_pool(genes_out)                                        # unique genes in the .genes.out GENE column
    n_pool = len(pool)                                                      # number of genes to sample from
    if n_pool < 2:                                                          # a set needs at least 2 genes
        raise ValueError("[!] ERROR: gene pool has fewer than 2 usable genes.")  # stop: pool too small
    logging.info(f"[+] INFO: Sampling pool: {n_pool} genes from {os.path.basename(genes_out)}")  # report pool size

    # ---- sizes of the real sets ----
    if isinstance(set_files, (str, os.PathLike)):                           # a single file was given
        set_files = [set_files]                                             # wrap it in a list
    sets = read_set_files(set_files)                                        # (name, genes) for every line of every file
    real = pd.DataFrame([(n, len(set(g))) for n, g in sets if n != cond_name],
                        columns=["set", "size"])                            # name and unique gene count of each real set
    max_size = int(max_frac * n_pool)                                       # largest set size allowed
    small, big = real["size"] < min_size, real["size"] > max_size           # sets outside the allowed size range
    if big.any():                                                           # some real sets are too large
        logging.warning(f"[!] WARNING: {int(big.sum())} real set(s) larger than {max_size} genes "
                        f"({100 * max_frac:g}% of the pool) left out of the size table: "
                        f"{', '.join(real.loc[big, 'set'].head(5))}{' ...' if big.sum() > 5 else ''}")  # name up to 5 of them
    real = real[~small & ~big]                                              # keep sets within the allowed range
    if real.empty:                                                          # nothing left to take sizes from
        raise ValueError(f"[!] ERROR: no real sets between {min_size} and {max_size} genes.")  # stop
    names, sizes = real["set"].to_numpy(object), real["size"].to_numpy(np.int64)  # arrays for fast indexing
    logging.info(f"[+] INFO: Real sets: {len(sizes)} ({int(small.sum())} below {min_size} genes dropped) "
                 f"| size range: {sizes.min()}-{sizes.max()}")              # report the size table

    # ---- size table: every real set once, then resampled up to n_sets ----
    n_real = len(sizes)                                                     # number of real sets
    pick = (rng.choice(n_real, n_sets, replace=False) if n_sets <= n_real   # fewer nulls than real sets: no replacement
            else np.r_[np.arange(n_real), rng.choice(n_real, n_sets - n_real)])  # else each real set once, rest at random
    src_names, src_sizes = names[pick], sizes[pick]                         # source set and size of each null set

    # ---- background (conditioning) line ----
    cond_genes = {g for n, gs in sets if n == cond_name for g in gs}        # genes in the background row(s)
    cond = [g for g in pool if g in cond_genes] if write_cond else []       # those in the pool, in pool order
    if write_cond and not cond:                                             # line requested but no genes to write
        logging.warning(f"[!] WARNING: no '{cond_name}' genes found in the pool; "
                        "no background line written.")                      # warn

    # ---- draw and write, one size at a time (different sizes can never be duplicates) ----
    logging.info(f"[+] INFO: Drawing {n_sets:,} null sets...")              # progress message
    n_done, n_dup_left = 0, 0                                               # sets written so far; duplicates kept
    with open(out_file, "w") as fh:                                         # open the output set file
        for s in np.unique(src_sizes):                                      # each distinct size, smallest first
            k = int((src_sizes == s).sum())                                 # number of null sets of this size
            idx = _draw_rows(rng, n_pool, int(s), k)                        # k random sets as pool indices, no gene twice
            if unique_sets and k > 1:                                       # optionally remove duplicated sets
                idx, left = _make_unique(rng, idx, n_pool, int(s), max_tries)  # redraw duplicates
                n_dup_left += left                                          # tally any that could not be fixed
            if s > 1 and np.any(np.diff(np.sort(idx, axis=1), axis=1) == 0):  # any gene twice within a set?
                raise RuntimeError("[!] Internal error: a null set contains a repeated gene.")  # should never happen
            fh.writelines(f"null.s{n_done + i + 1}\t" + "\t".join(row) + "\n"
                          for i, row in enumerate(pool[idx]))               # one line per set: name, then its genes
            n_done += k                                                     # advance the set counter
        if cond:                                                            # background genes to write
            fh.write("\t".join([cond_name] + cond) + "\n")                  # background line at the end

    if unique_sets:                                                         # summarise duplicate handling
        logging.info(f"[+] INFO: Duplicate sets remaining: {n_dup_left} "
                     f"({100 * n_dup_left / n_sets:.3g}%)")                 # count and percentage
    logging.info("[+] INFO: Check passed: no repeated genes within any null set.")  # every size group passed
    if cond:                                                                # background line was written
        logging.info(f"[+] INFO: Background line '{cond_name}' written ({len(cond)} genes).")  # report its size

    # ---- key: which real set each null set took its size from ----
    order = np.argsort(src_sizes, kind="stable")                            # file order: by size, then table order
    key = pd.DataFrame({"null_set": [f"null.s{i + 1}" for i in range(n_sets)],  # null set names, as in the file
                        "source_set": src_names[order],                     # real set whose size was used
                        "size": src_sizes[order]})                          # size of the null set
    key_file = (out_file[:-4] if out_file.endswith(".txt") else out_file) + ".key.tsv"  # key path next to the set file
    key.to_csv(key_file, sep="\t", index=False)                             # save the key
    logging.info(f"[+] INFO: Saved {out_file} and {key_file}")              # report the outputs
    return key                                                              # key as a DataFrame


# ---------------------------------------------------------------------------
# after MAGMA has been run on the null sets
# ---------------------------------------------------------------------------

def null_distribution(null_gsa_out, out_file, cond_name="druggable"):
    """STAT = BETA / SE for every null set in a MAGMA .gsa.out; saved for process_gsa().

    Written as .parquet (needs pyarrow) when out_file ends in .parquet, otherwise as a
    tab-delimited text file (.tsv / .tsv.gz). Both can be read by drugsets_post.load_null().
    """
    df = pd.read_csv(null_gsa_out, sep=r"\s+", comment="#")    # MAGMA results for the null sets, '#' header lines skipped
    name = df["FULL_NAME"] if "FULL_NAME" in df.columns else df["VARIABLE"]  # untruncated set names if MAGMA wrote them
    keep = name.astype(str) != cond_name                       # drop the background (conditioning) row
    if "TYPE" in df.columns:                                   # results include a TYPE column
        keep &= df["TYPE"].astype(str).str.upper() != "COVAR"  # also drop gene covariates
    df = df[keep]                                              # null sets only
    null = pd.DataFrame({"SET": name[keep].astype(str).to_numpy(),  # null set name
                         "NGENES": df["NGENES"].to_numpy(),    # genes MAGMA used for that set
                         "STAT": (df["BETA"] / df["SE"]).to_numpy()})  # test statistic, same as for the real sets
    if out_file.endswith(".parquet"):                          # parquet requested
        null.to_parquet(out_file, index=False)                 # save as parquet (needs pyarrow)
    else:                                                      # any other extension
        null.to_csv(out_file, sep="\t", index=False)           # save as tab-delimited text (gzipped if it ends in .gz)
    logging.info(f"[+] INFO: Null distribution: {len(null):,} statistics "
                 f"({int(null['STAT'].isna().sum())} NA) saved to {out_file}")  # report size, NAs and path
    return null                                                # the null distribution as a DataFrame


# ---------------------------------------------------------------------------
# standalone use
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Generate size-matched null gene sets for DRUGSETS.")  # command-line parser
    ap.add_argument("--set_file", nargs="+", required=True,
                    help="Real MAGMA set file(s); used only for the size distribution")  # one or more real set files
    ap.add_argument("--genes_out", required=True, help="MAGMA .genes.out defining the gene pool")  # sampling pool
    ap.add_argument("--out", required=True, help="Output null set file (.txt)")  # where to write the null sets
    ap.add_argument("--n_sets", type=int, required=True, help="Number of null sets")  # how many null sets
    ap.add_argument("--cond_name", default="druggable", help="Background row name (default: druggable)")  # background set name
    ap.add_argument("--no_cond", action="store_true", help="Do not write a background line")  # unconditioned runs
    ap.add_argument("--min_size", type=int, default=2, help="Minimum set size (default 2)")  # smallest real set size used
    ap.add_argument("--max_frac", type=float, default=0.25,
                    help="Drop real sets larger than this fraction of the pool (default 0.25)")  # size limit
    ap.add_argument("--unique_sets", action="store_true", help="Redraw duplicated null sets")  # remove duplicate null sets
    ap.add_argument("--max_tries", type=int, default=50)       # redraw rounds per size
    ap.add_argument("--seed", type=int, default=20260823)      # random seed
    a = ap.parse_args()                                        # read the arguments

    generate_null_sets(a.set_file, a.genes_out, a.out, a.n_sets, cond_name=a.cond_name,
                       write_cond=not a.no_cond, min_size=a.min_size, max_frac=a.max_frac,
                       unique_sets=a.unique_sets,
                       max_tries=a.max_tries, seed=a.seed)     # make the null sets


if __name__ == "__main__":                                     # only when run as a script, not when imported
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")           # show INFO messages with a timestamp
    main()                                                     # run
