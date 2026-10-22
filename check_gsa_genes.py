#!/usr/bin/env python3
"""
Check that the genes MAGMA reports per set in a .gsa.sets.genes.out file match the
genes in the gene-set file that was tested.

Used to catch the MAGMA quirk where overriding the gene statistic with a replacement
.genes.out file (e.g. PoPS-based gene-set analysis) gives wrong set members when the
gene lists don't overlap perfectly.

Per set it reports:
  extra    genes MAGMA lists that are NOT in the set file          -> always an error
  missing  genes in the set file that MAGMA does not list
  n_header NGENES in the set's header vs the number of gene rows   -> mismatch is an error

A missing gene is expected if it simply wasn't in the gene analysis (e.g. chrX genes, or
IDs not in the gene annotation). Pass gene_results (the .genes.out or .genes.raw used) to
separate those from unexplained misses, which are then also counted as errors.

Usage (command line):
  python check_gsa_genes.py results.gsa.sets.genes.out sets.txt [gene_results.genes.out]
"""
import sys

import pandas as pd


def _read_gsa_set_genes(path):
    """Parse a MAGMA .gsa.sets.genes.out file -> {set_name: (ngenes_header, [genes])}."""
    names, ngenes, genes = {}, {}, {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "#" and len(p) >= 4 and p[1].startswith("_SET"):
                key, field = p[1], p[2]
                value = line.split("=", 1)[1].strip()
                if field == "VARIABLE":
                    names[key] = value[:-len("(set)")].strip() if value.endswith("(set)") else value
                elif field == "NGENES":
                    ngenes[key] = int(value)
            elif p[0].startswith("_SET") and len(p) > 1 and p[1] != "GENE":
                genes.setdefault(p[0], []).append(p[1])
    return {names[k]: (ngenes.get(k), genes.get(k, [])) for k in names}


def _read_set_file(path):
    """Parse a MAGMA gene-set file (name<TAB>gene<TAB>gene...) -> {set_name: set(genes)}."""
    sets = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 2 or line.startswith("#"):
                continue
            sets.setdefault(p[0], set()).update(g for g in p[1:] if g)
    return sets


def _read_gene_results(path):
    """Gene IDs from a MAGMA .genes.out (GENE column) or .genes.raw (first field)."""
    ids = set()
    with open(path) as f:
        first = f.readline().split()
        if "GENE" in first:                      # .genes.out with header
            col = first.index("GENE")
            for line in f:
                p = line.split()
                if p:
                    ids.add(p[col])
        else:                                    # .genes.raw: '#' comment lines, ID first
            for line in [" ".join(first)] + list(f):
                p = line.split()
                if p and not p[0].startswith("#"):
                    ids.add(p[0])
    return ids


def check_gsa_set_genes(gsa_genes_out, set_file, gene_results=None, verbose=True):
    """
    Compare per-set genes in a MAGMA .gsa.sets.genes.out file with the gene-set file.

    Returns (ok, table): ok is True when no set has errors; table is a DataFrame with
    one row per set reported in the .gsa.sets.genes.out file.
    """
    reported = _read_gsa_set_genes(gsa_genes_out)
    expected = _read_set_file(set_file)
    analysed = _read_gene_results(gene_results) if gene_results else None

    rows = []
    for name, (n_header, genes) in reported.items():
        got = set(genes)
        if name not in expected:
            rows.append(dict(set=name, n_set_file=None, n_reported=len(got), n_header=n_header,
                             extra=sorted(got), missing=[], missing_unexplained=None,
                             status="FAIL: set not in set file"))
            continue
        exp = expected[name]
        extra = sorted(got - exp)
        missing = sorted(exp - got)
        unexplained = sorted(g for g in missing if g in analysed) if analysed is not None else None

        problems = []
        if extra:
            problems.append(f"{len(extra)} gene(s) not in set file")
        if n_header is not None and n_header != len(genes):
            problems.append(f"NGENES={n_header} but {len(genes)} rows")
        if len(genes) != len(got):
            problems.append("duplicate gene rows")
        if unexplained:
            problems.append(f"{len(unexplained)} analysed gene(s) missing")

        rows.append(dict(set=name, n_set_file=len(exp), n_reported=len(got), n_header=n_header,
                         extra=extra, missing=missing, missing_unexplained=unexplained,
                         status="FAIL: " + "; ".join(problems) if problems else "OK"))

    table = pd.DataFrame(rows)
    ok = bool(len(table)) and not table["status"].str.startswith("FAIL").any()

    if verbose:
        n_fail = int(table["status"].str.startswith("FAIL").sum()) if len(table) else 0
        print(f"Checked {len(table)} set(s) from {gsa_genes_out}")
        print(f"  {len(table) - n_fail} OK, {n_fail} FAIL")
        for r in table.itertuples():
            if r.status != "OK":
                print(f"  {r.set}: {r.status}")
                if r.extra:
                    print(f"    not in set file: {', '.join(r.extra)}")
                if r.missing_unexplained:
                    print(f"    analysed but missing: {', '.join(r.missing_unexplained)}")
        if analysed is None and len(table) and table["missing"].str.len().sum():
            print(f"  note: {int(table['missing'].str.len().sum())} set-file gene(s) are absent from the output; "
                  "pass gene_results to check whether they were analysed")
        if len(table) == 0:
            print("  no sets found in the .gsa.sets.genes.out file")
    return ok, table


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    ok, _ = check_gsa_set_genes(*sys.argv[1:])
    sys.exit(0 if ok else 1)
