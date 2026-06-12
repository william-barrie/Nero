#!/usr/bin/env python3
"""
make_ind_labels.py

Build a labelled eigenstrat .ind_ file from several sources, applying a fixed
label priority per sample.

Inputs
------
  --ind       Raw .ind from vcf_to_eigenstrat.sh (sample, sex, pop='?').
  --name2id   Whitespace table: <sample> <ref_group> [detail]. Sample-level
              ancient reference groups. 'allAncients' and 'African' are dropped
              so those samples keep a granular label (and African ancients can
              be tested individually as candidate sources).
  --ref-pop   Comma table from the ANCIENT reference panel, e.g.
                  RISE240,Yamnaya6,Yamnaya,1
              i.e. <sample>,<detailed>,<grouped>,<n>. The GROUPED label
              (0-based column 2 by default) is used: WHG / EHG / CHG /
              FarmerAnatolian / FarmerEarly / FarmerMiddle / FarmerLate /
              Yamnaya. This file is NOT the 1000G panel.
  --igsr      IGSR / 1000G sample TSV (igsr_samples.tsv) with columns
              'Sample name' and 'Population code'. Source of the MODERN 1000G
              population labels: CEU, GBR, TSI, IBS, FIN, YRI, ...
  --meta      sampleInfo TSV; 'groupLabel' keyed on 'sampleId'. Fallback for
              everything else, including the ancient African candidate sources.
  --out       Output .ind_ path.

Label priority (highest first)
------------------------------
  1. name2id ref_group        (curated ancient reference groups)
  2. ref-pop grouped label     (ancient reference panel, grouped)
  3. igsr Population code       (modern 1000G)
  4. metadata groupLabel        (everything else, incl. ancient African sources)
  5. '?'                        (genuinely unresolved; warned about)

All joins are on the ID prefix before the first '.', making them robust to
study-suffix drift between callsets (e.g. allentoft_2015 vs allentoft_2024).
"""

import argparse
import os
import sys
import pandas as pd

# Populations the downstream f4 analysis depends on. Used only for a presence
# warning; override with --expect if the analysis config changes.
EXPECTED_POPS = [
    "CEU", "GBR", "TSI", "IBS", "FIN", "YRI",                # 1000G targets + sister
    "WHG", "EHG", "CHG", "FarmerAnatolian", "FarmerEarly",   # ancient references
    "FarmerMiddle", "FarmerLate", "Yamnaya",
]


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def base_id(s):
    return s.split(".", 1)[0] if isinstance(s, str) else s


def require_file(path, what):
    if not path or not os.path.exists(path):
        die(f"{what} file not found: {path!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ind", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--name2id", required=True)
    ap.add_argument("--ref-pop", required=True)
    ap.add_argument("--igsr", required=True,
                    help="IGSR/1000G sample TSV (igsr_samples.tsv)")
    ap.add_argument("--ref-pop-grouped-col", type=int, default=2,
                    help="0-based column in --ref-pop holding the grouped "
                         "label (default 2)")
    ap.add_argument("--expect", default=None,
                    help="comma-separated populations to check for presence "
                         "(default: the f4 analysis set)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    for path, what in [(args.ind, "--ind"), (args.meta, "--meta"),
                       (args.name2id, "--name2id"),
                       (args.ref_pop, "--ref-pop"), (args.igsr, "--igsr")]:
        require_file(path, what)

    # --- input .ind ----------------------------------------------------------
    ind = pd.read_csv(args.ind, sep=r"\s+", header=None,
                      names=["sample_id", "sex", "pop"], dtype=str)
    ind["base_id"] = ind["sample_id"].map(base_id)

    # --- (1) name2id: ancient ref groups; drop allAncients + African ---------
    n2i = pd.read_csv(args.name2id, sep=r"\s+", header=None,
                      names=["n2i_id", "ref_group", "detail"], dtype=str)
    n2i = n2i[~n2i["ref_group"].isin(["allAncients", "African"])]
    n2i["base_id"] = n2i["n2i_id"].map(base_id)
    n2i = (n2i.dropna(subset=["ref_group"])
              .drop_duplicates(subset="base_id", keep="first")
              [["base_id", "ref_group"]])

    # --- (2) ref-pop: ancient reference panel, GROUPED label -----------------
    refpop_raw = pd.read_csv(args.ref_pop, header=None, dtype=str)
    gcol = args.ref_pop_grouped_col
    if refpop_raw.shape[1] <= gcol:
        die(f"--ref-pop has {refpop_raw.shape[1]} column(s); expected the "
            f"grouped label at 0-based column {gcol}. First row: "
            f"{refpop_raw.iloc[0].tolist()!r}. Is the file comma-separated "
            f"like 'RISE240,Yamnaya6,Yamnaya,1'?")
    refpop = (pd.DataFrame({
                  "base_id": refpop_raw.iloc[:, 0].map(base_id),
                  "ref_pop_grouped": refpop_raw.iloc[:, gcol]})
              .replace("", pd.NA).dropna()
              .drop_duplicates(subset="base_id", keep="first"))

    # --- (3) igsr: modern 1000G population codes -----------------------------
    igsr_raw = pd.read_csv(args.igsr, sep="\t", dtype=str)
    need = {"Sample name", "Population code"}
    have = set(igsr_raw.columns)
    if need - have:
        die(f"--igsr missing expected column(s): {sorted(need - have)}. "
            f"Found: {list(igsr_raw.columns)}. Expected the IGSR sample TSV "
            f"(igsr_samples.tsv) with 'Sample name' and 'Population code'.")
    igsr = (pd.DataFrame({
                "base_id": igsr_raw["Sample name"].map(base_id),
                "pop_1000g": igsr_raw["Population code"]})
            .replace("", pd.NA).dropna()
            .drop_duplicates(subset="base_id", keep="first"))

    # --- (4) metadata groupLabel ---------------------------------------------
    meta = pd.read_csv(args.meta, sep="\t", dtype=str, na_filter=False)
    if {"sampleId", "groupLabel"} - set(meta.columns):
        die("--meta must contain 'sampleId' and 'groupLabel'; found: "
            f"{list(meta.columns)}")
    meta["base_id"] = meta["sampleId"].map(base_id)
    meta = (meta[["base_id", "groupLabel"]]
            .replace("", pd.NA).dropna(subset=["groupLabel"])
            .drop_duplicates(subset="base_id", keep="first"))

    # --- merge & resolve label by priority -----------------------------------
    df = (ind
          .merge(n2i,    on="base_id", how="left")
          .merge(refpop, on="base_id", how="left")
          .merge(igsr,   on="base_id", how="left")
          .merge(meta,   on="base_id", how="left"))

    # The .ind "pop" is a placeholder ('?'); treat it as a real label only if
    # it is something other than '?', so '?' stays a true last resort.
    ind_pop = df["pop"].where(df["pop"].ne("?"))
    df["final_pop"] = (df["ref_group"]
                       .fillna(df["ref_pop_grouped"])
                       .fillna(df["pop_1000g"])
                       .fillna(df["groupLabel"])
                       .fillna(ind_pop))

    # Which rule actually resolved each sample (highest priority that hit).
    src = pd.Series("unresolved (-> '?')", index=df.index)
    src = src.mask(ind_pop.notna(),               ".ind pop column")
    src = src.mask(df["groupLabel"].notna(),       "metadata groupLabel")
    src = src.mask(df["pop_1000g"].notna(),        "1000G ref pop (igsr)")
    src = src.mask(df["ref_pop_grouped"].notna(),  "ref-pop grouped")
    src = src.mask(df["ref_group"].notna(),        "name2id ref_group")

    unresolved = df["final_pop"].isna()

    # --- write (fall back to '?' only for genuinely unresolved samples) ------
    out = df[["sample_id", "sex", "final_pop"]].copy()
    out["final_pop"] = out["final_pop"].fillna("?")
    out.to_csv(args.out, sep=" ", header=False, index=False)

    # --- summary -------------------------------------------------------------
    print(f"  Total samples:              {len(df)}", file=sys.stderr)
    for rule, n in src.value_counts().items():
        print(f"    resolved by {rule:<24} {n}", file=sys.stderr)
    print(f"  Unresolved (written '?'):   {int(unresolved.sum())}",
          file=sys.stderr)
    print(f"  Distinct populations:       {out['final_pop'].nunique()}",
          file=sys.stderr)

    # --- warnings ------------------------------------------------------------
    looks_1000g = df["sample_id"].astype(str).str.match(r"^(HG|NA)\d")
    n_1000g     = int(looks_1000g.sum())
    n_1000g_lab = int((looks_1000g & df["pop_1000g"].notna()).sum())
    if n_1000g and n_1000g_lab == 0:
        warn(f"{n_1000g} samples look like 1000G ids (HG*/NA*) but NONE matched "
             f"--igsr. Check igsr_samples.tsv covers them and that 'Sample "
             f"name' matches the .ind ids.")
    elif n_1000g and n_1000g_lab < n_1000g:
        warn(f"{n_1000g - n_1000g_lab}/{n_1000g} 1000G-looking samples were "
             f"not matched in --igsr.")

    if unresolved.any():
        un = df.loc[unresolved, "sample_id"].astype(str)
        n_un_1000g = int(un.str.match(r"^(HG|NA)\d").sum())
        examples = ", ".join(un.head(5))
        warn(f"{int(unresolved.sum())} samples unresolved (written '?'); "
             f"{n_un_1000g} look like 1000G ids. examples: {examples}")

    expect = ([p.strip() for p in args.expect.split(",") if p.strip()]
              if args.expect else EXPECTED_POPS)
    present = set(out["final_pop"].unique())
    miss = [p for p in expect if p not in present]
    if miss:
        warn("expected analysis populations absent from labels: "
             + ", ".join(miss))
    else:
        print("  All expected analysis populations present.", file=sys.stderr)


if __name__ == "__main__":
    main()
