#!/usr/bin/env python3
"""
make_ind_labels.py

Build a labelled eigenstrat .ind_ file from:
  - the raw .ind file produced by convertVCFtoEigenstrat.sh
  - name2id (sample-level reference-group mapping for ancient samples)
  - ref_pop_ids_mapped_grouped (1000G modern population labels)
  - the sample metadata TSV (groupLabel for everything else)

Label priority per sample:
  1. name2id ref_group  (excluding 'allAncients' and 'African' --
                         we drop 'African' so those samples retain their
                         granular groupLabel and can be tested individually
                         as candidate ancient African source populations)
  2. 1000G ref pop      (CEU, GBR, TSI, IBS, FIN, YRI, ...)
  3. metadata groupLabel
  4. existing pop column from the .ind file (fallback)

Sample IDs are matched on the prefix before the first '.'. This makes the
join robust to study-suffix drift between callsets (e.g. allentoft_2015 vs
allentoft_2024) without needing to know which suffix variant lives where.
"""

import argparse
import sys
import pandas as pd


def base_id(s):
    return s.split(".", 1)[0] if isinstance(s, str) else s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ind", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--name2id", required=True)
    ap.add_argument("--ref-pop", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # --- input .ind ----------------------------------------------------------
    ind = pd.read_csv(
        args.ind, sep=r"\s+", header=None,
        names=["sample_id", "sex", "pop"], dtype=str)
    ind["base_id"] = ind["sample_id"].map(base_id)

    # --- metadata: groupLabel keyed on base_id -------------------------------
    meta = pd.read_csv(args.meta, sep="\t", dtype=str, na_filter=False)
    meta["base_id"] = meta["sampleId"].map(base_id)
    meta = (meta[["base_id", "groupLabel"]]
            .replace("", pd.NA)
            .dropna(subset=["groupLabel"])
            .drop_duplicates(subset="base_id", keep="first"))

    # --- name2id: drop allAncients + African, keep first hit per base_id -----
    n2i = pd.read_csv(
        args.name2id, sep=r"\s+", header=None,
        names=["n2i_id", "ref_group", "detail"], dtype=str)
    n2i = n2i[~n2i["ref_group"].isin(["allAncients", "African"])].copy()
    n2i["base_id"] = n2i["n2i_id"].map(base_id)
    n2i = n2i.drop_duplicates(subset="base_id", keep="first")[["base_id", "ref_group"]]

    # --- 1000G modern ref pop ------------------------------------------------
    refpop = pd.read_csv(args.ref_pop, header=None, dtype=str)
    refpop.columns = ["sample"] + [f"c{i}" for i in range(1, refpop.shape[1])]
    # last non-NA column is taken as the population label
    refpop["pop_1000g"] = refpop.iloc[:, 1:].bfill(axis=1).iloc[:, 0]
    refpop["base_id"] = refpop["sample"].map(base_id)
    refpop = (refpop[["base_id", "pop_1000g"]]
              .dropna()
              .drop_duplicates(subset="base_id", keep="first"))

    # --- merge & resolve label by priority -----------------------------------
    df = (ind
          .merge(n2i,    on="base_id", how="left")
          .merge(refpop, on="base_id", how="left")
          .merge(meta,   on="base_id", how="left"))

    # The .ind "pop" column is a placeholder ('?') for these inputs, so treat
    # it as a genuine label only if it is something other than '?'. This keeps
    # '?' as a true last resort rather than letting it mask join failures.
    ind_pop = df["pop"].where(df["pop"].ne("?"))
    df["final_pop"] = (df["ref_group"]
                       .fillna(df["pop_1000g"])
                       .fillna(df["groupLabel"])
                       .fillna(ind_pop))

    # Which rule actually resolved each sample (highest priority that hit).
    src = pd.Series("unresolved (-> '?')", index=df.index)
    src = src.mask(ind_pop.notna(),         ".ind pop column")
    src = src.mask(df["groupLabel"].notna(), "metadata groupLabel")
    src = src.mask(df["pop_1000g"].notna(),  "1000G ref pop")
    src = src.mask(df["ref_group"].notna(),  "name2id ref_group")

    unresolved = df["final_pop"].isna()

    # --- write (fall back to '?' only for genuinely unresolved samples) ------
    out = df[["sample_id", "sex", "final_pop"]].copy()
    out["final_pop"] = out["final_pop"].fillna("?")
    out.to_csv(args.out, sep=" ", header=False, index=False)

    # --- summary -------------------------------------------------------------
    print(f"  Total samples:              {len(df)}", file=sys.stderr)
    for rule, n in src.value_counts().items():
        print(f"    resolved by {rule:<22} {n}", file=sys.stderr)
    print(f"  Unresolved (written '?'):   {unresolved.sum()}", file=sys.stderr)
    print(f"  Distinct populations:       {out['final_pop'].nunique()}", file=sys.stderr)

    # Break the unresolved set down by ID pattern -- 1000G samples look like
    # HG#####/NA##### and should have been caught by the 1000G ref pop join.
    if unresolved.any():
        un = df.loc[unresolved, "sample_id"].astype(str)
        n_1000g = un.str.match(r"^(HG|NA)\d").sum()
        print(f"    of which look like 1000G ids (HG*/NA*): {n_1000g}",
              file=sys.stderr)
        print(f"    examples: {', '.join(un.head(5))}", file=sys.stderr)


if __name__ == "__main__":
    main()
