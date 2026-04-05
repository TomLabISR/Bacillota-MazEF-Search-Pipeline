"""
Step 7 — Add NCBI taxonomy to the operon table.

Sources (in priority order for each field):
  1. data_summary.tsv                   → Organism name  (all 3,362 genomes)
  2. MazF_operons_analysis_taxonomy.csv → Phylum/Class/Order/Family/Genus
                                          (pre-existing, 2,849 unique refseqs)
  3. data/16s_cutsite.tsv               → fallback lineage from 16S analysis
                                          (2,211 unique refseqs)

Output: pipeline_output/MazF_operons_analysis_taxonomy.csv

Usage:
    python pipeline/step7_add_taxonomy.py [--force]
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    COMPLETE_CSV, TAXONOMY_CSV, OUTPUT_DIR, PROJECT_ROOT,
    CUTSITE_TSV, PREEXISTING_TAX,
)

# Path to NCBI assembly summary (organism names)
DATA_SUMMARY = os.path.join(
    PROJECT_ROOT, "Bacillota_ref_dataset", "ncbi_dataset",
    "data", "data_summary.tsv"
)

LINEAGE_COLS = ["Phylum", "Class", "Order", "Family", "Genus"]


def _load_organism_names() -> dict[str, str]:
    """
    Load {Assembly Accession → Organism Scientific Name} from data_summary.tsv.
    Covers all 3,362 genomes.
    """
    if not os.path.exists(DATA_SUMMARY):
        print(f"  WARNING: {DATA_SUMMARY} not found — organism names will be empty.")
        return {}

    df = pd.read_csv(DATA_SUMMARY, sep="\t", dtype=str,
                     usecols=["Assembly Accession", "Organism Scientific Name"])
    df = df.dropna(subset=["Assembly Accession"])
    mapping = dict(zip(df["Assembly Accession"].str.strip(),
                       df["Organism Scientific Name"].str.strip()))
    print(f"  Organism names loaded: {len(mapping)} from data_summary.tsv")
    return mapping


def _load_lineage_lookup() -> dict[str, dict]:
    """
    Build {refseq → {Phylum, Class, Order, Family, Genus}} from available sources.
    Tries pre-existing taxonomy CSV first, then 16s_cutsite.tsv as fallback.
    """
    lookup: dict[str, dict] = {}

    # Source 1: pre-existing MazF_operons_analysis_taxonomy.csv
    if os.path.exists(PREEXISTING_TAX):
        df = pd.read_csv(PREEXISTING_TAX, low_memory=False,
                         usecols=["refseq"] + LINEAGE_COLS)
        df = df.dropna(subset=["refseq"]).drop_duplicates("refseq")
        for _, row in df.iterrows():
            refseq = str(row["refseq"]).strip()
            if not any(str(row.get(c, "")).strip() for c in LINEAGE_COLS):
                continue  # skip empty rows
            lookup[refseq] = {c: str(row.get(c, "")).strip() for c in LINEAGE_COLS}
        print(f"  Lineage from pre-existing taxonomy CSV: {len(lookup)} refseqs")
    else:
        print(f"  WARNING: Pre-existing taxonomy CSV not found: {PREEXISTING_TAX}")

    # Source 2: 16s_cutsite.tsv (fills gaps)
    if os.path.exists(CUTSITE_TSV):
        cs_df = pd.read_csv(CUTSITE_TSV, sep="\t", dtype=str,
                             low_memory=False)
        cs_df.columns = [c.strip() for c in cs_df.columns]
        refseq_col = next(
            (c for c in cs_df.columns if c.lower() in ("refseq", "assembly")), None
        )
        if refseq_col:
            before = len(lookup)
            for _, row in cs_df.iterrows():
                refseq = str(row.get(refseq_col, "")).strip()
                if not refseq or refseq in lookup:
                    continue
                entry = {c: str(row.get(c, "")).strip() for c in LINEAGE_COLS
                         if c in cs_df.columns}
                if any(entry.values()):
                    lookup[refseq] = entry
            print(f"  Lineage from 16s_cutsite.tsv: +{len(lookup)-before} additional refseqs")
    else:
        print(f"  NOTE: 16s_cutsite.tsv not found at {CUTSITE_TSV} — skipping.")

    print(f"  Total lineage lookup: {len(lookup)} refseqs")
    return lookup


def run(force: bool = False) -> None:
    if os.path.exists(TAXONOMY_CSV) and not force:
        print(f"  Step 7 output already exists ({TAXONOMY_CSV}). Skipping.")
        return

    if not os.path.exists(COMPLETE_CSV):
        print(f"ERROR: {COMPLETE_CSV} not found. Run Step 6 first.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    operon_df = pd.read_csv(COMPLETE_CSV, low_memory=False)
    print(f"\nStep 7: Adding taxonomy to {len(operon_df)} operon records ...")

    org_names  = _load_organism_names()
    lin_lookup = _load_lineage_lookup()

    # Build new columns
    organisms  = []
    phyla      = []
    classes    = []
    orders     = []
    families   = []
    genera     = []

    for refseq in operon_df["refseq"].astype(str):
        refseq = refseq.strip()
        organism = org_names.get(refseq, "")
        organisms.append(organism)
        lin = lin_lookup.get(refseq, {})
        phyla.append(lin.get("Phylum", ""))
        classes.append(lin.get("Class", ""))
        orders.append(lin.get("Order", ""))
        families.append(lin.get("Family", ""))

        genus = lin.get("Genus", "")
        if not genus and organism:
            # Fallback: extract genus from the first word of the organism name.
            # Strip leading/trailing brackets (e.g. "[Eubacterium] cellulosolvens").
            first_word = organism.split()[0].strip("[]")
            if first_word[0].isupper():
                genus = first_word
        genera.append(genus)

    operon_df["Organism"] = organisms
    operon_df["Phylum"]   = phyla
    operon_df["Class"]    = classes
    operon_df["Order"]    = orders
    operon_df["Family"]   = families
    operon_df["Genus"]    = genera

    operon_df.to_csv(TAXONOMY_CSV, index=False)

    filled = sum(1 for p in phyla if p)
    print(f"\n  Operons with taxonomy: {filled} / {len(operon_df)}")
    if filled > 0:
        print(f"\n  Top phyla:")
        vc = pd.Series(phyla).value_counts()
        for phylum, cnt in vc[vc.index != ""].head(8).items():
            print(f"    {phylum:<35} {cnt}")

    print(f"\n  Output: {TAXONOMY_CSV}")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Add organism name and lineage to operon table."
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    run(force=args.force)


if __name__ == "__main__":
    main()
