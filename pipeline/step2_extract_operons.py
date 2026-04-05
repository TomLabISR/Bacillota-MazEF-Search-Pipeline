"""
Step 2 — Extract MazF operon context.

Reads every DIAMOND matches TSV from results/matches/ and records one row
per genome:
  - Genomes WITH a MazF hit: the FIRST MazF hit encountered (matching the
    original mazF_search.ipynb logic where `if ref not in operon:` meant
    only the first hit per genome was stored).
  - Genomes WITHOUT any MazF hit: a single empty row (refseq only) so that
    the final Operon_Summary includes all genomes in the dataset.

Output: pipeline_output/mazF_operons_analysis.csv

Usage:
    python pipeline/step2_extract_operons.py [--force]
"""

import glob
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    FAA_DIR, MATCHES_DIR, OPERON_CSV, OUTPUT_DIR,
    DIAMOND_COLS, MAZF_SUBJECT_ID, MAZE_SUBJECT_ID, ALR_SUBJECT_ID,
    MIN_MAZF_LENGTH,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_query_id(query_id: str) -> tuple[str, int]:
    """Split 'NC_003869.1_2082' → ('NC_003869.1', 2082)."""
    idx = query_id.rfind("_")
    fragment  = query_id[:idx]
    orf_index = int(query_id[idx + 1:])
    return fragment, orf_index


def _load_faa_sequences(refseq: str) -> dict[str, str]:
    """
    Load all protein sequences from results/faa/{refseq}.faa.

    Returns dict: { 'NC_003869.1_2082': 'MIVKRG...' }
    """
    faa_path = os.path.join(FAA_DIR, f"{refseq}.faa")
    if not os.path.exists(faa_path):
        return {}

    seqs: dict[str, str] = {}
    current_id = None
    current_seq: list[str] = []

    with open(faa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)

    return seqs


def _load_matches(tsv_path: str) -> pd.DataFrame | None:
    """Load one matches TSV into a DataFrame. Returns None if empty."""
    try:
        df = pd.read_csv(tsv_path, sep="\t", header=None, names=DIAMOND_COLS)
        if df.empty or df["Query_ID"].isna().all():
            return None
        return df
    except Exception:
        return None


def _classify_operon(adjacent: list[dict]) -> str:
    """Derive operon_structure from adjacent gene hits."""
    upstream_tokens   = []
    downstream_tokens = []

    for pos in sorted(set(g["position"] for g in adjacent if g["position"] < 0)):
        for g in adjacent:
            if g["position"] == pos:
                upstream_tokens.append(g["gene_name"])

    for pos in sorted(set(g["position"] for g in adjacent if g["position"] > 0)):
        for g in adjacent:
            if g["position"] == pos:
                downstream_tokens.append(g["gene_name"])

    if upstream_tokens and downstream_tokens:
        return "-".join(upstream_tokens) + "-mazF-" + "-".join(downstream_tokens)
    elif upstream_tokens:
        return "-".join(upstream_tokens) + "-mazF"
    elif downstream_tokens:
        return "mazF-" + "-".join(downstream_tokens)
    else:
        return "mazF_only"


# ---------------------------------------------------------------------------
# Core extraction — one record per genome
# ---------------------------------------------------------------------------

def _candidate_completeness(orf_idx: int, fragment: str, other_df: pd.DataFrame) -> int:
    """
    Score a candidate MazF hit by operon completeness.

    Returns:
      3  — both mazE (abs Δ==1) AND alr (abs Δ==2) adjacent
      2  — mazE adjacent only
      1  — alr adjacent only
      0  — lone mazF
    """
    frag_hits = other_df[other_df["Query_fragment"] == fragment].copy()
    frag_hits["delta"] = frag_hits["Query_orf_index"] - orf_idx
    nearby = frag_hits[frag_hits["delta"].abs() <= 2]

    has_mazE = any(
        r["Subject_ID"] == MAZE_SUBJECT_ID and abs(int(r["delta"])) == 1
        for _, r in nearby.iterrows()
    )
    has_alr = any(
        r["Subject_ID"] == ALR_SUBJECT_ID and abs(int(r["delta"])) <= 2
        for _, r in nearby.iterrows()
    )

    if has_mazE and has_alr:
        return 3
    if has_mazE:
        return 2
    if has_alr:
        return 1
    return 0


def _process_refseq(
    refseq: str,
    df: pd.DataFrame,
    faa_seqs: dict[str, str],
) -> dict | None:
    """
    Extract the operon record for one genome.

    Selection rule (reverse-engineered from original analysis):
      For each valid MazF hit, compute an operon-completeness score
      (alr+mazE+mazF=3 > mazE-mazF=2 > alr-mazF=1 > mazF_only=0).
      Pick the hit with the highest completeness score.
      Ties are broken by file order (first occurrence in the TSV).

    Returns a single record dict, or None if no valid MazF hit found.
    """
    df = df.copy()
    df[["Query_fragment", "Query_orf_index"]] = df["Query_ID"].apply(
        lambda qid: pd.Series(_parse_query_id(qid))
    )
    df["Identity"] = pd.to_numeric(df["Identity"], errors="coerce")
    df["E_value"]  = pd.to_numeric(df["E_value"],  errors="coerce")

    mazF_df  = df[df["Subject_ID"] == MAZF_SUBJECT_ID].copy()
    other_df = df[df["Subject_ID"] != MAZF_SUBJECT_ID].copy()

    if mazF_df.empty:
        return None

    # ---- Select best MazF hit: most complete operon, tie-break by file order ----
    best_score = -1
    mrow       = None
    for _, row in mazF_df.iterrows():
        qid = row["Query_ID"]
        seq = faa_seqs.get(qid, "")
        length = len(seq.rstrip("*"))
        if length > 0 and length < MIN_MAZF_LENGTH:
            continue  # skip fragments that are too short
        score = _candidate_completeness(
            int(row["Query_orf_index"]), row["Query_fragment"], other_df
        )
        if score > best_score:
            best_score = score
            mrow = row   # file order: first occurrence of each score wins

    if mrow is None:
        return None  # all hits too short

    mazf_qid   = mrow["Query_ID"]
    mazf_frag  = mrow["Query_fragment"]
    mazf_idx   = int(mrow["Query_orf_index"])
    mazf_ident = float(mrow["Identity"])

    mazf_seq = faa_seqs.get(mazf_qid, "")
    mazf_len = len(mazf_seq.rstrip("*"))

    operon_id = f"{refseq}_{mazf_frag}"

    # ---- Collect adjacent hits on the same fragment (|Δindex| ≤ 2) ---------
    same_frag = other_df[other_df["Query_fragment"] == mazf_frag].copy()
    same_frag["delta"] = same_frag["Query_orf_index"] - mazf_idx

    adjacent_genes: list[dict] = []
    has_mazE = False
    has_alr  = False

    for _, arow in same_frag[same_frag["delta"].abs() <= 2].iterrows():
        delta    = int(arow["delta"])
        gene_qid = arow["Query_ID"]
        gene_seq = faa_seqs.get(gene_qid, "")
        gene_len = len(gene_seq.rstrip("*"))
        subject  = arow["Subject_ID"]

        raw_name  = subject.split("_")[0]
        gene_name = "mazE" if raw_name.lower() == "ndoai" else raw_name

        adjacent_genes.append({
            "gene_name":  gene_name,
            "orf_index":  int(arow["Query_orf_index"]),
            "query_id":   gene_qid,
            "position":   delta,
            "identity":   str(round(float(arow["Identity"]), 1)),
            "length":     gene_len,
            "sequence":   gene_seq,
        })

        if subject == MAZE_SUBJECT_ID and abs(delta) == 1:
            has_mazE = True
        if subject == ALR_SUBJECT_ID and abs(delta) <= 2:
            has_alr = True

    operon_structure = _classify_operon(adjacent_genes)

    # ---- Prev / next ORF sequences (position ±1) for InterProScan ----------
    prev_qid = f"{mazf_frag}_{mazf_idx - 1}"
    next_qid = f"{mazf_frag}_{mazf_idx + 1}"
    prev_seq = faa_seqs.get(prev_qid, "")
    next_seq = faa_seqs.get(next_qid, "")

    return {
        "operon_id":             operon_id,
        "refseq":                refseq,
        "fragment":              mazf_frag,
        "mazF_orf_index":        mazf_idx,
        "mazF_identity":         mazf_ident,
        "mazF_query_id":         mazf_qid,
        "mazF_length":           mazf_len,
        "mazF_sequence":         mazf_seq,
        "adjacent_genes_count":  len(adjacent_genes),
        "adjacent_genes_data":   json.dumps(adjacent_genes),
        "prev_orf_query_id":     prev_qid if prev_seq else "",
        "prev_orf_length":       len(prev_seq.rstrip("*")) if prev_seq else 0,
        "prev_orf_sequence":     prev_seq,
        "next_orf_query_id":     next_qid if next_seq else "",
        "next_orf_length":       len(next_seq.rstrip("*")) if next_seq else 0,
        "next_orf_sequence":     next_seq,
        "has_mazE":              has_mazE,
        "has_alr":               has_alr,
        "operon_structure":      operon_structure,
    }


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def run(force: bool = False) -> None:
    if os.path.exists(OPERON_CSV) and not force:
        print(f"  Step 2 output already exists ({OPERON_CSV}). Skipping.")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(MATCHES_DIR, "*_matches.tsv")))
    total = len(tsv_files)
    if total == 0:
        print(f"ERROR: No matches TSVs found in {MATCHES_DIR}. Run Step 1 first.")
        sys.exit(1)

    print(f"\nStep 2: Extracting operons from {total} genomes "
          f"(one row per genome, first MazF hit only) ...")

    all_records: list[dict] = []
    mazf_found   = 0
    mazf_empty   = 0

    for i, tsv_path in enumerate(tsv_files, 1):
        basename = os.path.basename(tsv_path)
        refseq   = basename.replace("_matches.tsv", "")

        df = _load_matches(tsv_path)
        if df is None:
            # No DIAMOND hits at all — record empty row so genome stays in output
            all_records.append({"refseq": refseq})
            mazf_empty += 1
        else:
            faa_seqs = _load_faa_sequences(refseq)
            record   = _process_refseq(refseq, df, faa_seqs)
            if record is None:
                # Hits present but no valid MazF (all too short, etc.)
                all_records.append({"refseq": refseq})
                mazf_empty += 1
            else:
                all_records.append(record)
                mazf_found += 1

        if i % 200 == 0 or i == total:
            print(f"  Progress: {i}/{total} genomes  "
                  f"({mazf_found} with MazF, {mazf_empty} empty)")

    operon_df = pd.DataFrame(all_records)
    operon_df.to_csv(OPERON_CSV, index=False)

    with_mazf = operon_df["mazF_sequence"].notna().sum()
    print(f"\n  Total rows written:     {len(operon_df)}")
    print(f"  Genomes with MazF:      {with_mazf}")
    print(f"  Genomes without MazF:   {mazf_empty}")

    if with_mazf > 0:
        has_maze_count = operon_df["has_mazE"].sum()
        has_alr_count  = operon_df["has_alr"].sum()
        print(f"  Operons with direct MazE: {has_maze_count} "
              f"({100 * has_maze_count / with_mazf:.1f}%)")
        print(f"  Operons with Alr:         {has_alr_count} "
              f"({100 * has_alr_count  / with_mazf:.1f}%)")
        print(f"\n  Operon structure counts (MazF-positive genomes):")
        for structure, count in operon_df["operon_structure"].value_counts().items():
            print(f"    {structure:<30} {count}")

    print(f"\n  Output: {OPERON_CSV}")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Extract MazF operon context.")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    run(force=args.force)


if __name__ == "__main__":
    main()
