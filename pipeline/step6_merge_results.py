"""
Step 6 — Merge InterProScan results and enforce canonical operon logic.

Detection precedence:
  1. Direct DIAMOND homology (has_mazE == True from Step 2)
  2. InterProScan RHH domain (interpro_rhh_detected) — structural rescue

After merging, re-evaluates each operon's structure and sets is_canonical:

    is_canonical = True  iff
        combined_mazE_detected
        AND  MazE is at position -1 relative to MazF
        AND  alr is at position -2 relative to MazF

Empty genome rows (no MazF) are passed through unchanged.

Output: pipeline_output/mazF_operons_analysis_complete.csv

Usage:
    python pipeline/step6_merge_results.py [--force]
"""

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    OPERON_CSV, COMPLETE_CSV, OUTPUT_DIR,
    MAZE_CANONICAL_POSITION, ALR_CANONICAL_POSITION,
)
from pipeline.step5_interproscan import InterproResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _short_seq_id_from_interpro(interpro_seq_id: str) -> str:
    """
    Strip the leading refseq prefix from an InterProScan sequence_id.

    'GCF_014287895.1_NZ_JACOPP010000019.1_30' → 'NZ_JACOPP010000019.1_30'
    """
    for prefix in ("_NC_", "_NZ_", "_CP_", "_AP_", "_LN_"):
        idx = interpro_seq_id.find(prefix)
        if idx != -1:
            return interpro_seq_id[idx + 1:]
    parts = interpro_seq_id.split("_", 2)
    return "_".join(parts[2:]) if len(parts) >= 3 else interpro_seq_id


def _get_alr_position(adjacent_genes_data: str) -> int | None:
    """Return the position of alr relative to MazF, or None if absent."""
    try:
        genes = json.loads(adjacent_genes_data)
    except (json.JSONDecodeError, TypeError):
        return None
    for g in genes:
        if g.get("gene_name", "").lower() == "alr":
            return int(g.get("position", 0))
    return None


def _reclassify_operon(
    has_mazE_direct: bool,
    interpro_rhh: bool,
    interpro_position: int | None,
    alr_position: int | None,
) -> tuple[str, str, bool]:
    """
    Return (operon_structure, mazE_detection_method, is_canonical).

    Detection precedence: direct DIAMOND hit > InterProScan RHH.
    """
    maze_detected = has_mazE_direct or interpro_rhh
    maze_method   = ""
    maze_position = None

    if has_mazE_direct:
        maze_method   = "homology"
        maze_position = MAZE_CANONICAL_POSITION   # DIAMOND hit is always at -1
    elif interpro_rhh and interpro_position is not None:
        maze_method   = "InterPro"
        maze_position = interpro_position

    # Build operon structure tokens
    upstream_tokens   = []
    downstream_tokens = []

    if alr_position is not None:
        token = "alr"
        if alr_position < 0:
            upstream_tokens.append((alr_position, token))
        else:
            downstream_tokens.append((alr_position, token))

    if maze_detected and maze_position is not None:
        token = "mazE"
        if maze_position < 0:
            upstream_tokens.append((maze_position, token))
        else:
            downstream_tokens.append((maze_position, token))

    upstream_tokens.sort(key=lambda x: x[0])
    downstream_tokens.sort(key=lambda x: x[0])

    left  = [t for _, t in upstream_tokens]
    right = [t for _, t in downstream_tokens]

    if left and right:
        structure = "-".join(left) + "-mazF-" + "-".join(right)
    elif left:
        structure = "-".join(left) + "-mazF"
    elif right:
        structure = "mazF-" + "-".join(right)
    else:
        structure = "mazF_only"

    is_canonical = (
        maze_detected
        and maze_position == MAZE_CANONICAL_POSITION
        and alr_position  == ALR_CANONICAL_POSITION
    )

    return structure, maze_method, is_canonical


# ---------------------------------------------------------------------------
# Core merge
# ---------------------------------------------------------------------------

def run(
    interpro_results: dict[str, InterproResult],
    force: bool = False,
) -> None:
    if os.path.exists(COMPLETE_CSV) and not force:
        print(f"  Step 6 output already exists ({COMPLETE_CSV}). Skipping.")
        return

    if not os.path.exists(OPERON_CSV):
        print(f"ERROR: {OPERON_CSV} not found. Run Step 2 first.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(OPERON_CSV, low_memory=False)
    total = len(df)
    print(f"\nStep 6: Merging results for {total} rows ...")

    # Build lookup: short_seq_id → InterproResult
    interpro_lookup: dict[str, InterproResult] = {}
    for full_id, res in interpro_results.items():
        short = _short_seq_id_from_interpro(full_id)
        interpro_lookup[short] = res

    interpro_detected_col = []
    interpro_position_col = []
    interpro_seq_id_col   = []
    interpro_domains_col  = []
    combined_col          = []
    maze_method_col       = []
    structure_col         = []
    is_canonical_col      = []

    for _, row in df.iterrows():
        # Empty rows (no MazF) — pass through with all-False/empty values
        mazf_seq = str(row.get("mazF_sequence", ""))
        if not mazf_seq or mazf_seq == "nan":
            interpro_detected_col.append(False)
            interpro_position_col.append(None)
            interpro_seq_id_col.append("")
            interpro_domains_col.append("")
            combined_col.append(False)
            maze_method_col.append("")
            structure_col.append("")
            is_canonical_col.append(False)
            continue

        refseq     = str(row.get("refseq", ""))
        has_maze   = row.get("has_mazE") == True
        prev_qid   = str(row.get("prev_orf_query_id", ""))
        next_qid   = str(row.get("next_orf_query_id", ""))
        alr_pos    = _get_alr_position(str(row.get("adjacent_genes_data", "[]")))
        cur_struct = str(row.get("operon_structure", "mazF_only"))

        # InterProScan lookup (only for operons without a direct DIAMOND MazE hit)
        ipro_hit     = None
        ipro_pos     = None
        ipro_full_id = ""

        if not has_maze:
            for qid, pos in [(prev_qid, -1), (next_qid, 1)]:
                if not qid:
                    continue
                candidate = interpro_lookup.get(qid)
                if candidate is None:
                    full_key  = f"{refseq}_{qid}"
                    candidate = interpro_lookup.get(full_key)
                if candidate and candidate.is_valid_maze:
                    ipro_hit     = candidate
                    ipro_pos     = pos
                    ipro_full_id = candidate.seq_id
                    break

        ipro_ok = ipro_hit is not None and ipro_hit.is_valid_maze

        new_structure, maze_method, canonical = _reclassify_operon(
            has_mazE_direct   = has_maze,
            interpro_rhh      = ipro_ok,
            interpro_position = ipro_pos if ipro_ok else None,
            alr_position      = alr_pos,
        )

        combined = has_maze or ipro_ok

        interpro_detected_col.append(ipro_ok)
        interpro_position_col.append(ipro_pos if ipro_ok else None)
        interpro_seq_id_col.append(ipro_full_id if ipro_ok else "")
        interpro_domains_col.append(ipro_hit.domains_summary() if ipro_ok else "")
        combined_col.append(combined)
        maze_method_col.append(maze_method)
        structure_col.append(new_structure)
        is_canonical_col.append(canonical)

    df["interpro_rhh_detected"]    = interpro_detected_col
    df["interpro_rhh_position"]    = interpro_position_col
    df["interpro_rhh_sequence_id"] = interpro_seq_id_col
    df["interpro_rhh_domains"]     = interpro_domains_col
    df["combined_mazE_detected"]   = combined_col
    df["mazE_detection_method"]    = maze_method_col
    df["operon_structure"]         = structure_col
    df["is_canonical"]             = is_canonical_col

    df.to_csv(COMPLETE_CSV, index=False)

    with_mazf = df["mazF_sequence"].notna().sum()
    print(f"\n  Total rows:              {total}")
    print(f"  Rows with MazF:          {with_mazf}")
    print(f"  Combined MazE detection: {df['combined_mazE_detected'].sum()} operons")
    has_maze_bool = df["has_mazE"].fillna(False).astype(bool)
    print(f"    via DIAMOND homology:    {has_maze_bool.sum()}")
    print(f"    via InterProScan (new):  {(df['interpro_rhh_detected'] & ~has_maze_bool).sum()}")
    print(f"\n  Canonical operons [alr-mazE-mazF]: {df['is_canonical'].sum()}")
    print(f"\n  Operon structure counts (MazF-positive rows):")
    for structure, count in df["operon_structure"].value_counts().items():
        if structure:
            print(f"    {structure:<35} {count}")
    print(f"\n  Output: {COMPLETE_CSV}")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Merge InterProScan results with operon CSV."
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    from pipeline.step5_interproscan import run as ipro_run
    interpro_results = ipro_run()
    run(interpro_results, force=args.force)


if __name__ == "__main__":
    main()
