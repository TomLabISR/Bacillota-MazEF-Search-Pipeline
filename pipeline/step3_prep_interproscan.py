"""
Step 3 — Prepare InterProScan input FASTA.

Writes a FASTA file containing the adjacent ORF sequences (position -1 and
+1 relative to MazF) for all operons where no MazE was detected by direct
DIAMOND homology. This FASTA is submitted to InterProScan (Step 5) to look
for Ribbon-Helix-Helix (RHH) antitoxin domains.

Header format used in the output FASTA:
    >{refseq}_{fragment}_{orf_index}

This format matches what was used when generating the pre-computed
ORF_detection.tsv, ensuring compatibility during the merge step.

Output: pipeline_output/ORFs_adj_to_MazF.fasta

Usage:
    python pipeline/step3_prep_interproscan.py [--force]
"""

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import OPERON_CSV, ADJ_FASTA, OUTPUT_DIR


def run(force: bool = False) -> None:
    if os.path.exists(ADJ_FASTA) and not force:
        print(f"  Step 3 output already exists ({ADJ_FASTA}). Skipping.")
        return

    if not os.path.exists(OPERON_CSV):
        print(f"ERROR: {OPERON_CSV} not found. Run Step 2 first.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(OPERON_CSV, low_memory=False)
    total = len(df)
    print(f"\nStep 3: Preparing InterProScan input FASTA from {total} operons ...")

    # Only process operons without a direct DIAMOND MazE hit
    no_maze = df[df["has_mazE"] != True].copy()
    print(f"  Operons without direct MazE hit: {len(no_maze)}")

    written = 0
    seen: set[str] = set()  # deduplicate sequences

    with open(ADJ_FASTA, "w") as out:
        for _, row in no_maze.iterrows():
            refseq = str(row.get("refseq", ""))

            # Position -1 (prev_orf)
            prev_qid = str(row.get("prev_orf_query_id", ""))
            prev_seq = str(row.get("prev_orf_sequence", "")).strip()

            if prev_qid and prev_seq and prev_seq not in ("", "nan"):
                # Reconstruct full header: {refseq}_{fragment}_{orf_index}
                header = f"{refseq}_{prev_qid}" if not prev_qid.startswith(refseq) else prev_qid
                if header not in seen:
                    seen.add(header)
                    # Strip trailing stop codon '*' for InterProScan
                    clean_seq = prev_seq.rstrip("*")
                    out.write(f">{header}\n")
                    # Write in 60-char lines
                    for i in range(0, len(clean_seq), 60):
                        out.write(clean_seq[i:i+60] + "\n")
                    written += 1

            # Position +1 (next_orf) — also included in case MazE is downstream
            next_qid = str(row.get("next_orf_query_id", ""))
            next_seq = str(row.get("next_orf_sequence", "")).strip()

            if next_qid and next_seq and next_seq not in ("", "nan"):
                header = f"{refseq}_{next_qid}" if not next_qid.startswith(refseq) else next_qid
                if header not in seen:
                    seen.add(header)
                    clean_seq = next_seq.rstrip("*")
                    out.write(f">{header}\n")
                    for i in range(0, len(clean_seq), 60):
                        out.write(clean_seq[i:i+60] + "\n")
                    written += 1

    print(f"  Sequences written: {written} (unique ORFs)")
    print(f"\n  Output: {ADJ_FASTA}")
    print(
        "\n  NOTE: To run InterProScan from scratch (optional):\n"
        "    interproscan.sh -i pipeline_output/ORFs_adj_to_MazF.fasta \\\n"
        "        -f tsv -o detecting_mazE/ORF_detection.tsv \\\n"
        "        -appl Pfam,SUPERFAMILY,SMART,ProSiteProfiles \\\n"
        "        --goterms --pathways\n"
        "  The pipeline defaults to the pre-computed TSV on Zenodo."
    )


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Write adjacent-ORF FASTA for InterProScan."
    )
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output FASTA already exists.")
    args = parser.parse_args()
    run(force=args.force)


if __name__ == "__main__":
    main()
