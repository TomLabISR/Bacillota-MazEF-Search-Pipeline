"""
Step 5 — InterProScan RHH domain analysis.

For ORFs adjacent to MazF that did NOT get a direct MazE DIAMOND hit,
parse the pre-computed ORF_detection.tsv (or run InterProScan locally with
--run-interproscan) to find Ribbon-Helix-Helix (RHH) antitoxin domains.

Two-tier detection strategy: ORFs that did not yield a direct MazE DIAMOND
hit are re-evaluated for the Ribbon-Helix-Helix fold via InterProScan,
capturing divergent antitoxins that lack sequence similarity to the reference.

False-positive filters:
  - Reject hits containing IPR010994 (RuvA domain — unrelated protein)
  - Reject proteins longer than MAX_MAZE_LENGTH (150 aa)

Output: dict  { interpro_seq_id: InterproResult }
        Also saves pipeline_output/interpro_results.tsv as checkpoint.

Usage:
    python pipeline/step5_interproscan.py [--run-interproscan] [--force]
"""

import os
import subprocess
import sys
from dataclasses import dataclass, field

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    INTERPROSCAN_TSV, INTERPRO_COLS, RHH_ACCESSIONS, RHH_FALSE_POSITIVES,
    MAX_MAZE_LENGTH, ADJ_FASTA, OUTPUT_DIR,
)

INTERPRO_CHECKPOINT = os.path.join(OUTPUT_DIR, "interpro_results.tsv")


@dataclass
class InterproResult:
    """Aggregated InterProScan result for one sequence."""
    seq_id:         str
    seq_length:     int
    has_rhh:        bool
    is_valid_maze:  bool
    domains:        list[dict] = field(default_factory=list)

    def domains_summary(self) -> str:
        """Compact string representation of all RHH domains detected."""
        parts = []
        for d in self.domains:
            analysis = d.get("analysis", "")
            sig_acc  = d.get("signature_accession", "")
            interpro = d.get("interpro_accession",  "")
            desc     = d.get("interpro_description", "") or d.get("signature_description", "")
            if interpro and interpro != "-":
                parts.append(f"{analysis}:{interpro}({desc})")
            elif sig_acc:
                parts.append(f"{analysis}:{sig_acc}({desc})")
        return "; ".join(parts)


def _is_rhh_accession(accession: str) -> bool:
    return (
        accession in RHH_ACCESSIONS
        or any(accession.startswith(pfx) for pfx in ["G3DSA:1.10.1220"])
    )


def _short_seq_id(interpro_seq_id: str) -> str:
    """
    Convert InterProScan seq_id back to the Prodigal fragment_orfindex form.

    'GCF_014287895.1_NZ_JACOPP010000019.1_30' → 'NZ_JACOPP010000019.1_30'
    """
    for prefix in ("_NC_", "_NZ_", "_CP_", "_AP_", "_LN_"):
        idx = interpro_seq_id.find(prefix)
        if idx != -1:
            return interpro_seq_id[idx + 1:]
    parts = interpro_seq_id.split("_", 2)
    return "_".join(parts[2:]) if len(parts) >= 3 else interpro_seq_id


def _parse_interproscan_tsv(path: str) -> dict[str, InterproResult]:
    """
    Parse the InterProScan TSV and return RHH-positive sequences.

    The TSV has no header row — column names are assigned from INTERPRO_COLS.
    """
    if not os.path.exists(path):
        print(f"  WARNING: InterProScan TSV not found: {path}")
        print("  Run with --run-interproscan or download from Zenodo (step0).")
        return {}

    print(f"  Parsing InterProScan TSV: {path}")

    results: dict[str, InterproResult] = {}
    chunk_size = 50_000
    reader = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=INTERPRO_COLS,
        chunksize=chunk_size,
        dtype=str,
        low_memory=False,
        on_bad_lines="skip",
    )

    total_rows = 0
    for chunk in reader:
        total_rows += len(chunk)

        for _, row in chunk.iterrows():
            seq_id       = str(row.get("sequence_id", "")).strip()
            sig_acc      = str(row.get("signature_accession", "")).strip()
            interpro_acc = str(row.get("interpro_accession", "")).strip()
            seq_len_raw  = str(row.get("sequence_length", "0"))

            if not seq_id or seq_id == "nan":
                continue

            try:
                seq_len = int(float(seq_len_raw))
            except (ValueError, TypeError):
                seq_len = 0

            if seq_id not in results:
                results[seq_id] = InterproResult(
                    seq_id        = seq_id,
                    seq_length    = seq_len,
                    has_rhh       = False,
                    is_valid_maze = True,
                )

            result = results[seq_id]
            result.seq_length = max(result.seq_length, seq_len)

            if interpro_acc in RHH_FALSE_POSITIVES or sig_acc in RHH_FALSE_POSITIVES:
                result.is_valid_maze = False

            if _is_rhh_accession(sig_acc) or _is_rhh_accession(interpro_acc):
                result.has_rhh = True
                result.domains.append({
                    "analysis":              str(row.get("analysis", "")),
                    "signature_accession":   sig_acc,
                    "signature_description": str(row.get("signature_description", "")),
                    "interpro_accession":    interpro_acc,
                    "interpro_description":  str(row.get("interpro_description", "")),
                })

    print(f"  Processed {total_rows:,} rows")

    valid_rhh: dict[str, InterproResult] = {}
    for seq_id, res in results.items():
        if not res.has_rhh:
            continue
        if res.seq_length > MAX_MAZE_LENGTH:
            res.is_valid_maze = False
        if res.is_valid_maze:
            valid_rhh[seq_id] = res

    print(f"  RHH-positive sequences: {len(valid_rhh)}")
    return valid_rhh


def _run_interproscan(fasta: str, out_tsv: str) -> None:
    cmd = [
        "interproscan.sh",
        "-i",    fasta,
        "-f",    "TSV",
        "-o",    out_tsv,
        "-appl", "Pfam,SUPERFAMILY,SMART,ProSiteProfiles",
        "--goterms",
        "--pathways",
    ]
    print("  Running InterProScan (this may take several days for large datasets) ...")
    subprocess.run(cmd, check=True)


def _save_checkpoint(results: dict[str, InterproResult]) -> None:
    rows = [
        {
            "seq_id":        r.seq_id,
            "seq_length":    r.seq_length,
            "has_rhh":       r.has_rhh,
            "is_valid_maze": r.is_valid_maze,
            "domains":       r.domains_summary(),
        }
        for r in results.values()
    ]
    pd.DataFrame(rows).to_csv(INTERPRO_CHECKPOINT, sep="\t", index=False)
    print(f"  InterProScan checkpoint written: {INTERPRO_CHECKPOINT}")


def run(
    run_interproscan: bool = False,
    force: bool = False,
) -> dict[str, InterproResult]:
    """
    Return InterproResult dict keyed by the full InterProScan seq_id.
    """
    if os.path.exists(INTERPRO_CHECKPOINT) and not force and not run_interproscan:
        print(f"  Loading InterProScan results from checkpoint: {INTERPRO_CHECKPOINT}")
        df = pd.read_csv(INTERPRO_CHECKPOINT, sep="\t", dtype=str)
        results = {}
        for _, row in df.iterrows():
            seq_id = str(row["seq_id"])
            results[seq_id] = InterproResult(
                seq_id        = seq_id,
                seq_length    = int(float(str(row.get("seq_length", 0)))),
                has_rhh       = str(row.get("has_rhh", "False")).lower() == "true",
                is_valid_maze = str(row.get("is_valid_maze", "True")).lower() == "true",
            )
        print(f"  Loaded {len(results)} validated RHH hits from checkpoint.")
        return results

    print(f"\nStep 5: InterProScan RHH analysis")

    if run_interproscan:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_tsv = os.path.join(OUTPUT_DIR, "ORF_detection_new.tsv")
        _run_interproscan(ADJ_FASTA, out_tsv)
        tsv_path = out_tsv
    else:
        tsv_path = INTERPROSCAN_TSV

    results = _parse_interproscan_tsv(tsv_path)
    _save_checkpoint(results)
    return results


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Parse InterProScan RHH domain results."
    )
    parser.add_argument("--run-interproscan", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    results = run(run_interproscan=args.run_interproscan, force=args.force)
    print(f"\nStep 5 complete — {len(results)} validated RHH hits.")


if __name__ == "__main__":
    main()
