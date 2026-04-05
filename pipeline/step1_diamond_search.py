"""
Step 1 — DIAMOND blastp search.

For each genome, runs DIAMOND blastp using the three reference sequences
(MazF, MazE, alr) in ref_protein_mazEF.faa as subjects and the genome's Prodigal
protein FAA as the query.

Output: results/matches/{GCF_accession}_matches.tsv  (outfmt 6, 12 columns)

Usage:
    python pipeline/step1_diamond_search.py [--threads N] [--force]
"""

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    FAA_DIR, NCBI_DATA_DIR, MATCHES_DIR, REFERENCE_FASTA,
    DIAMOND_DB, DIAMOND_EVALUE, DIAMOND_THREADS, OUTPUT_DIR,
)


def _make_diamond_db(force: bool = False) -> None:
    """Build the DIAMOND database from the reference FASTA if needed."""
    if os.path.exists(DIAMOND_DB) and not force:
        print(f"  DIAMOND DB already exists: {DIAMOND_DB}")
        return

    os.makedirs(os.path.dirname(DIAMOND_DB), exist_ok=True)
    cmd = ["diamond", "makedb", "--in", REFERENCE_FASTA, "--db", DIAMOND_DB, "--quiet"]
    print(f"  Building DIAMOND database from {REFERENCE_FASTA} ...")
    subprocess.run(cmd, check=True)
    print(f"  Database written to {DIAMOND_DB}")


def _get_faa_path(refseq: str) -> str | None:
    """Return the FAA file path for *refseq*, or None if not found.

    Tries two locations:
      1. results/faa/{refseq}.faa  (Prodigal output from original analysis)
      2. Bacillota_ref_dataset/.../data/{refseq}/extracted_proteins.faa
    """
    primary = os.path.join(FAA_DIR, f"{refseq}.faa")
    if os.path.exists(primary):
        return primary

    # Fallback: NCBI data tree
    pattern = os.path.join(NCBI_DATA_DIR, f"{refseq}*", "extracted_proteins.faa")
    hits = glob.glob(pattern)
    if hits:
        return hits[0]

    return None


def _collect_refseqs() -> list[str]:
    """Gather all GCF accessions from available FAA files."""
    refseqs = []

    # Primary: results/faa/*.faa
    if os.path.isdir(FAA_DIR):
        for fname in os.listdir(FAA_DIR):
            if fname.endswith(".faa") and not fname.startswith("."):
                refseqs.append(fname[:-4])  # strip .faa

    if not refseqs:
        # Fallback: scan NCBI data tree
        for entry in os.scandir(NCBI_DATA_DIR):
            if entry.is_dir() and entry.name.startswith("GCF_"):
                faa = os.path.join(entry.path, "extracted_proteins.faa")
                if os.path.exists(faa):
                    refseqs.append(entry.name)

    return sorted(set(refseqs))


def _run_diamond(refseq: str, faa_path: str, threads: int, force: bool) -> str:
    """Run diamond blastp for one genome. Returns a status string."""
    out_tsv = os.path.join(MATCHES_DIR, f"{refseq}_matches.tsv")

    if os.path.exists(out_tsv) and not force:
        return f"SKIP  {refseq}"

    cmd = [
        "diamond", "blastp",
        "--db", DIAMOND_DB,
        "--query", faa_path,
        "--out", out_tsv,
        "--outfmt", "6",
        "--evalue", str(DIAMOND_EVALUE),
        "--threads", str(threads),
        "--quiet",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return f"OK    {refseq}"
    except subprocess.CalledProcessError as exc:
        return f"ERROR {refseq}: {exc.stderr.decode()[:120]}"


def run(threads: int = DIAMOND_THREADS, force: bool = False) -> None:
    os.makedirs(MATCHES_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    _make_diamond_db(force=force)

    refseqs = _collect_refseqs()
    total = len(refseqs)
    if total == 0:
        print("ERROR: No FAA files found. Run step0 first or check FAA_DIR in config.")
        sys.exit(1)

    print(f"\nStep 1: DIAMOND blastp — {total} genomes, {threads} threads")

    # Determine how many threads to give each DIAMOND call.
    # We parallelise at the genome level; each call gets 1 thread.
    genome_threads = max(1, threads // min(threads, 8))
    workers = min(threads, 8)

    done = 0
    errors = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {}
        for refseq in refseqs:
            faa = _get_faa_path(refseq)
            if faa is None:
                errors.append(f"  No FAA for {refseq}")
                continue
            fut = pool.submit(_run_diamond, refseq, faa, genome_threads, force)
            futures[fut] = refseq

        for fut in as_completed(futures):
            result = fut.result()
            done += 1
            if result.startswith("ERROR"):
                errors.append(result)
            if done % 100 == 0 or done == total:
                print(f"  Progress: {done}/{total} ({100*done//total}%)")

    print(f"\n  Done. Results in {MATCHES_DIR}/")
    if errors:
        print(f"\n  {len(errors)} errors:")
        for e in errors[:10]:
            print(f"    {e}")
        if len(errors) > 10:
            print(f"    ... and {len(errors)-10} more.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DIAMOND blastp for all genomes.")
    parser.add_argument("--threads", type=int, default=DIAMOND_THREADS,
                        help="Number of parallel threads (default: %(default)s)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output TSV already exists.")
    args = parser.parse_args()
    run(threads=args.threads, force=args.force)


if __name__ == "__main__":
    main()
