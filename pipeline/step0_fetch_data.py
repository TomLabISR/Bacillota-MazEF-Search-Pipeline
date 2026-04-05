"""
Step 0 — Fetch pre-computed data from Zenodo.

Heavy files (genome FAA archive, InterProScan TSV, tree, 16S data) are hosted
on Zenodo and are NOT committed to the GitHub repository.

Usage:
    python pipeline/step0_fetch_data.py            # download all missing files
    python pipeline/step0_fetch_data.py --force    # re-download even if present
"""

import argparse
import hashlib
import os
import sys
import tarfile
import urllib.request

# Add project root to path so we can import config
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    ZENODO_BASE_URL, ZENODO_RECORD_ID, ZENODO_FILES,
    FAA_DIR, NCBI_DATA_DIR, OUTPUT_DIR, PROJECT_ROOT,
    TREE_FILE, CUTSITE_TSV,
)


# ---------------------------------------------------------------------------
# File manifest: filename on Zenodo -> (local destination, description)
# ---------------------------------------------------------------------------
MANIFEST = {
    "extracted_proteins.tar.gz": (
        os.path.join(OUTPUT_DIR, "extracted_proteins.tar.gz"),
        "Prodigal protein FAA files for all 3,362 genomes",
    ),
    "ORF_detection.tsv": (
        os.path.join(PROJECT_ROOT, "detecting_mazE", "ORF_detection.tsv"),
        "Pre-computed InterProScan results for adjacent ORFs",
    ),
    "MazEFTree_min15.txt": (
        TREE_FILE,
        "Pre-computed Bacillota taxonomy tree pruned to genera with ≥15 genomes",
    ),
    "16s_cutsite.tsv": (
        CUTSITE_TSV,
        "16S rRNA sequences with TACATA cutsite counts (taxonomy fallback)",
    ),
}


def _download(url: str, dest: str) -> None:
    """Download *url* to *dest* with a simple progress indicator."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print(f"  Downloading {url}")
    print(f"  → {dest}")

    def _reporthook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            pct = min(100.0, downloaded * 100.0 / total_size)
            mb = downloaded / 1_048_576
            print(f"\r    {pct:5.1f}%  ({mb:.1f} MB)", end="", flush=True)

    urllib.request.urlretrieve(url, dest, reporthook=_reporthook)
    print()  # newline after progress


def _extract_tar(archive: str, dest_dir: str) -> None:
    """Extract a .tar.gz archive, placing contents under *dest_dir*."""
    print(f"  Extracting {os.path.basename(archive)} → {dest_dir}/")
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        for i, member in enumerate(members, 1):
            tar.extract(member, dest_dir)
            if i % 500 == 0 or i == len(members):
                print(f"\r    Extracted {i}/{len(members)} files", end="", flush=True)
    print()


def fetch_all(force: bool = False) -> None:
    if ZENODO_RECORD_ID == "XXXXXXX":
        print(
            "ERROR: Zenodo record ID has not been set yet.\n"
            "       Edit pipeline/config.py and replace ZENODO_RECORD_ID = 'XXXXXXX'\n"
            "       with the real record ID once the dataset is uploaded to Zenodo."
        )
        sys.exit(1)

    for filename, (local_path, description) in MANIFEST.items():
        print(f"\n[{filename}] — {description}")

        if os.path.exists(local_path) and not force:
            print(f"  Already present, skipping (use --force to re-download).")
            continue

        url = f"{ZENODO_BASE_URL}/{filename}"
        _download(url, local_path)

        # Post-processing for archive
        if filename.endswith(".tar.gz"):
            _extract_tar(local_path, FAA_DIR)
            print(f"  FAA files extracted to {FAA_DIR}/")

    print("\nStep 0 complete. All data files are in place.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download pre-computed data files from Zenodo."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist locally.",
    )
    args = parser.parse_args()
    fetch_all(force=args.force)


if __name__ == "__main__":
    main()
