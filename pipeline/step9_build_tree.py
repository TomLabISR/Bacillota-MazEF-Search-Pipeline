"""
Step 9 — Prune the full Bacillota topology tree to genera with ≥ MIN_GENUS_SIZE genomes.

The topology source (`data/MazEFTree_topology.tree`) covers all 140 Bacillota genera
with leaf names in the format `GenusName__count_`.  This step:
  1. Strips the embedded counts from leaf names to extract clean genus names.
  2. Re-counts genomes per genus from the taxonomy CSV (Step 7 output) so the
     pruning reflects the current dataset's genus-level representation.
  3. Prunes leaves where count < MIN_GENUS_SIZE (15).
  4. Writes the pruned tree with clean leaf names to pipeline_output/MazEFTree_min15.txt.

Output: pipeline_output/MazEFTree_min15.txt

Usage:
    python pipeline/step9_build_tree.py [--force]
"""

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    TAXONOMY_CSV, TREE_TOPOLOGY, TREE_FILE, OUTPUT_DIR, MIN_GENUS_SIZE,
)

_COUNT_RE = re.compile(r"__(\d+)_$")


def _strip_count(leaf_name: str) -> str:
    """Remove __count_ suffix from a leaf name, returning the clean genus name."""
    return _COUNT_RE.sub("", leaf_name)


def run(force: bool = False) -> None:
    if os.path.exists(TREE_FILE) and not force:
        print(f"  Step 9 output already exists ({TREE_FILE}). Skipping.")
        return

    if not os.path.exists(TAXONOMY_CSV):
        print(f"ERROR: {TAXONOMY_CSV} not found. Run Step 7 first.")
        sys.exit(1)

    if not os.path.exists(TREE_TOPOLOGY):
        print(f"ERROR: Topology tree not found at {TREE_TOPOLOGY}.")
        sys.exit(1)

    try:
        from ete3 import Tree
    except ImportError:
        print("ERROR: ete3 is required for Step 9. Install with: conda install -c conda-forge ete3")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Count genomes per genus from taxonomy table
    tax_df = pd.read_csv(TAXONOMY_CSV, low_memory=False)
    genus_counts = (
        tax_df["Genus"]
        .dropna()
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .value_counts()
    )

    print(f"\nStep 9: Pruning topology tree (min genus size = {MIN_GENUS_SIZE}) ...")
    print(f"  Taxonomy CSV: {len(tax_df)} genomes, {len(genus_counts)} named genera")

    # Load topology tree (leaves have __count_ suffixes from original dataset)
    tree = Tree(TREE_TOPOLOGY, format=1)
    initial_count = len(tree.get_leaves())
    print(f"  Topology tree: {initial_count} leaves")

    # Determine which leaves to keep based on our current genus counts
    leaves_to_keep = []
    for leaf in tree.get_leaves():
        genus = _strip_count(leaf.name)
        if genus_counts.get(genus, 0) >= MIN_GENUS_SIZE:
            leaves_to_keep.append(leaf)

    if not leaves_to_keep:
        print(f"ERROR: No genera meet the threshold of {MIN_GENUS_SIZE} genomes.")
        sys.exit(1)

    tree.prune(leaves_to_keep, preserve_branch_length=True)

    # Strip __count_ suffixes from leaf names in the pruned tree
    for leaf in tree.get_leaves():
        leaf.name = _strip_count(leaf.name)

    tree.write(outfile=TREE_FILE, format=1)

    final_count = len(tree.get_leaves())
    print(f"  Retained: {final_count} genera (removed {initial_count - final_count})")
    print(f"\n  Output: {TREE_FILE}")

    # Print top genera for verification
    top = sorted(
        [(leaf.name, genus_counts.get(leaf.name, 0)) for leaf in tree.get_leaves()],
        key=lambda x: x[1], reverse=True
    )
    print(f"\n  Top genera in pruned tree:")
    for g, n in top[:10]:
        print(f"    {g:<30} {n} genomes")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Prune Bacillota topology tree to genera with sufficient genomes."
    )
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output exists.")
    args = parser.parse_args()
    run(force=args.force)


if __name__ == "__main__":
    main()
