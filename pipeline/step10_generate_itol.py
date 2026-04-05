"""
Step 10 — Generate iTOL annotation files from pipeline output.

Reads the taxonomy CSV (Step 7 output) and the pruned tree (Step 9 output),
aggregates per-genus statistics, and writes 7 iTOL annotation text files to
pipeline_output/itol_annotations/.

Generated files:
  genus_mazE_prevalence.txt   — MazE antitoxin prevalence (Spectral_r colormap)
  genus_mazF_prevalence.txt   — MazF toxin prevalence (Spectral_r colormap)
  genus_tacata_prevalence.txt — TACATA site prevalence in 16S (Spectral_r)
  genus_primary_operon.txt    — Primary operon structure (discrete colors)
  genus_secondary_operon.txt  — Secondary operon structure (discrete colors)
  genus_mazF_identity.txt     — Average MazF DIAMOND identity (Spectral_r)
  genus_labels.txt            — Genus name + strain count labels

Internal column mapping (no changes to CSV files):
  combined_mazE_detected → total_mazE_detected
  operon_structure       → updatedOperon

Usage:
    python pipeline/step10_generate_itol.py [--force] [--palette PALETTE]
"""

import os
import sys
from collections import defaultdict

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    TAXONOMY_CSV, TREE_FILE, ITOL_DIR, OUTPUT_DIR, CUTSITE_TSV,
)

# Fixed discrete colors for canonical operon structures
OPERON_COLORS = {
    "alr-mazE-mazF": "#2166ac",   # canonical — blue
    "mazE-mazF":     "#92c5de",   # partial, no alr — light blue
    "alr-mazF":      "#f4a582",   # partial, no MazE — salmon
    "mazF_only":     "#d6604d",   # lone toxin — red
    "No MazF":       "#cccccc",   # no MazF detected — grey
}


# ---------------------------------------------------------------------------
# Tree parsing
# ---------------------------------------------------------------------------

def _parse_newick(newick_file: str) -> list[str]:
    """Extract leaf names from a Newick file using ete3."""
    try:
        from ete3 import Tree
        tree = Tree(newick_file, format=1)
        return [leaf.name for leaf in tree.get_leaves()]
    except ImportError:
        print("WARNING: ete3 not available; falling back to simple Newick parser.")
        with open(newick_file) as f:
            newick = f.read().strip()
        taxa = []
        current = ""
        for char in newick:
            if char in "(),":
                t = current.strip()
                if t:
                    colon = t.find(":")
                    taxa.append(t[:colon].strip() if colon != -1 else t)
                current = ""
            elif char != ":":
                current += char
        return [t for t in taxa if t]


# ---------------------------------------------------------------------------
# Per-genus aggregation
# ---------------------------------------------------------------------------

def _aggregate_by_genus(data: pd.DataFrame, genera_list: list[str]) -> dict:
    """Aggregate per-row data into per-genus summary statistics."""
    summary = {
        g: {
            "total": 0,
            "mazE_present": 0,
            "mazF_present": 0,
            "mazF_identity_sum": 0.0,
            "mazF_identity_count": 0,
            "tacata_present": 0,
            "has_16s": 0,
            "operon_structures": defaultdict(int),
        }
        for g in genera_list
    }

    for _, row in data.iterrows():
        genus = row.get("Genus")
        if genus not in summary:
            continue

        s = summary[genus]
        s["total"] += 1

        # MazE presence (uses internal alias total_mazE_detected)
        if row.get("total_mazE_detected") is True:
            s["mazE_present"] += 1

        # MazF presence + identity
        identity = row.get("mazF_identity")
        if identity is not None and not pd.isna(identity):
            s["mazF_present"] += 1
            s["mazF_identity_sum"] += float(identity)
            s["mazF_identity_count"] += 1

        # TACATA presence
        tacata = row.get("TACATA_count")
        if tacata is not None and not pd.isna(tacata):
            s["has_16s"] += 1
            if float(tacata) > 0:
                s["tacata_present"] += 1

        # Operon structure (uses internal alias updatedOperon)
        operon = row.get("updatedOperon")
        if operon is None or (isinstance(operon, float) and pd.isna(operon)):
            operon = "No MazF"
        s["operon_structures"][str(operon)] += 1

    # Calculate derived statistics
    for genus, s in summary.items():
        n = s["total"]
        if n == 0:
            continue
        s["mazE_prevalence"] = s["mazE_present"] / n
        s["mazF_prevalence"] = s["mazF_present"] / n
        s["tacata_prevalence"] = (
            s["tacata_present"] / s["has_16s"] if s["has_16s"] > 0 else 0.0
        )
        s["avg_mazF_identity"] = (
            s["mazF_identity_sum"] / s["mazF_identity_count"]
            if s["mazF_identity_count"] > 0 else 0.0
        )

        sorted_structs = sorted(
            s["operon_structures"].items(), key=lambda x: x[1], reverse=True
        )
        s["top_operon"] = sorted_structs[0][0] if sorted_structs else "unknown"
        s["top_operon_count"] = sorted_structs[0][1] if sorted_structs else 0
        s["top_operon_pct"] = s["top_operon_count"] / n

        if len(sorted_structs) >= 2:
            s["second_operon"] = sorted_structs[1][0]
            s["second_operon_count"] = sorted_structs[1][1]
            s["second_operon_pct"] = sorted_structs[1][1] / n
            s["top_to_second_ratio"] = (
                sorted_structs[0][1] / sorted_structs[1][1]
                if sorted_structs[1][1] > 0 else float("inf")
            )
        else:
            s["second_operon"] = None
            s["second_operon_count"] = 0
            s["second_operon_pct"] = 0.0
            s["top_to_second_ratio"] = float("inf")

    return summary


# ---------------------------------------------------------------------------
# iTOL file generators
# ---------------------------------------------------------------------------

def _spectral_color(value: float, palette: list) -> str:
    """Map a 0–1 float to a hex color from the 256-entry Spectral_r palette."""
    idx = max(0, min(255, int(value * 255)))
    return mcolors.rgb2hex(palette[idx])


def _prevalence_strip(
    genera: list[str],
    summary: dict,
    stat_key: str,
    label: str,
    legend_title: str,
    palette: list,
    no_data_label: str = "No data",
) -> str:
    """Generic DATASET_COLORSTRIP for a 0–1 prevalence metric."""
    out = [
        "DATASET_COLORSTRIP",
        "SEPARATOR COMMA",
        f"DATASET_LABEL,{label}",
        f"COLOR,{mcolors.rgb2hex(palette[-1])}",
        "STRIP_WIDTH,50",
        "MARGIN,5",
        "BORDER_WIDTH,0.5",
        "BORDER_COLOR,#000000",
        f"LEGEND_TITLE,{legend_title}",
        "LEGEND_SHAPES," + ",".join(["1"] * 11),
    ]
    legend_colors = [_spectral_color(i / 100.0, palette) for i in range(0, 101, 10)]
    legend_labels = [str(i) for i in range(0, 101, 10)]
    out.append("LEGEND_COLORS," + ",".join(legend_colors))
    out.append("LEGEND_LABELS," + ",".join(legend_labels))
    out.append("")
    out.append("DATA")

    for genus in genera:
        s = summary.get(genus, {})
        if s.get("total", 0) > 0 and stat_key in s:
            color = _spectral_color(s[stat_key], palette)
            out.append(f"{genus},{color},{s[stat_key]*100:.1f}%")
        else:
            out.append(f"{genus},#CCCCCC,{no_data_label}")

    return "\n".join(out) + "\n"


def _identity_strip(
    genera: list[str],
    summary: dict,
    palette: list,
) -> str:
    """DATASET_COLORSTRIP for average MazF DIAMOND identity (0–100 scale)."""
    out = [
        "DATASET_COLORSTRIP",
        "SEPARATOR COMMA",
        "DATASET_LABEL,Average mazF Identity",
        "COLOR,#ff0000",
        "STRIP_WIDTH,50",
        "MARGIN,5",
        "BORDER_WIDTH,0.5",
        "BORDER_COLOR,#000000",
        "LEGEND_TITLE,mazF Identity (%)",
        "LEGEND_SHAPES," + ",".join(["1"] * 11),
    ]
    legend_colors = [_spectral_color(i / 100.0, palette) for i in range(0, 101, 10)]
    legend_labels = [str(i) for i in range(0, 101, 10)]
    out.append("LEGEND_COLORS," + ",".join(legend_colors))
    out.append("LEGEND_LABELS," + ",".join(legend_labels))
    out.append("")
    out.append("DATA")

    for genus in genera:
        s = summary.get(genus, {})
        if s.get("mazF_identity_count", 0) > 0:
            identity = s["avg_mazF_identity"]
            color = _spectral_color(identity / 100.0, palette)
            out.append(f"{genus},{color},{identity:.1f}%")
        else:
            out.append(f"{genus},#CCCCCC,No mazF")

    return "\n".join(out) + "\n"


def _operon_strips(genera: list[str], summary: dict) -> tuple[str, str]:
    """Two DATASET_COLORSTRIP files: primary and secondary operon structures."""
    # Collect all structures seen; assign colors (fixed map + hash for unknowns)
    all_structs = set()
    for s in summary.values():
        if s.get("top_operon"):
            all_structs.add(s["top_operon"])
        if s.get("second_operon"):
            all_structs.add(s["second_operon"])
    all_structs.discard(None)

    color_map = dict(OPERON_COLORS)
    for op in all_structs:
        if op not in color_map:
            color_map[op] = "#%06X" % (hash(op) % 0xFFFFFF)

    sorted_structs = sorted(color_map.keys())
    legend_shapes  = ",".join(["1"] * len(sorted_structs))
    legend_colors  = ",".join(color_map[op] for op in sorted_structs)
    legend_labels  = ",".join(sorted_structs)

    def _header(label: str, default_color: str) -> list[str]:
        return [
            "DATASET_COLORSTRIP",
            "SEPARATOR COMMA",
            f"DATASET_LABEL,{label}",
            f"COLOR,{default_color}",
            "STRIP_WIDTH,50",
            "MARGIN,5",
            "BORDER_WIDTH,0.5",
            "BORDER_COLOR,#000000",
            "LEGEND_TITLE,Operon Structure",
            f"LEGEND_SHAPES,{legend_shapes}",
            f"LEGEND_COLORS,{legend_colors}",
            f"LEGEND_LABELS,{legend_labels}",
            "",
            "DATA",
        ]

    p_out = _header("Primary Operon Structure", "#ff0000")
    s_out = _header("Secondary Operon Structure", "#00ff00")

    for genus in genera:
        s = summary.get(genus, {})
        if s.get("total", 0) > 0:
            top_op  = s.get("top_operon", "unknown")
            top_pct = int(round(s.get("top_operon_pct", 0) * 100))
            p_out.append(f"{genus},{color_map.get(top_op, '#A9A9A9')},{top_op} ({top_pct}%)")

            second_op = s.get("second_operon")
            if second_op and s.get("second_operon_count", 0) > 0:
                sec_pct = int(round(s.get("second_operon_pct", 0) * 100))
                ratio   = s.get("top_to_second_ratio", 0)
                if ratio < float("inf"):
                    s_out.append(f"{genus},{color_map.get(second_op, '#A9A9A9')},{second_op} ({sec_pct}%) [1:{ratio:.1f}]")
                else:
                    s_out.append(f"{genus},{color_map.get(second_op, '#A9A9A9')},{second_op} ({sec_pct}%)")
            else:
                s_out.append(f"{genus},#EEEEEE,None")
        else:
            p_out.append(f"{genus},#CCCCCC,No data")
            s_out.append(f"{genus},#EEEEEE,None")

    return "\n".join(p_out) + "\n", "\n".join(s_out) + "\n"


def _label_file(genera: list[str], summary: dict) -> str:
    out = ["LABELS", "SEPARATOR COMMA", "DATA"]
    for genus in genera:
        n = summary.get(genus, {}).get("total", 0)
        out.append(f"{genus},{genus} ({n})")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(force: bool = False, palette: str = "Spectral_r") -> None:
    # Check outputs
    if os.path.isdir(ITOL_DIR) and not force:
        existing = [f for f in os.listdir(ITOL_DIR) if f.endswith(".txt")]
        if len(existing) >= 7:
            print(f"  Step 10 output already exists ({ITOL_DIR}/, {len(existing)} files). Skipping.")
            return

    if not os.path.exists(TAXONOMY_CSV):
        print(f"ERROR: {TAXONOMY_CSV} not found. Run Step 7 first.")
        sys.exit(1)

    if not os.path.exists(TREE_FILE):
        print(f"ERROR: {TREE_FILE} not found. Run Step 9 first.")
        sys.exit(1)

    os.makedirs(ITOL_DIR, exist_ok=True)
    print(f"\nStep 10: Generating iTOL annotations (palette: {palette}) ...")

    # Load taxonomy data
    data = pd.read_csv(TAXONOMY_CSV, low_memory=False)
    print(f"  Taxonomy CSV: {len(data)} genomes")

    # Internal column mapping — do NOT modify the CSV
    data = data.rename(columns={
        "combined_mazE_detected": "total_mazE_detected",
        "operon_structure": "updatedOperon",
    })

    # Join TACATA_count from 16s_cutsite.tsv if available
    has_tacata = False
    if os.path.exists(CUTSITE_TSV):
        cs = pd.read_csv(CUTSITE_TSV, sep="\t", dtype=str, low_memory=False)
        cs.columns = [c.strip() for c in cs.columns]
        refseq_col = next(
            (c for c in cs.columns if c.lower() in ("refseq", "assembly")), None
        )
        if refseq_col and "TACATA_count" in cs.columns:
            cs = cs[[refseq_col, "TACATA_count"]].rename(columns={refseq_col: "refseq"})
            cs["TACATA_count"] = pd.to_numeric(cs["TACATA_count"], errors="coerce")
            data = data.merge(cs, on="refseq", how="left")
            has_tacata = True
            print(f"  TACATA_count joined from 16s_cutsite.tsv")
    if not has_tacata:
        print("  WARNING: 16s_cutsite.tsv not available — tacata_prevalence.txt will be skipped.")
        data["TACATA_count"] = float("nan")

    # Parse genera from pruned tree
    genera_list = _parse_newick(TREE_FILE)
    print(f"  Pruned tree: {len(genera_list)} genera")

    # Aggregate per genus
    summary = _aggregate_by_genus(data, genera_list)
    covered = sum(1 for g in genera_list if summary.get(g, {}).get("total", 0) > 0)
    print(f"  Genera with data: {covered}/{len(genera_list)}")

    # Build Spectral_r palette once
    sp = sns.color_palette(palette, as_cmap=False, n_colors=256)

    # Generate files
    files = {}
    files["genus_mazE_prevalence.txt"] = _prevalence_strip(
        genera_list, summary, "mazE_prevalence",
        "MazE Antitoxin Prevalence", "MazE Prevalence (%)", sp,
    )
    files["genus_mazF_prevalence.txt"] = _prevalence_strip(
        genera_list, summary, "mazF_prevalence",
        "MazF Toxin Prevalence", "MazF Prevalence (%)", sp,
    )
    if has_tacata:
        files["genus_tacata_prevalence.txt"] = _prevalence_strip(
            genera_list, summary, "tacata_prevalence",
            "TACATA Site Prevalence (16S)", "TACATA Prevalence (%)", sp,
            no_data_label="No 16S data",
        )
    else:
        print("  Skipping genus_tacata_prevalence.txt (no TACATA data).")

    files["genus_mazF_identity.txt"] = _identity_strip(genera_list, summary, sp)

    primary, secondary = _operon_strips(genera_list, summary)
    files["genus_primary_operon.txt"]   = primary
    files["genus_secondary_operon.txt"] = secondary
    files["genus_labels.txt"]           = _label_file(genera_list, summary)

    # Write all files
    for fname, content in files.items():
        fpath = os.path.join(ITOL_DIR, fname)
        with open(fpath, "w") as fh:
            fh.write(content)

    print(f"\n  Output: {ITOL_DIR}/")
    for fname in files:
        print(f"    {fname}")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate iTOL annotation files from pipeline output."
    )
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output files already exist.")
    parser.add_argument("--palette", default="Spectral_r",
                        help="Seaborn/matplotlib colormap name (default: Spectral_r).")
    args = parser.parse_args()
    run(force=args.force, palette=args.palette)


if __name__ == "__main__":
    main()
