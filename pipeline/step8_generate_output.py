"""
Step 8 — Generate the 6-sheet supplementary Excel file (sup_new_fig5.xlsx).

Sheet structure (matching the reference file):
  1. README          — column descriptions
  2. Operon_Summary  — one row per genome; includes all genomes (empty rows for
                       genomes with no MazF hit)
  3. MazF_Sequences  — protein sequences with coordinates
  4. MazE_Sequences  — protein sequences + detection method
  5. Alr_Sequences   — protein sequences with coordinates
  6. 16S_rRNA        — 16S sequences with TACATA counts

Column naming convention (matching reference exactly):
  - All sheets except 16S_rRNA use '  mazF_operon_id' (2 leading spaces)
  - 16S_rRNA uses 'operon_id' (no spaces)

Usage:
    python pipeline/step8_generate_output.py [--force]
    python pipeline/step8_generate_output.py --validate
"""

import argparse
import json
import os
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils.dataframe import dataframe_to_rows

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.config import (
    TAXONOMY_CSV, OUTPUT_EXCEL, OUTPUT_DIR,
    FAA_DIR, REFERENCE_SUPFIG, PROJECT_ROOT, CUTSITE_TSV,
)


# ---------------------------------------------------------------------------
# Coordinate lookup (Prodigal FAA headers)
# ---------------------------------------------------------------------------

def _load_faa_coords(faa_dir: str) -> dict[str, dict[str, tuple]]:
    coords: dict[str, dict] = {}
    if not os.path.isdir(faa_dir):
        print(f"  WARNING: FAA directory not found: {faa_dir}. Coordinates will be empty.")
        return coords

    faa_files = [f for f in os.listdir(faa_dir)
                 if f.endswith(".faa") and not f.startswith(".")]
    print(f"  Loading coordinates from {len(faa_files)} FAA files ...")

    for faa_file in faa_files:
        refseq = faa_file[:-4]
        coords[refseq] = {}
        with open(os.path.join(faa_dir, faa_file)) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                parts = line[1:].strip().split(" # ")
                if len(parts) >= 4:
                    coords[refseq][parts[0]] = (parts[1], parts[2], parts[3])
    return coords


def _get_coords(coords, refseq, fragment, orf_index) -> tuple[str, str, str]:
    if not refseq or not fragment or pd.isna(orf_index):
        return "", "", ""
    try:
        orf_str = str(int(float(orf_index)))
    except (ValueError, TypeError):
        orf_str = str(orf_index)
    seq_id = f"{fragment}_{orf_str}"
    return coords.get(refseq, {}).get(seq_id, ("", "", ""))


def _safe_str(val) -> str:
    """Return '' for NaN/None, otherwise str(val)."""
    if val is None:
        return ""
    try:
        if pd.isna(val):
            return ""
    except (TypeError, ValueError):
        pass
    return str(val)


# ---------------------------------------------------------------------------
# Sheet builders
# ---------------------------------------------------------------------------

def _build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        operon_id   = _safe_str(row.get("operon_id"))
        mazf_seq    = _safe_str(row.get("mazF_sequence"))
        mazf_detect = bool(mazf_seq and mazf_seq != "nan")
        has_maze    = row.get("combined_mazE_detected")
        maze_detect = bool(has_maze) if pd.notna(has_maze) else False
        has_alr     = row.get("has_alr")
        alr_detect  = bool(has_alr) if pd.notna(has_alr) else False

        rows.append({
            "  mazF_operon_id":      operon_id,
            "refseq":                _safe_str(row.get("refseq")),
            "fragment":              _safe_str(row.get("fragment")),
            "organism":              _safe_str(row.get("Organism")),
            "phylum":                _safe_str(row.get("Phylum")),
            "class":                 _safe_str(row.get("Class")),
            "order":                 _safe_str(row.get("Order")),
            "family":                _safe_str(row.get("Family")),
            "genus":                 _safe_str(row.get("Genus")),
            "operon_structure":      _safe_str(row.get("operon_structure")),
            "mazF_detected":         mazf_detect,
            "mazE_detected":         maze_detect,
            "mazE_detection_method": _safe_str(row.get("mazE_detection_method")),
            "alr_detected":          alr_detect,
            "is_canonical":          bool(row.get("is_canonical")) if pd.notna(row.get("is_canonical")) else False,
        })
    return pd.DataFrame(rows)


def _build_mazf_sequences(df: pd.DataFrame, coords: dict) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        seq = _safe_str(row.get("mazF_sequence"))
        if not seq or seq == "nan":
            continue
        start, end, strand = _get_coords(
            coords,
            _safe_str(row.get("refseq")),
            _safe_str(row.get("fragment")),
            row.get("mazF_orf_index"),
        )
        rows.append({
            "  mazF_operon_id": _safe_str(row.get("operon_id")),
            "orf_index":        row.get("mazF_orf_index", ""),
            "start":            start,
            "end":              end,
            "direction":        strand,
            "identity_pct":     row.get("mazF_identity", ""),
            "length":           row.get("mazF_length", ""),
            "sequence":         seq,
        })
    return pd.DataFrame(rows)


def _build_maze_sequences(df: pd.DataFrame, coords: dict) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        method = _safe_str(row.get("mazE_detection_method"))
        if not method:
            continue

        adjacent_raw = _safe_str(row.get("adjacent_genes_data")) or "[]"
        try:
            adjacent = json.loads(adjacent_raw)
        except (json.JSONDecodeError, ValueError):
            adjacent = []

        refseq   = _safe_str(row.get("refseq"))
        fragment = _safe_str(row.get("fragment"))

        if method == "homology":
            # Find mazE (or legacy ndoAI) at position -1
            maze_gene = next(
                (g for g in adjacent
                 if g.get("gene_name", "").lower() in ("maze", "ndoai")
                 and int(g.get("position", 0)) == -1),
                None,
            )
            if maze_gene is None:
                maze_gene = next(
                    (g for g in adjacent
                     if g.get("gene_name", "").lower() in ("maze", "ndoai")),
                    None,
                )
            if maze_gene is None:
                continue
            orf_idx  = maze_gene.get("orf_index", "")
            seq      = _safe_str(maze_gene.get("sequence"))
            identity = maze_gene.get("identity", "")
            length   = maze_gene.get("length", "")
            start, end, strand = _get_coords(coords, refseq, fragment, orf_idx)
            rows.append({
                "  mazF_operon_id":  _safe_str(row.get("operon_id")),
                "orf_index":         orf_idx,
                "start":             start,
                "end":               end,
                "direction":         strand,
                "identity_pct":      identity,
                "length":            length,
                "sequence":          seq,
                "detection_method":  "homology",
                "interpro_domains":  "",
            })

        elif method == "InterPro":
            pos_val = row.get("interpro_rhh_position")
            try:
                pos = int(float(pos_val))
            except (TypeError, ValueError):
                pos = -1

            if pos == -1:
                qid    = _safe_str(row.get("prev_orf_query_id"))
                seq    = _safe_str(row.get("prev_orf_sequence"))
                length = row.get("prev_orf_length", "")
            else:
                qid    = _safe_str(row.get("next_orf_query_id"))
                seq    = _safe_str(row.get("next_orf_sequence"))
                length = row.get("next_orf_length", "")

            if not qid or not seq or seq == "nan":
                continue

            orf_idx = qid.rsplit("_", 1)[-1] if qid else ""
            start, end, strand = _get_coords(coords, refseq, fragment, orf_idx)

            # Format InterPro domains
            domains_raw = _safe_str(row.get("interpro_rhh_domains"))
            try:
                domain_list = json.loads(domains_raw) if domains_raw.startswith("[") else []
            except (json.JSONDecodeError, ValueError):
                domain_list = []
            domains_str = _format_interpro_domains(domain_list) if domain_list else domains_raw

            rows.append({
                "  mazF_operon_id":  _safe_str(row.get("operon_id")),
                "orf_index":         orf_idx,
                "start":             start,
                "end":               end,
                "direction":         strand,
                "identity_pct":      "",
                "length":            length,
                "sequence":          seq,
                "detection_method":  "InterPro",
                "interpro_domains":  domains_str,
            })

    return pd.DataFrame(rows)


def _format_interpro_domains(domains: list[dict]) -> str:
    """Format domain list to match reorganize_mazF_supplementary.py output."""
    parts = []
    for d in domains:
        analysis    = d.get("analysis", "")
        sig_desc    = d.get("signature_description", "")
        interpro_d  = d.get("interpro_description", "")
        interpro_a  = d.get("interpro_accession", "")
        desc = interpro_d if interpro_d and interpro_d != "-" else sig_desc
        if desc == "-":
            desc = ""
        if analysis and (desc or interpro_a):
            if desc:
                entry = f"{analysis}: {desc} ({interpro_a})" if interpro_a else f"{analysis}: {desc}"
            else:
                entry = f"{analysis}: {interpro_a}"
            parts.append(entry)
    return "; ".join(parts)


def _build_alr_sequences(df: pd.DataFrame, coords: dict) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        adjacent_raw = _safe_str(row.get("adjacent_genes_data")) or "[]"
        try:
            adjacent = json.loads(adjacent_raw)
        except (json.JSONDecodeError, ValueError):
            continue

        alr = next((g for g in adjacent if g.get("gene_name", "").lower() == "alr"), None)
        if alr is None:
            continue

        refseq   = _safe_str(row.get("refseq"))
        fragment = _safe_str(row.get("fragment"))
        orf_idx  = alr.get("orf_index", "")
        seq      = _safe_str(alr.get("sequence"))
        if not seq or seq == "nan":
            continue

        start, end, strand = _get_coords(coords, refseq, fragment, orf_idx)
        rows.append({
            "  mazF_operon_id": _safe_str(row.get("operon_id")),
            "orf_index":        orf_idx,
            "start":            start,
            "end":              end,
            "direction":        strand,
            "identity_pct":     alr.get("identity", ""),
            "length":           alr.get("length", ""),
            "sequence":         seq,
        })
    return pd.DataFrame(rows)


def _build_16s_rrna(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build 16S rRNA sheet from 16s_cutsite.tsv (keyed by refseq).
    One row per operon for refseqs that have a 16S record.
    Column name 'operon_id' (no spaces) matches the reference exactly.
    """
    if not os.path.exists(CUTSITE_TSV):
        print(f"  NOTE: {CUTSITE_TSV} not found — 16S_rRNA sheet will be empty.")
        return pd.DataFrame(columns=["operon_id", "refseq",
                                     "TACATA_count", "sequence_length", "sequence"])

    cs = pd.read_csv(CUTSITE_TSV, sep="\t", dtype=str, low_memory=False)
    cs.columns = [c.strip() for c in cs.columns]

    refseq_col = next(
        (c for c in cs.columns if c.lower() in ("refseq", "assembly")), None
    )
    if refseq_col is None:
        print(f"  WARNING: No refseq column in {CUTSITE_TSV}. 16S_rRNA sheet will be empty.")
        return pd.DataFrame(columns=["operon_id", "refseq",
                                     "TACATA_count", "sequence_length", "sequence"])

    cs = cs.rename(columns={refseq_col: "refseq"})
    needed = [c for c in ["refseq", "TACATA_count", "sequence_length", "sequence"]
              if c in cs.columns]
    cs_lookup = cs[needed].drop_duplicates("refseq").set_index("refseq").to_dict("index")

    rows = []
    for _, row in df.iterrows():
        refseq = _safe_str(row.get("refseq"))
        if refseq not in cs_lookup:
            continue
        rec = cs_lookup[refseq]
        rows.append({
            "operon_id":       _safe_str(row.get("operon_id")),
            "refseq":          refseq,
            "TACATA_count":    rec.get("TACATA_count", ""),
            "sequence_length": rec.get("sequence_length", ""),
            "sequence":        rec.get("sequence", ""),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Excel writer
# ---------------------------------------------------------------------------

def _write_sheet(ws, df: pd.DataFrame, bold_header: bool = True) -> None:
    for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), 1):
        for c_idx, value in enumerate(row, 1):
            cell = ws.cell(row=r_idx, column=c_idx, value=value)
            if r_idx == 1 and bold_header:
                cell.font = Font(bold=True)
                cell.fill = PatternFill("solid", fgColor="D9E1F2")


def _auto_width(ws, df: pd.DataFrame, max_width: int = 50) -> None:
    for idx, col in enumerate(df.columns, 1):
        col_letter = ws.cell(row=1, column=idx).column_letter
        max_len = max(
            len(str(col)),
            max((len(str(v)) for v in df[col].astype(str).head(100)), default=0),
        )
        ws.column_dimensions[col_letter].width = min(max_len + 2, max_width)


def _readme_sheet(ws) -> None:
    ws.title = "README"
    content = [
        ["Supplementary Data: MazF/MazE Operon Analysis"],
        [""],
        ["This workbook contains 6 sheets:"],
        ["Sheet",           "Description"],
        ["README",          "This sheet — describes workbook contents"],
        ["Operon_Summary",  "One row per genome: taxonomy, structure, detection flags"],
        ["MazF_Sequences",  "MazF protein sequences with genomic coordinates"],
        ["MazE_Sequences",  "MazE/NdoAI sequences with detection method"],
        ["Alr_Sequences",   "Alr (alanine racemase) protein sequences"],
        ["16S_rRNA",        "16S rRNA sequences with TACATA cutsite counts"],
        [""],
        ["is_canonical",    "TRUE only when the exact [alr]-2 → [mazE]-1 → [mazF]0 architecture is confirmed"],
        ["mazE_detection_method", "'homology' = direct BLASTp hit; 'InterPro' = RHH domain annotation"],
        ["interpro_domains","Semicolon-separated InterPro domain annotations (MazE_Sequences only)"],
        ["TACATA_count",    "Number of MazF cutsites (TACATA) in the 16S rRNA sequence"],
    ]
    for row_data in content:
        ws.append(row_data)
    ws["A1"].font = Font(bold=True, size=14)
    ws.column_dimensions["A"].width = 25
    ws.column_dimensions["B"].width = 80


def _write_excel(out_path, summary_df, mazf_df, maze_df, alr_df, rrna_df) -> None:
    wb = Workbook()
    ws_readme = wb.active
    _readme_sheet(ws_readme)

    for sheet_name, df in [
        ("Operon_Summary", summary_df),
        ("MazF_Sequences", mazf_df),
        ("MazE_Sequences", maze_df),
        ("Alr_Sequences",  alr_df),
        ("16S_rRNA",       rrna_df),
    ]:
        ws = wb.create_sheet(sheet_name)
        if len(df) > 0:
            _write_sheet(ws, df)
            _auto_width(ws, df)
        else:
            ws.append([f"No {sheet_name} data found"])

    wb.save(out_path)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _normalize_for_validation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the reference dataframe for comparison against our output.

    The reference uses 'ndoAI' in operon_structure strings; we use 'mazE'.
    Replace ndoAI → mazE in operon_structure so assert_frame_equal passes.
    """
    df = df.copy()
    if "operon_structure" in df.columns:
        df["operon_structure"] = (
            df["operon_structure"]
            .astype(str)
            .str.replace("ndoAI", "mazE", regex=False)
            .replace("nan", "")
        )
    return df


def validate(ref_path: str, new_path: str) -> bool:
    """
    Sheet-by-sheet comparison using pandas.testing.assert_frame_equal.

    - ndoAI is normalised to mazE in the reference before comparison.
    - Columns present only in the pipeline output (e.g. is_canonical) are
      excluded from the comparison.
    - README sheet is skipped.

    Returns True if all shared columns in all sheets match.
    """
    print(f"\n=== Validation ===")
    print(f"  Reference: {ref_path}")
    print(f"  Pipeline:  {new_path}")

    if not os.path.exists(ref_path):
        print(f"  WARNING: Reference file not found — cannot validate.")
        return False
    if not os.path.exists(new_path):
        print(f"  ERROR: Pipeline output not found — run Step 7 first.")
        return False

    ref_sheets = pd.read_excel(ref_path, sheet_name=None, dtype=str)
    new_sheets = pd.read_excel(new_path, sheet_name=None, dtype=str)

    all_ok = True
    for sheet_name, ref_df in ref_sheets.items():
        if sheet_name == "README":
            continue
        if sheet_name not in new_sheets:
            print(f"  MISSING sheet: '{sheet_name}'")
            all_ok = False
            continue

        new_df = new_sheets[sheet_name]

        # Columns only in reference → report as missing
        # Columns only in pipeline → excluded (e.g. is_canonical)
        shared_cols = [c for c in ref_df.columns if c in new_df.columns]
        extra_ref   = [c for c in ref_df.columns if c not in new_df.columns]
        extra_new   = [c for c in new_df.columns if c not in ref_df.columns]
        if extra_ref:
            print(f"    INFO: cols in reference but missing from pipeline: {extra_ref}")
        if extra_new:
            print(f"    INFO: pipeline-only cols (not compared): {extra_new}")

        # Normalize ndoAI → mazE in reference for operon_structure comparison
        ref_norm = _normalize_for_validation(ref_df[shared_cols])
        new_norm = new_df[shared_cols].copy()
        if "operon_structure" in new_norm.columns:
            new_norm["operon_structure"] = (
                new_norm["operon_structure"].astype(str).replace("nan", "")
            )

        try:
            pd.testing.assert_frame_equal(
                ref_norm.reset_index(drop=True).fillna(""),
                new_norm.reset_index(drop=True).fillna(""),
                check_like=True,
                check_dtype=False,
            )
            print(f"  ✓ Sheet '{sheet_name}' matches "
                  f"({len(ref_df)} rows × {len(shared_cols)} shared columns)")
        except AssertionError as exc:
            print(f"  ✗ Sheet '{sheet_name}' DIFFERS:")
            for line in str(exc).split("\n")[:12]:
                print(f"      {line}")
            all_ok = False

    if all_ok:
        print("\n  All sheets match — pipeline output validated successfully.")
    else:
        print("\n  Validation FAILED — see differences above.")
    return all_ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(force: bool = False) -> None:
    if os.path.exists(OUTPUT_EXCEL) and not force:
        print(f"  Step 8 output already exists ({OUTPUT_EXCEL}). Skipping.")
        return

    if not os.path.exists(TAXONOMY_CSV):
        print(f"ERROR: {TAXONOMY_CSV} not found. Run Step 7 first.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(TAXONOMY_CSV, low_memory=False)
    total = len(df)
    print(f"\nStep 8: Generating supplementary Excel from {total} rows ...")

    coords = _load_faa_coords(FAA_DIR)

    summary_df = _build_summary(df)
    mazf_df    = _build_mazf_sequences(df, coords)
    maze_df    = _build_maze_sequences(df, coords)
    alr_df     = _build_alr_sequences(df, coords)
    rrna_df    = _build_16s_rrna(df)

    print(f"\n  Sheet row counts:")
    print(f"    Operon_Summary:  {len(summary_df)}")
    print(f"    MazF_Sequences:  {len(mazf_df)}")
    print(f"    MazE_Sequences:  {len(maze_df)}")
    print(f"    Alr_Sequences:   {len(alr_df)}")
    print(f"    16S_rRNA:        {len(rrna_df)}")
    print(f"\n  Canonical operons: {summary_df['is_canonical'].sum()}")

    _write_excel(OUTPUT_EXCEL, summary_df, mazf_df, maze_df, alr_df, rrna_df)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate 6-sheet supplementary Excel and/or validate output."
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--validate", action="store_true",
                        help="Compare output against reference sup_new_fig5.xlsx.")
    args = parser.parse_args()

    if args.validate:
        validate(REFERENCE_SUPFIG, OUTPUT_EXCEL)
    else:
        run(force=args.force)
        if os.path.exists(REFERENCE_SUPFIG):
            print()
            validate(REFERENCE_SUPFIG, OUTPUT_EXCEL)


if __name__ == "__main__":
    main()
