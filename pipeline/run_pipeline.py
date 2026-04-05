"""
Master orchestrator for the MazF/MazE operon detection pipeline.

Runs Steps 0–10 in sequence, with checkpointing (each step is skipped if its
output file already exists, unless --force is specified).

Detection method:
  Step 2: DIAMOND blastp → first MazF hit per genome + empty rows for no-MazF genomes
  Step 5: InterProScan RHH domain rescue (for ORFs that failed DIAMOND MazE BLASTp)
  Step 6: Merge InterProScan results; set combined_mazE_detected and is_canonical
  Step 9: Prune Bacillota taxonomy tree to genera with ≥ 15 genomes
  Step 10: Generate 7 iTOL annotation files

Usage examples:
    # Default run (uses pre-computed DIAMOND/InterProScan results)
    python pipeline/run_pipeline.py --threads 8

    # Resume from Step 9 (tree + iTOL only)
    python pipeline/run_pipeline.py --start-from 9

    # Force tree rebuild (overrides checkpointing for Step 9 only)
    python pipeline/run_pipeline.py --start-from 9 --force-tree

    # Full from-scratch run (re-runs all tools)
    python pipeline/run_pipeline.py --force --run-interproscan --threads 8

    # Validate pipeline output against reference sup_new_fig5.xlsx
    python pipeline/run_pipeline.py --validate
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline.config import OUTPUT_DIR, REFERENCE_SUPFIG, OUTPUT_EXCEL


def _banner(step: int, name: str) -> None:
    bar = "─" * 60
    print(f"\n{bar}")
    print(f"  Step {step}: {name}")
    print(f"{bar}")


def run(
    start_from: int        = 0,
    threads: int           = 4,
    force: bool            = False,
    force_tree: bool       = False,
    run_interproscan: bool = False,
) -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t0 = time.time()

    # ------------------------------------------------------------------
    # Step 0 — Fetch data from Zenodo
    # ------------------------------------------------------------------
    if start_from <= 0:
        _banner(0, "Fetch data from Zenodo")
        from pipeline.step0_fetch_data import fetch_all
        fetch_all(force=force)

    # ------------------------------------------------------------------
    # Step 1 — DIAMOND blastp
    # ------------------------------------------------------------------
    if start_from <= 1:
        _banner(1, "DIAMOND blastp search")
        from pipeline.step1_diamond_search import run as step1
        step1(threads=threads, force=force)

    # ------------------------------------------------------------------
    # Step 2 — Extract operons (one row per genome; first MazF hit only)
    # ------------------------------------------------------------------
    if start_from <= 2:
        _banner(2, "Extract MazF operon context")
        from pipeline.step2_extract_operons import run as step2
        step2(force=force)

    # ------------------------------------------------------------------
    # Step 3 — Prepare InterProScan input
    # ------------------------------------------------------------------
    if start_from <= 3:
        _banner(3, "Prepare InterProScan input FASTA")
        from pipeline.step3_prep_interproscan import run as step3
        step3(force=force)

    # ------------------------------------------------------------------
    # Step 4 — InterProScan RHH domain analysis
    # ------------------------------------------------------------------
    if start_from <= 4:
        _banner(4, "InterProScan RHH domain analysis")
    from pipeline.step5_interproscan import run as step4_interpro
    interpro_results = step4_interpro(
        run_interproscan = run_interproscan,
        force            = force if start_from <= 4 else False,
    )

    # ------------------------------------------------------------------
    # Step 5 — Merge results
    # ------------------------------------------------------------------
    if start_from <= 5:
        _banner(5, "Merge all detection results")
        from pipeline.step6_merge_results import run as step5_merge
        step5_merge(interpro_results, force=force)

    # ------------------------------------------------------------------
    # Step 6 — Add taxonomy
    # ------------------------------------------------------------------
    if start_from <= 6:
        _banner(6, "Add NCBI taxonomy")
        from pipeline.step7_add_taxonomy import run as step6_tax
        step6_tax(force=force)

    # ------------------------------------------------------------------
    # Step 7 — Generate output Excel
    # ------------------------------------------------------------------
    if start_from <= 7:
        _banner(7, "Generate supplementary Excel")
        from pipeline.step8_generate_output import run as step7_excel
        step7_excel(force=force)

    # ------------------------------------------------------------------
    # Step 8 — Prune Bacillota taxonomy tree
    # ------------------------------------------------------------------
    if start_from <= 8:
        _banner(8, "Prune Bacillota taxonomy tree (min 15 genomes/genus)")
        from pipeline.step9_build_tree import run as step8_tree
        step8_tree(force=force or force_tree)

    # ------------------------------------------------------------------
    # Step 9 — Generate iTOL annotation files
    # ------------------------------------------------------------------
    if start_from <= 9:
        _banner(9, "Generate iTOL annotation files")
        from pipeline.step10_generate_itol import run as step9_itol
        step9_itol(force=force)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    mins, secs = divmod(int(elapsed), 60)
    print(f"\n{'═' * 60}")
    print(f"  Pipeline complete in {mins}m {secs}s")
    print(f"  Output: {OUTPUT_EXCEL}")
    print(f"{'═' * 60}\n")


def validate_only() -> None:
    from pipeline.step8_generate_output import validate
    ok = validate(REFERENCE_SUPFIG, OUTPUT_EXCEL)
    sys.exit(0 if ok else 1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MazF/MazE operon detection pipeline (Steps 0–7).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--start-from", type=int, default=0, metavar="N",
        help="Resume pipeline from step N (0–9). Earlier steps are skipped.",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for DIAMOND (default: 4).",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run all steps even if output files already exist.",
    )
    parser.add_argument(
        "--run-interproscan", action="store_true",
        help="Run InterProScan locally instead of using pre-computed TSV.",
    )
    parser.add_argument(
        "--force-tree", action="store_true",
        help="Re-run Step 9 (tree pruning) even if MazEFTree_min15.txt already exists.",
    )
    parser.add_argument(
        "--validate", action="store_true",
        help="Compare pipeline output against reference sup_new_fig5.xlsx.",
    )
    args = parser.parse_args()

    if args.validate:
        validate_only()
    else:
        run(
            start_from       = args.start_from,
            threads          = args.threads,
            force            = args.force,
            force_tree       = args.force_tree,
            run_interproscan = args.run_interproscan,
        )


if __name__ == "__main__":
    main()
