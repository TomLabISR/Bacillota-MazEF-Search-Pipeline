# MazEF Operon Detection Pipeline

Reproducible pipeline for identifying MazF toxin-antitoxin operon systems across 3,362 *Bacillota* genomes.

## Overview

The pipeline detects the canonical `alr → mazE → mazF` operon architecture using:

1. **DIAMOND blastp** — identifies MazF/MazE/alr homologues in each genome
2. **InterProScan** — rescues MazE antitoxins missed by DIAMOND via Ribbon-Helix-Helix (RHH) domain detection

Identified systems are then mapped onto the NCBI taxonomy reference tree to analyse the prevalence and phylogenetic distribution of the canonical MazEF system across the *Bacillota* phylum.

A system is **canonical** (`is_canonical = True`) if and only if alr is at position −2 and MazE is at position −1 relative to MazF (all contiguous, any strand orientation).

## Quick Start

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate mazef_pipeline

# 2. Download heavy data from Zenodo (genomes + pre-computed results)
python pipeline/run_pipeline.py --start-from 0

# 3. Run full pipeline (Steps 1–10) using pre-computed data
python pipeline/run_pipeline.py --threads 8

# 4. Outputs
#   pipeline_output/sup_new_fig5.xlsx        ← 6-sheet supplementary table
#   pipeline_output/MazEFTree_min15.txt      ← pruned Bacillota genus tree
#   pipeline_output/itol_annotations/        ← 7 iTOL annotation files
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 0 | `step0_fetch_data.py` | Download genome FAAs + pre-computed results from Zenodo |
| 1 | `step1_diamond_search.py` | DIAMOND blastp of each genome vs. 3 reference sequences |
| 2 | `step2_extract_operons.py` | Extract MazF operon context; one row per genome |
| 3 | `step3_prep_interproscan.py` | Write ORFs adjacent to MazF into a FASTA for InterProScan |
| 4–5 | `step5_interproscan.py` | Parse InterProScan RHH domain results (or run locally) |
| 6 | `step6_merge_results.py` | Merge all detections; assign `is_canonical` |
| 7 | `step7_add_taxonomy.py` | Add NCBI taxonomy |
| 8 | `step8_generate_output.py` | Generate `sup_new_fig5.xlsx` |
| 9 | `step9_build_tree.py` | Prune Bacillota taxonomy tree to genera with ≥15 genomes |
| 10 | `step10_generate_itol.py` | Generate iTOL annotation files from pipeline output |

## Running InterProScan from Scratch (Optional)

To re-run InterProScan from raw sequences (this process might take a while):

```bash
python pipeline/run_pipeline.py --force --run-interproscan --threads 8
```

## Reference Sequences (`ref_protein_mazEF.faa`)

| Sequence ID | Gene | Organism |
|-------------|------|----------|
| `ndoA_NP_388347.1` | MazF toxin | *B. subtilis* str. 168 |
| `ndoAI_WP_144530614.1` | MazE antitoxin | *B. subtilis* str. 168 |
| `alr_WP_003234284.1` | Alr (alanine racemase) | *B. subtilis* str. 168 |

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| DIAMOND e-value | ≤ 1×10⁻⁵ | DIAMOND blastp cutoff |
| Min MazF length | 85 aa | Filter fragment hits |
| Max MazE length | 150 aa | InterProScan false-positive filter |
| Adjacency window | ≤ 2 ORFs | alr / MazE must be within 2 ORF positions of MazF |

## Output (`pipeline_output/sup_new_fig5.xlsx`)

| Sheet | Contents |
|-------|----------|
| README | Column descriptions |
| Operon_Summary | One row per genome: taxonomy, structure, detection flags |
| MazF_Sequences | MazF protein sequences with genomic coordinates |
| MazE_Sequences | MazE sequences with detection method (homology / InterPro) |
| Alr_Sequences | Alr sequences with genomic coordinates |
| 16S_rRNA | 16S rRNA sequences with TACATA cutsite counts |

## Data Availability

Pre-computed genome FAAs (Prodigal), DIAMOND results, and InterProScan TSV are hosted on Zenodo:  
**DOI: [10.5281/zenodo.19420781](https://doi.org/10.5281/zenodo.19420781)**

## Requirements

- Python ≥ 3.11
- DIAMOND ≥ 2.1
- ETE 3 (required for Step 9 tree pruning)
- InterProScan (optional; only needed with `--run-interproscan`)

See `environment.yml` for the full conda environment.

## References / Tools Used

- **DIAMOND** — Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods*, 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x
- **InterProScan** — Blum M, et al. The InterPro protein families and domains database: 20 years on. *Nucleic Acids Research*, 49(D1), D344–D354 (2021). https://doi.org/10.1093/nar/gkaa977
- **iTOL** — Letunic I, Bork P. Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. *Nucleic Acids Research*, 52(W1), W78–W82 (2024). https://doi.org/10.1093/nar/gkae268
- **Prodigal** — Hyatt D, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*, 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119
- **NCBI Taxonomy** — Schoch CL, et al. NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database*, 2020, baaa062 (2020). https://doi.org/10.1093/database/baaa062
- **ETE 3** — Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. *Molecular Biology and Evolution*, 33(6), 1635–1638 (2016). https://doi.org/10.1093/molbev/msw046
