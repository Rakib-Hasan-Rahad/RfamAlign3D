# Methods in detail

The pipeline runs in five steps. Each step reads only the frozen outputs of the previous one, writes to its own
folder under `outputs/`, and is resumable. Script names are given so every statement can be checked against code.

## Step 1 — Inventory of Rfam 15.1 seed alignments (`scripts/step1_inventory.py`)
Parses all 4,227 families of `Rfam.seed.gz` (count checked against the release README) and records, per family and
seed row, the ungapped sequence, its sha256, the SS_cons and RF lines, and every explicit seed-row → PDB-chain link
(`#=GR <PDB>_<chain>_SS`). Free-text PDB mentions are kept as leads. Output: `outputs/01_inventory/`
(`validation.json` lists the integrity checks).

## Step 2 — Structures, verification and pairs (`scripts/step2/`)
| Stage | Script | Does |
|---|---|---|
| s01 | `s01_discover.py` | Candidate chains from explicit tags, Step 1 leads and Rfam's family→PDB mapping (`Rfam.pdb.gz` ∪ `pdb_full_region`) |
| s02 | `s02_metadata_yield.py` | Exact seed-row search in deposited sequences, species from the RNA entity's source records |
| s03 | `s03_sequence_evidence.py` | Natural-origin evidence: genome assembly search (prokaryotes) or NCBI BLAST, each hit confirmed by exact string search of the fetched record |
| s04 | `s04_coordinates.py` | mmCIF parsing (gemmi): first model, first altloc, verified interval, modified-residue parents, backbone completeness |
| s05 | `s05_select_prepare.py` | Origin decision per chain, eligible pairs, ranking, staging of STAR3D-readable PDB files |
| s06 | `s06_export_validate.py` | Export and 10 acceptance gates (round-trip of every staged file through STAR3D's fixed-column reader rules) |
| s07 | `s07_snapshot.py` | Frozen snapshot with sha256 manifest |
Independent cross-check: `exact_local_verify.py` searched versioned NCBI genome assemblies and transcripts locally for exact
whole-chain matches (`outputs/verification_archive/`); `import_exact_local.py` re-confirmed every hit against the NCBI record
before it was used.

## Step 3 — STAR3D (`scripts/step3/`)
`build_jobs.py` writes one job per pair; `step3_runner.py` runs STAR3D v1.2 (`Preprocess`, then
`STAR3D.jar -o -p` at default parameters) in a fresh copy of the STAR3D folder per run, either natively on Linux x86_64 or in the
`star3d-runtime:1` Docker image (`env/Dockerfile.star3d`, Ubuntu 22.04 + OpenJDK 8, linux/amd64). 540 runs ran natively and 30 in Docker; `run_provenance.tsv`
records which. When the cohort was extended (R8), only new pairs were run; earlier runs were reused only if their
recorded identity (direction, input checksums, STAR3D jar, parameters, command) equalled the new job's (`assemble_returned.py`). A QC rerun of one pair (RF00522) in Docker reproduced the native result: the same 32 mappings and RMSD, with byte-identical preprocessing files. Every returned file is listed in `RETURN.sha256`.

## Step 4 — Comparison (`scripts/step4/`)
- `s41_compare.py`: correspondences, outcome per nucleotide (both sides), pair scores, disagreement regions.
- `s42_motifs.py`: BGSU RNA 3D Hub loops, Motif Atlas 4.16 groups, Rfam 15.1 motif hits at each disagreement.
- `s43_diagnostics.py`: offset along RNA B, nearness to seed gaps, mutually exclusive patterns.
- `s44_bgsu_context.py`: hairpin / internal / junction loop context of every assessable nucleotide.
- `s45_overall.py`: headline metric and internal consistency checks (all must be true).
- `s46_statistics.py`: pair-bootstrap intervals, family-balanced rates, leave-one-family-out, Mantel–Haenszel.
- `audit_interim.py`: independent audit questions Q1–Q8 (`outputs/04_comparison/<cohort>/audit/`).

## Step 5 — Evidence, database and browser (`scripts/step5/`)
- `build_step5.py`: SQLite database of every table, summary report.
- `build_app.py`: per-pair records (C1′ distances in STAR3D's superposition, MC-Annotate base pairs of both chains)
  and the interactive browser (`app/`, served by `serve.py`).
- `s52_basepair_preservation.py`: for each base pair of RNA A, what each alignment maps it to.
- `s53_candidate_corrections.py`: verdict per disagreement region from interactions and distances.

## Definitions
- **Correspondence (Rfam)**: two residues whose seed rows share an alignment column.
- **Correspondence (STAR3D)**: a mapping line `A:<res> <-> B:<res>` in STAR3D's output.
- **Outcome of an RNA A nucleotide**: same partner, different partner, Rfam only, STAR3D only, neither.
- **Assessable**: the nucleotide and each proposed counterpart have coordinates.
- **Offset**: position of the STAR3D counterpart minus position of the Rfam counterpart along ungapped RNA B.
- **Kept base pair**: both mapped nucleotides form a base pair of the same Leontis–Westhof class, compared in the
  same nucleotide order.
- **STAR3D RMSD**: STAR3D's own value, Kabsch superposition over per-residue centroids of C3′, C4′, C5′, O3′, O5′, P
  (checked in STAR3D's `PDB.java`/`STAR3D.java`). **C1′ distances** are measured in that superposition without refitting.
