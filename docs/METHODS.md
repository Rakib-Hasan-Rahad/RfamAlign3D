# Methods in detail

The pipeline runs in six steps. Each step reads only the frozen outputs of the previous one, writes to its own
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
- `s52_basepair_preservation.py`, `s53_candidate_corrections.py`: first versions of the interaction analysis, kept for
  provenance and **superseded by Step 6** (s53 counted unmapped endpoints as gains and let a distance measured in STAR3D's
  frame veto interaction evidence; see `docs/DECISIONS.md`, A4–A5).

## Definitions
- **Correspondence (Rfam)**: two residues whose seed rows share an alignment column.
- **Correspondence (STAR3D)**: a mapping line `A:<res> <-> B:<res>` in STAR3D's output.
- **Outcome of an RNA A nucleotide**: same partner, different partner, Rfam only, STAR3D only, neither.
- **Assessable**: the nucleotide and each proposed counterpart have coordinates.
- **Offset**: position of the STAR3D counterpart minus position of the Rfam counterpart along ungapped RNA B.
- **Preserved base pair**: see Step 6 (three definitions; the class-preserving one is used for headline numbers).
- **STAR3D RMSD**: STAR3D's own value, Kabsch superposition over per-residue centroids of C3′, C4′, C5′, O3′, O5′, P
  (checked in STAR3D's `PDB.java`/`STAR3D.java`). **C1′ distances** are measured in that superposition without refitting.

## Step 6 — Origin reconciliation and interaction evidence (`scripts/step6/`)
| Script | Does | Output (`outputs/06_audit/step2_r8/`) |
|---|---|---|
| `s61_origin_reconciliation.py` | Every origin check ever cached, its state at R7 and R8, its disposition, and the kind of record behind each cohort chain | `origin_check_reconciliation.tsv`, `origin_summary.json` |
| `s62_fr3d_annotate.py` | FR3D (fr3d-python, commit 288f98c, as pinned in `config.yaml`) base pairs of every cohort chain from the deposited mmCIF, first model | `fr3d/` (+ sha256 manifest) |
| `s63_interactions_regions.py` | Interaction preservation and region support (rules in `rules.py`) | `interactions_long.tsv.gz`, `interaction_preservation_by_family.tsv`, `interaction_summary.json`, `region_support_audit.tsv`, `region_summary.json` |
| `s64_correction_case.py` | Row-level correction candidate for RF00001 and its tests | `correction_case/` |
| `s65_figures.py` | Manuscript figures | `figures/` |
| `s66_claims.py` | Every headline claim with count, denominator, definition, source and status | `claim_evidence.tsv` |

**Annotators.** MC-Annotate output is the file STAR3D's own preprocessing writes; STAR3D's loop alignment rewards
matching these W/H/S-edge pairs (`LoopAlign.multi_pairing_align`), so an evaluation that uses only MC-Annotate is not
independent of STAR3D. FR3D is run separately on the deposited structures and is the independent annotator.

**Endpoint states** of each mapped source nucleotide, kept separate: `src_no_coords`, `unmapped`, `target_no_coords`
(also: target outside the analysed chain/model), `mapped`. Residues are identified by PDB, model, chain, number and
insertion code.

**Preservation** of a source pair (u, v) whose counterparts are t, w:
A *any pair*: t–w carries any annotated base pair. B *same class*: t–w carries the same Leontis–Westhof class compared
in mapped order (a class stored as w–t is reversed: cWH ↔ cHW). C *class + identity*: B and identical nucleotides.
Headline numbers use B.

**Types** (from the source pair; parent bases for modified nucleotides): canonical = cWW A–U/G–C; G·U wobble = cWW G–U;
noncanonical = every other LW class and cWW between other bases.

**Comparison sets.** Paired comparisons use the common eligible set (both nucleotides mapped under both methods);
coverage and intention-to-treat rates (unmapped = not preserved) are reported next to them. Scopes: whole compared
interval, and pairs touching a disagreement (either nucleotide different partner, Rfam only or STAR3D only).
Direction A→B (as STAR3D was run) is primary; B→A is a check and is never pooled with A→B.

**Uncertainty.** 95% percentile intervals, 2,000 replicates, seed 20261125, for the paired difference of pooled rates,
(i) resampling RNA pairs and (ii) resampling sequences (seed rows), a pair weighted by the product of its two
sequences' multiplicities. Family-specific results, family-balanced means and leave-one-family-out ranges are given;
with three informative families they are descriptive.

**Regions.** A disagreement region (Step 4) favours a method when, on the common eligible set, pairs preserved only
under that method outnumber those preserved only under the other under both annotators; *strong* when the net gain is
at least 2 under both and the other method gains none; *weak* when only one annotator shows a difference. C1′
distances (STAR3D's superposition) are reported and never decide.

**Correction candidate.** Seed-row segments that STAR3D shifts the same way against most partners (≥60% of ≥5
partners, ≥3 consecutive nucleotides) point to the row rather than to a pair. For the strongest (RF00001
mycobacterial rows) each residue is moved to the seed column chosen by ≥2/3 of partners; widths, sequences and all
other rows are unchanged. The proposal is tested without STAR3D: each FR3D pair of the mycobacterial chains is placed
at its column pair, and the number of other cohort 5S rows whose own FR3D pair of the same class sits at the same
columns (their seed placement unchanged) is compared before and after. Canonical pairing at SS_cons columns is
reported as a second check.
