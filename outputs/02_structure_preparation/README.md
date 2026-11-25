# Step 2 — Verified structure pairs and frozen Rfam references for original STAR3D

*Automated pipeline (Claude); no human review has occurred. STAR3D has **not** been run in this step.*
Current figures: `professor_summary.md` (one page) and `validation_report.md` (acceptance gates).
Rules: `step2_policy.yaml`, including the dated revisions R1–R4. Each revision changes execution order or clarifies evidence handling; none relaxes the seed-row match.

## Open first
| File | Content |
|---|---|
| `SQARMA_Step_2_STAR3D_inputs.xlsx` | Browsing copy with sheets `Family_status` (all 4,227 families), `Structures`, `Chain_review`, `Pairs_to_run` and `Yield_by_family`. Labelled **INTERIM** until every family has an outcome |
| `pairs_to_run.tsv` | Approved primary jobs only. Each row has both seed rows, species, PDB/chain, staged files and SHA-256, flags, the fixed query→target order, and the Step 3 commands |
| `inputs_star3d/` | One single-chain legacy-PDB file per representative, plus `<stem>.conversion.tsv` (reversible residue and chain mapping) and `file_manifest.tsv` |
| `residue_crosswalk.tsv.gz` | One row per seed nucleotide: Rfam column → ungapped position → entity/label chain/label position → author chain/number/insertion code → staged residue. Unobserved residues have no staged endpoint |
| `rfam_pair_reference.tsv.gz` | The original Rfam correspondence for every pair, column by column: RF class, both characters, residues, coordinate flags, and residue_pair / query_only / target_only / both_gap |

## Stage outputs (machine-readable, authoritative)
| Stage | Files |
|---|---|
| s01 discovery | `input_reconciliation.json` (Step 1 figures reproduced), `discovery_evidence.tsv` (every evidence record), `candidate_chains.tsv`, `candidate_intervals.tsv`, `family_discovery.tsv` |
| s02 metadata / yield | `entry_metadata.tsv` (current/obsolete; replacements recorded, never substituted), `interval_screening.tsv`, `yield_by_family.tsv`, `possible_partner.json` |
| s03 sequence evidence | `sequence_evidence.tsv`; per-key evidence in `cache/origin/`, genome assemblies in `cache/genomes/`, BLAST XML and job IDs in `cache/blast/`, record regions in `cache/ncbi/regions/` |
| s04 coordinates | parsed mmCIF summaries in `cache/mmcif_parsed/`, `wildcard_rescreen.tsv`, `validation_sample.tsv` (10 mapping-only candidates re-checked) |
| s05 selection | `chain_review.tsv`, `pair_candidates.tsv` (every combination: selected / eligible_alternate / pruned / not eligible, with reason), `origin_queue_extra.tsv`, `selection_state.json` |
| s06 export | `family_status.tsv`, `structures.tsv`, workbook, `validation.json`, `validation_report.md`, `professor_summary.md`, `pytest_results.txt` |

## Rerun / resume
```
.venv/bin/python scripts/step2/run_all.py             # full run or resume after any interruption
.venv/bin/python scripts/step2/run_all.py --from s05  # re-decide and export from cached evidence
```
Every download, BLAST job ID, genome search, mmCIF parse and chain review is cached and written atomically. Failed network work is never cached, so a rerun repeats only unfinished work. `checkpoint.json` records each stage's last completion.

## Step 3 (not executed here)
For each row of `pairs_to_run.tsv`: copy a fresh `STAR3D_source`, put both staged files in `PDB/` as `<stem>.pdb`, then run
`step3_preprocess` and `step3_align` (defaults, `-o -p`, one forward run) in `star3d-runtime:1`.
Compare the output against `rfam_pair_reference.tsv.gz` via `residue_crosswalk.tsv.gz`. Rfam's correspondence is never given to STAR3D.
