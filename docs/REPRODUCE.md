# Reproducing and checking the results

Three levels, from a five-minute check to a full rerun.

## 1. Check the published results (no installation)
```bash
git clone https://github.com/Rakib-Hasan-Rahad/RfamAlign3D.git && cd RfamAlign3D
python3 scripts/verify.py
```
Needs only Python 3 (standard library). It verifies the sha256 of every file in the Step 2 snapshot and recounts the
headline numbers from the raw comparison table. Checks that need the release archives report `SKIPPED` until you
unpack them (section 2); a skipped check is never counted as passed.

## 2. Add the large files (release archives)
Download the archives from the GitHub release of the version you are checking and unpack them in the repository root:

| Archive | Contents | Unpacks to |
|---|---|---|
| `star3d_runs.tar.gz` | Raw STAR3D alignments, superpositions and preprocessing for all 570 runs | `outputs/03_star3d/step2_r8/returned/` |
| `verification_archive.tar.gz` | Exact local search evidence and STAR3D run outputs as collected | `outputs/verification_archive/` |
| `browser_and_database.tar.gz` | Interactive browser and SQLite database | `outputs/05_results/step2_r8/` |
| `pinned_inputs.tar.gz` | The exact Rfam 15.1, wwPDB and NCBI taxonomy files used (STAR3D excluded, see below) | `inputs/` |
| `evidence_blast.tar.gz` | NCBI BLAST XML behind every natural-origin decision | `cache/blast/` |

Then rerun `python3 scripts/verify.py`; every check should pass.

Browse the results: `python3 scripts/step5/serve.py outputs/05_results/step2_r8/app 8770` and open
<http://127.0.0.1:8770>.

## 3. Rerun the analysis
```bash
bash setup.sh            # Python 3.13 venv, Docker image for STAR3D, pinned downloads (sha256-checked), environment check
```
STAR3D v1.2 is downloaded from its authors' site (genome.ucf.edu) and checked against the sha256 in `config.yaml`;
it is not redistributed here.

Recompute everything after STAR3D from the published files (minutes):
```bash
.venv/bin/python scripts/step4/s41_compare.py step2_r8
.venv/bin/python scripts/step4/s42_motifs.py step2_r8
.venv/bin/python scripts/step4/s43_diagnostics.py step2_r8
.venv/bin/python scripts/step4/s44_bgsu_context.py step2_r8
.venv/bin/python scripts/step4/audit_interim.py step2_r8
.venv/bin/python scripts/step4/s45_overall.py step2_r8
.venv/bin/python scripts/step4/s46_statistics.py step2_r8
.venv/bin/python scripts/step5/build_step5.py step2_r8
.venv/bin/python scripts/step5/build_app.py step2_r8
.venv/bin/python scripts/step5/s52_basepair_preservation.py step2_r8
.venv/bin/python scripts/step5/s53_candidate_corrections.py step2_r8
.venv/bin/python -m pytest -q tests
```

Full rerun from Rfam (hours to days; depends on NCBI):
```bash
.venv/bin/python scripts/step1_inventory.py
.venv/bin/python scripts/step2/run_all.py           # resumable; NCBI searches are cached per (sequence, species)
.venv/bin/python scripts/step2/s07_snapshot.py snapshot    # freeze the selection as a cohort snapshot
.venv/bin/python scripts/reconcile.py <old_cohort> <new_cohort>   # which pairs can reuse an earlier STAR3D run
.venv/bin/python scripts/step3/build_jobs.py step2_r8
# in the job bundle written by build_jobs.py (cluster), or locally with STAR3D_DOCKER_IMAGE=star3d-runtime:1:
python3 step3_runner.py setup && python3 step3_runner.py run [--only <pair_id> ...]
.venv/bin/python scripts/step3/assemble_returned.py step2_r8 step2_final <copied_runs_dir>   # merge reused + new runs, identity-checked
```
Online sources change over time (PDB revisions, NCBI records, BGSU releases). The caches and archives above hold the
exact copies used, so a rerun against them reproduces the published numbers; a rerun against live services may
differ and should be compared with `docs/results.json`.

## Tests
`pytest tests` (29 tests) covers the Step 1 parser, Step 2 matching/origin/taxonomy rules and the Step 4 comparison
logic, including mutation tests of known edge cases. The taxonomy test needs `inputs/` (section 2 or `setup.sh`).
