# Decisions and their justification

Every rule that shapes the dataset or the analysis, when it was fixed, and why. Policy files are frozen; any
change is appended here and in `outputs/02_structure_preparation/step2_policy.yaml` (`revisions`) with the
stage that was rerun. Dates are UTC.

## Fixed before any data step (config v1, 2026-11-13)

| Rule | Value | Why |
|---|---|---|
| Reference alignment | Rfam 15.1 `Rfam.seed.gz` only (sha256 `41f014f4…`) | Seed alignments are the curated alignments; full-family and cmalign-added rows are model output, not curation |
| Sequence normalisation | uppercase, T→U, nothing else | Any further change would test something other than the published rows |
| Seed row ↔ structure | the whole ungapped seed row must occur exactly in the deposited chain | Removes any doubt that the structure is the sequence Rfam aligned |
| Natural sequences only | engineered constructs excluded; masking does not rescue them | Engineered substitutions are not part of the family's evolution and would penalise Rfam unfairly |
| Coordinates | first model, first alternate location, author numbering | Standard, reproducible choice; insertion codes kept |
| STAR3D | v1.2, default parameters, one forward run per pair, no replicates or retries | Avoids tuning STAR3D towards or away from Rfam; the query is the lexicographically smaller chain ID |

## Step 2 policy v2.0 (frozen 2026-11-13, before candidate selection)

Full text: `outputs/02_structure_preparation/step2_policy.yaml`. Key points:
- **Natural origin** requires the full deposited sequence to occur exactly, in the stated species, in an
  independent sequence record (not PDB-derived; older than the structure or ≥50 nt longer). BLAST scores never verify
  on their own; every hit is confirmed by exact string search of the fetched record.
- **Not found ≠ engineered.** Chains without a confirmed record are `provisional_origin` and are not used; they are
  never labelled engineered.
- **Pairs**: same family; both chains verified natural and technically usable; different sequences, species and
  PDB entries. One pair per unordered pair of distinct sequences, chosen by a fixed ranking (observed coverage,
  backbone completeness, reported resolution, then stable ID).
- **Coverage**: no hard cutoff (`coverage_cutoff: none`); coverage ranks candidates. This supersedes the
  `min_observed_fraction_of_row: 0.90` value in config v1. In the final cohort (`step2_r8`) 137 of 138 chains have
  ≥90% of the seed row observed (minimum 89.5%).

## Revisions R1–R8

| ID | Date | What | Why | Effect on decisions |
|---|---|---|---|---|
| R1 | 11-13 | Evidence searched once per (sequence, species); prokaryotes in one genome assembly per taxon, others by one batched BLAST per species | Execution order only | None |
| R2 | 11-13 | Parent base of modified residues: entry declaration, then entity sequence, gemmi table, then wwPDB Chemical Component Dictionary | Many entries omit the parent field (e.g. tRNA modifications) | Correct exact matching of modified residues |
| R3 | 11-13 | Organism given by scientific name without taxon ID resolved by exact NCBI name match | e.g. 3D2G (*A. thaliana*) | Species resolved for such entries |
| R4 | 11-13 | Non-native-flank exclusion applies only to RNA not isolated from the organism | 3 natural rRNA chains had been labelled engineered | 3 chains no longer excluded |
| R5 | 11-13 | Placeholder taxa (metagenome, unidentified, synthetic) never resolve species | They carry species rank in NCBI | Such chains stay unresolved |
| R6 | 11-13 | Prokaryotic keys searched locally first; BLAST only for misses; stall watchdog | Workers idle behind NCBI queue | None (execution) |
| R7 | 11-13 | Step 2 closed by the researcher; checks still lacking evidence become `SEARCH_INCOMPLETE` (not exclusion, not engineering) | NCBI BLAST stalled for hours | 39 checks left open; cohort frozen as `step2_final` |
| R8 | 11-25 | The 39 checks left open by R7 re-queried at NCBI with the unchanged rules; still-queued BLAST jobs keep their job ID (`poll_blast.py`); five species queued for >2 h searched with megablast on the same database and taxon filter (`megablast_fallback.py`; every whole-chain query ≥28 nt, every hit still confirmed by exact string search) | R7 left checks open that could change the cohort | All 234 checks resolved. 535 pairs unchanged; 35 RF00001 pairs added (4V8P_B3 *T. thermophila*, 6XU8_A7 *D. melanogaster* newly verified); only those 35 run with STAR3D, the 535 earlier runs reused after an identity check (`assemble_returned.py`). Cohort `step2_r8` |

## Analysis decisions (Steps 4–5)

| Decision | Choice | Justification |
|---|---|---|
| Correspondence under Rfam | two residues share a seed column | Definition of an MSA |
| Correspondence under STAR3D | STAR3D's raw mapping lines, resolved through the residue crosswalk | No reinterpretation of STAR3D output |
| Assessable | the nucleotide and every proposed counterpart have coordinates | Unobserved residues cannot be compared structurally |
| Headline metric | (different + Rfam only + STAR3D only) / (assessable mapped nucleotides), RNA A view | Counts every kind of difference once; `neither` is not a difference |
| Base-pair annotation | MC-Annotate, the annotator bundled with and run by STAR3D v1.2 | Same annotation STAR3D's preprocessing produces; Leontis–Westhof classes for W/H/S edges. Config v1 listed FR3D; it was not run, so all interaction results use MC-Annotate (stated in the paper) |
| Class comparison | same class compared in the same nucleotide order (cWH listed reversed = cHW) | MC-Annotate writes edges in listing order |
| Missing annotation | "no annotated interaction", never "absent" | Annotators can miss weak contacts |
| Rfam motif-scan coverage | family level only | Rfam publishes no row-level scan record |
| Uncertainty | 95% percentile intervals from 2,000 bootstrap resamples of whole pairs (seed 20261125); family-balanced rates; leave-one-family-out | Nucleotides of one pair are not independent; tRNA and 5S rRNA dominate |
| Candidate corrections (`s53`) | region verdict from base pairs kept only by one method, with distance as a veto (margin 1 Å) | Interaction evidence does not depend on STAR3D's superposition; distances do. **Rule defined on 2026-11-25 after the aggregate preservation results were seen**; the margin is reported and can be varied (`s53 <cohort> <margin>`) |
| One STAR3D direction | query = smaller stable ID | Fixed before running; a reverse-direction run was not part of the design and is listed as a limitation |

## Recorded but not used as filters
- Resolution, experimental method and revision of each entry (`structures.tsv`).
- Backbone completeness flags (`ready_with_coverage_flags`).
- Pair sequence identity over shared columns (descriptive only).

## Audit revisions (2026-11-25, step 6)
An audit of the evidence behind the manuscript claims (`outputs/06_audit/step2_r8/claim_evidence.tsv`). No input,
snapshot, STAR3D run or earlier output was changed; step 6 adds analyses and replaces two earlier claims.

| ID | What | Why | Effect |
|---|---|---|---|
| A1 | Origin checks reconciled key by key (236 ever cached, 234 current) | "All checks resolved" meant "all searches finished" | 103 exact independent support, 85 not found, 36 flank differs, 10 unknown residues; 121 without a whole-chain match are not used. Cohort chains: 72 genome assembly, 38 longer record, 28 only an older record |
| A2 | R8 megablast wording | "Finds the same hits" overstated it | Every 100%-identical full-length hit of a query ≥28 nt contains megablast's seed; every hit re-confirmed by exact string search; hitlists (50) saturated for 3 of 5 queries |
| A3 | FR3D added as independent annotator | STAR3D's loop score uses the MC-Annotate pairs that s52 evaluated | Direction of all main contrasts replicated with FR3D |
| A4 | Preservation measured on the common eligible set, three definitions, by family, with a sequence-level bootstrap | Pair resampling ignores that sequences recur in up to 18 pairs | The pooled noncanonical difference (CI 8.7–20.2 by pairs) spans zero by sequences (−4.2 to 32.2; non-cWW −7.8 to 36.7). Effect is family specific: 5S rRNA favours STAR3D, tRNA favours Rfam |
| A5 | Region verdicts of s53 withdrawn; interaction-only tiers from two annotators | s53 counted unmapped Rfam endpoints as STAR3D gains and vetoed with a distance measured in STAR3D's frame | 221/111 replaced by 218/208 (strong 74/55) |
| A6 | Row-level correction candidate RF00001 (mycobacterial helix IV) with Stockholm patch | Pairwise regions cannot justify a seed edit | Candidate for curator review; conflicting SS_cons evidence documented |
