# Step 2 INTERIM snapshot 2026-11-13T12:03:14Z

Frozen copy of the Step 2 export validated at 2026-11-13T12:03:06Z (all 10 acceptance gates: PASS). Evidence checks were still running; STAR3D has not been run.

| Group | Pairs | Staged files | Families |
|---|---:|---:|---|
| frozen_family_complete | 535 | 131 | RF00001, RF00002, RF00005, RF00015, RF00023, RF00026, RF00169, RF00177, RF00522, RF01852, RF03087 |
| eligible_family_incomplete | 0 | 0 |  |

- **frozen_family_complete**: the family has an outcome for every sequence pair; these selections are final.
- **eligible_family_incomplete**: both members are verified and usable now, and every combination that could outrank the selection was already evaluated; but the family still has pending sequence pairs, so the final export re-derives it and reports any change against this snapshot.
- Chains with no confirmed whole-chain record (not_found) are provisional_origin and are not in any pair.

Files: `pairs_frozen.tsv`, `pairs_selection_may_change.tsv`, `rfam_pair_reference__<group>.tsv.gz`, `residue_crosswalk.tsv.gz`, `inputs_star3d/` + `input_file_hashes.tsv`, workbook, `validation_report.md`, `snapshot_info.json`, `snapshot_manifest.tsv` (SHA-256 of every file here).
