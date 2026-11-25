# Step 1 — Master inventory of explicit seed-row ↔ PDB-chain links (Rfam 15.1)

*Generated 2026-11-13 by an AI agent (Claude); not yet reviewed by a person.* This step records what each seed
record says. It does not judge whether a structure is suitable for STAR3D, and it excludes no family.

## Source
`Rfam.seed.gz`, release 15.1, https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.seed.gz,
SHA-256 `41f014f4ab5628620935f089b1bb77987ebc162fda7c72bff368a32c9cccabac`. This is the same archive as the earlier
inventory, reused after a checksum match. The parser version, run time and decompressed hash are in `provenance.json`.

## Results
| | Count |
|---|---|
| Families parsed / unique accessions / release README count | 4,227 / 4,227 / 4,227 |
| Seed rows (all stored, gapped + ungapped) | 109,390 |
| Explicit row–chain links | 383 (172 rows, 306 PDB entries, 381 distinct chains) |
| Families with ≥ 1 explicit link | 74 (of these, 31 have ≥ 2 distinct linked sequences) |
| Families with no explicit link in this seed record | 4,153 (kept) |
| Unresolved PDB leads | 195 in 56 families. 1 is a row tag without a chain (`RF00051 9OBM`). 194 are free-text `#=GF CC` mentions, 16 of which name a PDB entry not explicitly linked in that family |
| Chains linked in > 1 family | 6MJ0 A and B (RF00233, RF00390), both kept |

Links are identical to the previous inventory (383 links, none added or removed). The independent RF00522 check
matches the expected counts: 43 rows, 8 linked rows, 7 distinct linked sequences, 16 chains.

## Rules
- **Explicit link:** `#=GR <row_id> <PDB>_<chain>_SS`, where the tag matches `^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]{1,4}_SS$`
  and `<row_id>` is a sequence row of the same family. The PDB ID is uppercased; the chain is kept exactly as written.
  The sequence comes from the row's own line. The `_SS` symbols are kept separately as `structure_annotation`.
- **Lead:** either (a) a row tag that looks like a PDB ID but is not `<PDB>_<chain>_SS`, or (b) a PDB-ID-like token in
  the `#=GF CC` text of a family whose text mentions "PDB". Tokens like `22nt` are excluded as lengths. Leads are
  never converted into links.
- **Feature formats in this release** (`feature_census.tsv`): there are no `#=GS` lines and no PDB references in `#=GF DR`.
  `#=GR` tags are `PP`, `<PDB>_<chain>_SS`, `9OBM` and a generic `SS` on one row (RF00668, not a PDB link).
  All `#=GC` tracks are preserved, including `RF`, `SS_cons`, `RNA_ligand_*` and `RNA_structural_elements`.
- **Ungapped** = the gapped row with only `-` and `.` removed.
  **`sequence_group_hash`** = SHA-256 of the ASCII upper-case ungapped sequence, with no T→U conversion.
  Flags per row: 0 rows contain T, 474 contain non-ACGU symbols (e.g. `N`), 0 contain lowercase residues.
- `/start-end` in the row ID gives the source-sequence interval and strand. It is not PDB numbering; that bridge is checked in Step 2.

## Files
| File | Content |
|---|---|
| `rfam15.1_step1_inventory.xlsx` | Browsing copy with sheets `About`, `Families`, `Row_PDB_links` and `PDB_leads`. All identifiers are stored as text. No cell exceeds Excel's limit; the `SEE_COMPANION_FILE` redirect exists but was not needed |
| `families.tsv`, `row_pdb_links.tsv`, `pdb_leads.tsv` | Lossless copies of the three sheets |
| `seed_rows.jsonl.gz` | One line per seed row in every family: ID, range, gapped and ungapped sequences, hash, flags, and all `#=GR` tracks |
| `family_tracks.jsonl.gz` | One line per family: width, all `#=GC` tracks, all `#=GF` fields |
| `feature_census.tsv`, `anomalies.tsv` | Markup tags found, and parser notes. There is 1 note: one `#=GF RA` line is Latin-1, not UTF-8 |
| `validation.json` | Every check (19/19 gates pass) and the release-wide counts |
| `RF00522_raw_record.sto`, `rf00522_independent_check.json` | Raw record and the separate check |

## Rerun
```
.venv/bin/python scripts/step1_inventory.py      # ~14 s; exit 0 only if all gates pass
.venv/bin/python scripts/step1_rf00522_check.py  # independent example check
.venv/bin/python -m pytest -q tests              # parser tests on invented data
```
