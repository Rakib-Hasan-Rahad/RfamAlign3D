# Data dictionary

Columns are read from the files themselves (`docs/make_data_dictionary` step in `docs/REPRODUCE.md`). Cohort: `step2_r8`.
Unit IDs follow the BGSU convention `PDB|model|chain|nucleotide|number|||insertion_code` (unused fields empty).
Outcome categories (RNA A view): `same_partner`, `different_partner`, `rfam_only` (STAR3D leaves it unmapped), `star3d_only` (it faces a seed gap), `neither`.

## `outputs/01_inventory/families.tsv`
Step 1. One row per Rfam 15.1 family: seed size, explicit PDB links found in the seed record.
*4,227 rows.* Columns: `family`, `family_id`, `description`, `rfam_type`, `alignment_width`, `total_seed_rows`, `distinct_row_sequences`, `linked_rows`, `distinct_linked_sequences`, `distinct_pdb_entries`, `distinct_pdb_chains`, `unresolved_leads`, `link_status`

## `outputs/01_inventory/row_pdb_links.tsv`
Step 1. Every explicit seed-row to PDB-chain link (#=GR <PDB>_<chain>_SS).
*383 rows.* Columns: `family`, `family_id`, `seed_row_id`, `source_accession`, `source_start`, `source_end`, `source_strand`, `pdb_id`, `chain_id`, `evidence_type`, `raw_tag`, `sequence_group_hash`, `ungapped_length`, `ungapped_sequence`, `original_gapped_sequence`, `structure_annotation`

## `outputs/02_structure_preparation/snapshots/step2_r8/family_status.tsv`
Step 2. Outcome for every family (why it did or did not yield pairs).
*4,227 rows.* Columns: `family`, `family_id`, `discovery_status`, `outcome`, `family_complete`, `evidence_records`, `candidate_entries`, `candidate_chains`, `exact_unambiguous_intervals`, `distinct_matched_sequences`, `screening_seq_pairs_diff_species`, `screening_seq_pairs_species_unresolved`, `chains_reviewed`, `verified_natural`, `verified_but_technical_unusable`, `pairs_to_run`, `seq_pairs_pending`, `principal_reasons`, `obsolete_entries_not_substituted`

## `outputs/02_structure_preparation/snapshots/step2_r8/chain_review.tsv`
Step 2. One row per candidate chain interval: sequence match, species, origin decision and reason, technical checks.
*2,333 rows.* Columns: `candidate_id`, `family`, `pdb_id`, `auth_chain`, `label_asym_id`, `entity_id`, `seed_row_id`, `matching_rows`, `screening_status`, `seqres_interval`, `mmcif_interval_confirmed`, `species_taxid`, `species_name`, `species_basis`, `decision`, `reason`, `origin_evidence`, `technical_status`, `technical_problems`, `observed_fraction`, `backbone_complete_fraction`, `interval_length`, `cropped_5prime`, `cropped_3prime`, `author_number_range`, `author_numbers_monotonic`, `backbone_complete`, `decision_mode`, `decision_utc`, `deposited_equals_seqres`, `deposited_length`, `deposited_parent_sha256`, `doc_description`, `doc_details`, `doc_documented_engineering`, `doc_engineering_words`, `doc_mutation_documented`, `doc_pdbx_mutation`, `doc_src_method`, `evidence_key`, `expression_hosts`, `insertion_codes`, `longest_internal_missing`, `method`, `mmcif_row_occurrences`, `mmcif_sha256`, `mmcif_status`, `modified_residues`, `observed`, `other_rows_inside_interval`, `residues_missing_backbone_atoms`, `resolution`, `rna_sources`, `screen_species`, `species_source_taxa`, `unobserved_residues`

## `outputs/02_structure_preparation/snapshots/step2_r8/structures.tsv`
Step 2. Every PDB entry used: method, resolution, revision, mmCIF checksum.
*1,035 rows.* Columns: `pdb_id`, `rcsb_page`, `download_url`, `mmcif`, `mmcif_sha256`, `mmcif_bytes`, `method`, `resolution`, `revision`, `revision_date`, `deposit_date`, `models`, `title`, `status`, `retrieved`

## `outputs/02_structure_preparation/snapshots/step2_r8/pairs_to_run.tsv`
Step 2. The 535 selected pairs with both chains, species and staged files.
*570 rows.* Columns: `pair_id`, `family`, `family_id`, `query_candidate`, `target_candidate`, `query_seed_row`, `target_seed_row`, `query_species`, `target_species`, `query_pdb_chain`, `target_pdb_chain`, `query_label_interval`, `target_label_interval`, `query_stem`, `query_chain`, `query_file`, `query_sha256`, `target_stem`, `target_chain`, `target_file`, `target_sha256`, `query_technical`, `target_technical`, `query_origin`, `target_origin`, `ranking_key`, `pair_identity`, `flags`, `step3_preprocess`, `step3_align`, `star3d_run`

## `outputs/02_structure_preparation/snapshots/step2_r8/residue_crosswalk.tsv.gz`
Step 2. Residue crosswalk: seed-row position and column <-> PDB unit ID (author numbering, insertion code) <-> staged STAR3D residue.
*17,514 rows.* Columns: `candidate_id`, `family`, `seed_row_id`, `rfam_column`, `row_char`, `ungapped_position`, `pdb_id`, `entity_id`, `label_asym_id`, `label_seq_id`, `auth_chain`, `auth_seq_num`, `auth_ins_code`, `comp_id`, `parent`, `parent_source`, `modified`, `model`, `observed`, `missing_backbone_atoms`, `staged_file`, `staged_chain`, `staged_resid`, `transformation`

## `outputs/03_star3d/step2_r8/returned/run_provenance.tsv`
Step 3. Where and how each STAR3D run was executed (natively on Linux or in Docker), tool versions.
*570 rows.* Columns: `pair_id`, `source`, `execution_environment`, `status`, `identity_matches_final_job`

## `outputs/04_comparison/step2_r8/run_outcomes.tsv`
Step 4. One row per pair: STAR3D status (valid_alignment / no_alignment), RMSD, comparison status.
*570 rows.* Columns: `pair_id`, `family`, `step3_status`, `declared`, `rmsd`, `superposition`, `aln_sha256`, `step4_status`, `problems`

## `outputs/04_comparison/step2_r8/residue_comparisons.tsv.gz`
Step 4. One row per nucleotide per side per pair: Rfam and STAR3D counterparts, outcome category, coordinates flags. The core table.
*99,570 rows.* Columns: `pair_id`, `side`, `unit_id`, `ungapped_position`, `nucleotide`, `comp_id`, `rfam_column`, `rf_class`, `ss_char`, `ss_context`, `coords`, `category`, `rfam_partner`, `rfam_partner_column`, `rfam_partner_coords`, `star3d_partner`, `star3d_partner_column`, `raw_line`, `scored_different`

## `outputs/04_comparison/step2_r8/pair_scores.tsv`
Step 4. One row per aligned pair: shared / Rfam-only / STAR3D-only correspondences and reproduction rates.
*538 rows.* Columns: `pair_id`, `family`, `raw_rfam_edges`, `raw_star3d_edges`, `raw_shared`, `rfam_edges`, `star3d_edges`, `shared_edges`, `rfam_recovery`, `star3d_agreement`, `rfam_only_edges`, `star3d_only_edges`, `rf_available`, `ss_cons_parse`, `rfam_consensus_reproduced`, `rfam_insert_reproduced`, `rfam_unknown_reproduced`, `star3d_both_consensus_in_rfam`, `star3d_mixed_consensus_insert_in_rfam`, `star3d_both_insert_in_rfam`, `star3d_unknown_involved_in_rfam`, `consensus_only_rfam_recovery`, `consensus_only_star3d_agreement`

## `outputs/04_comparison/step2_r8/regions.tsv`
Step 4. Disagreement regions: runs of consecutive different-partner nucleotides.
*2,323 rows.* Columns: `region_id`, `pair_id`, `side`, `length`, `start_position`, `end_position`, `start_unit`, `end_unit`, `rf_composition`, `ss_composition`

## `outputs/04_comparison/step2_r8/diagnostics/disagreement_diagnostics.tsv.gz`
Step 4. Each different-partner nucleotide: offset, nearness to a seed gap, region, pattern.
*4,676 rows.* Columns: `pair_id`, `unit_id`, `position`, `rfam_column`, `rf_class`, `ss_context`, `offset`, `near_seed_gap`, `region_id`, `pair_disagreement`, `region_length`, `region_offsets`, `pattern`

## `outputs/04_comparison/step2_r8/motif/bgsu_loop_context.tsv.gz`
Step 4. BGSU loop membership (hairpin/internal/junction) of each assessable A-side nucleotide.
*44,079 rows.* Columns: `pair_id`, `family`, `unit_id`, `outcome`, `loops`, `loop_types`, `category`

## `outputs/04_comparison/step2_r8/motif/disagreement_observations.tsv`
Step 4. Motif evidence class for each different-partner nucleotide (BGSU loop, Motif Atlas group, Rfam motif hit).
*4,676 rows.* Columns: `observation_id`, `pair_id`, `family`, `a_residue`, `b_rfam_partner`, `b_star3d_partner`, `rfam_column`, `rf_class`, `ss_context`, `raw_line`, `classification`, `rfam_scan_coverage`, `rfam_coverage_basis`, `endpoints_with_named`, `endpoints_with_loop`, `coverage`

## `outputs/04_comparison/step2_r8/overall/pair_overall.tsv`
Step 4. Per-pair counts of the five outcomes, assessable and unassessable.
*538 rows.* Columns: `pair_id`, `family`, `same_partner`, `different_partner`, `rfam_only`, `star3d_only`, `neither`, `total`, `same_partner_unassessable`, `different_partner_unassessable`, `rfam_only_unassessable`, `star3d_only_unassessable`

## `outputs/04_comparison/step2_r8/statistics/statistics.tsv`
Step 4. Every headline rate with pair-bootstrap 95% CI, family-balanced rate and leave-one-family-out range.
*11 rows.* Columns: `metric`, `n`, `N`, `rate`, `ci95_low`, `ci95_high`, `pairs`, `families`, `family_balanced_rate`, `lofo_min`, `lofo_max`, `lofo_largest_effect`

## `outputs/05_results/step2_r8/basepair_preservation_pairs.tsv.gz`
Step 5. Every MC-Annotate base pair of RNA A and what each alignment maps it to (same class / different / not comparable / no annotation / unassessable).
*24,185 rows.* Columns: `pair_id`, `family`, `u1`, `u2`, `source_class`, `mc_annotate`, `rfam_outcome`, `rfam_target_class`, `star3d_outcome`, `star3d_target_class`, `touches_disagreement`

## `outputs/05_results/step2_r8/candidate_corrections.tsv`
Step 5. Each A-side disagreement region with its distance and interaction evidence and the verdict (supports STAR3D / Rfam / conflicting / ...).
*1,155 rows.* Columns: `region_id`, `pair_id`, `family`, `length`, `start_unit`, `end_unit`, `ss_composition`, `median_d_rfam`, `median_d_star3d`, `distance_favours`, `interactions_kept_only_star3d`, `interactions_kept_only_rfam`, `noncanonical_kept_only_star3d`, `noncanonical_kept_only_rfam`, `interactions_favour`, `verdict`

## JSON summaries
- `outputs/04_comparison/step2_r8/overall/overall_metric.json`: headline metric, definitions and internal consistency checks.
- `outputs/04_comparison/step2_r8/statistics/statistics.json`: rates with intervals; hairpin vs outside test.
- `outputs/05_results/step2_r8/basepair_preservation.json`: preservation overall, at disagreements, and by interaction type.
- `outputs/05_results/step2_r8/candidate_corrections.json`: verdict counts.
- `docs/results.json`: every number quoted in the README and paper, generated from the files above.
