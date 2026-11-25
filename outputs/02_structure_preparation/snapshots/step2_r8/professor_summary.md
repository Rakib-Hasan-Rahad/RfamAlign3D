# SQARMA Step 2 — structure pairs prepared for original STAR3D (FINAL)

*Automated pipeline, 2026-11-25; no human review has occurred. **STAR3D has not yet been run in this phase.***

**Scope.** All 4,227 Rfam 15.1 families were searched (explicit seed tags, Step 1 leads, and both Rfam 15.1 family-to-PDB mapping files). 105 families have a PDB chain matching an entire seed row exactly; 4,058 had no candidate in the searched sources and 64 had candidates without an exact seed match.

**Result.** 11 families yield eligible pairs; **570 jobs are ready** (RF00005 tRNA (381), RF00001 5S_rRNA (171), RF00002 5_8S_rRNA (10), RF00015 U4 (1), RF00023 tmRNA (1), RF00026 U6 (1), RF00169 Bacteria_small_SRP (1), RF00177 SSU_rRNA_bacteria (1), RF00522 PreQ1 (1), RF01852 tRNA-Sec (1), RF03087 ROOL (1)). Each job pairs two verified-natural RNAs from different species, different PDB entries and different seed-row sequences; inputs: 462 verified chains fully modelled, 381 with coverage flags (recorded per residue).

**Why others did not qualify** (chain decisions; one chain can block many pairs):

- GENOME_NOT_FOUND: 381
- SPECIES_UNRESOLVED: 114
- TRNA_3PRIME_CCA_POSSIBLY_NOT_ENCODED: 84
- NON_NATIVE_FLANK_AT_LOCUS: 24
- UNKNOWN_RESIDUES_IN_CHAIN: 23
- DOCUMENTED_CONSTRUCT_CHANGE: 10

Family outcomes: no_candidate_in_searched_sources 4058; exact_match_single_sequence_no_pair_possible 64; candidate_found_no_exact_seed_match 64; pair_possible_no_eligible_pair 24; eligible_pairs_selected 11; no_different_species_pair_possible 6.

**Traceable example.**

> Pair RF00522__rf00522_6vuh_Au_1_33__vs__rf00522_7rex_Au_1_34 (RF00522 PreQ1): query 6VUH_A [Caldanaerobacter subterraneus (911092)] seed row URS000080E32E_119072/1-33; target 7REX_A [Carnobacterium antarcticum (2126436)] seed row URS00023119CB_2126436/1-34.
> Rfam column 6 (consensus): query C = A:1 (staged rf00522_6vuh_Au_1_33 A:1), target U = A:1 (staged rf00522_7rex_Au_1_34 A:1).
>   crosswalk RF00522__6VUH_A__1-33: col 6 -> ungapped 1 -> A/1 -> auth A:1 (C) -> inputs_star3d/rf00522_6vuh_Au_1_33.pdb A:1
>   crosswalk RF00522__7REX_A__1-34: col 6 -> ungapped 1 -> A/1 -> auth A:1 (U) -> inputs_star3d/rf00522_7rex_Au_1_34.pdb A:1

**Limits.** Species and natural origin are automated decisions from deposited annotations and NCBI records; tRNA pairs dominate the job list; pairs within a family share structures (flagged); a pair being eligible says nothing about whether the methods will agree.

Validation: 10/10 gates pass (validation_report.md).
