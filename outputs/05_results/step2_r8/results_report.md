# RfamAlign3D — interim results (cohort `step2_r8`)

**INTERIM.** Generated from `rfamalign3d.sqlite`. Frozen interim cohort step2_r8; Step 2 checks were still pending for other pairs. Prepared by an AI agent; not yet reviewed by a person.

**1 · Question and cohort.** Do two original Rfam 15.1 seed rows and one forward run of original STAR3D v1.2 assign the same residue partners across the whole family interval, for verified natural RNAs of different species?

**2 · Coverage.** 4,227 families screened; 11 families have selected pairs in this cohort. 570 pairs attempted: **538 valid alignments in 10 families**, 32 “No similar stacks” (NA), 0 technical failures. Pairs reuse RNAs; tRNA and 5S rRNA supply 552 of 570 pairs.

**3 · Alignment result.** 39,403 exact residue pairs are shared (Rfam 47,313, STAR3D 45,230 assessable pairs). *Exclusion rule:* a residue pair is assessable only if both endpoints are observed in the staged coordinates (defined from Step 2, not from either method). Of 47,857 raw Rfam pairs, 544 are excluded (A endpoint unobserved: 123; A endpoint unobserved;B endpoint unobserved: 53; B endpoint unobserved: 368); 0 STAR3D pairs are excluded. Consensus and insert strata partition the 47,313 assessable Rfam pairs. Rfam consensus-column pairs reproduced by STAR3D: **39,059/46,857** (83.4%); insert-column pairs: 344/456. A-side partner disagreements: **4,676/44,079 (10.6%)** jointly assessable observations.

| Family | Pairs | Valid | No stacks | A-side different / assessable |
|---|---:|---:|---:|---:|
| RF00005 tRNA | 381 | 350 | 31 | 2,444/21,980 (11.1%) |
| RF00001 5S_rRNA | 171 | 171 | 0 | 1,821/18,734 (9.7%) |
| RF00002 5_8S_rRNA | 10 | 10 | 0 | 152/1,244 (12.2%) |
| RF03087 ROOL | 1 | 1 | 0 | 126/347 (36.3%) |
| RF01852 tRNA-Sec | 1 | 1 | 0 | 5/86 (5.8%) |
| RF00522 PreQ1 | 1 | 1 | 0 | 15/31 (48.4%) |
| RF00177 SSU_rRNA_bacteria | 1 | 1 | 0 | 78/1,480 (5.3%) |
| RF00169 Bacteria_small_SRP | 1 | 1 | 0 | 32/32 (100.0%) |
| RF00026 U6 | 1 | 1 | 0 | 3/86 (3.5%) |
| RF00023 tmRNA | 1 | 0 | 1 | NA |
| RF00015 U4 | 1 | 1 | 0 | 0/59 (0.0%) |

**Where and how the alignments differ** (descriptive; each disagreement in exactly one pattern): 3,896/4,676 disagreements are offsets of exactly one position along RNA B's ungapped sequence (insertion codes preserved). Constant-offset runs next to a seed gap: 3,003/4,676; pair-wide difference (≥50% of residues mapped by both differ): 785/4,676 in 26 pairs, not one constant offset (391/785 share their pair's most common offset). Seed gaps lie within two columns of 3695/4676 disagreements and of 27931/39403 agreements. Disagreement rate in consensus hairpin loops 1262/7761 vs 3414/36318 elsewhere.

**4 · Motif result** (A-side disagreements, any of 3 endpoints, counted once): Rfam computational motif only 1,217; BGSU Atlas group only 15; both 32; BGSU loop only 3,021; consensus loop context only 54; no documented record (BGSU loop lookup complete; Rfam absence not verified per row) 286; Rfam coverage unknown with no other record 51. Baseline, query endpoint only: Atlas group 11/128 disagreements vs 306/1465 agreements (Atlas covers few of these structures); Rfam motif 774/4661 vs 8345/39387. Descriptive, not an enrichment test. *Rfam coverage:* `motif_family_stats` reconciles with all seed rows of the covered families, so coverage is established **only at family level** for every family: no row-level scan record exists, and an absent Rfam hit is not a verified negative for any individual seed row; RF00522 is not listed, so its Rfam coverage is **unknown** for all 15 of its disagreement observations, including those with documented BGSU loop evidence. Rfam hits are position spans, shown with their source records.

**5 · Inspectable example.** tRNA-Sec 3RG5 A U47C ↔ 3W3S B G47B (Rfam column 66) vs U47C (STAR3D `A:47C<->B:47C`); both targets are in BGSU loop IL_3W3S_003, Atlas group IL_10892.5 — a disagreement within the same motif that does not establish which partner is correct. Site page: `site/evidence/trna_sec_example.html`.
