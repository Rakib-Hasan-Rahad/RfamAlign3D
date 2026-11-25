# Where RNA sequence and structural alignments agree—and where they differ

*SQARMA final results. Every number is a query on the study database. Prepared by an AI research agent; not yet reviewed by a person.*

**Question.** We compared the nucleotide correspondences implied by two original Rfam 15.1 seed rows with those computed by original STAR3D v1.2 from the experimental structures, over the whole family interval, for verified natural RNAs from different species.

**Study scope.** The study screened 4,227 Rfam families and obtained usable structural comparisons in 10 families (503 usable alignments of 535 attempted; 32 pairs gave “No similar stacks”). Comparisons reuse RNAs, and tRNA and 5S rRNA dominate the dataset (517 of 535 pairs). Alignment differences identify candidates for investigation; they do not establish which method is correct.

**Main results.**
- STAR3D reproduces **35,655/42,990 (82.9%)** assessable Rfam consensus correspondences.
- **All assessable alignment differences, RNA A view: 8,439/44,428 (19.0%)** — different partner 4,281, Rfam only 3,164, STAR3D only 994. Excluded: 771 neither-mapped and 536 not assessable (a nucleotide or its counterpart lacks coordinates).
- Of different-partner cases, the alternatives are one nucleotide apart in **3,528/4,281 (82.4%)**.
- Correspondence sets (direction-independent; unit = residue correspondences): shared 35,989; Rfam only 7,445; STAR3D only 5,339; set difference 12,784/48,773 (26.2%).

The percentages use different denominators and are not complements.

| Family | Pairs | Overall assessable difference rate | Different-partner rate |
|---|---:|---:|---:|
| tRNA (RF00005) | 350 | 5,433/24,969 (21.8%) | 2444/21980 |
| 5S_rRNA (RF00001) | 136 | 2,265/15,764 (14.4%) | 1426/14925 |
| 5_8S_rRNA (RF00002) | 10 | 312/1,404 (22.2%) | 152/1244 |
| U4 (RF00015) | 1 pair studied | 5/64 (7.8%) | 0/59 |
| U6 (RF00026) | 1 pair studied | 10/93 (10.8%) | 3/86 |
| Bacteria_small_SRP (RF00169) | 1 pair studied | 33/33 (100.0%) | 32/32 |
| SSU_rRNA_bacteria (RF00177) | 1 pair studied | 106/1,508 (7.0%) | 78/1480 |
| PreQ1 (RF00522) | 1 pair studied | 17/33 (51.5%) | 15/31 |
| tRNA-Sec (RF01852) | 1 pair studied | 5/86 (5.8%) | 5/86 |
| ROOL (RF03087) | 1 pair studied | 253/474 (53.4%) | 126/347 |

**Where and how the alignments differ** (descriptive; different-partner observations; each in exactly one pattern): 3,528/4,281 different-partner observations are offsets of exactly one position along RNA B's ungapped sequence (insertion codes preserved). Constant-offset runs next to a seed gap: 2,663/4,281; pair-wide difference (≥50% of residues mapped by both differ): 785/4,281 in 26 pairs, not one constant offset (391/785 share their pair's most common offset). Seed gaps lie within two columns of 3327/4281 disagreements and of 25193/35989 agreements. Disagreement rate in consensus hairpin loops 1226/7441 vs 3055/32829 elsewhere.

**Motif evidence** (different-partner observations, any of three endpoints, counted once): Rfam computational motif span 1,162; BGSU Atlas group 15; both 32; BGSU structural loop only 2,760; consensus loop context only 21; no documented record 286; Rfam coverage unknown with no other record 5. Rfam scan coverage is established per family, not per row; a missing record is not a confirmed absence. Descriptive only.

**Example.** tRNA-Sec 3RG5 A U47C corresponds to 3W3S B G47B in the seed alignment (column 66) and to 3W3S B U47C in STAR3D (`A:47C<->B:47C`); both 3W3S residues belong to BGSU loop IL_3W3S_003 (Motif Atlas 4.16 group IL_10892.5) — a correspondence difference within one motif that does not establish which partner is correct.

**Limitations.** Constructs without sufficient independent natural-sequence evidence were not included in the verified comparison cohort; this includes chains whose sequence search did not complete, which are unresolved rather than shown to be engineered. Methods, selection counts and provenance (execution environments, checksums, policy history) are in the website's Methods & limitations page and the downloadable tables.
