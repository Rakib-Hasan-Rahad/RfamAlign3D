# Where RNA sequence and structural alignments agree—and where they differ

*RfamAlign3D results. Every number is a query on the study database. Prepared by an AI research agent; not yet reviewed by a person.*

**Question.** We compared the nucleotide correspondences implied by two original Rfam 15.1 seed rows with those computed by original STAR3D v1.2 from the experimental structures, over the whole family interval, for verified natural RNAs from different species.

**Study scope.** The study screened 4,227 Rfam families and obtained usable structural comparisons in 10 families (538 usable alignments of 570 attempted; 32 pairs gave “No similar stacks”). Comparisons reuse RNAs, and tRNA and 5S rRNA dominate the dataset (552 of 570 pairs). Alignment differences identify candidates for investigation; they do not establish which method is correct.

**Main results.**
- STAR3D reproduces **39,059/46,857 (83.4%)** assessable Rfam consensus correspondences.
- **All assessable alignment differences, RNA A view: 8,996/48,399 (18.6%)** — different partner 4,676, Rfam only 3,234, STAR3D only 1,086. Excluded: 796 neither-mapped and 544 not assessable (a nucleotide or its counterpart lacks coordinates).
- Of different-partner cases, the alternatives are one nucleotide apart in **3,896/4,676 (83.3%)**.
- Correspondence sets (direction-independent; unit = residue correspondences): shared 39,403; Rfam only 7,910; STAR3D only 5,827; set difference 13,737/53,140 (25.9%).

The percentages use different denominators and are not complements.

| Family | Pairs | Overall assessable difference rate | Different-partner rate |
|---|---:|---:|---:|
| tRNA (RF00005) | 350 | 5,433/24,969 (21.8%) | 2444/21980 |
| 5S_rRNA (RF00001) | 171 | 2,822/19,735 (14.3%) | 1821/18734 |
| 5_8S_rRNA (RF00002) | 10 | 312/1,404 (22.2%) | 152/1244 |
| U4 (RF00015) | 1 pair studied | 5/64 (7.8%) | 0/59 |
| U6 (RF00026) | 1 pair studied | 10/93 (10.8%) | 3/86 |
| Bacteria_small_SRP (RF00169) | 1 pair studied | 33/33 (100.0%) | 32/32 |
| SSU_rRNA_bacteria (RF00177) | 1 pair studied | 106/1,508 (7.0%) | 78/1480 |
| PreQ1 (RF00522) | 1 pair studied | 17/33 (51.5%) | 15/31 |
| tRNA-Sec (RF01852) | 1 pair studied | 5/86 (5.8%) | 5/86 |
| ROOL (RF03087) | 1 pair studied | 253/474 (53.4%) | 126/347 |

**Where and how the alignments differ** (descriptive; different-partner observations; each in exactly one pattern): 3,896/4,676 different-partner observations are offsets of exactly one position along RNA B's ungapped sequence (insertion codes preserved). Constant-offset runs next to a seed gap: 3,003/4,676; pair-wide difference (≥50% of residues mapped by both differ): 785/4,676 in 26 pairs, not one constant offset (391/785 share their pair's most common offset). Seed gaps lie within two columns of 3695/4676 disagreements and of 27931/39403 agreements. Disagreement rate in consensus hairpin loops 1262/7761 vs 3414/36318 elsewhere.

**Motif evidence** (different-partner observations, any of three endpoints, counted once): Rfam computational motif span 1,217; BGSU Atlas group 15; both 32; BGSU structural loop only 3,021; consensus loop context only 54; no documented record 286; Rfam coverage unknown with no other record 51. Rfam scan coverage is established per family, not per row; a missing record is not a confirmed absence. Descriptive only.

**Structural evidence: base-pair preservation at disagreements** (step 6). For each base pair of RNA A that touches a disagreement and whose two nucleotides are mapped under both methods, does each method's counterpart pair carry the same Leontis-Westhof class (compared in mapped order)? Two annotators: MC-Annotate (the file STAR3D itself uses for scoring) and FR3D (independent). Rfam vs STAR3D:

| Interaction, family | MC-Annotate | FR3D |
|---|---|---|
| canonical A-U/G-C, 5S rRNA | 43.2% vs 49.2% (n=577) | 42.6% vs 46.7% (n=561) |
| canonical A-U/G-C, tRNA | 57.4% vs 48.6% (n=627) | 56.3% vs 44.2% (n=616) |
| canonical A-U/G-C, other families | 45.3% vs 54.7% (n=86) | 47.6% vs 53.6% (n=84) |
| canonical A-U/G-C, pooled | 50.2% vs 49.3% (n=1,290) | 49.6% vs 45.9% (n=1,261) |
| noncanonical, 5S rRNA | 11.9% vs 43.7% (n=293) | 10.8% vs 50.6% (n=259) |
| noncanonical, tRNA | 18.0% vs 8.8% (n=250) | 21.5% vs 13.1% (n=274) |
| noncanonical, other families | 2.2% vs 26.7% (n=45) | 11.1% vs 36.1% (n=36) |
| noncanonical, pooled | 13.8% vs 27.6% (n=588) | 16.0% vs 31.6% (n=569) |

The pooled noncanonical difference (+13.8 points) is not robust once shared sequences are resampled (95% interval -4.2 to 32.2); the effect is family specific: STAR3D's correspondences preserve noncanonical pairs far more often in 5S rRNA, and Rfam's do better in tRNA.

**Disagreement regions** (1,155): interaction evidence agreed between both annotators in favour of STAR3D in 218 regions (74 strong) and of Rfam in 208 (55 strong); 295 had no pair eligible under both methods and 271 showed no difference. C1' distances are measured in STAR3D's own superposition and never decide.

**Candidate correction (for curator review).** RF00001 seed rows X75601.1/91-203 (*M. tuberculosis*) and X55251.1/3-116 (*M. smegmatis*), helix IV: STAR3D places the same segment in the same alternative columns against most of the 18 partners. Moving 15 and 14 residues (sequences unchanged) places the six helix IV cWW pairs and the noncanonical pairs U76-A95 (tWH) and A77-A94 (tHS) of 7KGB chain B at columns where 11-16 of 17 other structures have the same FR3D pair, against 0-1 now. Conflicting evidence: fewer canonical pairs at SS_cons-paired columns (30 to 27), but SS_cons in this window departs from the 3D structures. Proposed Stockholm file: `outputs/06_audit/step2_r8/correction_case/`.

**Contrasting case.** Bovine tRNA 7O7Y AT vs yeast 8ASW X: STAR3D shifts the acceptor-stem 3' strand by one and loses four Watson-Crick pairs that the seed alignment keeps, although distances in STAR3D's frame favour STAR3D.

**Limitations.** STAR3D is one aligner run in one direction, and its loop score uses the MC-Annotate base pairs; FR3D is the independent check. One sequence enters many pairs and three families (tRNA, 5S rRNA, 5.8S rRNA) carry nearly all the interaction evidence, so family results are descriptive. Annotators can miss weak contacts, and structures differ in resolution and context. Chains without a whole-chain match in an independent record were not used; they are unresolved, not shown to be engineered. Methods, selection counts and provenance (execution environments, checksums, policy history) are in the website's Methods & limitations page and the downloadable tables.
