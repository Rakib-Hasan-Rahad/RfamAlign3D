# Limitations

Stated plainly, with what was done about each.

1. **One structural aligner, one direction.** All 3D correspondences come from STAR3D v1.2, run once per pair with
   the smaller chain ID as query. Disagreements show where Rfam and STAR3D differ, not which is right; neither is
   treated as ground truth. STAR3D can misregister short helices (e.g. tRNA acceptor stems).
2. **STAR3D chooses the superposition.** C1′ distances are measured in STAR3D's fitted frame, which favours STAR3D's
   correspondences. Distance-only verdicts are therefore reported separately as weaker evidence.
3. **Annotation can miss interactions, and MC-Annotate is not independent of STAR3D.** STAR3D's scoring uses the
   MC-Annotate pairs; FR3D (step 6) is the independent check, but any 3D aligner favours geometrically matched pairs.
   "Not preserved" means no annotated pair of the same class, not proof of absence.
4. **Small and uneven cohort.** Only families with natural structures from at least two species can be tested:
   538 aligned pairs in 10 families, dominated by tRNA (350) and 5S rRNA (171). Pooled and family-balanced rates
   are both reported, with leave-one-family-out ranges.
5. **Origin searches depend on NCBI records.** All 234 checks are resolved (R8). A chain whose sequence is in no
   searched record is not used; this can only remove pairs, never add questionable ones. Five searches used megablast,
   which can miss short interval queries (<28 nt) that only label flank differences.
6. **Rfam motif-scan coverage is family level.** A missing Rfam motif hit is not verified absence for a given row;
   RF00522's scan scope is undocumented.
7. **Candidate corrections are pairwise.** Revising a seed alignment requires checking a proposal against all rows
   of the family.
8. **Analysis choices made after seeing results.** The step 6 region tiers and the correction-case selection were
   defined after the aggregate results were known; both are reported with their rules and alternatives.
9. **Dependence and family imbalance.** One sequence enters many pairs and three families carry nearly all the
   information; results are family specific (5S rRNA vs tRNA) and family-level inference is weak.
10. **Conformational and contextual variation.** Structures differ in resolution, ligand and ribosome context, which
   can change annotated interactions independently of alignment.
11. **Origin evidence strength varies.** 28 cohort chains rest on early RNA-sequencing records rather than genome
   assemblies.
