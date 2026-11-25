# Limitations

Stated plainly, with what was done about each.

1. **One structural aligner, one direction.** All 3D correspondences come from STAR3D v1.2, run once per pair with
   the smaller chain ID as query. Disagreements show where Rfam and STAR3D differ, not which is right. The
   interaction-based verdicts (`s53`) are the main evidence for that question.
2. **STAR3D chooses the superposition.** C1′ distances are measured in STAR3D's fitted frame, which favours STAR3D's
   correspondences. Distance-only verdicts are therefore reported separately as weaker evidence.
3. **Annotation can miss interactions.** Base pairs come from MC-Annotate. "Not kept" means no annotated pair of
   the same class, not proof of absence. FR3D, listed in the original config, was not run.
4. **Small and uneven cohort.** Only families with natural structures from at least two species can be tested:
   503 aligned pairs in 10 families, dominated by tRNA (350) and 5S rRNA (136). Pooled and family-balanced rates
   are both reported, with leave-one-family-out ranges.
5. **Open origin checks.** 39 natural-origin checks were open when Step 2 closed (policy R7) and are being re-queried
   (R8). Chains with unconfirmed origin were not used, which can only remove pairs, never add questionable ones.
6. **Rfam motif-scan coverage is family level.** A missing Rfam motif hit is not verified absence for a given row;
   RF00522's scan scope is undocumented.
7. **Candidate corrections are pairwise.** Revising a seed alignment requires checking a proposal against all rows
   of the family.
8. **The candidate-correction rule was defined after the aggregate preservation results were seen.** Its margin is
   a parameter, and the verdict counts should be read with that in mind.
