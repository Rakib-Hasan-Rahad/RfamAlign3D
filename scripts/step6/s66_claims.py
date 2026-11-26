"""Step 6.6: claim_evidence.tsv - each headline claim, its count, denominator, definition, source and status."""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import csv, json, sys
ROOT, A = _ARGS[0] if len(_ARGS) > 0 else _REPO, _ARGS[1] if len(_ARGS) > 1 else _os.path.join(_REPO, "outputs", "06_audit", "step2_r8")
S = json.load(open(f"{A}/interaction_summary.json"))["results"]
FL = json.load(open(f"{A}/interaction_summary.json"))["family_level"]
O = json.load(open(f"{A}/origin_summary.json")); RS = json.load(open(f"{A}/region_summary.json"))
om = json.load(open(f"{ROOT}/outputs/04_comparison/step2_r8/overall/overall_metric.json"))["overall"]
st = json.load(open(f"{ROOT}/outputs/04_comparison/step2_r8/statistics/statistics.json"))
old = json.load(open(f"{ROOT}/outputs/05_results/step2_r8/candidate_corrections.json"))
pc = json.load(open(f"{A}/correction_case/patch_check.json"))["core"]
def k(a, d, s, t): return S[f"{a}|{d}|{s}|{t}"]
def rate(o, m): return f"{o[f'{m}_class']}/{o['common_eligible']} ({o[f'{m}_class_rate']*100:.1f}%)"
rows = []
def add(claim, value, denom, definition, source, status, note=""):
    rows.append(dict(claim=claim, value=value, denominator=denom, definition=definition, source=source, status=status, note=note))
add("Pairs attempted / aligned / no alignment / technical failures", "570 / 538 / 32 / 0", "570 jobs",
    "STAR3D v1.2 one forward run per pair; 'No similar stacks' in stdout = no alignment", "outputs/03_star3d/step2_r8/returned/ledger.json; runs/*/result.json", "CONFIRMED",
    "32 = 31 RF00005 + 1 RF00023; exit code 0 in all 32")
add("Cohort composition", "11 families; 138 chains; 126 PDB entries; 72 seed rows (distinct sequences); 37 species", "570 pairs",
    "from step3_jobs.tsv", "outputs/03_star3d/step2_r8/step3_jobs.tsv", "CONFIRMED")
add("535 earlier runs reused unchanged", "535 identical job fields; 570/570 result identity = job", "570",
    "direction, a/b input sha256, STAR3D jar sha256 e03ff7bc..., parameters, command", "audit recheck of result.json vs step3_jobs (final and r8)", "CONFIRMED",
    "505 reused runs native Linux, 30 in an x86-64 container on macOS, 35 new native; same jar and Java 8")
add("Origin checks resolved", f"{O['keys_in_current_queue']} current checks: " + ", ".join(f"{k_} {v}" for k_, v in O["current_outcomes"].items()),
    f"{O['keys_total']} keys ever cached", "search completed for every current key; 'resolved' = search finished, NOT natural origin established",
    "origin_check_reconciliation.tsv", "CORRECTED WORDING", "Only 103 keys have exact independent support; 85 not found and 36 flank-differs remain of unresolved origin and are not used. 2 cached keys left the queue (no partner / superseded)")
add("Cohort chains are natural", "138/138 verified_natural", "138 chains",
    "exact full-chain match in an independent non-PDB record of the stated species", "origin_summary.json", "CONFIRMED with caveat",
    f"evidence: genome assembly {O['cohort_chain_evidence']['genome_assembly']}, record >= chain+50 nt {O['cohort_chain_evidence']['long_record(>=chain+50nt)']}, only a record older than the structure {O['cohort_chain_evidence']['older_record_only']} (mostly 1980s direct RNA-sequencing entries of tRNA/5S)")
add("R8 megablast searches", "5 species; 2 whole-chain found (C. albicans X16634.1; T. thermophila J01893.1), 3 not found", "5 searches",
    "core_nt, txid[Organism:exp], megablast word 28, E=10, hitlist 50, no filter; hits confirmed by exact string search of efetched record", "cache/blast/*.megablast.xml; SCIENTIFIC_AUDIT.md table", "CORRECTED WORDING",
    "replace 'megablast finds the same hits' with: every 100%-identical full-length hit of a query >=28 nt contains a 28-nt seed; 3 of 5 hit lists were saturated at 50")
add("Overall difference rate", f"{om['different_partner']+om['rfam_only']+om['star3d_only']}/{om['same_partner']+om['different_partner']+om['rfam_only']+om['star3d_only']} (18.6%)", "assessable RNA A nucleotides",
    "(D+R+S)/(A+D+R+S), A->B view", "overall_metric.json; verify.py recount", "CONFIRMED", "pair-bootstrap CI 16.8-20.5 ignores shared molecules; describe as descriptive")
for t, lbl in (("canonical_WC", "Canonical A-U/G-C kept at disagreements"), ("GU_wobble", "G-U wobble kept at disagreements"),
               ("noncanonical_non_cWW", "Non-cWW noncanonical kept at disagreements (old '8.9 vs 23.5')"), ("noncanonical", "All noncanonical kept at disagreements")):
    o = k("MC-Annotate", "A->B", "disagreement", t); f = k("FR3D", "A->B", "disagreement", t)
    add(lbl, f"MC-Annotate Rfam {rate(o,'rfam')} vs STAR3D {rate(o,'star3d')}; FR3D Rfam {rate(f,'rfam')} vs STAR3D {rate(f,'star3d')}",
        f"common eligible source pairs touching a disagreement; {o['rna_pairs_common']} RNA pairs, {o['families_common']} families",
        "same oriented LW class between mapped counterparts; A->B; source identities define type", "interaction_summary.json; interaction_preservation_by_family.tsv",
        "CONFIRMED, REINTERPRETED" if t != "canonical_WC" else "CONFIRMED, WORDING",
        f"diff MC {o['diff_class']*100:+.1f} pts, pair CI {[round(x*100,1) for x in o['ci95_diff_class_pair_bootstrap']]}, sequence-vertex CI {[round(x*100,1) for x in o['ci95_diff_class_sequence_vertex_bootstrap']]}")
fl = FL["MC-Annotate|disagreement|noncanonical"]
add("Noncanonical effect by family", json.dumps(fl["per_family_diff_class"]), "families with >=10 common pairs",
    "STAR3D - Rfam same-class rate", "interaction_summary.json family_level", "CORRECTED", "positive in 5S rRNA, NEGATIVE in tRNA (Rfam higher); pooled effect driven by 5S rRNA")
add("14.6-point CI 8.7-20.2", "reproduced exactly (pair bootstrap)", "481 common non-cWW pairs",
    "2,000 resamples of whole RNA pairs, seed 20261125", "s52 / audit", "NOT DEFENSIBLE AS STATED",
    "under a sequence-vertex bootstrap that respects shared molecules: -7.8 to 36.7; report descriptively with family results")
o = k("MC-Annotate", "A->B", "whole", "noncanonical")
add("Whole compared interval, noncanonical", f"Rfam {rate(o,'rfam')} vs STAR3D {rate(o,'star3d')}", "all common eligible source pairs",
    "same class", "interaction_summary.json", "NEW", "effect ~2 points; disagreement-conditioned result must not be generalised")
add("Candidate corrections 221 vs 111", f"old: {old['verdicts']}", "1,155 A-side regions", "s53 verdict with STAR3D-frame distance veto", "candidate_corrections.json", "WITHDRAWN",
    "s53 counted unmapped Rfam endpoints as STAR3D gains and vetoed with a STAR3D-frame distance")
add("Region support (audit)", json.dumps(RS["interaction_tiers"]), "1,155 regions", "interaction evidence on common eligible set; both annotators; distance never decides",
    "region_support_audit.tsv", "REPLACES 221/111", "STAR3D strong+moderate 218 vs Rfam 208; strong 74 vs 55")
add("Proposed correction RF00001 mycobacterial helix IV", f"2 rows changed, {pc['residues_moved']} residues moved; FR3D pairs at proposed columns supported by 11-16 of {pc['reference_rows']} other rows (0-1 before)",
    "8 FR3D pairs per M. tuberculosis chain", "row-level, STAR3D majority over 18 partners; tested with FR3D column consistency", "correction_case/", "CANDIDATE FOR CURATOR REVIEW",
    "conflict: canonical/G-U pairs at SS_cons columns 30 -> 27; SS_cons itself disagrees with 3D in this window")
add("Fig.3 example C1' RMS 3.94 -> 1.21 A", "pair 7KGB_B vs 9QT5_2", "1 pair", "C1' RMS in STAR3D's fitted frame", "pair.json geom", "DOWNWEIGHTED", "favours STAR3D by construction; not used as evidence")
add("Independence of interaction evaluation", "MC-Annotate annotations are STAR3D's own preprocessing; STAR3D's loop score rewards matching annotated W/H/S pairs (LoopAlign.multi_pairing_align)",
    "", "", "work/star3d_src/STAR3D_source/src/LoopAlign.java, Lib.java", "CORRECTED", "FR3D (pinned 288f98c) added as independent annotator; direction replicated")
add("Rfam uses no 3D information", "false: Rfam 15.0 improved 65 families with 3D structures; RF00001 seed in 15.1 has no 3D-curated rows", "", "", "Rfam 15 paper; Rfam docs", "CORRECTED")
add("AlphaFold/evolutionary benefit", "not tested", "", "", "", "REWORDED as potential relevance")
with open(f"{A}/claim_evidence.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "claims")
