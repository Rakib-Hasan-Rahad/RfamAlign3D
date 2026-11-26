"""Step 6.3: interaction preservation recomputed from source annotations, with explicit definitions and
dependence-aware uncertainty.

  <project>/.venv/bin/python scripts/step6/s63_interactions_regions.py [<repo_root> <out_dir>] [--reps 2000] [--seed 20261125]

Inputs (read only): step2_r8 pair records (pair.json: residues, both alignments' counterparts, MC-Annotate base pairs
as produced by STAR3D v1.2 preprocessing), step3_jobs.tsv (seed rows), and the FR3D annotation written by
a2_fr3d_annotate.py (independent of STAR3D's inputs).

For every source base pair (u,v) annotated in the SOURCE RNA, each alignment maps u->t, v->w in the TARGET RNA.
Endpoint states (kept separate):  src_no_coords | unmapped | target_no_coords | mapped
Eligible under a method: both endpoints mapped to target residues with coordinates, in the target chain instance.
Common eligible set: eligible under BOTH methods (paired comparison); mapping coverage is reported separately.
Preservation definitions (never interchanged):
  A any_pair       the target annotator lists ANY base pair between t and w
  B same_class     t-w has the source pair's Leontis-Westhof class, compared in mapped endpoint order
                   (target class reversed when stored as w-t: cWH stored reversed = cHW)
  C class+identity same_class AND identity(t)=identity(u) and identity(w)=identity(v)
Interaction types (source pair; identities are parent bases):
  canonical_WC  cWW and {AU,UA,GC,CG};  GU_wobble  cWW and {GU,UG};
  noncanonical  every other LW class, plus cWW between any other bases (reported also as noncanonical_cWW_other)
Scopes: whole compared interval (all source pairs) and disagreement (either endpoint has category different_partner,
rfam_only or star3d_only); 'dp_only' reproduces the earlier definition (different_partner only).
Direction: A->B (primary; query->target, as STAR3D was run) and B->A (sensitivity). The two directions share the same
pairs and are never pooled as independent observations.
Uncertainty: 95% percentile intervals, 2,000 replicates, seed 20261125, for the paired difference STAR3D - Rfam in the
pooled ratio over the common eligible set, under (i) bootstrap of whole RNA pairs and (ii) a vertex bootstrap over
unique sequences (seed rows): sequences are resampled and a pair receives weight m_a * m_b, which carries the
dependence created by one sequence occurring in many pairs. Family-level results are descriptive.
Outputs (in audit_dir): interactions_long.tsv.gz, interaction_preservation_by_family.tsv, interaction_summary.json,
region_support_audit.tsv, region_summary.json
"""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import argparse
import collections
import csv
import glob
import gzip
import json
import os
import statistics

import numpy as np

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from rules import DIS, btype, endpoint_state, flip, key, preservation, region_tier  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=_REPO)
ap.add_argument("out", nargs="?", default=_os.path.join(_REPO, "outputs", "06_audit", "step2_r8"))
ap.add_argument("--reps", type=int, default=2000); ap.add_argument("--seed", type=int, default=20261125)
a = ap.parse_args()
ROOT, OUT = a.root, a.out
C5 = f"{ROOT}/outputs/05_results/step2_r8"
C4 = f"{ROOT}/outputs/04_comparison/step2_r8"


# ------------------------------------------------------------------ FR3D annotations, keyed by residue identity
fr3d = collections.defaultdict(dict)      # pdb -> {(k1,k2): class oriented k1->k2}
for f in glob.glob(f"{OUT}/fr3d/*_basepair.txt"):
    pdb = os.path.basename(f).split("_")[0]
    for line in open(f):
        x = line.rstrip("\n").split("\t")
        if len(x) < 3 or x[1].startswith("n"):
            continue
        c = x[1]
        c = c[0] + c[1:3].upper() if len(c) >= 3 else c
        fr3d[pdb][(key(x[0]), key(x[2]))] = c

jobs = {j["pair_id"]: j for j in csv.DictReader(open(f"{ROOT}/outputs/03_star3d/step2_r8/step3_jobs.tsv"), delimiter="\t")}
prov = {r["pair_id"]: r["source"] for r in csv.DictReader(open(f"{ROOT}/outputs/03_star3d/step2_r8/returned/run_provenance.tsv"), delimiter="\t")}


def src_pairs_mca(p, side):
    out = []
    for u1, u2, lw, raw in p["inter"][side]:
        out.append((u1, u2, lw or ""))
    return out


def src_pairs_fr3d(p, side, res):
    pdb = p[side.lower()]["chain"].split("_")[0]
    units = {key(u): u for u in res}
    seen, out = set(), []
    for (k1, k2), c in fr3d.get(pdb, {}).items():
        if k1 in units and k2 in units and (k2, k1) not in seen:
            seen.add((k1, k2)); out.append((units[k1], units[k2], c))
    return out


def tgt_lookup_mca(p, side):
    d = {}
    for x0, x1, c, raw in p["inter"][side]:
        d[(x0, x1)] = c or ""
        d[(x1, x0)] = flip(c) if c else ""
    return d


def tgt_lookup_fr3d(p, side, res):
    pdb = p[side.lower()]["chain"].split("_")[0]
    units = {key(u): u for u in res}
    return {(units[k1], units[k2]): c for (k1, k2), c in fr3d.get(pdb, {}).items() if k1 in units and k2 in units}


rows = []
pair_meta = {}
for f in sorted(glob.glob(f"{C5}/app/data/pairs/*/pair.json")):
    p = json.load(open(f))
    if not p.get("res") or not p.get("inter"):
        continue
    pid = p["id"]
    j = jobs[pid]
    pair_meta[pid] = dict(family=p["family"], seed_a=j["a_seed_row"], seed_b=j["b_seed_row"],
                          r8_new=prov[pid].startswith("new"))
    R = {s: {x["u"]: x for x in p["res"][s]} for s in "AB"}
    for src, tgt, direction in (("A", "B", "A->B"), ("B", "A", "B->A")):
        tchain = tuple(p[tgt.lower()]["chain"].split("_", 1))
        for ann in ("MC-Annotate", "FR3D"):
            if ann == "MC-Annotate":
                sp_, tl = src_pairs_mca(p, src), tgt_lookup_mca(p, tgt)
            else:
                sp_, tl = src_pairs_fr3d(p, src, R[src]), tgt_lookup_fr3d(p, tgt, R[tgt])
            for u1, u2, c in sp_:
                r1, r2 = R[src].get(u1), R[src].get(u2)
                if not (r1 and r2):
                    continue
                rec = dict(pair_id=pid, family=p["family"], direction=direction, annotator=ann, u1=u1, u2=u2,
                           source_class=c, n1=r1["n"], n2=r2["n"], type=btype(c, r1["n"], r2["n"]),
                           cat1=r1["c"], cat2=r2["c"],
                           disagreement="yes" if (r1["c"] in DIS or r2["c"] in DIS) else "no",
                           dp_only="yes" if "different_partner" in (r1["c"], r2["c"]) else "no")
                for m, via in (("rfam", "rp"), ("star3d", "sp")):
                    st, t1, t2 = endpoint_state(r1, r2, via, R[tgt], tchain)
                    rec[f"{m}_state"] = st
                    anyp = cls = ident = ""
                    tc = ""
                    if st == "mapped":
                        tc = tl.get((t1, t2))
                        anyp, cls, ident = preservation(c, tc, R[tgt][t1]["n"] == r1["n"] and R[tgt][t2]["n"] == r2["n"])
                    rec[f"{m}_target"] = f"{t1}~{t2}" if st == "mapped" else ""
                    rec[f"{m}_target_class"] = tc or ""
                    rec[f"{m}_any"], rec[f"{m}_class"], rec[f"{m}_class_ident"] = anyp, cls, ident
                rows.append(rec)

with gzip.open(f"{OUT}/interactions_long.tsv.gz", "wt", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ aggregation
pids = sorted(pair_meta)
pidx = {p: i for i, p in enumerate(pids)}
seqs = sorted({pair_meta[p][k] for p in pids for k in ("seed_a", "seed_b")})
sidx = {s: i for i, s in enumerate(seqs)}
rng = np.random.default_rng(a.seed)
BOOT_P = rng.integers(0, len(pids), size=(a.reps, len(pids)))
PW = np.zeros((a.reps, len(pids)))
for b in range(a.reps):
    PW[b] = np.bincount(BOOT_P[b], minlength=len(pids))
VM = np.zeros((a.reps, len(seqs)))
for b in range(a.reps):
    VM[b] = np.bincount(rng.integers(0, len(seqs), size=len(seqs)), minlength=len(seqs))
ia = np.array([sidx[pair_meta[p]["seed_a"]] for p in pids]); ib = np.array([sidx[pair_meta[p]["seed_b"]] for p in pids])
VW = VM[:, ia] * VM[:, ib]
pc = lambda v: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]  # noqa: E731

TYPES = ("all_LW", "canonical_WC", "GU_wobble", "noncanonical", "noncanonical_cWW_other", "noncanonical_non_cWW")


def in_type(r, t):
    if t == "all_LW":
        return r["type"] != "no_LW_class"
    if t == "noncanonical_cWW_other":
        return r["type"] == "noncanonical" and r["source_class"] == "cWW"
    if t == "noncanonical_non_cWW":
        return r["type"] == "noncanonical" and r["source_class"] != "cWW"
    return r["type"] == t


def scope_ok(r, s):
    return s == "whole" or (s == "disagreement" and r["disagreement"] == "yes") or (s == "dp_only" and r["dp_only"] == "yes")


def summarise(sel, with_ci=True):
    """sel: rows. Paired comparison on the common eligible set, plus coverage and intention-to-treat rates."""
    n = len(sel)
    elig = [r for r in sel if r["rfam_state"] == "mapped" and r["star3d_state"] == "mapped"]
    o = dict(source_interactions=n,
             unique_source_interactions=len({(r["direction"], r["u1"], r["u2"]) for r in sel}),
             rna_pairs=len({r["pair_id"] for r in sel}), families=len({r["family"] for r in sel}),
             eligible_rfam=sum(r["rfam_state"] == "mapped" for r in sel),
             eligible_star3d=sum(r["star3d_state"] == "mapped" for r in sel), common_eligible=len(elig),
             rna_pairs_common=len({r["pair_id"] for r in elig}), families_common=len({r["family"] for r in elig}))
    for m in ("rfam", "star3d"):
        o[f"{m}_states"] = dict(collections.Counter(r[f"{m}_state"] for r in sel))
        for d in ("any", "class", "class_ident"):
            k = sum(r[f"{m}_{d}"] == "yes" for r in elig)
            o[f"{m}_{d}"] = k
            o[f"{m}_{d}_rate"] = round(k / len(elig), 4) if elig else None
            o[f"{m}_{d}_itt_rate"] = round(sum(r[f"{m}_{d}"] == "yes" for r in sel) / n, 4) if n else None
    for d in ("any", "class", "class_ident"):
        o[f"diff_{d}"] = round((o[f"star3d_{d}"] - o[f"rfam_{d}"]) / len(elig), 4) if elig else None
    if with_ci and elig:
        N = np.zeros(len(pids)); Rr = np.zeros(len(pids)); Ss = np.zeros(len(pids))
        Ra = np.zeros(len(pids)); Sa = np.zeros(len(pids))
        for r in elig:
            i = pidx[r["pair_id"]]; N[i] += 1
            Rr[i] += r["rfam_class"] == "yes"; Ss[i] += r["star3d_class"] == "yes"
            Ra[i] += r["rfam_any"] == "yes"; Sa[i] += r["star3d_any"] == "yes"
        for nm, W in (("pair_bootstrap", PW), ("sequence_vertex_bootstrap", VW)):
            n_ = W @ N
            ok = n_ > 0
            o[f"ci95_diff_class_{nm}"] = pc(((W @ Ss) - (W @ Rr))[ok] / n_[ok])
            o[f"ci95_diff_any_{nm}"] = pc(((W @ Sa) - (W @ Ra))[ok] / n_[ok])
            o[f"ci95_rfam_class_{nm}"] = pc((W @ Rr)[ok] / n_[ok])
            o[f"ci95_star3d_class_{nm}"] = pc((W @ Ss)[ok] / n_[ok])
    return o


by_fam_rows = []
summary = dict(definitions=__doc__.split("Outputs")[0].strip(), cohort="step2_r8",
               bootstrap=dict(replicates=a.reps, seed=a.seed, level=0.95, units=["RNA pair", "sequence (seed row) vertex"]))
G = collections.defaultdict(list)
for r in rows:
    G[(r["annotator"], r["direction"])].append(r)
main = {}
for (ann, d), rr in sorted(G.items()):
    for s in ("whole", "disagreement", "dp_only"):
        for t in TYPES:
            sel = [r for r in rr if scope_ok(r, s) and in_type(r, t)]
            o = summarise(sel)
            main[f"{ann}|{d}|{s}|{t}"] = o
            fams = sorted({r["family"] for r in sel})
            for fam in fams + ["ALL", "ALL_excl_tRNA_5S", "ALL_excl_R8_new"]:
                if fam == "ALL":
                    ss = sel
                elif fam == "ALL_excl_tRNA_5S":
                    ss = [r for r in sel if r["family"] not in ("RF00005", "RF00001")]
                elif fam == "ALL_excl_R8_new":
                    ss = [r for r in sel if not pair_meta[r["pair_id"]]["r8_new"]]
                else:
                    ss = [r for r in sel if r["family"] == fam]
                of = summarise(ss, with_ci=fam.startswith("ALL") and s != "dp_only" and t in ("canonical_WC", "noncanonical", "all_LW", "GU_wobble"))
                by_fam_rows.append(dict(annotator=ann, direction=d, scope=s, type=t, family=fam,
                                        **{k: (json.dumps(v) if isinstance(v, (dict, list)) else v) for k, v in of.items()}))
summary["results"] = main

# family-balanced and leave-one-family-out (descriptive) for the key contrasts
fb = {}
for ann in ("MC-Annotate", "FR3D"):
    for s in ("whole", "disagreement"):
        for t in ("canonical_WC", "GU_wobble", "noncanonical"):
            fr = [x for x in by_fam_rows if x["annotator"] == ann and x["direction"] == "A->B" and x["scope"] == s
                  and x["type"] == t and not x["family"].startswith("ALL") and x["common_eligible"]]
            diffs = {x["family"]: x["diff_class"] for x in fr if x["common_eligible"] >= 10}
            full = main[f"{ann}|A->B|{s}|{t}"]["diff_class"]
            sel = [r for r in G[(ann, "A->B")] if scope_ok(r, s) and in_type(r, t)]
            lofo = {}
            for fam in sorted({r["family"] for r in sel}):
                o = summarise([r for r in sel if r["family"] != fam], with_ci=False)
                lofo[fam] = o["diff_class"]
            fb[f"{ann}|{s}|{t}"] = dict(
                families_with_ge10_common=len(diffs), per_family_diff_class=diffs,
                family_balanced_mean_diff=round(statistics.mean(diffs.values()), 4) if diffs else None,
                families_star3d_higher=sum(v > 0 for v in diffs.values()), families_rfam_higher=sum(v < 0 for v in diffs.values()),
                families_equal=sum(v == 0 for v in diffs.values()), pooled_diff=full,
                leave_one_family_out=lofo, lofo_range=[min(lofo.values()), max(lofo.values())] if lofo else None)
summary["family_level"] = fb
json.dump(summary, open(f"{OUT}/interaction_summary.json", "w"), indent=1)
with open(f"{OUT}/interaction_preservation_by_family.tsv", "w", newline="") as fh:
    fns = list(dict.fromkeys(k for r in by_fam_rows for k in r))
    w = csv.DictWriter(fh, fieldnames=fns, delimiter="\t", restval=""); w.writeheader(); w.writerows(by_fam_rows)

# ------------------------------------------------------------------ region audit (A side regions from Step 4)
members = collections.defaultdict(set)
with gzip.open(f"{C4}/region_members.tsv.gz", "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["side"] == "A":
            members[r["region_id"]].add(r["unit_id"])
old = {r["region_id"]: r for r in csv.DictReader(open(f"{C5}/candidate_corrections.tsv"), delimiter="\t")}
ix = collections.defaultdict(list)
for r in rows:
    if r["direction"] == "A->B" and r["type"] != "no_LW_class":
        ix[(r["pair_id"], r["annotator"])].append(r)
reg_rows = []
for reg in csv.DictReader(open(f"{C4}/regions.tsv"), delimiter="\t"):
    if reg["side"] != "A":
        continue
    pid, units = reg["pair_id"], members[reg["region_id"]]
    o = old.get(reg["region_id"], {})
    rec = dict(region_id=reg["region_id"], pair_id=pid, family=reg["pair_id"].split("__")[0], length=int(reg["length"]),
               ss_composition=reg["ss_composition"], old_verdict=o.get("verdict", ""),
               median_d_rfam=o.get("median_d_rfam", ""), median_d_star3d=o.get("median_d_star3d", ""),
               distance_favours=o.get("distance_favours", ""))
    for ann, tag in (("MC-Annotate", "mca"), ("FR3D", "fr3d")):
        sel = [r for r in ix[(pid, ann)] if r["u1"] in units or r["u2"] in units]
        el = [r for r in sel if r["rfam_state"] == "mapped" and r["star3d_state"] == "mapped"]
        g = lambda m, d: sum(r[f"{m}_{d}"] == "yes" for r in el)  # noqa: E731
        only_s = [r for r in el if r["star3d_class"] == "yes" and r["rfam_class"] != "yes"]
        only_r = [r for r in el if r["rfam_class"] == "yes" and r["star3d_class"] != "yes"]
        rec.update({f"{tag}_source_interactions": len(sel), f"{tag}_eligible_both": len(el),
                    f"{tag}_rfam_eligible": sum(r["rfam_state"] == "mapped" for r in sel),
                    f"{tag}_star3d_eligible": sum(r["star3d_state"] == "mapped" for r in sel),
                    f"{tag}_rfam_any": g("rfam", "any"), f"{tag}_star3d_any": g("star3d", "any"),
                    f"{tag}_rfam_class": g("rfam", "class"), f"{tag}_star3d_class": g("star3d", "class"),
                    f"{tag}_rfam_class_ident": g("rfam", "class_ident"), f"{tag}_star3d_class_ident": g("star3d", "class_ident"),
                    f"{tag}_gain_star3d": len(only_s), f"{tag}_gain_rfam": len(only_r),
                    f"{tag}_gain_star3d_noncanonical": sum(r["type"] == "noncanonical" for r in only_s),
                    f"{tag}_gain_rfam_noncanonical": sum(r["type"] == "noncanonical" for r in only_r)})
    m_s, m_r, f_s, f_r = (rec["mca_gain_star3d"], rec["mca_gain_rfam"], rec["fr3d_gain_star3d"], rec["fr3d_gain_rfam"])
    net_m, net_f = m_s - m_r, f_s - f_r
    dist = rec["distance_favours"]
    if rec["mca_eligible_both"] == 0 and rec["fr3d_eligible_both"] == 0:
        tier = "unassessable"
    elif m_s == m_r == f_s == f_r == 0:
        tier = "insufficient_no_interaction_difference"
    elif (net_m > 0 and net_f < 0) or (net_m < 0 and net_f > 0):
        tier = "mixed_annotators_disagree"
    elif net_m == 0 and net_f == 0:
        tier = "tied"
    else:
        side = "STAR3D" if (net_m > 0 or net_f > 0) else "Rfam"
        other = "Rfam" if side == "STAR3D" else "STAR3D"
        gains_other = (m_r + f_r) if side == "STAR3D" else (m_s + f_s)
        both_agree = (net_m > 0 and net_f > 0) if side == "STAR3D" else (net_m < 0 and net_f < 0)
        net_min = min(abs(net_m), abs(net_f))
        if dist == other:
            tier = "mixed_distance_conflicts"
        elif both_agree and net_min >= 2 and gains_other == 0:
            tier = f"{side}_supported_strong"
        elif both_agree:
            tier = f"{side}_supported_moderate"
        else:
            tier = f"{side}_supported_weak_single_annotator"
    rec["audit_tier"] = tier            # interaction tier with the (STAR3D-frame) distance veto, for comparison only
    # primary tier: interaction evidence only (distance is measured in STAR3D's own superposition, so it is reported
    # but never decides); strong = both annotators net >= 2 the same way and no class gain for the other method
    itier = region_tier(rec["mca_eligible_both"], rec["fr3d_eligible_both"], m_s, m_r, f_s, f_r)
    rec["interaction_tier"] = itier
    reg_rows.append(rec)
with open(f"{OUT}/region_support_audit.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(reg_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(reg_rows)
rs = dict(regions=len(reg_rows),
          interaction_tiers=dict(collections.Counter(r["interaction_tier"] for r in reg_rows)),
          interaction_tiers_by_family={f: dict(collections.Counter(r["interaction_tier"] for r in reg_rows if r["family"] == f))
                                       for f in sorted({r["family"] for r in reg_rows})},
          interaction_tier_vs_distance={t: dict(collections.Counter(r["distance_favours"] for r in reg_rows if r["interaction_tier"] == t))
                                        for t in sorted({r["interaction_tier"] for r in reg_rows})},
          tiers_with_distance_veto=dict(collections.Counter(r["audit_tier"] for r in reg_rows)),
          old_verdicts=dict(collections.Counter(r["old_verdict"] for r in reg_rows)),
          old_vs_tier={k: dict(v) for k, v in
                       {ov: collections.Counter(r["audit_tier"] for r in reg_rows if r["old_verdict"] == ov)
                        for ov in sorted({r["old_verdict"] for r in reg_rows})}.items()},
          old_supports_star3d_by_mca_gain=dict(collections.Counter(min(r["mca_gain_star3d"] - r["mca_gain_rfam"], 5)
                                                                    for r in reg_rows if r["old_verdict"] == "supports_STAR3D")),
          length_by_tier={t: round(statistics.mean([r["length"] for r in reg_rows if r["audit_tier"] == t]), 2)
                          for t in sorted({r["audit_tier"] for r in reg_rows})},
          strong_star3d_regions=sorted({r["region_id"] for r in reg_rows if r["audit_tier"] == "STAR3D_supported_strong"}),
          strong_rfam_regions=sorted({r["region_id"] for r in reg_rows if r["audit_tier"] == "Rfam_supported_strong"}))
json.dump(rs, open(f"{OUT}/region_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in rs.items() if "regions" not in k or k == "regions"}, indent=1))
