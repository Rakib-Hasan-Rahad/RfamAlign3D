"""Step 5 / part 2: cohort-wide base-pair preservation under each alignment (Rfam seed vs STAR3D).

  python scripts/step5/s52_basepair_preservation.py <cohort> [--reps 2000] [--seed 20261125]

Input: the per-pair records written by build_app.py (outputs/05_results/<cohort>/app/data/pairs/*/pair.json), which
carry, for both RNAs, MC-Annotate base pairs (the binary bundled with STAR3D v1.2, run by STAR3D's own preprocessing
on the staged files) and both alignments' counterparts. Nothing is recomputed from structures here.

For every base pair (u1, u2) MC-Annotate annotates in RNA A, each alignment maps u1 -> t1 and u2 -> t2 in RNA B:
  same_class      MC-Annotate annotates t1-t2 with the same Leontis-Westhof class, compared in the SAME nucleotide
                  order (the target class is flipped when MC-Annotate lists t2 before t1: cWH listed reversed = cHW)
  different_class t1-t2 annotated with another class
  not_comparable  t1-t2 annotated, but one of the two classes is not W/H/S edge-based (no LW class)
  no_annotation   MC-Annotate lists no base pair t1-t2 (not proof of absence)
  unassessable    u1 or u2 unmapped, or t1/t2 without coordinates
Comparison population: base pairs assessable under BOTH alignments; reported for all such pairs and for the subset
touching a disagreement (u1 or u2 has category different_partner). Paired difference STAR3D - Rfam in the
same-class rate, with a 95% cluster-bootstrap interval resampling whole structure pairs.
Outputs: outputs/05_results/<cohort>/basepair_preservation.json and basepair_preservation_pairs.tsv.gz
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTCOMES = ("same_class", "different_class", "not_comparable", "no_annotation", "unassessable")


def classify(res_a, res_b, inter_b, u1, u2, lw, via):
    r1, r2 = res_a.get(u1), res_a.get(u2)
    t1, t2 = (r1 or {}).get(via), (r2 or {}).get(via)
    if not (r1 and r2 and t1 and t2 and r1.get("xyz") and r2.get("xyz")):
        return "unassessable", None
    x1, x2 = res_b.get(t1), res_b.get(t2)
    if not (x1 and x2 and x1.get("xyz") and x2.get("xyz")):
        return "unassessable", None
    hit = inter_b.get((t1, t2)) or inter_b.get((t2, t1))
    if hit is None:
        return "no_annotation", None
    cls = hit[2]
    tc = cls[0] + cls[2] + cls[1] if cls and hit[0] != t1 else cls
    if not lw or not tc:
        return "not_comparable", tc
    return ("same_class" if tc == lw else "different_class"), tc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261125)
    a = ap.parse_args()
    base = os.path.join(ROOT, "outputs", "05_results", a.cohort)
    rows, per_pair = [], {}
    for f in sorted(glob.glob(os.path.join(base, "app", "data", "pairs", "*", "pair.json"))):
        p = json.load(open(f))
        if not p.get("res") or not p.get("inter"):
            continue
        ra = {x["u"]: x for x in p["res"]["A"]}
        rb = {x["u"]: x for x in p["res"]["B"]}
        ib = {(x[0], x[1]): x for x in p["inter"]["B"]}
        cnt = collections.Counter()
        for u1, u2, lw, raw in p["inter"]["A"]:
            o_r, c_r = classify(ra, rb, ib, u1, u2, lw, "rp")
            o_s, c_s = classify(ra, rb, ib, u1, u2, lw, "sp")
            touch = any(ra.get(u, {}).get("c") == "different_partner" for u in (u1, u2))
            rows.append(dict(pair_id=p["id"], family=p["family"], u1=u1, u2=u2, source_class=lw or "", mc_annotate=raw,
                             rfam_outcome=o_r, rfam_target_class=c_r or "", star3d_outcome=o_s,
                             star3d_target_class=c_s or "", touches_disagreement="yes" if touch else "no"))
            both = o_r != "unassessable" and o_s != "unassessable"
            for scope, ok in (("all", both), ("touching_disagreement", both and touch)):
                if ok:
                    cnt[(scope, "n")] += 1
                    cnt[(scope, "rfam_same")] += o_r == "same_class"
                    cnt[(scope, "star_same")] += o_s == "same_class"
        per_pair[p["id"]] = (p["family"], cnt)
    with gzip.open(os.path.join(base, "basepair_preservation_pairs.tsv.gz"), "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    pids = sorted(per_pair)
    rng = np.random.default_rng(a.seed)
    boot = rng.integers(0, len(pids), size=(a.reps, len(pids)))
    out = dict(cohort=a.cohort, source="MC-Annotate (bundled with STAR3D v1.2), base pairs only",
               rule="class compared in the same nucleotide order; missing annotation is not proof of absence",
               base_pairs_in_A=len(rows), structure_pairs=len(pids))
    for meth, key in (("rfam", "rfam_outcome"), ("star3d", "star3d_outcome")):
        out[f"{meth}_outcomes_all_A_pairs"] = dict(collections.Counter(r[key] for r in rows))
        out[f"{meth}_outcomes_touching_disagreement"] = dict(collections.Counter(
            r[key] for r in rows if r["touches_disagreement"] == "yes"))
    for scope in ("all", "touching_disagreement"):
        n = np.array([per_pair[p][1][(scope, "n")] for p in pids], float)
        rs = np.array([per_pair[p][1][(scope, "rfam_same")] for p in pids], float)
        ss = np.array([per_pair[p][1][(scope, "star_same")] for p in pids], float)
        N = n[boot].sum(1)
        ok = N > 0
        d = (ss[boot].sum(1)[ok] - rs[boot].sum(1)[ok]) / N[ok]
        pc = lambda v: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]  # noqa: E731
        out[scope] = dict(assessable_under_both=int(n.sum()),
                          rfam_same_class=f"{int(rs.sum())}/{int(n.sum())}",
                          star3d_same_class=f"{int(ss.sum())}/{int(n.sum())}",
                          difference_star3d_minus_rfam=round(float((ss.sum() - rs.sum()) / n.sum()), 4) if n.sum() else None,
                          ci95_pair_bootstrap=pc(d) if len(d) else None,
                          rfam_ci95=pc(rs[boot].sum(1)[ok] / N[ok]), star3d_ci95=pc(ss[boot].sum(1)[ok] / N[ok]))
    # ---- by interaction type, at base pairs touching a disagreement and assessable under both alignments
    def btype(r):
        c, m = r["source_class"], r["mc_annotate"]
        if not c:
            return None                          # no LW class: cannot be 'same class' under either method
        if c == "cWW":
            return "canonical_watson_crick" if m[:3] in ("A-U", "U-A", "G-C", "C-G") else "cww_non_watson_crick"
        return "noncanonical"
    idx = {p: i for i, p in enumerate(pids)}
    by = {}
    for t in ("canonical_watson_crick", "cww_non_watson_crick", "noncanonical"):
        n = np.zeros(len(pids)); rs = np.zeros(len(pids)); ss = np.zeros(len(pids))
        for r in rows:
            if (r["touches_disagreement"] == "yes" and btype(r) == t
                    and "unassessable" not in (r["rfam_outcome"], r["star3d_outcome"])):
                i = idx[r["pair_id"]]
                n[i] += 1
                rs[i] += r["rfam_outcome"] == "same_class"
                ss[i] += r["star3d_outcome"] == "same_class"
        N = n[boot].sum(1)
        ok = N > 0
        pc = lambda v: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]  # noqa: E731
        by[t] = dict(n=int(n.sum()), rfam_same=int(rs.sum()), star3d_same=int(ss.sum()),
                     rfam_rate=round(rs.sum() / n.sum(), 4), star3d_rate=round(ss.sum() / n.sum(), 4),
                     rfam_ci95=pc(rs[boot].sum(1)[ok] / N[ok]), star3d_ci95=pc(ss[boot].sum(1)[ok] / N[ok]),
                     difference_ci95=pc((ss[boot].sum(1)[ok] - rs[boot].sum(1)[ok]) / N[ok]),
                     pairs_contributing=int((n > 0).sum()))
    out["touching_disagreement_by_type"] = by
    out["type_rule"] = ("canonical_watson_crick = cWW A-U/U-A/G-C/C-G; cww_non_watson_crick = other cWW (e.g. G-U wobble); "
                        "noncanonical = any other Leontis-Westhof class; pairs without an LW class are excluded")
    json.dump(out, open(os.path.join(base, "basepair_preservation.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    main()
