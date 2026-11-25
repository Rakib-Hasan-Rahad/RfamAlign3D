"""Step 4 / part 6: uncertainty and robustness for every headline rate (reads frozen Step 4 tables only).

  python scripts/step4/s46_statistics.py <cohort> [--reps 2000] [--seed 20261125]

Why: nucleotides of one pair are not independent (one STAR3D run decides them all) and tRNA contributes most
observations. Each rate is therefore reported with
  * a 95% cluster-bootstrap interval: whole PAIRS resampled with replacement, ratio of sums recomputed
    (percentile interval, fixed seed, `reps` replicates);
  * a family-balanced rate: the unweighted mean of per-family rates (each family counts once);
  * leave-one-family-out (LOFO): the pooled rate recomputed with each family removed; min and max reported.
Hairpin vs outside-hairpin (BGSU loop membership) is additionally tested with a Cochran-Mantel-Haenszel test
stratified by family (common odds ratio, Mantel-Haenszel estimate), alongside the cluster-bootstrap interval of the
rate difference.
Populations and assessability rules are exactly those of s43/s44/s45 (not redefined here).
Outputs: outputs/04_comparison/<cohort>/statistics/statistics.json and statistics.tsv
"""
import argparse
import collections
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402


def frac(s):
    a, b = s.split("/")
    return int(a), int(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261125)
    a = ap.parse_args()
    base = os.path.join(K.ROOT, "outputs", "04_comparison", a.cohort)
    out = os.path.join(base, "statistics")
    os.makedirs(out, exist_ok=True)

    # ---- per-pair numerator / denominator for each metric (pair_id -> {metric: [num, den]})
    M = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    fam = {}
    for r in K.tsv(os.path.join(base, "overall", "pair_overall.tsv")):
        p = r["pair_id"]
        fam[p] = r["family"]
        A, D, R, S = (int(r[k]) for k in ("same_partner", "different_partner", "rfam_only", "star3d_only"))
        M[p]["different_partner_rate"] = [D, A + D]
        M[p]["comprehensive_difference_rate"] = [D + R + S, A + D + R + S]
    for r in K.tsv(os.path.join(base, "pair_scores.tsv")):
        p = r["pair_id"]
        if p not in fam:
            continue
        sh, ro, so = (int(r[k]) for k in ("shared_edges", "rfam_only_edges", "star3d_only_edges"))
        M[p]["correspondence_set_difference"] = [ro + so, sh + ro + so]
        if r.get("rfam_consensus_reproduced") and "/" in r["rfam_consensus_reproduced"]:
            M[p]["rfam_consensus_reproduced"] = list(frac(r["rfam_consensus_reproduced"]))
    for r in K.tsv(os.path.join(base, "diagnostics", "disagreement_diagnostics.tsv.gz")):
        m = M[r["pair_id"]]["one_position_offset_share"]
        m[0] += abs(int(r["offset"])) == 1
        m[1] += 1
    loops = K.tsv(os.path.join(base, "motif", "bgsu_loop_context.tsv.gz"))
    for r in loops:
        if r["category"] == "coverage_unknown":
            continue
        t = r["loop_types"].split(";") if r["loop_types"] else []
        d = r["outcome"] == "different"
        for key, inside in (("hairpin_loop", "HL" in t), ("internal_loop", "IL" in t), ("junction", "junction" in t),
                            ("any_bgsu_loop", bool(t)), ("no_bgsu_loop", not t), ("outside_hairpin", "HL" not in t)):
            if inside:
                m = M[r["pair_id"]][key + "_difference_rate"]
                m[0] += d
                m[1] += 1
    pairs = sorted(fam)
    fams = sorted(set(fam.values()))
    metrics = ["different_partner_rate", "comprehensive_difference_rate", "correspondence_set_difference",
               "rfam_consensus_reproduced", "one_position_offset_share", "hairpin_loop_difference_rate",
               "outside_hairpin_difference_rate", "internal_loop_difference_rate", "junction_difference_rate",
               "any_bgsu_loop_difference_rate", "no_bgsu_loop_difference_rate"]
    arr = {k: np.array([M[p][k] if k in M[p] else [0, 0] for p in pairs], dtype=float) for k in metrics}
    rng = np.random.default_rng(a.seed)
    boot = rng.integers(0, len(pairs), size=(a.reps, len(pairs)))   # same resamples for every metric
    famidx = np.array([fams.index(fam[p]) for p in pairs])

    def ci(x):
        n, d = x[boot, 0].sum(1), x[boot, 1].sum(1)
        v = n[d > 0] / d[d > 0]
        return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]

    res, rows = {}, []
    for k in metrics:
        x = arr[k]
        n, d = int(x[:, 0].sum()), int(x[:, 1].sum())
        if not d:
            continue
        fr = {f: (int(x[famidx == i, 0].sum()), int(x[famidx == i, 1].sum())) for i, f in enumerate(fams)}
        fr = {f: v for f, v in fr.items() if v[1]}
        lofo = {f: round((n - v[0]) / (d - v[1]), 4) for f, v in fr.items() if d - v[1]}
        res[k] = dict(n=n, N=d, rate=round(n / d, 4), ci95_pair_bootstrap=ci(x),
                      pairs_contributing=int((x[:, 1] > 0).sum()),
                      family_balanced_rate=round(float(np.mean([v[0] / v[1] for v in fr.values()])), 4),
                      families_contributing=len(fr),
                      lofo_min=min(lofo.values()) if lofo else None, lofo_max=max(lofo.values()) if lofo else None,
                      lofo_family_with_largest_effect=max(lofo, key=lambda f: abs(lofo[f] - n / d)) if lofo else None,
                      by_family={f: f"{v[0]}/{v[1]}" for f, v in fr.items()})
        r = res[k]
        rows.append(dict(metric=k, n=n, N=d, rate=r["rate"], ci95_low=r["ci95_pair_bootstrap"][0],
                         ci95_high=r["ci95_pair_bootstrap"][1], pairs=r["pairs_contributing"],
                         families=r["families_contributing"], family_balanced_rate=r["family_balanced_rate"],
                         lofo_min=r["lofo_min"], lofo_max=r["lofo_max"],
                         lofo_largest_effect=r["lofo_family_with_largest_effect"]))

    # ---- hairpin vs outside: bootstrap difference + family-stratified CMH
    h, o = arr["hairpin_loop_difference_rate"], arr["outside_hairpin_difference_rate"]
    diff = (h[boot, 0].sum(1) / h[boot, 1].sum(1)) - (o[boot, 0].sum(1) / o[boot, 1].sum(1))
    tables, num, den = [], 0.0, 0.0
    for i, f in enumerate(fams):
        a1, n1 = h[famidx == i, 0].sum(), h[famidx == i, 1].sum()
        c1, n2 = o[famidx == i, 0].sum(), o[famidx == i, 1].sum()
        if n1 and n2:
            t = np.array([[a1, n1 - a1], [c1, n2 - c1]])
            tables.append(t)
            T = t.sum()
            num += t[0, 0] * t[1, 1] / T
            den += t[0, 1] * t[1, 0] / T
    # Mantel-Haenszel chi-square with continuity correction, computed directly
    A = sum(t[0, 0] for t in tables)
    E = sum(t[0].sum() * t[:, 0].sum() / t.sum() for t in tables)
    V = sum(t[0].sum() * t[1].sum() * t[:, 0].sum() * t[:, 1].sum() / (t.sum() ** 2 * (t.sum() - 1)) for t in tables)
    chi2 = (abs(A - E) - 0.5) ** 2 / V
    p = float(stats.chi2.sf(chi2, 1))
    res["hairpin_vs_outside"] = dict(
        hairpin=f"{int(h[:, 0].sum())}/{int(h[:, 1].sum())}", outside=f"{int(o[:, 0].sum())}/{int(o[:, 1].sum())}",
        rate_difference=round(float(h[:, 0].sum() / h[:, 1].sum() - o[:, 0].sum() / o[:, 1].sum()), 4),
        ci95_pair_bootstrap=[round(float(np.percentile(diff, 2.5)), 4), round(float(np.percentile(diff, 97.5)), 4)],
        cmh_common_odds_ratio_by_family=round(num / den, 3), cmh_chi2=round(float(chi2), 2), cmh_p=p,
        cmh_strata=len(tables))
    res["_method"] = dict(reps=a.reps, seed=a.seed, unit_resampled="pair", pairs=len(pairs), families=fams,
                          interval="percentile of the ratio of sums", populations="as defined in s43/s44/s45")
    K.wtsv(os.path.join(out, "statistics.tsv"), rows)
    json.dump(res, open(os.path.join(out, "statistics.json"), "w"), indent=1)
    for r in rows:
        print(f"{r['metric']:36s} {r['n']:>6}/{r['N']:<6} {r['rate']:.3f} [{r['ci95_low']:.3f}, {r['ci95_high']:.3f}] "
              f"fam-bal {r['family_balanced_rate']:.3f} LOFO {r['lofo_min']}-{r['lofo_max']} ({r['lofo_largest_effect']})")
    print("hairpin vs outside:", json.dumps(res["hairpin_vs_outside"]))


if __name__ == "__main__":
    main()
