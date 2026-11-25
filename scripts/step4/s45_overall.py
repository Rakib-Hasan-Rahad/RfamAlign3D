"""Step 4 / part 5: comprehensive alignment-difference metric (RNA A view) + direction-independent set summary.

  python scripts/step4/s45_overall.py <cohort>

RNA A view, one outcome per A-side nucleotide (from residue_comparisons, unchanged frozen tables):
  same_partner      Rfam and STAR3D give the same counterpart        assessable if A and the counterpart are observed
  different_partner both map it, to different counterparts           assessable if A and BOTH counterparts are observed
  rfam_only         Rfam maps it; STAR3D leaves it unmapped           assessable if A and the Rfam counterpart are observed
  star3d_only       STAR3D maps it; Rfam leaves it unmapped           assessable if A and the STAR3D counterpart are observed
  neither           neither maps it (not a difference; excluded from the denominator)
Overall difference rate = (D + R + S) / (A + D + R + S), assessable observations only; unassessable and neither are
reported separately. Correspondence-set summary counts EDGES (residue pairs), not nucleotides, using the same
coordinate rule (both endpoints observed) as pair_scores.
Outputs: outputs/04_comparison/<cohort>/overall/overall_metric.json, family_overall.tsv, pair_overall.tsv
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402

CATS = ("same_partner", "different_partner", "rfam_only", "star3d_only")


def main():
    cohort = sys.argv[1]
    base = os.path.join(K.ROOT, "outputs", "04_comparison", cohort)
    out = os.path.join(base, "overall")
    os.makedirs(out, exist_ok=True)
    snap = os.path.join(K.ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    fam = {p["pair_id"]: (p["family"], p["family_id"]) for p in K.tsv(os.path.join(snap, "pairs_to_run.tsv"))}
    res = K.tsv(os.path.join(base, "residue_comparisons.tsv.gz"))
    bco = {(r["pair_id"], r["unit_id"]): r["coords"] == "yes" for r in res if r["side"] == "B"}
    per_pair = collections.defaultdict(collections.Counter)
    seen = collections.Counter()
    for r in res:
        if r["side"] != "A":
            continue
        pid, c = r["pair_id"], r["category"]
        seen[(pid, r["unit_id"])] += 1
        a = r["coords"] == "yes"
        rp = bco.get((pid, r["rfam_partner"]), False) if r["rfam_partner"] else None
        sp = bco.get((pid, r["star3d_partner"]), False) if r["star3d_partner"] else None
        if c == "neither":
            key = "neither"
        elif c == "same_partner":
            key = c if a and rp else c + "_unassessable"
        elif c == "different_partner":
            key = c if (a and rp and sp) else c + "_unassessable"
        elif c == "rfam_only":
            key = c if a and rp else c + "_unassessable"
        else:  # star3d_only
            key = c if a and sp else c + "_unassessable"
        per_pair[pid][key] += 1
        per_pair[pid]["total"] += 1
    checks = {}
    checks["every_A_observation_one_outcome"] = all(v == 1 for v in seen.values())

    def summarize(counters):
        t = collections.Counter()
        for c in counters:
            t.update(c)
        A, D, R, S = (t[k] for k in CATS)
        return dict(total_observations=t["total"], same_partner=A, different_partner=D, rfam_only=R, star3d_only=S,
                    differences=D + R + S, denominator=A + D + R + S,
                    rate=round((D + R + S) / (A + D + R + S), 4) if A + D + R + S else None,
                    different_partner_rate=f"{D}/{A + D}",
                    neither=t["neither"], unassessable={k: t[k + "_unassessable"] for k in CATS},
                    unassessable_total=sum(t[k + "_unassessable"] for k in CATS))
    overall = summarize(per_pair.values())
    checks["partition_sums_to_total"] = overall["denominator"] + overall["neither"] + overall["unassessable_total"] == overall["total_observations"]
    checks["numerator_equals_categories"] = overall["differences"] == overall["different_partner"] + overall["rfam_only"] + overall["star3d_only"]
    # old metric reproducible
    old = json.load(open(os.path.join(base, "audit", "audit.json")))["Q1"]
    checks["different_partner_matches_audit"] = (overall["different_partner"], overall["same_partner"] + overall["different_partner"]) == (old["assessable_different"], old["jointly_assessable"])
    # direction-independent set summary (edges), from pair_scores (assessable edges)
    ps = K.tsv(os.path.join(base, "pair_scores.tsv"))
    sh = sum(int(p["shared_edges"]) for p in ps)
    ro = sum(int(p["rfam_only_edges"]) for p in ps)
    so = sum(int(p["star3d_only_edges"]) for p in ps)
    sets = dict(unit="residue correspondences (edges), both endpoints observed", shared=sh, rfam_only=ro, star3d_only=so,
                union=sh + ro + so, symmetric_difference=ro + so, set_difference_rate=round((ro + so) / (sh + ro + so), 4))
    checks["edge_sets_consistent"] = sh + ro == sum(int(p["rfam_edges"]) for p in ps) and sh + so == sum(int(p["star3d_edges"]) for p in ps)
    # offsets restricted to different-partner observations
    dg = json.load(open(os.path.join(base, "diagnostics", "summary.json")))
    offs = dict(one_position=dg["abs_offset"].get("1", 0), different_partner_observations=dg["observations"])
    checks["offset_population_is_different_partner"] = dg["observations"] == overall["different_partner"]
    fams = collections.defaultdict(list)
    for pid, c in per_pair.items():
        fams[fam[pid]].append(c)
    rows = []
    for (f, fid), cs in sorted(fams.items()):
        s_ = summarize(cs)
        rows.append(dict(family=f, family_id=fid, pairs=len(cs), **{k: v for k, v in s_.items() if not isinstance(v, dict)},
                         **{f"unassessable_{k}": v for k, v in s_["unassessable"].items()}))
    K.wtsv(os.path.join(out, "family_overall.tsv"), rows)
    checks["families_sum_to_total"] = sum(r["differences"] for r in rows) == overall["differences"] and sum(r["denominator"] for r in rows) == overall["denominator"]
    K.wtsv(os.path.join(out, "pair_overall.tsv"), [dict(pair_id=p, family=fam[p][0], **{k: c[k] for k in list(CATS) + ["neither", "total"]},
                                                        **{k + "_unassessable": c[k + "_unassessable"] for k in CATS}) for p, c in sorted(per_pair.items())])
    res_out = dict(cohort=cohort, view="RNA A view (one outcome per A-side nucleotide)", overall=overall, sets=sets, offsets=offs,
                   definitions=dict(rate="(different + Rfam only + STAR3D only) / (same + different + Rfam only + STAR3D only), assessable only",
                                    assessable="the A nucleotide and every counterpart it is mapped to have coordinates in the staged structures",
                                    excluded="neither-mapped (not a difference) and unassessable observations; both reported separately"),
                   checks=checks)
    json.dump(res_out, open(os.path.join(out, "overall_metric.json"), "w"), indent=1)
    print(json.dumps(dict(overall=overall, sets=sets, offsets=offs, checks=checks), indent=1))


if __name__ == "__main__":
    main()
