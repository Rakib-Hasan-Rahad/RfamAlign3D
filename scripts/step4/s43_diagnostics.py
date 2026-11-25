"""Step 4 / part 3: descriptive diagnostics of WHY two correspondences differ (patterns, not proven causes).

  python scripts/step4/s43_diagnostics.py <cohort>

Reads frozen Step 4 tables + the frozen Rfam pair reference; writes outputs/04_comparison/<cohort>/diagnostics/.
For every A-side scored different_partner observation:
  offset            B ungapped position of STAR3D's partner minus that of Rfam's partner (signed)
  region_offsets    the set of offsets inside its A-side disagreement region
  near_seed_gap     a gap ('-' or '.') in either seed row within +/-2 columns of the A residue's Rfam column
  context           RF class (consensus/insert) and SS_cons context of the A residue's column
  pair_disagreement fraction of both-mapped A-side residues that differ in this pair
Primary pattern (first rule that applies):
  global_misregistration     pair-level disagreement >= 50% of both-mapped residues
  register_shift_at_gap      region (>=2 nt) with one constant offset, |offset| <= 3, next to a seed gap
  register_shift             region (>=2 nt) with one constant offset, not next to a seed gap
  isolated_shift_at_gap      single-residue region next to a seed gap
  loop_rearrangement         offsets vary within the region and the residue is in a consensus loop/insert column
  other                      everything else
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402

LOOPS = {"hairpin_loop", "internal_loop", "bulge", "multiloop", "exterior"}


def main():
    cohort = sys.argv[1]
    base = os.path.join(K.ROOT, "outputs", "04_comparison", cohort)
    snap = os.path.join(K.ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    out = os.path.join(base, "diagnostics")
    os.makedirs(out, exist_ok=True)
    res = K.tsv(os.path.join(base, "residue_comparisons.tsv.gz"))
    bpos = {(r["pair_id"], r["unit_id"]): int(r["ungapped_position"]) for r in res if r["side"] == "B"}
    side = {(s["pair_id"], s["side"]): s for s in K.tsv(os.path.join(base, "side_scores.tsv"))}
    ref = {}
    for r in K.tsv(os.path.join(snap, "rfam_pair_reference.tsv.gz")):
        ref.setdefault(r["pair_id"], {})[int(r["rfam_column"])] = (r["query_char"], r["target_char"])
    mem = {}
    for m in K.tsv(os.path.join(base, "region_members.tsv.gz")):
        if m["side"] == "A":
            mem[(m["pair_id"], m["unit_id"])] = m["region_id"]
    obs = []
    for r in res:
        if r["side"] != "A" or r["scored_different"] != "yes":
            continue
        pid = r["pair_id"]
        off = bpos[(pid, r["star3d_partner"])] - bpos[(pid, r["rfam_partner"])]
        col = int(r["rfam_column"])
        cols = ref.get(pid, {})
        near = any(c in cols and ("-" in cols[c][0] + cols[c][1] or "." in cols[c][0] + cols[c][1])
                   for c in range(col - 2, col + 3))
        s = side[(pid, "A")]
        both = int(s["both_mapped"]) or 1
        obs.append(dict(pair_id=pid, unit_id=r["unit_id"], position=r["ungapped_position"], rfam_column=col,
                        rf_class=r["rf_class"], ss_context=r["ss_context"], offset=off,
                        near_seed_gap="yes" if near else "no", region_id=mem.get((pid, r["unit_id"]), ""),
                        pair_disagreement=round(int(s["scored_different"]) / both, 3)))
    # baseline: same rule for scored AGREEMENT observations (A side), so 'near a gap' can be compared fairly
    agree_near = agree_n = 0
    for r in res:
        if r["side"] == "A" and r["category"] == "same_partner" and r["coords"] == "yes" and r["rfam_partner_coords"] == "yes":
            col, cols = int(r["rfam_column"]), ref.get(r["pair_id"], {})
            agree_n += 1
            agree_near += any(c in cols and ("-" in cols[c][0] + cols[c][1] or "." in cols[c][0] + cols[c][1]) for c in range(col - 2, col + 3))
    by_region = collections.defaultdict(list)
    for o in obs:
        by_region[o["region_id"]].append(o)
    for rid, os_ in by_region.items():
        offs = {o["offset"] for o in os_}
        for o in os_:
            o["region_length"] = len(os_)
            o["region_offsets"] = ";".join(str(x) for x in sorted(offs))
            if o["pair_disagreement"] >= 0.5:
                pat = "global_misregistration"
            elif len(os_) >= 2 and len(offs) == 1 and abs(o["offset"]) <= 3 and any(x["near_seed_gap"] == "yes" for x in os_):
                pat = "register_shift_at_gap"
            elif len(os_) >= 2 and len(offs) == 1:
                pat = "register_shift"
            elif len(os_) == 1 and o["near_seed_gap"] == "yes":
                pat = "isolated_shift_at_gap"
            elif len(offs) > 1 and (o["ss_context"] in LOOPS or o["rf_class"] == "insert"):
                pat = "loop_rearrangement"
            else:
                pat = "other"
            o["pattern"] = pat
    K.wtsv(os.path.join(out, "disagreement_diagnostics.tsv.gz"), obs)
    summ = dict(observations=len(obs), pattern=dict(collections.Counter(o["pattern"] for o in obs)),
                abs_offset=dict(sorted(collections.Counter(min(abs(o["offset"]), 10) for o in obs).items())),
                near_seed_gap=dict(collections.Counter(o["near_seed_gap"] for o in obs)),
                near_seed_gap_baseline=dict(disagreements=f"{sum(o['near_seed_gap'] == 'yes' for o in obs)}/{len(obs)}",
                                            agreements=f"{agree_near}/{agree_n}"),
                ss_context=dict(collections.Counter(o["ss_context"] for o in obs)),
                rf_class=dict(collections.Counter(o["rf_class"] for o in obs)),
                by_family_pattern={})
    pf = collections.defaultdict(collections.Counter)
    for o in obs:
        pf[o["pair_id"].split("__")[0]][o["pattern"]] += 1
    summ["by_family_pattern"] = {k: dict(v) for k, v in pf.items()}
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
