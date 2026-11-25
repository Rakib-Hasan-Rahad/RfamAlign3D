"""Step 4 / part 4: BGSU structural-loop context of every jointly assessable A-side observation (query residue only).

  python scripts/step4/s44_bgsu_context.py <cohort>

Population (identical to the motif baseline): A-side observations that are scored different_partner, or same_partner
with both endpoints observed. For each, the query residue's BGSU RNA 3D Hub loops (exact PDB, model 1, chain,
unit ID incl. insertion code) are read from the cached loops files.
  loop type  HL = hairpin loop, IL = internal loop, J3/J4/... = multi-helix junction (from the loop ID prefix)
  membership (non-exclusive): a residue can lie in several loops (e.g. a hairpin's closing pair also in an IL)
  exclusive category: HL only / IL only / junction only / several loop types / no BGSU loop / coverage unknown
  coverage unknown: the chain has no residue in any loop of the structure's loop file (cannot tell 'no loop')
Outputs: outputs/04_comparison/<cohort>/motif/bgsu_loop_context.tsv.gz and bgsu_loop_context.json
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402
import s42_motifs as M  # noqa: E402


def ltype(lid):
    p = lid.split("_")[0]
    return "junction" if p.startswith("J") else p


def frac(n, d):
    return f"{n}/{d}"


def main():
    cohort = sys.argv[1]
    base = os.path.join(K.ROOT, "outputs", "04_comparison", cohort)
    snap = os.path.join(K.ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    fam = {p["pair_id"]: (p["family"], p["family_id"]) for p in K.tsv(os.path.join(snap, "pairs_to_run.tsv"))}
    res = K.tsv(os.path.join(base, "residue_comparisons.tsv.gz"))
    loops, index, chains = {}, {}, {}

    def lookup(unit):
        pdb, ch = unit.split("|")[0], unit.split("|")[2]
        if pdb not in loops:
            L = M.load_loops(pdb) or {}
            loops[pdb] = L
            idx = collections.defaultdict(set)
            cs = set()
            for lid, us in L.items():
                for u in us:
                    idx[M.ukey(u)].add(lid)
                    cs.add(u.split("|")[2])
            index[pdb], chains[pdb] = idx, cs
        if ch not in chains[pdb]:
            return None
        return sorted(index[pdb].get(M.ukey(unit), set()))

    rows = []
    for r in res:
        if r["side"] != "A":
            continue
        ok_same = r["category"] == "same_partner" and r["coords"] == "yes" and r["rfam_partner_coords"] == "yes"
        if not (ok_same or r["scored_different"] == "yes"):
            continue
        lids = lookup(r["unit_id"])
        types = sorted({ltype(l) for l in lids}) if lids is not None else []
        if lids is None:
            cat = "coverage_unknown"
        elif not lids:
            cat = "no_bgsu_loop"
        elif len(types) > 1:
            cat = "several_loop_types"
        else:
            cat = {"HL": "hairpin_loop_only", "IL": "internal_loop_only", "junction": "junction_only"}[types[0]]
        rows.append(dict(pair_id=r["pair_id"], family=fam[r["pair_id"]][0], unit_id=r["unit_id"],
                         outcome="different" if r["scored_different"] == "yes" else "same",
                         loops=";".join(lids or []), loop_types=";".join(types), category=cat))
    K.wtsv(os.path.join(base, "motif", "bgsu_loop_context.tsv.gz"), rows)

    def rate(sel):
        xs = [x for x in rows if sel(x)]
        return frac(sum(x["outcome"] == "different" for x in xs), len(xs))
    covered = lambda x: x["category"] != "coverage_unknown"  # noqa: E731
    member = {t: rate(lambda x, t=t: t in x["loop_types"].split(";")) for t in ("HL", "IL", "junction")}
    member["any_bgsu_loop"] = rate(lambda x: covered(x) and x["category"] != "no_bgsu_loop")
    member["no_bgsu_loop"] = rate(lambda x: x["category"] == "no_bgsu_loop")
    cats = ["hairpin_loop_only", "internal_loop_only", "junction_only", "several_loop_types", "no_bgsu_loop", "coverage_unknown"]
    exclusive = {c: rate(lambda x, c=c: x["category"] == c) for c in cats}
    nd = sum(x["outcome"] == "different" for x in rows)
    na = len(rows) - nd
    share = {c: dict(disagreements=frac(sum(x["category"] == c and x["outcome"] == "different" for x in rows), nd),
                     agreements=frac(sum(x["category"] == c and x["outcome"] == "same" for x in rows), na)) for c in cats}
    hl_out = rate(lambda x: covered(x) and "HL" not in x["loop_types"].split(";"))
    by_fam = []
    for f in sorted({x["family"] for x in rows}):
        fx = [x for x in rows if x["family"] == f]
        g = lambda sel: frac(sum(x["outcome"] == "different" for x in fx if sel(x)), sum(1 for x in fx if sel(x)))  # noqa: E731
        by_fam.append(dict(family=f, name=next(v[1] for v in fam.values() if v[0] == f),
                           hairpin=g(lambda x: "HL" in x["loop_types"].split(";")),
                           internal=g(lambda x: "IL" in x["loop_types"].split(";")),
                           junction=g(lambda x: "junction" in x["loop_types"].split(";")),
                           no_loop=g(lambda x: x["category"] == "no_bgsu_loop"),
                           unknown=sum(1 for x in fx if x["category"] == "coverage_unknown"),
                           total=g(lambda x: True)))
    out = dict(cohort=cohort, population="A-side query residues: scored different_partner + same_partner with both endpoints observed",
               observations=len(rows), disagreements=nd, agreements=na,
               membership_rate=member, hairpin_vs_outside=dict(hairpin=member["HL"], outside_hairpin=hl_out),
               exclusive_rate=exclusive, exclusive_share=share, by_family=by_fam,
               rules=dict(membership="non-exclusive: a residue in several loops counts in each loop type",
                          exclusive="each observation in exactly one category; shares sum to 100%",
                          coverage_unknown="chain absent from every loop of the structure's BGSU loop file; excluded from no-loop",
                          source="BGSU RNA 3D Hub loops/download/<PDB>, model 1, exact unit IDs (retrieved 2026-11-13)"))
    json.dump(out, open(os.path.join(base, "motif", "bgsu_loop_context.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("observations", "membership_rate", "hairpin_vs_outside", "exclusive_rate")}, indent=1))


if __name__ == "__main__":
    main()
