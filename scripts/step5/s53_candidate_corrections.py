"""Step 5 / part 3: which disagreement regions does the structural evidence resolve, and in whose favour?

  python scripts/step5/s53_candidate_corrections.py <cohort>

Unit: an A-side disagreement region (consecutive different-partner nucleotides, Step 4 regions.tsv).
Two independent lines of evidence per region, both already computed upstream:
  distance     median C1'-C1' distance (STAR3D superposition) from the region's nucleotides to their Rfam
               counterparts vs to their STAR3D counterparts (pair.json dR/dS). Favours the method whose median
               is smaller by > `margin` Angstrom (default 1.0); otherwise neutral.
  interactions MC-Annotate base pairs of RNA A with at least one nucleotide in the region
               (basepair_preservation_pairs). gain_STAR3D = kept (same class) only under STAR3D;
               gain_Rfam = kept only under Rfam. Favours the method with more gains; neutral if equal.
Verdict (pre-stated, symmetric):
  supports_STAR3D  interactions favour STAR3D and distance does not favour Rfam      -> candidate correction of Rfam
  supports_Rfam    interactions favour Rfam and distance does not favour STAR3D
  conflicting      interactions and distance favour different methods
  distance_only_*  no interaction difference; distance favours one method (weaker evidence)
  unresolved       neither line of evidence separates the two
A verdict is evidence for review, not proof: STAR3D chose the superposition, so distances favour it by design;
the interaction test does not depend on the superposition.
Outputs: outputs/05_results/<cohort>/candidate_corrections.tsv and candidate_corrections.json
"""
import collections
import csv
import gzip
import json
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    coh = sys.argv[1]
    margin = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    c4 = os.path.join(ROOT, "outputs", "04_comparison", coh)
    c5 = os.path.join(ROOT, "outputs", "05_results", coh)
    members = collections.defaultdict(list)
    for r in tsv(os.path.join(c4, "region_members.tsv.gz")):
        if r["side"] == "A":
            members[r["region_id"]].append(r["unit_id"])
    bp = collections.defaultdict(list)
    for r in tsv(os.path.join(c5, "basepair_preservation_pairs.tsv.gz")):
        bp[r["pair_id"]].append(r)
    cache = {}
    rows = []
    for reg in tsv(os.path.join(c4, "regions.tsv")):
        if reg["side"] != "A":
            continue
        pid = reg["pair_id"]
        if pid not in cache:
            cache[pid] = {x["u"]: x for x in json.load(open(os.path.join(c5, "app", "data", "pairs", pid, "pair.json")))["res"]["A"]}
        res = cache[pid]
        units = set(members[reg["region_id"]])
        dR = [res[u]["dR"] for u in units if res.get(u, {}).get("dR") is not None]
        dS = [res[u]["dS"] for u in units if res.get(u, {}).get("dS") is not None]
        mR, mS = (statistics.median(dR) if dR else None), (statistics.median(dS) if dS else None)
        dist = "neutral"
        if mR is not None and mS is not None:
            dist = "STAR3D" if mR - mS > margin else ("Rfam" if mS - mR > margin else "neutral")
        gs = gr = ncs = ncr = 0
        for b in bp[pid]:
            if b["u1"] in units or b["u2"] in units:
                s, f = b["star3d_outcome"] == "same_class", b["rfam_outcome"] == "same_class"
                nc = b["source_class"] not in ("", "cWW")
                if s and not f:
                    gs += 1; ncs += nc
                elif f and not s:
                    gr += 1; ncr += nc
        inter = "STAR3D" if gs > gr else ("Rfam" if gr > gs else "neutral")
        if inter != "neutral":
            other = "Rfam" if inter == "STAR3D" else "STAR3D"
            verdict = f"supports_{inter}" if dist != other else "conflicting"
        else:
            verdict = f"distance_only_{dist}" if dist != "neutral" else "unresolved"
        rows.append(dict(region_id=reg["region_id"], pair_id=pid, family=pid.split("__")[0], length=reg["length"],
                         start_unit=reg["start_unit"], end_unit=reg["end_unit"], ss_composition=reg["ss_composition"],
                         median_d_rfam=mR, median_d_star3d=mS, distance_favours=dist,
                         interactions_kept_only_star3d=gs, interactions_kept_only_rfam=gr,
                         noncanonical_kept_only_star3d=ncs, noncanonical_kept_only_rfam=ncr,
                         interactions_favour=inter, verdict=verdict))
    with open(os.path.join(c5, "candidate_corrections.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    v = collections.Counter(r["verdict"] for r in rows)
    nt = collections.Counter()
    for r in rows:
        nt[r["verdict"]] += int(r["length"])
    sup = [r for r in rows if r["verdict"] == "supports_STAR3D"]
    out = dict(cohort=coh, margin_angstrom=margin, regions=len(rows), verdicts=dict(v), nucleotides_by_verdict=dict(nt),
               supports_star3d=dict(regions=len(sup), pairs=len({r["pair_id"] for r in sup}),
                                    families=sorted({r["family"] for r in sup}),
                                    with_noncanonical_gain=sum(1 for r in sup if r["noncanonical_kept_only_star3d"] > 0)),
               supports_rfam=dict(regions=v["supports_Rfam"],
                                  pairs=len({r["pair_id"] for r in rows if r["verdict"] == "supports_Rfam"})))
    json.dump(out, open(os.path.join(c5, "candidate_corrections.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
