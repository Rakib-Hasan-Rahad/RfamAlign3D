"""Step 4 / part 2: motif evidence at disagreement endpoints, AFTER the comparison tables were frozen.

  python scripts/step4/s42_motifs.py <cohort>

Refuses to run unless every table in freeze_before_motifs.json still has its frozen sha256 (motif data cannot
influence the comparison). Evidence types are kept separate:
  rfam_motif   Rfam 15.1 motif_matches hit on this exact seed row (rfamseq_acc/start/stop equal the row ID);
               positions query_start..query_stop = 1-based ungapped seed-row positions. Coverage: every seed row of
               release 15.1 (Rfam scans all seed alignments with its motif CMs; absence = assessed, no hit).
  bgsu_loop    BGSU RNA 3D Hub loop instance of this exact PDB, model 1, chain, containing the unit ID.
  bgsu_group   that exact loop's unit set equals an instance of a Motif Atlas 4.16 group (HL/IL/J3).
  consensus    Rfam SS_cons context of the endpoint's column (hairpin/internal/bulge/multiloop) - context only.
Per A-side scored different_partner observation, three endpoints are tested (A residue, B residue chosen by Rfam,
B residue chosen by STAR3D) and the observation is counted once.
Outputs in outputs/04_comparison/<cohort>/motif/.
"""
import collections
import csv
import gzip
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402

ROOT = K.ROOT
CACHE = os.path.join(ROOT, "cache")
LOOP_CTX = {"hairpin_loop", "internal_loop", "bulge", "multiloop"}


def ukey(u):
    p = u.split("|")
    return (p[0].upper(), p[1], p[2], p[4], p[7] if len(p) > 7 else "")


def load_rfam():
    hits = collections.defaultdict(list)
    with gzip.open(os.path.join(CACHE, "rfam15.1_motif", "motif_matches.txt.gz"), "rt") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            hits[(x[1], x[2], int(x[3]), int(x[4]))].append(dict(motif_acc=x[0], q1=int(x[5]), q2=int(x[6]),
                                                                  e_value=x[9], bits=x[10]))
    names = {}
    with gzip.open(os.path.join(CACHE, "rfam15.1_motif", "motif.txt.gz"), "rt") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            names[x[0]] = x[1] if len(x) > 1 else ""
    return hits, names


def rfam_coverage():
    """Family-level Rfam motif-scan coverage (verified, not assumed): a family is covered when it appears in
    motif_family_stats AND every stats row reconciles (num_hits = exact seed rows hit; frac_hits = hits / seed rows
    of the 15.1 seed, +/-0.0015). Otherwise coverage_unknown. Returns {family: (status, detail)}."""
    nrows = {}
    with gzip.open(os.path.join(ROOT, "outputs", "01_inventory", "family_tracks.jsonl.gz"), "rt") as f:
        for line in f:
            t = json.loads(line)
            nrows[t["family"]] = t["n_rows"]
    hit_rows = collections.defaultdict(set)
    with gzip.open(os.path.join(CACHE, "rfam15.1_motif", "motif_matches.txt.gz"), "rt") as f:
        for line in f:
            x = line.split("\t")
            hit_rows[(x[1], x[0])].add((x[2], x[3], x[4]))
    st = collections.defaultdict(list)
    with gzip.open(os.path.join(CACHE, "rfam15.1_motif", "motif_family_stats.txt.gz"), "rt") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            st[x[0]].append((x[1], int(x[2]), float(x[3])))
    cov = {}
    for fam, n in nrows.items():
        if fam not in st:
            cov[fam] = ("coverage_unknown", "family absent from motif_family_stats (scan scope not documented)")
            continue
        bad = [m for m, k, fr in st[fam] if len(hit_rows[(fam, m)]) != k or abs(len(hit_rows[(fam, m)]) / n - fr) > 0.0015]
        cov[fam] = ("covered", f"stats reconcile with {n} seed rows") if not bad else \
                   ("coverage_unknown", f"stats do not reconcile for {','.join(bad)}")
    return cov


def load_atlas():
    inst = {}
    for k in ("hl", "il", "j3"):
        grp = None
        for line in open(os.path.join(CACHE, "bgsu", f"atlas_{k}_4.16.csv")):
            line = line.strip()
            if line.startswith(">"):
                grp = line[1:]
            elif line:
                units = [u.strip('"') for u in line.split(",")]
                inst[frozenset(ukey(u) for u in units)] = grp
    pdbch = {(next(iter(s))[0], next(iter(s))[2]) for s in inst}
    return inst, pdbch


def load_loops(pdb):
    p = os.path.join(CACHE, "bgsu", "loops", f"loops_{pdb}.csv")
    if not os.path.exists(p):
        return None
    loops = {}
    for line in open(p):
        if '","' not in line:
            continue
        lid, units = line.strip().split('","', 1)
        us = [u for u in units.strip('"').split(",") if u.split("|")[1] == "1"]
        if us:
            loops[lid.strip('"')] = us
    return loops


def main():
    cohort = sys.argv[1]
    base = os.path.join(ROOT, "outputs", "04_comparison", cohort)
    fr = json.load(open(os.path.join(base, "freeze_before_motifs.json")))
    changed = [n for n, h in fr["sha256"].items() if K.sha(os.path.join(base, n)) != h]
    if changed:
        sys.exit(f"frozen comparison tables changed: {changed}")
    out = os.path.join(base, "motif")
    os.makedirs(out, exist_ok=True)
    snap = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    pairs = {p["pair_id"]: p for p in K.tsv(os.path.join(snap, "pairs_to_run.tsv"))}
    res = K.tsv(os.path.join(base, "residue_comparisons.tsv.gz"))
    by_unit = {(r["pair_id"], r["side"], r["unit_id"]): r for r in res}
    rhits, rnames = load_rfam()
    rcov = rfam_coverage()
    atlas, atlas_pdbch = load_atlas()
    loops_cache, cover = {}, []

    def loops_for(pdb):
        if pdb not in loops_cache:
            loops_cache[pdb] = load_loops(pdb)
        return loops_cache[pdb]

    def seed_key(fam, row):
        acc, rng = row.rsplit("/", 1)
        a, b = rng.split("-")
        return (fam, acc, int(a), int(b))

    def endpoint(pid, side, unit, role):
        """All evidence for one residue endpoint."""
        p = pairs[pid]
        r = by_unit.get((pid, side, unit))
        row = p["query_seed_row"] if side == "A" else p["target_seed_row"]
        pos = int(r["ungapped_position"])
        ev = []
        # Rfam motif hits on this exact seed row
        mh = [h for h in rhits.get(seed_key(p["family"], row), []) if h["q1"] <= pos <= h["q2"]]
        rf_status = "hit" if mh else ("assessed_no_hit" if rcov[p["family"]][0] == "covered" else "coverage_unknown")
        for h in mh:
            ev.append(dict(source="Rfam 15.1 motif_matches", evidence_type="Rfam computational motif hit",
                           record=h["motif_acc"], name=rnames.get(h["motif_acc"], ""), membership="seed-row position span",
                           detail=f"seed row {row} positions {h['q1']}-{h['q2']} (E={h['e_value']}, {h['bits']} bits)",
                           url=f"https://rfam.org/motif/{h['motif_acc']}"))
        # BGSU loops / Atlas
        pdb, ch = unit.split("|")[0], unit.split("|")[2]
        L = loops_for(pdb)
        if L is None:
            bg_status, at_status = "lookup_failed", "lookup_failed"
        else:
            in_chain = any(u.split("|")[2] == ch for us in L.values() for u in us)
            k = ukey(unit)
            hit = [(lid, us) for lid, us in L.items() if k in {ukey(u) for u in us}]
            bg_status = "hit" if hit else ("assessed_no_hit" if in_chain else "coverage_unknown")
            at_cov = (pdb.upper(), ch) in atlas_pdbch
            at_status = "not_covered"
            for lid, us in hit:
                grp = atlas.get(frozenset(ukey(u) for u in us))
                ev.append(dict(source="BGSU RNA 3D Hub loops", evidence_type="BGSU structural loop", record=lid, name="",
                               membership="exact unit-ID member", detail=f"{len(us)} members (model 1)",
                               url=f"http://rna.bgsu.edu/rna3dhub/loops/view/{lid}"))
                if grp:
                    ev.append(dict(source="BGSU Motif Atlas 4.16", evidence_type="BGSU motif group", record=grp, name="",
                                   membership="loop unit set equals Atlas instance", detail=lid,
                                   url=f"http://rna.bgsu.edu/rna3dhub/motif/view/{grp}"))
                    at_status = "hit"
            if at_status != "hit" and at_cov:
                at_status = "assessed_no_group"
        ctx = r["ss_context"]
        return dict(pair_id=pid, side=side, role=role, unit_id=unit, seed_row=row, ungapped_position=pos,
                    rfam_column=r["rfam_column"], rf_class=r["rf_class"], ss_context=ctx,
                    rfam_motif=rf_status, bgsu_loop=bg_status, bgsu_group=at_status,
                    consensus_loop_context="yes" if ctx in LOOP_CTX else ("unknown" if ctx == "unknown" else "no"),
                    evidence=ev)

    obs_rows, ev_rows, base_rows = [], [], []
    for r in res:
        if r["side"] != "A" or r["category"] not in ("different_partner", "same_partner"):
            continue
        q = endpoint(r["pair_id"], "A", r["unit_id"], "A_residue")
        rt = endpoint(r["pair_id"], "B", r["rfam_partner"], "B_rfam_partner")
        if r["category"] == "same_partner":
            if r["coords"] == "yes" and r["rfam_partner_coords"] == "yes":
                base_rows.append(dict(pair_id=r["pair_id"], family=pairs[r["pair_id"]]["family"], group="agreement",
                                      q=q, t=rt))
            continue
        if r["scored_different"] != "yes":
            continue
        st = endpoint(r["pair_id"], "B", r["star3d_partner"], "B_star3d_partner")
        eps = [q, rt, st]
        base_rows.append(dict(pair_id=r["pair_id"], family=pairs[r["pair_id"]]["family"], group="disagreement", q=q, t=rt))
        grp = any(e["bgsu_group"] == "hit" for e in eps)
        rfm = any(e["rfam_motif"] == "hit" for e in eps)
        named = grp or rfm
        loop = any(e["bgsu_loop"] == "hit" for e in eps)
        ctx = any(e["consensus_loop_context"] == "yes" for e in eps)
        complete = all(e["rfam_motif"] in ("hit", "assessed_no_hit") and e["bgsu_loop"] in ("hit", "assessed_no_hit")
                       for e in eps)
        cls = ("bgsu_group_and_rfam_motif" if grp and rfm else "bgsu_motif_group" if grp else
               "rfam_computational_motif" if rfm else "structural_loop_only" if loop else
               "consensus_loop_context_only" if ctx and complete else
               "no_documented_hit" if complete else "incomplete")
        oid = f"{r['pair_id']}:A:{r['unit_id']}"
        obs_rows.append(dict(observation_id=oid, pair_id=r["pair_id"], family=pairs[r["pair_id"]]["family"],
                             a_residue=r["unit_id"], b_rfam_partner=r["rfam_partner"], b_star3d_partner=r["star3d_partner"],
                             rfam_column=r["rfam_column"], rf_class=r["rf_class"], ss_context=r["ss_context"],
                             raw_line=r["raw_line"], classification=cls,
                             rfam_scan_coverage=("family_inferred" if rcov[pairs[r["pair_id"]]["family"]][0] == "covered" else "unknown"),
                             rfam_coverage_basis=("motif_family_stats reconcile with all seed rows of the family: inference that every seed "
                                                  "row (incl. these) was scanned; no row-level scan record exists"
                                                  if rcov[pairs[r["pair_id"]]["family"]][0] == "covered" else rcov[pairs[r["pair_id"]]["family"]][1]),
                             endpoints_with_named=";".join(e["role"] for e in eps if e["rfam_motif"] == "hit" or e["bgsu_group"] == "hit"),
                             endpoints_with_loop=";".join(e["role"] for e in eps if e["bgsu_loop"] == "hit"),
                             coverage=";".join(f"{e['role']}:rfam={e['rfam_motif']},loop={e['bgsu_loop']},group={e['bgsu_group']}" for e in eps)))
        for e in eps:
            for x in e["evidence"] or [dict(source="", evidence_type="none", record="", name="", membership="", detail="", url="")]:
                ev_rows.append(dict(observation_id=oid, role=e["role"], unit_id=e["unit_id"], seed_row=e["seed_row"],
                                    ungapped_position=e["ungapped_position"], rfam_motif=e["rfam_motif"],
                                    bgsu_loop=e["bgsu_loop"], bgsu_group=e["bgsu_group"],
                                    consensus_loop_context=e["consensus_loop_context"], **x))
    K.wtsv(os.path.join(out, "disagreement_observations.tsv"), obs_rows)
    K.wtsv(os.path.join(out, "overlap_evidence.tsv.gz"), ev_rows)

    # ---------------- equal-opportunity baseline (query-only and query-or-Rfam-target, identical coverage rules)
    def member(e, src):
        return {"bgsu_loop": e["bgsu_loop"] == "hit", "rfam_motif": e["rfam_motif"] == "hit",
                "bgsu_atlas_group": e["bgsu_group"] == "hit"}[src]

    def covered(e, src):
        return {"bgsu_loop": e["bgsu_loop"] in ("hit", "assessed_no_hit"),
                "rfam_motif": e["rfam_motif"] in ("hit", "assessed_no_hit"),
                "bgsu_atlas_group": e["bgsu_group"] in ("hit", "assessed_no_group")}[src]
    bl = []
    for scope, fams in [("pooled", None)] + [(f, {f}) for f in sorted({b["family"] for b in base_rows})]:
        for src in ("bgsu_loop", "rfam_motif", "bgsu_atlas_group"):
            for mode in ("query_only", "query_or_rfam_target"):
                row = dict(scope=scope, source=src, endpoints=mode)
                for g in ("disagreement", "agreement"):
                    xs = [b for b in base_rows if b["group"] == g and (fams is None or b["family"] in fams)]
                    if mode == "query_only":
                        ok = [b for b in xs if covered(b["q"], src)]
                        hit = sum(member(b["q"], src) for b in ok)
                    else:
                        ok = [b for b in xs if covered(b["q"], src) and covered(b["t"], src)]
                        hit = sum(member(b["q"], src) or member(b["t"], src) for b in ok)
                    row[f"{g}_hits"] = f"{hit}/{len(ok)}"
                    row[f"{g}_not_assessable"] = len(xs) - len(ok)
                bl.append(row)
    K.wtsv(os.path.join(out, "motif_baseline.tsv"), bl)

    # ---------------- region-level: members overlapping each source
    oc = {o["observation_id"]: o for o in obs_rows}
    rg = []
    mem = collections.defaultdict(list)
    for m in K.tsv(os.path.join(base, "region_members.tsv.gz")):
        if m["side"] == "A":
            mem[m["region_id"]].append(m)
    for rid, ms in sorted(mem.items()):
        cls = collections.Counter(oc.get(f"{m['pair_id']}:A:{m['unit_id']}", {}).get("classification", "not_observed") for m in ms)
        rg.append(dict(region_id=rid, pair_id=ms[0]["pair_id"], members=len(ms),
                       members_bgsu_motif_group=cls["bgsu_motif_group"] + cls["bgsu_group_and_rfam_motif"],
                       members_rfam_computational_motif=cls["rfam_computational_motif"] + cls["bgsu_group_and_rfam_motif"], members_structural_loop_only=cls["structural_loop_only"],
                       members_consensus_context_only=cls["consensus_loop_context_only"],
                       members_no_documented_hit=cls["no_documented_hit"], members_incomplete=cls["incomplete"]))
    K.wtsv(os.path.join(out, "region_motif_overlap.tsv"), rg)
    src = dict(rfam_motif_matches=dict(release="15.1", url="https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/motif_matches.txt.gz",
                                       sha256=K.sha(os.path.join(CACHE, "rfam15.1_motif", "motif_matches.txt.gz"))),
               bgsu_atlas=dict(release="4.16 (current at retrieval, checked 2026-11-13)",
                               sha256={k: K.sha(os.path.join(CACHE, "bgsu", f"atlas_{k}_4.16.csv")) for k in ("hl", "il", "j3")}),
               bgsu_loops=dict(url="http://rna.bgsu.edu/rna3dhub/loops/download/<PDB>", retrieved="2026-11-13",
                               files={p: K.sha(os.path.join(CACHE, "bgsu", "loops", f"loops_{p}.csv"))
                                      for p in sorted(loops_cache) if loops_cache[p] is not None}),
               rfam_scan_coverage={f: rcov[f] for f in sorted({p["family"] for p in pairs.values()})},
               frozen_comparison_verified=True, motif_predictor_run=False)
    json.dump(src, open(os.path.join(out, "sources.json"), "w"), indent=1)
    c = collections.Counter(o["classification"] for o in obs_rows)
    print("observations", len(obs_rows), dict(c))
    print([b for b in bl if b["scope"] == "pooled"])


if __name__ == "__main__":
    main()
