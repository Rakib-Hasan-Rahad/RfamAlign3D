"""Step 5 interactive results app (single page + per-pair JSON), generated from sqarma_interim.sqlite and Step 4 diagnostics.

  python scripts/step5/build_app.py <cohort>
  python3 -m http.server 8766 --bind 127.0.0.1 --directory "outputs/05_results/step2_interim_440/app"   # http://127.0.0.1:8766/
"""
import collections
import glob
import gzip
import json
import os
import shutil
import sqlite3
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_step5 as B  # noqa: E402

ROOT = B.ROOT
PATTERN_TEXT = {
    "register_shift_at_gap": ("Constant-offset run next to a seed gap",
                              "Rule: a disagreement region of ≥2 consecutive residues with one shared offset of 1–3 positions, with at least one member within two columns of a gap in either seed row. The two alignments place an insertion/deletion at different positions."),
    "register_shift": ("Constant-offset run away from seed gaps",
                       "Rule: a region of ≥2 residues with one shared offset, not meeting the gap condition above."),
    "isolated_shift_at_gap": ("Single residue next to a seed gap",
                              "Rule: a one-residue region within two columns of a seed gap."),
    "global_misregistration": ("Pair-wide difference (≥50% of residues)",
                               "Rule: in this pair, at least half of the residues mapped by both methods have different partners. This is not necessarily one constant offset (see the audit figure shown below)."),
    "loop_rearrangement": ("Varying offsets in a loop or insert column",
                           "Rule: offsets differ within the region and the residue lies in a consensus loop (SS_cons) or an RF insert column."),
    "other": ("Other", "Rule: none of the above applies."),
}


def parse_c1(path):
    """C1' coordinates from STAR3D's superposition PDB: model 1 = RNA A (chain A), model 2 = RNA B (chain B)."""
    out, m = {}, 0
    for l in open(path):
        if l.startswith("MODEL"):
            m += 1
        elif l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "C1'":
            out[("A" if m == 1 else "B", l[22:26].strip() + l[26:27].strip())] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out


def dist(p, q):
    return round(sum((a - b) ** 2 for a, b in zip(p, q)) ** 0.5, 2) if p and q else None


MCA_RES = __import__("re").compile(r"^'?([A-Za-z0-9])'?(-?\d+)(?:\.([A-Za-z]))?$")


def parse_mca_pairs(path):
    """Base pairs from MC-Annotate (the binary bundled with STAR3D v1.2, run during STAR3D preprocessing).
    Returns [((chain, resid), (chain, resid), raw_label, lw)]; lw (e.g. cWW) only when both edges are W/H/S."""
    out, on = [], False
    if not path or not os.path.exists(path):
        return out
    for l in open(path):
        if l.startswith("Base-pairs"):
            on = True
            continue
        if on:
            if not l.strip() or l.endswith("-----\n"):
                break
            left, _, right = l.partition(" : ")
            a, b = left.strip().split("-", 1)
            ma, mb = MCA_RES.match(a), MCA_RES.match(b)
            if not (ma and mb):
                continue
            toks = right.split()
            edges = next((t for t in toks if "/" in t), "")
            cis = "cis" if " cis" in right else ("trans" if " trans" in right else "")
            lw = ""
            if edges and cis:
                e1, e2 = edges.split("/")
                if e1[:1] in "WHS" and e2[:1] in "WHS":
                    lw = ("c" if cis == "cis" else "t") + e1[0] + e2[0]
            out.append(((ma.group(1), ma.group(2) + (ma.group(3) or "")), (mb.group(1), mb.group(2) + (mb.group(3) or "")),
                        right.strip(), lw))
    return out


VIEWDIR = [""]


def enrich(det, pid, p, ret, db, staged):
    """Gapped columns, staged IDs, C1' distances in STAR3D's frame, RMS values and MC-Annotate base pairs."""
    det["cols"] = [dict(c=int(c), a=qc, b=tc, ap=int(qp) if qp else None, bp=int(tp) if tp else None, rf=rf)
                   for c, qc, tc, qp, tp, rf in q(db, """SELECT rfam_column, query_char, target_char, query_ungapped_position,
                       target_ungapped_position, rf_class FROM rfam_reference WHERE pair_id=? ORDER BY CAST(rfam_column AS INT)""", pid)]
    cand = {"A": p["query_candidate"], "B": p["target_candidate"]}
    for side in ("A", "B"):
        for x in det["res"][side]:
            sr = staged.get((cand[side], x["p"]), ("", ""))
            x["st"], x["sch"] = sr
    pdbp = os.path.join(ret, pid, "out", f"{pid}__forward.aln.pdb")
    c1 = parse_c1(pdbp) if os.path.exists(pdbp) else {}
    # viewer copy: residues renumbered 1..n per model (insertion codes make author numbers ambiguous for selection)
    vidx = {}
    if os.path.exists(pdbp):
        lines, m, cur, n = [], 0, None, 0
        for l in open(pdbp):
            if l.startswith("MODEL"):
                m, cur, n = m + 1, None, 0
                lines.append(l)
            elif l.startswith(("ATOM", "HETATM")):
                key = l[22:26].strip() + l[26:27].strip()
                if key != cur:
                    cur, n = key, n + 1
                    vidx[("A" if m == 1 else "B", key)] = n
                lines.append(l[:22] + f"{n:>4}" + " " + l[27:])
            else:
                lines.append(l)
        open(os.path.join(VIEWDIR[0], pid, "view.pdb"), "w").write("".join(lines))
        det["files"]["view"] = f"data/pairs/{pid}/view.pdb"
    byu = {s: {x["u"]: x for x in det["res"][s]} for s in ("A", "B")}
    coord = lambda side, u: c1.get((side, byu[side][u]["st"])) if u and u in byu[side] and byu[side][u]["st"] else None  # noqa: E731
    for side, other in (("A", "B"), ("B", "A")):
        for x in det["res"][side]:
            me = coord(side, x["u"])
            x["dR"] = dist(me, coord(other, x["rp"])) if x["rp"] else None
            x["dS"] = dist(me, coord(other, x["sp"])) if x["sp"] else None
            x["c1"] = me is not None
            x["vi"] = vidx.get((side, x["st"])) if x["st"] else None
    def rms(pairs_):
        ds = [dist(coord("A", a), coord("B", b)) for a, b in pairs_]
        ds = [d for d in ds if d is not None]
        return dict(value=round((sum(d * d for d in ds) / len(ds)) ** 0.5, 2) if ds else None, n=len(ds), of=len(pairs_))
    rp = [(x["u"], x["rp"]) for x in det["res"]["A"] if x["rp"]]
    sp = [(x["u"], x["sp"]) for x in det["res"]["A"] if x["sp"]]
    det["geom"] = dict(star3d_rmsd=det.get("rmsd"), c1_rms_star3d=rms(sp), c1_rms_rfam=rms(rp), c1_available=bool(c1),
                       frame="STAR3D superposition (.aln.pdb); distances measured without refitting")
    # MC-Annotate base pairs per RNA, mapped to unit IDs through staged chain + residue ID
    jb = one_job(db, pid)
    inter = {}
    for side, sid, ch, cnd in (("A", jb[0], jb[1], cand["A"]), ("B", jb[2], jb[3], cand["B"])):
        files = glob.glob(os.path.join(os.path.dirname(ret), "preprocess", f"{sid}__{ch}__*", "STAR3D_struct_info", f"{sid}.mca"))
        st2u = {(x["sch"], x["st"]): x["u"] for x in det["res"][side] if x["st"]}
        lst = []
        for (c1_, r1), (c2_, r2), raw, lw in parse_mca_pairs(files[0] if files else None):
            u1, u2 = st2u.get((c1_, r1)), st2u.get((c2_, r2))
            if u1 and u2:
                lst.append([u1, u2, lw, raw])
        inter[side] = lst
    def has(side, u1, u2):
        return any({a, b} == {u1, u2} for a, b, *_ in inter[side])
    comp = []
    for a1, a2, lw, raw in inter["A"]:
        r1, r2 = byu["A"][a1]["rp"], byu["A"][a2]["rp"]
        s1, s2 = byu["A"][a1]["sp"], byu["A"][a2]["sp"]
        comp.append(dict(a1=a1, a2=a2, lw=lw, raw=raw,
                         rfam=[r1, r2, ("present" if has("B", r1, r2) else "absent") if r1 and r2 else "unmapped"],
                         star=[s1, s2, ("present" if has("B", s1, s2) else "absent") if s1 and s2 else "unmapped"]))
    det["inter"] = dict(A=inter["A"], B=inter["B"], compare=comp,
                        source="MC-Annotate binary bundled with STAR3D v1.2 (tools/MC-Annotate), run by STAR3D Preprocess on each staged file; "
                               "base pairs only; classification as written by MC-Annotate (edge/edge, cis/trans); LW shorthand only where both edges are W/H/S")


def one_job(db, pid):
    return db.execute("SELECT a_star3d_id, a_chain, b_star3d_id, b_chain FROM jobs WHERE pair_id=?", (pid,)).fetchone()


def q(db, sql, *a):
    return db.execute(sql, a).fetchall()


def build(cohort):
    out = os.path.join(ROOT, "outputs", "05_results", cohort)
    db = sqlite3.connect(os.path.join(out, "sqarma_interim.sqlite"))
    H = B.headline(db)
    app = os.path.join(out, "app")
    if os.path.exists(app):
        shutil.rmtree(app)
    os.makedirs(os.path.join(app, "data", "pairs"))
    VIEWDIR[0] = os.path.join(app, "data", "pairs")
    diag_p = os.path.join(ROOT, "outputs", "04_comparison", cohort, "diagnostics")
    diag = {(d["pair_id"], d["unit_id"]): d for d in B.tsv(os.path.join(diag_p, "disagreement_diagnostics.tsv.gz"))}
    dsum = json.load(open(os.path.join(diag_p, "summary.json")))
    ret = os.path.join(ROOT, "outputs", "03_star3d", cohort, "returned", "runs")

    # ---------------- per-family A-side category totals
    famcat = collections.defaultdict(collections.Counter)
    for f, c, n in q(db, """SELECT p.family, c.category, COUNT(*) FROM residue_comparisons c JOIN pairs p USING(pair_id)
                          WHERE c.side='A' GROUP BY p.family, c.category"""):
        famcat[f][c] = n
    ctx_agree = collections.Counter(dict(q(db, """SELECT ss_context, COUNT(*) FROM residue_comparisons WHERE side='A' AND
                                   category='same_partner' AND coords='yes' AND rfam_partner_coords='yes' GROUP BY ss_context""")))
    fam_pattern = collections.defaultdict(collections.Counter)
    for d in diag.values():
        fam_pattern[d["pair_id"].split("__")[0]][d["pattern"]] += 1
    fams = []
    for f, fid, n, v, ns in H["fam"]:
        a, d = H["fam_ad"].get(f, (0, 0))
        cn, cd = H["fam_cons"].get(f, (0, 0))
        fams.append(dict(id=f, name=fid, pairs=n, valid=v, no_stacks=ns, assessable=a, different=d, cons=[cn, cd],
                         cats=dict(famcat[f]), patterns=dict(fam_pattern[f]), rfam_cov=H["rfam_cov"].get(f, "")))
    # ---------------- pairs index + per-pair files
    pairs = {r[0]: dict(zip([c[0] for c in db.execute("SELECT * FROM pairs").description], r)) for r in db.execute("SELECT * FROM pairs")}
    runs = {r[0]: r for r in q(db, "SELECT pair_id, step3_status, rmsd, declared, superposition FROM runs")}
    scores = {r[0]: dict(zip([c[0] for c in db.execute("SELECT * FROM pair_scores").description], r)) for r in db.execute("SELECT * FROM pair_scores")}
    sides = collections.defaultdict(dict)
    cols = [c[0] for c in db.execute("SELECT * FROM side_scores").description]
    for r in db.execute("SELECT * FROM side_scores"):
        d = dict(zip(cols, r))
        sides[d["pair_id"]][d["side"]] = d
    obs_cls = {r[0]: r for r in q(db, "SELECT pair_id||'|'||a_residue, classification, rfam_scan_coverage FROM motif_observations")}
    ev = collections.defaultdict(list)
    loopmem = {}

    def members(lid):
        pdb = lid.split("_")[1]
        if pdb not in loopmem:
            loopmem[pdb] = {}
            f = os.path.join(ROOT, "cache", "bgsu", "loops", f"loops_{pdb}.csv")
            if os.path.exists(f):
                for line in open(f):
                    if '","' in line:
                        a, b = line.strip().split('","', 1)
                        loopmem[pdb][a.strip('"')] = [u for u in b.strip('"').split(",") if u.split("|")[1] == "1"]
        return loopmem[pdb].get(lid, [])
    for oid, role, unit, et, rec, nm, url, det in q(db, """SELECT observation_id, role, unit_id, evidence_type, record, name, url, detail
                                                         FROM overlap_evidence WHERE evidence_type!='none'"""):
        e = dict(role=role, unit=unit, type=et, record=rec, name=nm, url=url, detail=det)
        if et == "BGSU structural loop":
            e["members"] = members(rec)
        elif et == "BGSU motif group":
            e["loop"] = det
        ev[oid].append(e)
    staged = {}
    for cid, pos, sr, sc in q(db, "SELECT candidate_id, ungapped_position, staged_resid, staged_chain FROM residues"):
        staged[(cid, int(pos))] = (sr, sc)
    index = []
    for pid, p in pairs.items():
        r, s = runs[pid], scores.get(pid)
        sa, sb = sides.get(pid, {}).get("A"), sides.get(pid, {}).get("B")
        pats = collections.Counter(d["pattern"] for (pp, _), d in diag.items() if pp == pid)
        index.append(dict(id=pid, family=p["family"], fname=p["family_id"], a=p["query_pdb_chain"], b=p["target_pdb_chain"],
                          a_sp=p["query_species"].split(" (")[0], b_sp=p["target_species"].split(" (")[0], status=r[1],
                          diff=int(sa["scored_different"]) if sa else None, both=int(sa["both_mapped"]) if sa else None,
                          shared=int(s["shared_edges"]) if s else None, rfam=int(s["rfam_edges"]) if s else None,
                          star=int(s["star3d_edges"]) if s else None, top=pats.most_common(1)[0][0] if pats else ""))
        det = dict(id=pid, family=p["family"], fname=p["family_id"], status=r[1], rmsd=r[2],
                   a=dict(chain=p["query_pdb_chain"], species=p["query_species"], seed=p["query_seed_row"]),
                   b=dict(chain=p["target_pdb_chain"], species=p["target_species"], seed=p["target_seed_row"]),
                   files={})
        dd = os.path.join(app, "data", "pairs", pid)
        os.makedirs(dd)
        for fn, key in ((f"{pid}__forward.aln", "aln"), (f"{pid}__forward.aln.pdb", "pdb")):
            src = os.path.join(ret, pid, "out", fn)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(dd, key + (".txt" if key == "aln" else ".pdb")))
                det["files"][key] = f"data/pairs/{pid}/{key}" + (".txt" if key == "aln" else ".pdb")
        if s:
            det["scores"] = {k: s[k] for k in ("rfam_edges", "star3d_edges", "shared_edges", "rfam_recovery", "star3d_agreement",
                                              "rfam_only_edges", "star3d_only_edges", "rfam_consensus_reproduced", "rfam_insert_reproduced")}
            det["sides"] = {k: {x: sides[pid][k][x] for x in ("nucleotides", "observed", "same_partner", "different_partner", "rfam_only",
                                                             "star3d_only", "neither", "both_mapped", "scored_different")} for k in ("A", "B")}
            res = {"A": [], "B": []}
            posB = {}
            for side in ("A", "B"):
                for u, pos, nt, cat, rp, sp, col, rf, ssc, coords, line, sc in q(db, """SELECT unit_id, ungapped_position, nucleotide, category,
                        rfam_partner, star3d_partner, rfam_column, rf_class, ss_context, coords, raw_line, scored_different FROM residue_comparisons
                        WHERE pair_id=? AND side=? ORDER BY CAST(ungapped_position AS INT)""", pid, side):
                    res[side].append(dict(u=u, p=int(pos), n=nt, c=cat, rp=rp, sp=sp, col=int(col), rf=rf, ss=ssc, xyz=coords == "yes",
                                          line=line, sc=sc == "yes"))
                    if side == "B":
                        posB[u] = int(pos)
            posA = {x["u"]: x["p"] for x in res["A"]}
            for x in res["A"]:
                x["rpp"] = posB.get(x["rp"]) if x["rp"] else None
                x["spp"] = posB.get(x["sp"]) if x["sp"] else None
                d = diag.get((pid, x["u"]))
                if d:
                    x["off"], x["pat"], x["gap"] = int(d["offset"]), d["pattern"], d["near_seed_gap"] == "yes"
                oc = obs_cls.get(f"{pid}|{x['u']}")
                if oc:
                    x["mcls"], x["rcov"] = oc[1], oc[2]
                    x["ev"] = ev.get(f"{pid}:A:{x['u']}", [])
            for x in res["B"]:
                x["rpp"] = posA.get(x["rp"]) if x["rp"] else None
                x["spp"] = posA.get(x["sp"]) if x["sp"] else None
            det["res"] = res
            det["regions"] = [dict(id=rid.split(":", 2)[2], side=sd, len=int(ln), start=su, end=eu, rf=rfc, ss=ssc)
                              for rid, sd, ln, su, eu, rfc, ssc in q(db, "SELECT region_id, side, length, start_unit, end_unit, rf_composition, ss_composition FROM regions WHERE pair_id=?", pid)]
        if s:
            enrich(det, pid, p, ret, db, staged)
        json.dump(det, open(os.path.join(dd, "pair.json"), "w"), separators=(",", ":"))
    index.sort(key=lambda x: (x["family"], -(x["diff"] or -1)))
    cons = H["rfam_consensus_reproduced"]
    aud = json.load(open(os.path.join(ROOT, "outputs", "04_comparison", cohort, "audit", "audit.json")))
    B.set_label(cohort)
    bctx_p = os.path.join(ROOT, "outputs", "04_comparison", cohort, "motif", "bgsu_loop_context.json")
    ov_p = os.path.join(ROOT, "outputs", "04_comparison", cohort, "overall")
    overall = json.load(open(os.path.join(ov_p, "overall_metric.json"))) if os.path.exists(os.path.join(ov_p, "overall_metric.json")) else None
    fam_overall = B.tsv(os.path.join(ov_p, "family_overall.tsv")) if overall else []
    sel = dict(q(db, "SELECT decision||':'||reason, COUNT(*) FROM chain_review WHERE technical_status!='' GROUP BY 1"))
    sel_other = dict(q(db, "SELECT decision, COUNT(*) FROM chain_review WHERE technical_status='' GROUP BY 1"))
    closed = B.json.load(open(os.path.join(ROOT, "outputs", "02_structure_preparation", "STEP2_CLOSED.json"))) if cohort.startswith("step2_final") else None
    # downloadable tables (portable copies of the authoritative TSVs)
    tdir = os.path.join(app, "data", "tables")
    os.makedirs(tdir, exist_ok=True)
    for src in [os.path.join(ov_p, "family_overall.tsv"), os.path.join(ov_p, "pair_overall.tsv"), os.path.join(ov_p, "overall_metric.json"),
                os.path.join(ROOT, "outputs", "04_comparison", cohort, "pair_scores.tsv"),
                os.path.join(ROOT, "outputs", "04_comparison", cohort, "side_scores.tsv"),
                os.path.join(ROOT, "outputs", "04_comparison", cohort, "motif", "disagreement_observations.tsv"),
                os.path.join(ROOT, "outputs", "04_comparison", cohort, "diagnostics", "summary.json"),
                os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort, "family_status.tsv"),
                os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort, "pairs_to_run.tsv")]:
        if os.path.exists(src):
            dst = os.path.basename(src) if "diagnostics" not in src else "difference_patterns_summary.json"
            shutil.copy2(src, os.path.join(tdir, dst))
    prov = os.path.join(ROOT, "outputs", "03_star3d", cohort, "returned", "run_provenance.tsv")
    if os.path.exists(prov):
        shutil.copy2(prov, os.path.join(tdir, "star3d_run_provenance.tsv"))
    shutil.make_archive(os.path.join(app, "data", "sqarma_tables"), "zip", tdir)
    summary = dict(overall=overall, fam_overall=fam_overall, selection=dict(reviewed=sel, other=sel_other,
                   search_incomplete_checks=len(closed["unresolved_checks"]) if closed else None), build=__import__("datetime").datetime.now().strftime("%Y%m%d%H%M%S"), bgsu_ctx=json.load(open(bctx_p)) if os.path.exists(bctx_p) else None,cohort=cohort, label=B.LABEL["v"], label_text=B.LABEL["text"], audit_q7=aud.get("Q7", {}), audit_batch=None, families_screened=H["families_screened"], families_with_pairs=H["families_with_pairs"],
                   pairs=H["pairs"], valid=H["valid"], families_valid=H["families_valid"], no_stacks=H["no_stacks"],
                   failures=H["failures"], assessable=H["assessable"], different=H["different"], shared=H["shared_edges"],
                   rfam_edges=H["rfam_edges"], star3d_edges=H["star3d_edges"], rfam_raw=H["rfam_raw"],
                   rfam_excluded=H["rfam_raw"] - H["rfam_edges"], cons=list(cons), ins=list(H["rfam_insert_reproduced"]),
                   obs_class=H["obs_class"], baseline={k: list(v) for k, v in H["baseline"].items()}, atlas_any=H["atlas_any"],
                   rfam_any=H["rfam_any"], diag=dsum, ctx_agree=dict(ctx_agree), patterns=PATTERN_TEXT, families=fams,
                   trna_5s=sum(f["pairs"] for f in fams if f["id"] in ("RF00005", "RF00001")))
    json.dump(summary, open(os.path.join(app, "data", "summary.json"), "w"), separators=(",", ":"))
    json.dump(index, open(os.path.join(app, "data", "pairs_index.json"), "w"), separators=(",", ":"))
    shutil.copy2(os.path.join(out, "sqarma_interim.sqlite"), os.path.join(app, "data", "sqarma.sqlite"))
    shutil.copy2(os.path.join(out, "professor_report.md"), os.path.join(app, "data", "professor_report.md"))
    shutil.copy2(os.path.join(os.path.dirname(os.path.abspath(__file__)), "app_index.html"), os.path.join(app, "index.html"))
    os.makedirs(os.path.join(app, "lib"), exist_ok=True)
    shutil.copy2(os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib", "3Dmol-min.js"), os.path.join(app, "lib", "3Dmol-min.js"))
    print("app:", os.path.relpath(app, ROOT), "| pairs", len(index))


if __name__ == "__main__":
    build(sys.argv[1])
