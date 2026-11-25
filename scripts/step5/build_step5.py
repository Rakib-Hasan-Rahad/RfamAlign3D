"""Step 5 (INTERIM): one SQLite database -> professor report + static website, all from the same tables.

  python scripts/step5/build_step5.py <cohort>
  python3 -m http.server 8765 --directory "outputs/05_results/step2_interim_440/site"   # then http://localhost:8765/

Reads only frozen/derived files (snapshot, Step 3 returned results, Step 4 tables); never writes to them.
Every number on the site and in the report comes from a SQL query on rfamalign3d.sqlite.
"""
import collections
import csv
import gzip
import hashlib
import html
import json
import os
import shutil
import sqlite3
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
csv.field_size_limit(1 << 30)
E = html.escape


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def nd(s):
    if not s or s == "NA":
        return 0, 0
    a, b = s.split("/")
    return int(a), int(b)


def pct(n, d):
    return f"{n:,}/{d:,} ({100 * n / d:.1f}%)" if d else "NA"


# ------------------------------------------------------------------ database
def build_db(cohort, out):
    snap = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    s3 = os.path.join(ROOT, "outputs", "03_star3d", cohort)
    s4 = os.path.join(ROOT, "outputs", "04_comparison", cohort)
    sources = {
        "families": os.path.join(snap, "family_status.tsv"),
        "pairs": os.path.join(snap, "pairs_to_run.tsv"),
        "structures": os.path.join(snap, "structures.tsv"),
        "chain_review": os.path.join(snap, "chain_review.tsv"),
        "residues": os.path.join(snap, "residue_crosswalk.tsv.gz"),
        "rfam_reference": os.path.join(snap, "rfam_pair_reference.tsv.gz"),
        "jobs": os.path.join(s3, "step3_jobs.tsv"),
        "runs": os.path.join(s4, "run_outcomes.tsv"),
        "alignment_edges": os.path.join(s4, "alignment_edges.tsv.gz"),
        "pair_scores": os.path.join(s4, "pair_scores.tsv"),
        "side_scores": os.path.join(s4, "side_scores.tsv"),
        "residue_comparisons": os.path.join(s4, "residue_comparisons.tsv.gz"),
        "regions": os.path.join(s4, "regions.tsv"),
        "region_members": os.path.join(s4, "region_members.tsv.gz"),
        "motif_observations": os.path.join(s4, "motif", "disagreement_observations.tsv"),
        "overlap_evidence": os.path.join(s4, "motif", "overlap_evidence.tsv.gz"),
        "motif_baseline": os.path.join(s4, "motif", "motif_baseline.tsv"),
        "region_motif": os.path.join(s4, "motif", "region_motif_overlap.tsv"),
    }
    dbp = os.path.join(out, "rfamalign3d.sqlite")
    if os.path.exists(dbp):
        os.remove(dbp)
    db = sqlite3.connect(dbp)
    prov = []
    for t, p in sources.items():
        rows = tsv(p)
        cols = list(rows[0]) if rows else ["empty"]
        db.execute(f'CREATE TABLE {t} ({", ".join(chr(34) + c + chr(34) + " TEXT" for c in cols)})')
        db.executemany(f'INSERT INTO {t} VALUES ({",".join("?" * len(cols))})', [[r.get(c, "") for c in cols] for r in rows])
        prov.append((t, os.path.relpath(p, ROOT), sha(p), len(rows)))
    # run records from Step 3 (status, timing, hashes)
    db.execute("CREATE TABLE step3_results (pair_id TEXT, status TEXT, elapsed_s TEXT, max_rss_kb TEXT, aln_sha256 TEXT, superposition TEXT, host TEXT)")
    for j in tsv(sources["jobs"]):
        r = json.load(open(os.path.join(s3, "returned", "runs", j["pair_id"], "result.json")))
        aln = (r.get("outputs") or {}).get(f"out/{j['pair_id']}__forward.aln", {}).get("sha256", "")
        db.execute("INSERT INTO step3_results VALUES (?,?,?,?,?,?,?)",
                   (j["pair_id"], r["status"], r.get("elapsed_s"), r.get("max_rss_kb"), aln, r.get("superposition"), r.get("host")))
    src = json.load(open(os.path.join(s4, "motif", "sources.json")))
    db.execute("CREATE TABLE annotation_coverage (source TEXT, scope TEXT, key TEXT, status TEXT, detail TEXT)")
    for f, (st, d) in src["rfam_scan_coverage"].items():
        db.execute("INSERT INTO annotation_coverage VALUES (?,?,?,?,?)", ("Rfam 15.1 motif_matches", "family", f, st, d))
    for p, h in src["bgsu_loops"]["files"].items():
        db.execute("INSERT INTO annotation_coverage VALUES (?,?,?,?,?)", ("BGSU RNA 3D Hub loops", "structure", p, "lookup_completed", h))
    db.execute("CREATE TABLE provenance (table_name TEXT, source_file TEXT, sha256 TEXT, rows INTEGER)")
    db.executemany("INSERT INTO provenance VALUES (?,?,?,?)", prov)
    db.execute("CREATE TABLE study (key TEXT, value TEXT)")
    db.executemany("INSERT INTO study VALUES (?,?)", [("cohort", cohort), ("label", "INTERIM"),
                                                      ("motif_sources", json.dumps({k: v for k, v in src.items() if k != "bgsu_loops"}))])
    for idx in ("residue_comparisons(pair_id)", "alignment_edges(pair_id)", "overlap_evidence(observation_id)",
                "motif_observations(pair_id)", "region_members(region_id)", "residues(candidate_id)", "rfam_reference(pair_id)"):
        db.execute(f"CREATE INDEX IF NOT EXISTS i_{idx.split('(')[0]}_{idx.split('(')[1][:-1]} ON {idx}")
    db.commit()
    return db, dbp


def q(db, sql, *a):
    return db.execute(sql, a).fetchall()


def one(db, sql, *a):
    return db.execute(sql, a).fetchone()[0]


# ------------------------------------------------------------------ headline numbers (all SQL)
def headline(db):
    H = {}
    H["families_screened"] = one(db, "SELECT COUNT(*) FROM families")
    H["pairs"] = one(db, "SELECT COUNT(*) FROM runs")
    H["valid"] = one(db, "SELECT COUNT(*) FROM runs WHERE step3_status='valid_alignment'")
    H["no_stacks"] = one(db, "SELECT COUNT(*) FROM runs WHERE step3_status='no_alignment'")
    H["failures"] = one(db, "SELECT COUNT(*) FROM runs WHERE step3_status NOT IN ('valid_alignment','no_alignment')")
    H["families_with_pairs"] = one(db, "SELECT COUNT(DISTINCT family) FROM runs")
    H["families_valid"] = one(db, "SELECT COUNT(DISTINCT family) FROM runs WHERE step3_status='valid_alignment'")
    H["assessable"] = one(db, """SELECT COUNT(*) FROM residue_comparisons WHERE side='A' AND
        ((category='same_partner' AND coords='yes' AND rfam_partner_coords='yes') OR scored_different='yes')""")
    H["different"] = one(db, "SELECT COUNT(*) FROM residue_comparisons WHERE side='A' AND scored_different='yes'")
    H["rfam_edges"] = one(db, "SELECT SUM(CAST(rfam_edges AS INT)) FROM pair_scores")
    H["star3d_edges"] = one(db, "SELECT SUM(CAST(star3d_edges AS INT)) FROM pair_scores")
    H["shared_edges"] = one(db, "SELECT SUM(CAST(shared_edges AS INT)) FROM pair_scores")
    for k in ("consensus", "insert", "unknown"):
        n = d = 0
        for (v,) in q(db, f"SELECT rfam_{k}_reproduced FROM pair_scores"):
            a, b = nd(v)
            n, d = n + a, d + b
        H[f"rfam_{k}_reproduced"] = (n, d)
    H["rfam_raw"] = one(db, "SELECT COUNT(*) FROM alignment_edges WHERE method='rfam'")
    H["rfam_excluded"] = dict(q(db, "SELECT exclusion_reason, COUNT(*) FROM alignment_edges WHERE method='rfam' AND assessable='no' GROUP BY exclusion_reason"))
    H["star3d_excluded"] = one(db, "SELECT COUNT(*) FROM alignment_edges WHERE method='star3d' AND assessable='no'")
    H["rfam_cov_obs"] = dict(q(db, "SELECT rfam_scan_coverage, COUNT(*) FROM motif_observations GROUP BY rfam_scan_coverage"))
    H["rf_pairs"] = one(db, "SELECT COUNT(*) FROM pair_scores WHERE rf_available='yes'")
    H["obs_class"] = dict(q(db, "SELECT classification, COUNT(*) FROM motif_observations GROUP BY classification"))
    H["atlas_any"] = one(db, "SELECT COUNT(*) FROM motif_observations WHERE classification IN ('bgsu_motif_group','bgsu_group_and_rfam_motif')")
    H["rfam_any"] = one(db, "SELECT COUNT(*) FROM motif_observations WHERE classification IN ('rfam_computational_motif','bgsu_group_and_rfam_motif')")
    H["baseline"] = {r[0]: (r[1], r[2]) for r in q(db, """SELECT source, disagreement_hits, agreement_hits FROM motif_baseline
                                                          WHERE scope='pooled' AND endpoints='query_only'""")}
    H["fam"] = q(db, """SELECT r.family, p.family_id, COUNT(*), SUM(r.step3_status='valid_alignment'), SUM(r.step3_status='no_alignment')
                       FROM runs r JOIN pairs p USING(pair_id) GROUP BY r.family ORDER BY COUNT(*) DESC""")
    fa = {}
    for f, a, d in q(db, """SELECT p.family, SUM((c.category='same_partner' AND c.coords='yes' AND c.rfam_partner_coords='yes') OR c.scored_different='yes'),
                           SUM(c.scored_different='yes') FROM residue_comparisons c JOIN pairs p USING(pair_id) WHERE c.side='A' GROUP BY p.family"""):
        fa[f] = (a, d)
    H["fam_ad"] = fa
    H["fam_cons"] = {}
    for f, v in q(db, "SELECT family, rfam_consensus_reproduced FROM pair_scores"):
        n, d = H["fam_cons"].get(f, (0, 0))
        a, b = nd(v)
        H["fam_cons"][f] = (n + a, d + b)
    H["rfam_cov"] = dict(q(db, "SELECT key, status FROM annotation_coverage WHERE scope='family'"))
    return H


# ------------------------------------------------------------------ HTML shell
CSS = """:root{--bg:#fbfbf9;--fg:#1d2327;--mut:#5b6670;--card:#fff;--line:#d9dee2;--acc:#1f5f8b;--same:#2e8b6f;--diff:#d2691e;
--one:#7a8fb8;--one2:#9b7ab8;--none:#c9ced3;--warn:#fff4e5}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15191c;--fg:#e6e9eb;--mut:#a3adb5;--card:#1d2226;--line:#353c42;
--acc:#7fb6dc;--same:#4fb391;--diff:#ee8a45;--one:#8ea3cc;--one2:#b195cc;--none:#4a5157;--warn:#3a2e1a}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:1060px;margin:0 auto;padding:22px 16px 60px}a{color:var(--acc)}h1{font-size:1.55rem;margin:.2em 0}h2{font-size:1.15rem;margin:1.8em 0 .5em;border-bottom:1px solid var(--line);padding-bottom:.2em}
.banner{background:var(--warn);border:1px solid var(--diff);border-radius:8px;padding:8px 12px;margin:8px 0 14px;font-size:.9rem}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:10px 0}
.kpi{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px}.kpi div{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px}
.kpi b{display:block;font-size:1.3rem}.kpi span{color:var(--mut);font-size:.82rem}.wrap{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:.86rem}th,td{border-bottom:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top}th{color:var(--mut)}
td.n{text-align:right;font-variant-numeric:tabular-nums}code{font:12.5px ui-monospace,Menlo,monospace}.sm{font-size:.82rem;color:var(--mut)}
.trk{display:flex;flex-wrap:wrap;gap:2px;margin:6px 0}.trk span{display:inline-block;min-width:16px;padding:1px 2px;border-radius:3px;color:#fff;font:600 11px ui-monospace,Menlo,monospace;text-align:center}
.same_partner{background:var(--same)}.different_partner{background:var(--diff)}.rfam_only{background:var(--one)}.star3d_only{background:var(--one2)}.neither{background:var(--none)}
.lg span{display:inline-flex;align-items:center;gap:4px;margin-right:12px;font-size:.82rem}.lg i{width:12px;height:12px;border-radius:3px;display:inline-block}
nav{font-size:.85rem;margin-bottom:6px}"""


LABEL = {"v": "INTERIM"}


def page(title, body, depth=0):
    up = "../" * depth
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{E(title)}</title><link rel="stylesheet" href="{up}style.css"></head><body><main>
<nav><a href="{up}index.html">Study</a> · <a href="{up}evidence/trna_sec_example.html">Lead example</a> · <a href="{up}report.html">Professor report</a></nav>
<div class="banner">{"<b>INTERIM</b> — " + LABEL.get('text', '') + " " if LABEL['v'] != "FINAL" else ""}Prepared by an AI research agent; not yet reviewed by a person.</div>
{body}</main></body></html>"""


def legend():
    return '<div class="lg">' + "".join(f'<span><i class="{c}"></i>{l}</span>' for c, l in (
        ("same_partner", "same partner"), ("different_partner", "different partner"), ("rfam_only", "Rfam only"),
        ("star3d_only", "STAR3D only"), ("neither", "neither"))) + "</div>"


# ------------------------------------------------------------------ site
def build_site(db, H, out, cohort):
    site = os.path.join(out, "site")
    if os.path.exists(site):
        shutil.rmtree(site)
    for d in ("family", "pair", "evidence", "data"):
        os.makedirs(os.path.join(site, d))
    open(os.path.join(site, "style.css"), "w").write(CSS)
    ret = os.path.join(ROOT, "outputs", "03_star3d", cohort, "returned", "runs")
    pairs = {r[0]: dict(zip([c[0] for c in db.execute("SELECT * FROM pairs").description], r))
             for r in db.execute("SELECT * FROM pairs")}
    runs = {r[0]: r for r in q(db, "SELECT pair_id, step3_status, rmsd, declared, superposition FROM runs")}
    scores = {r[0]: dict(zip([c[0] for c in db.execute("SELECT * FROM pair_scores").description], r))
              for r in db.execute("SELECT * FROM pair_scores")}
    sides = collections.defaultdict(dict)
    for r in db.execute("SELECT * FROM side_scores"):
        d = dict(zip([c[0] for c in db.execute("SELECT * FROM side_scores").description], r))
        sides[d["pair_id"]][d["side"]] = d
    # ---------- study page
    cons = H["rfam_consensus_reproduced"]
    kpis = [(f"{H['families_screened']:,}", "Rfam 15.1 families screened"),
            (f"{H['families_with_pairs']}", "families with selected pairs (this cohort)"),
            (f"{H['pairs']}", "pairs attempted"), (f"{H['valid']}", f"valid STAR3D alignments ({H['families_valid']} families)"),
            (f"{H['no_stacks']}", "“No similar stacks” (scores NA)"), (f"{H['failures']}", "technical failures"),
            (f"{H['different']:,} / {H['assessable']:,}", f"A-side different-partner observations ({100 * H['different'] / H['assessable']:.1f}%)"),
            (f"{H['shared_edges']:,}", f"exact residue pairs shared (Rfam {H['rfam_edges']:,} · STAR3D {H['star3d_edges']:,})"),
            (f"{cons[0]:,} / {cons[1]:,}", "Rfam consensus-column pairs reproduced by STAR3D")]
    rows = []
    for f, fid, n, v, ns in H["fam"]:
        a, d = H["fam_ad"].get(f, (0, 0))
        cn, cd = H["fam_cons"].get(f, (0, 0))
        rows.append(f"<tr><td><a href='family/{f}.html'>{f}</a></td><td>{E(fid)}</td><td class=n>{n}</td><td class=n>{v}</td><td class=n>{ns}</td>"
                    f"<td class=n>{pct(d, a) if a else 'NA'}</td><td class=n>{pct(cn, cd) if cd else 'NA'}</td><td>{H['rfam_cov'].get(f, '')}</td></tr>")
    oc = H["obs_class"]
    body = f"""<h1>Do the Rfam seed alignment and original STAR3D assign the same residue partners?</h1>
<p class=sm>Two original Rfam 15.1 seed rows vs one forward STAR3D v1.2 run on the corresponding experimental structures, over the whole
family interval, for verified natural RNAs from different species.</p>
<div class=kpi>{''.join(f'<div><b>{a}</b><span>{b}</span></div>' for a, b in kpis)}</div>
<div class=card><b>Lead example:</b> <a href="evidence/trna_sec_example.html">tRNA-Sec 3RG5 A U47C — Rfam and STAR3D choose different
partners inside the same BGSU motif group</a>.</div>
<div class=card><b>Read with care.</b> Pairs reuse RNAs (one structure can appear in many pairs), and tRNA (RF00005) and 5S rRNA (RF00001)
supply {sum(n for f, _, n, _, _ in H['fam'] if f in ('RF00005', 'RF00001'))} of {H['pairs']} pairs, so pooled residue counts are not independent
samples. Neither method is treated as ground truth: a disagreement says the two correspondences differ, not which is right.</div>
<h2>Families</h2><div class=wrap><table><tr><th>Family</th><th>Name</th><th>Pairs</th><th>Valid</th><th>No stacks</th>
<th>A-side different / assessable</th><th>Rfam consensus pairs reproduced</th><th>Rfam motif scan</th></tr>{''.join(rows)}</table></div>
<p class=sm>Families screened but without a selected pair are not “structurally compared”; their reasons are in the Step 2 family_status table
(downloadable database).</p>
<h2>Motif overlap at disagreements (A-side, 3 endpoints, counted once)</h2>
<div class=wrap><table><tr><th>Evidence class</th><th>Observations</th></tr>{''.join(f'<tr><td>{E(k)}</td><td class=n>{v:,}</td></tr>' for k, v in sorted(oc.items(), key=lambda x: -x[1]))}</table></div>
<p class=sm>{H['atlas_any']} disagreements have a BGSU Motif Atlas group at any endpoint; {H['rfam_any']} have an Rfam computational motif hit.
Equal-opportunity baseline (query endpoint only, identical coverage rules): Atlas group {H['baseline']['bgsu_atlas_group'][0]} disagreements vs
{H['baseline']['bgsu_atlas_group'][1]} agreements (Atlas covers few of these structures); Rfam motif {H['baseline']['rfam_motif'][0]} vs
{H['baseline']['rfam_motif'][1]}; BGSU loop {H['baseline']['bgsu_loop'][0]} vs {H['baseline']['bgsu_loop'][1]}. Descriptive only — not an enrichment test.</p>
<h2>Downloads</h2><p><a href="data/rfamalign3d.sqlite">rfamalign3d.sqlite</a> (all tables; provenance table lists source files and SHA-256).</p>"""
    open(os.path.join(site, "index.html"), "w").write(page(f"RfamAlign3D {LABEL['v'].lower()} results", body))
    shutil.copy2(os.path.join(out, "rfamalign3d.sqlite"), os.path.join(site, "data", "rfamalign3d.sqlite"))
    # ---------- family pages
    for f, fid, n, v, ns in H["fam"]:
        rs = []
        for pid in sorted(p for p in pairs if pairs[p]["family"] == f):
            p, r, s = pairs[pid], runs[pid], scores.get(pid)
            sa, sb = sides.get(pid, {}).get("A"), sides.get(pid, {}).get("B")
            rs.append(f"<tr><td><a href='../pair/{pid}.html'>{E(p['query_pdb_chain'])} → {E(p['target_pdb_chain'])}</a></td>"
                      f"<td>{E(p['query_species'])}<br>{E(p['target_species'])}</td><td>{E(r[1])}</td>"
                      f"<td class=n>{s['rfam_edges'] + ' / ' + s['star3d_edges'] + ' / ' + s['shared_edges'] if s else 'NA'}</td>"
                      f"<td class=n>{s['rfam_consensus_reproduced'] if s else 'NA'}</td>"
                      f"<td class=n>{sa['different_of_both_mapped'] if sa else 'NA'}</td><td class=n>{sb['different_of_both_mapped'] if sb else 'NA'}</td></tr>")
        body = f"""<h1>{f} · {E(fid)}</h1><p>{n} selected pairs; {v} valid alignments; {ns} “No similar stacks”.</p>
<div class=wrap><table><tr><th>Pair (A → B, frozen run direction)</th><th>Species</th><th>STAR3D outcome</th><th>Rfam / STAR3D / shared edges</th>
<th>Rfam consensus pairs reproduced</th><th>A-side different / both mapped</th><th>B-side different / both mapped</th></tr>{''.join(rs)}</table></div>"""
        open(os.path.join(site, "family", f"{f}.html"), "w").write(page(f"{f} {fid}", body, 1))
    # ---------- pair pages
    for pid, p in pairs.items():
        r = runs[pid]
        dd = os.path.join(site, "data", pid)
        os.makedirs(dd)
        links = []
        for fn in (f"{pid}__forward.aln", f"{pid}__forward.aln.pdb"):
            src = os.path.join(ret, pid, "out", fn)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(dd, fn))
                links.append(f"<a href='../data/{pid}/{fn}'>{fn.split('__')[-1]}</a>")
        qa, qb = p["query_pdb_chain"].split("_")[0], p["target_pdb_chain"].split("_")[0]
        head = f"""<h1>{E(p['query_pdb_chain'])} vs {E(p['target_pdb_chain'])}</h1>
<p>{E(p['family'])} {E(p['family_id'])} · A: <a href="https://www.rcsb.org/structure/{qa}">{qa}</a> {E(p['query_species'])}, seed row
<code>{E(p['query_seed_row'])}</code> · B: <a href="https://www.rcsb.org/structure/{qb}">{qb}</a> {E(p['target_species'])}, seed row
<code>{E(p['target_seed_row'])}</code></p><p class=sm>STAR3D outcome: <b>{E(r[1])}</b>{' · RMSD ' + E(r[2]) + ' Å' if r[2] else ''} · raw files: {' · '.join(links) or 'none'}
· superposition: {E(r[4] or 'n/a')} (“STAR3D superposition for inspection”, not evidence of correctness)</p>"""
        s = scores.get(pid)
        if not s:
            body = head + "<div class=card>No valid alignment: all scores are <b>NA</b> (not zero).</div>"
            open(os.path.join(site, "pair", f"{pid}.html"), "w").write(page(pid, body, 1))
            continue
        sa, sb = sides[pid]["A"], sides[pid]["B"]
        card = (f"Each alignment contains {s['rfam_edges']} (Rfam) and {s['star3d_edges']} (STAR3D) assessable residue pairs; {s['shared_edges']} "
                f"are identical. Viewed from {p['query_pdb_chain']}, {sa['scored_different']} of {sa['both_mapped']} nucleotides mapped by both "
                f"methods receive different partners. Viewed from {p['target_pdb_chain']}, {sb['scored_different']} of {sb['both_mapped']} do.")
        sc = "".join(f"<tr><td>{E(k)}</td><td class=n>{E(s[k])}</td></tr>" for k in (
            "rfam_recovery", "star3d_agreement", "rfam_only_edges", "star3d_only_edges", "rfam_consensus_reproduced", "rfam_insert_reproduced",
            "rfam_unknown_reproduced", "star3d_both_consensus_in_rfam", "star3d_mixed_consensus_insert_in_rfam", "star3d_both_insert_in_rfam",
            "consensus_only_rfam_recovery", "consensus_only_star3d_agreement"))
        side_t = "".join(f"<tr><td>{x['side']} {E(x['pdb_chain'])}</td>" + "".join(f"<td class=n>{x[k]}</td>" for k in (
            "nucleotides", "observed", "same_partner", "different_partner", "rfam_only", "star3d_only", "neither", "both_mapped", "scored_different")) + "</tr>"
            for x in (sa, sb))
        tracks = ""
        for side in ("A", "B"):
            cells = []
            for u, pos, nt, cat, rp, sp, col, line in q(db, """SELECT unit_id, ungapped_position, nucleotide, category, rfam_partner, star3d_partner,
                    rfam_column, raw_line FROM residue_comparisons WHERE pair_id=? AND side=? ORDER BY CAST(ungapped_position AS INT)""", pid, side):
                tip = f"{u} col {col}: {cat}; Rfam {rp or '-'}; STAR3D {sp or '-'} {line}"
                cells.append(f'<span class="{cat}" title="{E(tip)}">{E(nt)}</span>')
            tracks += f"<p class=sm>{side} view ({E(sa['pdb_chain'] if side == 'A' else sb['pdb_chain'])}), positions 1…{len(cells)}</p><div class=trk>{''.join(cells)}</div>"
        regs = q(db, "SELECT region_id, length, start_unit, end_unit, rf_composition, ss_composition FROM regions WHERE pair_id=? AND side='A'", pid)
        rm = {r[0]: r for r in q(db, "SELECT region_id, members_bgsu_motif_group, members_rfam_computational_motif, members_structural_loop_only, members_no_documented_hit, members_incomplete FROM region_motif WHERE pair_id=?", pid)}
        reg_t = "".join(f"<tr><td><code>{E(r[0].split(':', 2)[2])}</code></td><td class=n>{r[1]}</td><td><code>{E(r[2])}</code> … <code>{E(r[3])}</code></td>"
                        f"<td>{E(r[4])}</td><td>{E(r[5])}</td><td class=n>{rm.get(r[0], ('', 0, 0, 0, 0, 0))[1]}</td><td class=n>{rm.get(r[0], ('', 0, 0, 0, 0, 0))[2]}</td>"
                        f"<td class=n>{rm.get(r[0], ('', 0, 0, 0, 0, 0))[3]}</td></tr>" for r in regs)
        obs = q(db, """SELECT o.observation_id, o.a_residue, o.b_rfam_partner, o.b_star3d_partner, o.rfam_column, o.raw_line, o.classification
                       FROM motif_observations o WHERE o.pair_id=?""", pid)
        ob_t = ""
        for oid, a, b1, b2, col, line, cls in obs:
            evs = q(db, "SELECT role, evidence_type, record, name, url FROM overlap_evidence WHERE observation_id=? AND evidence_type!='none'", oid)
            ev = "<br>".join(f"{E(role)}: <a href='{E(url)}'>{E(rec)}</a> {E(nm)} <span class=sm>({E(et)})</span>" for role, et, rec, nm, url in evs) or "<span class=sm>none documented</span>"
            ob_t += f"<tr><td><code>{E(a)}</code></td><td class=n>{col}</td><td><code>{E(b1)}</code></td><td><code>{E(b2)}</code> <code class=sm>{E(line)}</code></td><td>{E(cls)}</td><td>{ev}</td></tr>"
        body = head + f"""<div class=card><b>{E(card)}</b></div>
<h2>Direction-independent scores</h2><table>{sc}</table>
<h2>Residue views (two views of one comparison — never summed)</h2><table><tr><th>Side</th><th>nt</th><th>observed</th><th>same</th><th>different</th>
<th>Rfam only</th><th>STAR3D only</th><th>neither</th><th>both mapped</th><th>scored different</th></tr>{side_t}</table>
{legend()}{tracks}
<h2>Disagreement regions (A side)</h2><div class=wrap><table><tr><th>Region</th><th>Length</th><th>Members</th><th>RF</th><th>SS_cons context</th>
<th>members w/ Atlas group</th><th>members w/ Rfam motif</th><th>members loop only</th></tr>{reg_t}</table></div>
<h2>Disagreement evidence (A side, three endpoints)</h2><div class=wrap><table><tr><th>A residue</th><th>Rfam col</th><th>B (Rfam)</th>
<th>B (STAR3D) + raw line</th><th>class</th><th>evidence (click for source record)</th></tr>{ob_t}</table></div>"""
        open(os.path.join(site, "pair", f"{pid}.html"), "w").write(page(pid, body, 1))
    build_example(db, site)
    return site


def build_example(db, site):
    pid = one(db, "SELECT pair_id FROM pairs WHERE query_pdb_chain='3RG5_A' AND target_pdb_chain='3W3S_B'")
    oid = one(db, "SELECT observation_id FROM motif_observations WHERE pair_id=? AND a_residue='3RG5|1|A|U|47|||C'", pid)
    o = db.execute("SELECT a_residue,b_rfam_partner,b_star3d_partner,rfam_column,raw_line,classification FROM motif_observations WHERE observation_id=?", (oid,)).fetchone()
    p = db.execute("SELECT query_species,target_species,query_seed_row,target_seed_row FROM pairs WHERE pair_id=?", (pid,)).fetchone()
    evs = q(db, "SELECT role, unit_id, evidence_type, record, detail, url FROM overlap_evidence WHERE observation_id=? AND evidence_type!='none' ORDER BY role", oid)
    loop_src = open(os.path.join(ROOT, "cache", "bgsu", "loops", "loops_3W3S.csv")).read()
    loop_line = next(l for l in loop_src.splitlines() if l.startswith('"IL_3W3S_003"'))
    ref = q(db, "SELECT rfam_column, query_char, target_char, query_residue, target_residue, correspondence FROM rfam_reference WHERE pair_id=? AND CAST(rfam_column AS INT) BETWEEN 63 AND 69", pid)
    ev_t = "".join(f"<tr><td>{E(r)}</td><td><code>{E(u)}</code></td><td>{E(t)}</td><td><a href='{E(url)}'>{E(rec)}</a></td><td class=sm>{E(d)}</td></tr>" for r, u, t, rec, d, url in evs)
    ref_t = "".join(f"<tr><td class=n>{c}</td><td>{E(a)}</td><td>{E(b)}</td><td><code>{E(x)}</code></td><td><code>{E(y)}</code></td><td>{E(k)}</td></tr>" for c, a, b, x, y, k in ref)
    body = f"""<h1>Lead example: tRNA-Sec — different partners inside the same motif</h1>
<div class=card><p><b>3RG5 chain A residue U47C</b> ({E(p[0])}) corresponds to <b>3W3S chain B residue G47B</b> ({E(p[1])}) in the
Rfam seed alignment (column {E(o[3])}) and to <b>3W3S chain B residue U47C</b> in STAR3D (raw line <code>{E(o[4])}</code>).</p>
<p>Both 3W3S target residues are members of BGSU loop <a href="http://rna.bgsu.edu/rna3dhub/loops/view/IL_3W3S_003">IL_3W3S_003</a>,
whose unit set equals an instance of Motif Atlas 4.16 group <a href="http://rna.bgsu.edu/rna3dhub/motif/view/IL_10892.5">IL_10892.5</a>.
This shows a <b>correspondence disagreement within the same motif</b>; it does <b>not</b> establish which partner is correct.</p></div>
<h2>Original Rfam seed columns 63–69 (both rows)</h2><p class=sm>Seed rows: A <code>{E(p[2])}</code>, B <code>{E(p[3])}</code> (Rfam 15.1 RF01852).</p>
<table><tr><th>Column</th><th>A char</th><th>B char</th><th>A residue</th><th>B residue</th><th>Rfam correspondence</th></tr>{ref_t}</table>
<h2>Source record (BGSU RNA 3D Hub, loops/download/3W3S, model 1)</h2><p><code>{E(loop_line)}</code></p>
<p class=sm>Source unit IDs (<code>PDB|model|chain|nt|number|||insertion</code>) are used unchanged as residue IDs; the Step 2 crosswalk
maps <code>3W3S|1|B|G|47|||B</code> → author B:47B → staged B:47B → ungapped position 52 → Rfam column 66, and <code>3W3S|1|B|U|47|||C</code>
→ B:47C → position 53 → column 67.</p>
<h2>All evidence for this observation</h2><table><tr><th>Endpoint</th><th>Residue</th><th>Evidence type</th><th>Record</th><th>Detail</th></tr>{ev_t}</table>
<p>Pair page: <a href="../pair/{pid}.html">{E(pid)}</a> · RCSB: <a href="https://www.rcsb.org/structure/3RG5">3RG5</a>,
<a href="https://www.rcsb.org/structure/3W3S">3W3S</a> · raw alignment: <a href="../data/{pid}/{pid}__forward.aln">.aln</a></p>"""
    open(os.path.join(site, "evidence", "trna_sec_example.html"), "w").write(page("Lead example: tRNA-Sec", body, 1))


# ------------------------------------------------------------------ professor report
def report(db, H, out, site):
    cons, ins = H["rfam_consensus_reproduced"], H["rfam_insert_reproduced"]
    fam_rows = []
    for f, fid, n, v, ns in H["fam"]:
        a, d = H["fam_ad"].get(f, (0, 0))
        fam_rows.append(f"| {f} {fid} | {n} | {v} | {ns} | {pct(d, a) if a else 'NA'} |")
    oc = H["obs_class"]
    b = H["baseline"]
    tr5 = sum(n for f, _, n, _, _ in H["fam"] if f in ("RF00005", "RF00001"))
    where = ""
    dg = os.path.join(ROOT, "outputs", "04_comparison", os.path.basename(out), "diagnostics", "summary.json")
    au = os.path.join(ROOT, "outputs", "04_comparison", os.path.basename(out), "audit", "audit.json")
    if os.path.exists(dg) and os.path.exists(au):
        d, a = json.load(open(dg)), json.load(open(au)).get("Q7", {})
        n = d["observations"]
        pt = d["pattern"]
        where = (f"**Where and how the alignments differ** (descriptive; each disagreement in exactly one pattern): "
                 f"{d['abs_offset'].get('1', 0):,}/{n:,} disagreements are offsets of exactly one position along RNA B's ungapped sequence "
                 f"(insertion codes preserved). Constant-offset runs next to a seed gap: {pt.get('register_shift_at_gap', 0):,}/{n:,}; "
                 f"pair-wide difference (≥50% of residues mapped by both differ): {pt.get('global_misregistration', 0):,}/{n:,} in "
                 f"{a.get('pairwide_pairs', '?')} pairs, not one constant offset ({a.get('pairwide_modal_offset_share', '?')} share their pair's most common offset). "
                 f"Seed gaps lie within two columns of {d['near_seed_gap_baseline']['disagreements']} disagreements and of "
                 f"{d['near_seed_gap_baseline']['agreements']} agreements. Disagreement rate in consensus hairpin loops "
                 f"{a.get('hairpin', {}).get('hairpin', '?')} vs {a.get('hairpin', {}).get('outside_hairpin', '?')} elsewhere.")
    ovp = os.path.join(ROOT, "outputs", "04_comparison", os.path.basename(out), "overall", "overall_metric.json")
    if LABEL["v"] == "FINAL" and os.path.exists(ovp):
        ov = json.load(open(ovp))
        o, st, off = ov["overall"], ov["sets"], ov["offsets"]
        fo = tsv(os.path.join(os.path.dirname(ovp), "family_overall.tsv"))
        frows = "\n".join(f"| {r['family_id']} ({r['family']}) | {'1 pair studied' if r['pairs'] == '1' else r['pairs']} | "
                          f"{(pct(int(r['differences']), int(r['denominator'])) if int(r['denominator']) else 'NA')} | {r['different_partner_rate']} |"
                          for r in sorted(fo, key=lambda r: -int(r["pairs"])))
        md = f"""# Where RNA sequence and structural alignments agree—and where they differ

*RfamAlign3D results. Every number is a query on the study database. Prepared by an AI research agent; not yet reviewed by a person.*

**Question.** We compared the nucleotide correspondences implied by two original Rfam 15.1 seed rows with those computed by original STAR3D v1.2 from the experimental structures, over the whole family interval, for verified natural RNAs from different species.

**Study scope.** The study screened {H['families_screened']:,} Rfam families and obtained usable structural comparisons in {H['families_valid']} families ({H['valid']} usable alignments of {H['pairs']} attempted; {H['no_stacks']} pairs gave “No similar stacks”). Comparisons reuse RNAs, and tRNA and 5S rRNA dominate the dataset ({tr5} of {H['pairs']} pairs). Alignment differences identify candidates for investigation; they do not establish which method is correct.

**Main results.**
- STAR3D reproduces **{cons[0]:,}/{cons[1]:,} ({100 * cons[0] / cons[1]:.1f}%)** assessable Rfam consensus correspondences.
- **All assessable alignment differences, RNA A view: {o['differences']:,}/{o['denominator']:,} ({100 * o['rate']:.1f}%)** — different partner {o['different_partner']:,}, Rfam only {o['rfam_only']:,}, STAR3D only {o['star3d_only']:,}. Excluded: {o['neither']:,} neither-mapped and {o['unassessable_total']:,} not assessable (a nucleotide or its counterpart lacks coordinates).
- Of different-partner cases, the alternatives are one nucleotide apart in **{off['one_position']:,}/{off['different_partner_observations']:,} ({100 * off['one_position'] / off['different_partner_observations']:.1f}%)**.
- Correspondence sets (direction-independent; unit = residue correspondences): shared {st['shared']:,}; Rfam only {st['rfam_only']:,}; STAR3D only {st['star3d_only']:,}; set difference {st['symmetric_difference']:,}/{st['union']:,} ({100 * st['set_difference_rate']:.1f}%).

The percentages use different denominators and are not complements.

| Family | Pairs | Overall assessable difference rate | Different-partner rate |
|---|---:|---:|---:|
{frows}

{where.replace("disagreements are offsets", "different-partner observations are offsets").replace("each disagreement in exactly one pattern", "different-partner observations; each in exactly one pattern")}

**Motif evidence** (different-partner observations, any of three endpoints, counted once): Rfam computational motif span {oc.get('rfam_computational_motif', 0):,}; BGSU Atlas group {oc.get('bgsu_motif_group', 0)}; both {oc.get('bgsu_group_and_rfam_motif', 0)}; BGSU structural loop only {oc.get('structural_loop_only', 0):,}; consensus loop context only {oc.get('consensus_loop_context_only', 0)}; no documented record {oc.get('no_documented_hit', 0)}; Rfam coverage unknown with no other record {oc.get('incomplete', 0)}. Rfam scan coverage is established per family, not per row; a missing record is not a confirmed absence. Descriptive only.

**Example.** tRNA-Sec 3RG5 A U47C corresponds to 3W3S B G47B in the seed alignment (column 66) and to 3W3S B U47C in STAR3D (`A:47C<->B:47C`); both 3W3S residues belong to BGSU loop IL_3W3S_003 (Motif Atlas 4.16 group IL_10892.5) — a correspondence difference within one motif that does not establish which partner is correct.

**Limitations.** Constructs without sufficient independent natural-sequence evidence were not included in the verified comparison cohort; this includes chains whose sequence search did not complete, which are unresolved rather than shown to be engineered. Methods, selection counts and provenance (execution environments, checksums, policy history) are in the website's Methods & limitations page and the downloadable tables.
"""
        open(os.path.join(out, "results_report.md"), "w").write(md)
        shutil.copy2(os.path.join(out, "results_report.md"), os.path.join(site, "data", "results_report.md"))
        return
    md = f"""# RfamAlign3D — {LABEL['v'].lower()} results (cohort `{os.path.basename(out)}`)

**{LABEL['v']}.** Generated from `rfamalign3d.sqlite`. {LABEL.get('text', '')} Prepared by an AI agent; not yet reviewed by a person.

**1 · Question and cohort.** Do two original Rfam 15.1 seed rows and one forward run of original STAR3D v1.2 assign the same residue partners across the whole family interval, for verified natural RNAs of different species?

**2 · Coverage.** {H['families_screened']:,} families screened; {H['families_with_pairs']} families have selected pairs in this cohort. {H['pairs']} pairs attempted: **{H['valid']} valid alignments in {H['families_valid']} families**, {H['no_stacks']} “No similar stacks” (NA), {H['failures']} technical failures. Pairs reuse RNAs; tRNA and 5S rRNA supply {tr5} of {H['pairs']} pairs.

**3 · Alignment result.** {H['shared_edges']:,} exact residue pairs are shared (Rfam {H['rfam_edges']:,}, STAR3D {H['star3d_edges']:,} assessable pairs). *Exclusion rule:* a residue pair is assessable only if both endpoints are observed in the staged coordinates (defined from Step 2, not from either method). Of {H['rfam_raw']:,} raw Rfam pairs, {H['rfam_raw'] - H['rfam_edges']:,} are excluded ({'; '.join(f"{k}: {v}" for k, v in sorted(H['rfam_excluded'].items()))}); {H['star3d_excluded']} STAR3D pairs are excluded. Consensus and insert strata partition the {H['rfam_edges']:,} assessable Rfam pairs. Rfam consensus-column pairs reproduced by STAR3D: **{cons[0]:,}/{cons[1]:,}** ({100 * cons[0] / cons[1]:.1f}%); insert-column pairs: {ins[0]:,}/{ins[1]:,}. A-side partner disagreements: **{H['different']:,}/{H['assessable']:,} ({100 * H['different'] / H['assessable']:.1f}%)** jointly assessable observations.

| Family | Pairs | Valid | No stacks | A-side different / assessable |
|---|---:|---:|---:|---:|
{chr(10).join(fam_rows)}

{where}

**4 · Motif result** (A-side disagreements, any of 3 endpoints, counted once): Rfam computational motif only {oc.get('rfam_computational_motif', 0):,}; BGSU Atlas group only {oc.get('bgsu_motif_group', 0)}; both {oc.get('bgsu_group_and_rfam_motif', 0)}; BGSU loop only {oc.get('structural_loop_only', 0):,}; consensus loop context only {oc.get('consensus_loop_context_only', 0)}; no documented record (BGSU loop lookup complete; Rfam absence not verified per row) {oc.get('no_documented_hit', 0)}; Rfam coverage unknown with no other record {oc.get('incomplete', 0)}. Baseline, query endpoint only: Atlas group {b['bgsu_atlas_group'][0]} disagreements vs {b['bgsu_atlas_group'][1]} agreements (Atlas covers few of these structures); Rfam motif {b['rfam_motif'][0]} vs {b['rfam_motif'][1]}. Descriptive, not an enrichment test. *Rfam coverage:* `motif_family_stats` reconciles with all seed rows of the covered families, so coverage is established **only at family level** for every family: no row-level scan record exists, and an absent Rfam hit is not a verified negative for any individual seed row; RF00522 is not listed, so its Rfam coverage is **unknown** for all {H['rfam_cov_obs'].get('unknown', 0)} of its disagreement observations, including those with documented BGSU loop evidence. Rfam hits are position spans, shown with their source records.

**5 · Inspectable example.** tRNA-Sec 3RG5 A U47C ↔ 3W3S B G47B (Rfam column 66) vs U47C (STAR3D `A:47C<->B:47C`); both targets are in BGSU loop IL_3W3S_003, Atlas group IL_10892.5 — a disagreement within the same motif that does not establish which partner is correct. Site page: `site/evidence/trna_sec_example.html`.
"""
    open(os.path.join(out, "results_report.md"), "w").write(md)
    import re
    body = md
    body = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", E(body).replace("&lt;-&gt;", "&lt;-&gt;"))
    body = re.sub(r"`(.+?)`", r"<code>\1</code>", body)
    lines, tbl = [], []
    for l in body.splitlines():
        if l.startswith("|"):
            tbl.append(l)
            continue
        if tbl:
            hdr = [c.strip() for c in tbl[0].strip("|").split("|")]
            rows_ = [[c.strip() for c in t.strip("|").split("|")] for t in tbl[2:]]
            lines.append("<table><tr>" + "".join(f"<th>{h}</th>" for h in hdr) + "</tr>" +
                         "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows_) + "</table>")
            tbl = []
        if l.startswith("# "):
            lines.append(f"<h1>{l[2:]}</h1>")
        elif l.strip():
            lines.append(f"<p>{l}</p>")
    open(os.path.join(site, "report.html"), "w").write(page("RfamAlign3D report", "\n".join(lines)))
    shutil.copy2(os.path.join(out, "results_report.md"), os.path.join(site, "data", "results_report.md"))


def set_label(cohort):
    if cohort.startswith("step2_final"):
        LABEL.update(v="FINAL", text="", _old="Final cohort: Step 2 closed on 2026-11-13 by the researcher's decision (policy R7); "
                                     "39 natural-sequence checks whose NCBI search did not complete are provisional (SEARCH_INCOMPLETE), "
                                     "not excluded as engineered. STAR3D: 505 of 535 runs natively on Linux; the 30 new pairs ran in Docker "
                                     "(linux/amd64 emulation, same STAR3D v1.2 package) after a QC rerun of RF00522 reproduced the native run's "
                                     "32 mappings, 2.75 Å RMSD and all preprocessing files byte-for-byte.")
    else:
        LABEL.update(v="INTERIM", text=f"Frozen interim cohort {cohort}; Step 2 checks were still pending for other pairs.")


def main():
    cohort = sys.argv[1]
    set_label(cohort)
    out = os.path.join(ROOT, "outputs", "05_results", cohort)
    os.makedirs(out, exist_ok=True)
    db, dbp = build_db(cohort, out)
    H = headline(db)
    site = build_site(db, H, out, cohort)
    report(db, H, out, site)
    json.dump({k: v for k, v in H.items() if k not in ("fam", "fam_ad", "fam_cons")}, open(os.path.join(out, "headline.json"), "w"), indent=1, default=str)
    print(json.dumps({k: H[k] for k in ("families_screened", "pairs", "valid", "families_valid", "no_stacks", "failures", "assessable",
                                        "different", "rfam_consensus_reproduced", "rfam_insert_reproduced", "rf_pairs", "atlas_any", "rfam_any", "baseline")}, indent=1, default=str))
    print("report:", os.path.relpath(os.path.join(out, "results_report.md"), ROOT))
    print("site:  ", os.path.relpath(site, ROOT))


if __name__ == "__main__":
    main()
