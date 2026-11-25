"""Step 2 / stage 6: Family_status for all families, Structures table, workbook (INTERIM until every family has an
outcome), acceptance gates (validation_report.md / validation.json) and professor_summary.md.
Reads only Step 2 outputs; writes only under outputs/02_structure_preparation/."""
import collections
import glob
import gzip
import json
import os
import subprocess
import sys

import gemmi
import openpyxl
from openpyxl.styles import Font

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

OK_DEC = "verified_natural"
TECH = {"ready", "ready_with_coverage_flags"}


def family_status():
    disc = {r["family"]: r for r in C.read_tsv(C.P("family_discovery.tsv"))}
    yld = {r["family"]: r for r in C.read_tsv(C.P("yield_by_family.tsv"))}
    chains = C.read_tsv(C.P("candidate_chains.tsv"))
    ents = {r["pdb_id"]: r for r in C.read_tsv(C.P("entry_metadata.tsv"))}
    wild = C.read_tsv(C.P("wildcard_rescreen.tsv")) if os.path.exists(C.P("wildcard_rescreen.tsv")) else []
    pre = C.jload(C.P("prefetch_status.json"), {})
    state = C.jload(C.P("selection_state.json"), {}).get("families", {})
    cr = C.read_tsv(C.P("chain_review.tsv"))
    pairs = C.read_tsv(C.P("pairs_to_run.tsv"))
    by_cr = collections.defaultdict(list)
    for r in cr:
        by_cr[r["family"]].append(r)
    by_ch = collections.defaultdict(list)
    for r in chains:
        by_ch[r["family"]].append(r)
    wild_fam = collections.Counter(r["family"] for r in wild)
    npairs = collections.Counter(p["family"] for p in pairs)
    out = []
    for fam, d in sorted(disc.items()):
        cs = by_ch.get(fam, [])
        obsolete = sorted({r["pdb_id"] for r in cs if ents.get(r["pdb_id"], {}).get("status") == "obsolete"})
        wl = [r for r in cs if r["screening"] == "wildcard_only_needs_mmcif"]
        wl_unchecked = [r for r in wl if pre.get(r["pdb_id"]) != "ok" and "deposited_unknown=-" not in r["diagnostic"]]
        other_incomplete = [r for r in cs if r["screening"] == "search_incomplete" and r["pdb_id"] not in obsolete]
        st = d["discovery_status"]
        if st == "search_incomplete":
            if other_incomplete or wl_unchecked:
                st = "search_incomplete"
            elif int(d["chains_exact_seed_match"]) or wild_fam.get(fam):
                st = "exact_seed_match_found"
            else:
                st = "candidate_found_no_exact_seed_match"
        y = yld.get(fam, {})
        s = state.get(fam, {})
        reasons = collections.Counter(r["reason"] or r["decision"] for r in by_cr.get(fam, [])
                                      if r["decision"] not in (OK_DEC,) and r["decision"])
        tech_bad = sum(1 for r in by_cr.get(fam, []) if r["decision"] == OK_DEC and r.get("technical_status") not in TECH)
        if st == "no_candidate_in_searched_sources":
            outcome = "no_candidate_in_searched_sources"
        elif st == "candidate_found_no_exact_seed_match":
            outcome = "candidate_found_no_exact_seed_match"
        elif st == "search_incomplete":
            outcome = "search_incomplete"
        elif s.get("seq_pairs_pending"):
            outcome = "in_progress_evidence_pending"
        elif npairs.get(fam):
            outcome = "eligible_pairs_selected"
        elif y.get("different_sequence_pair_possible") != "yes" and not wild_fam.get(fam):
            outcome = ("exact_match_single_sequence_no_pair_possible" if int(y.get("distinct_matched_sequences") or 0) <= 1
                       else "no_different_species_pair_possible")
        else:
            outcome = "pair_possible_no_eligible_pair"
        out.append(dict(family=fam, family_id=d["family_id"], discovery_status=st, outcome=outcome,
                        family_complete="no" if outcome in ("in_progress_evidence_pending", "search_incomplete") else "yes",
                        evidence_records=d["evidence_records"], candidate_entries=d["candidate_entries"],
                        candidate_chains=d["candidate_chains"], exact_unambiguous_intervals=y.get("exact_unambiguous_intervals", 0),
                        distinct_matched_sequences=y.get("distinct_matched_sequences", 0),
                        screening_seq_pairs_diff_species=y.get("seq_pairs_known_different_species", 0),
                        screening_seq_pairs_species_unresolved=y.get("seq_pairs_species_unresolved", 0),
                        chains_reviewed=sum(1 for r in by_cr.get(fam, []) if r.get("technical_status")),
                        verified_natural=sum(1 for r in by_cr.get(fam, []) if r["decision"] == OK_DEC),
                        verified_but_technical_unusable=tech_bad,
                        pairs_to_run=npairs.get(fam, 0), seq_pairs_pending=s.get("seq_pairs_pending", 0),
                        principal_reasons=";".join(f"{k}:{v}" for k, v in reasons.most_common(4)),
                        obsolete_entries_not_substituted=";".join(obsolete)))
    C.write_tsv(C.P("family_status.tsv"), out)
    return out


def structures():
    cr = C.read_tsv(C.P("chain_review.tsv"))
    ids = sorted({r["pdb_id"] for r in cr if r.get("technical_status")})
    meta = {r["pdb_id"]: r for r in C.read_tsv(C.P("entry_metadata.tsv"))}
    out = []
    import s04_coordinates as X
    for p in ids:
        d = X.entry_data(p)["provenance"]
        out.append(dict(pdb_id=p, rcsb_page=d["rcsb_page"], download_url=d["download_url"], mmcif=d["mmcif"],
                        mmcif_sha256=d["mmcif_sha256"], mmcif_bytes=d["mmcif_bytes"], method=d["method"],
                        resolution=d["resolution"], revision=d.get("revision", ""), revision_date=d.get("revision_date", ""),
                        deposit_date=d.get("deposit_date", ""), models=d.get("models", ""), title=d["title"],
                        status=meta.get(p, {}).get("status", "current"),
                        retrieved=next((r["retrieved_utc"] for r in C.read_tsv(C.MANIFEST) if r["key"] == f"mmcif:{p}"), "")))
    C.write_tsv(C.P("structures.tsv"), out)
    return out


SHEETS = [("Family_status", "family_status.tsv"), ("Structures", "structures.tsv"), ("Chain_review", "chain_review.tsv"),
          ("Pairs_to_run", "pairs_to_run.tsv"), ("Yield_by_family", "yield_by_family.tsv")]


def workbook(interim):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "About"
    for line in [("Status", "INTERIM - not every family has an outcome yet" if interim else "FINAL - every family has an outcome"),
                 ("Generated", C.now()), ("Policy", "step2_policy.yaml"), ("STAR3D run in this step", "no"),
                 ("Note", "TSV files and companion tables are authoritative; this workbook is a browsing copy. "
                          "Residue-level tables (residue_crosswalk.tsv.gz, rfam_pair_reference.tsv.gz) are not in Excel.")]:
        ws.append(list(line))
    for name, f in SHEETS:
        rows = C.read_tsv(C.P(f))
        w = wb.create_sheet(name)
        cols = list(rows[0]) if rows else ["empty"]
        w.append(cols)
        for c in w[1]:
            c.font = Font(bold=True)
        for r in rows:
            vals = [(r.get(c) or "")[:32000] for c in cols]
            w.append(vals)
            for cell in w[w.max_row]:
                cell.data_type = "s"
                cell.number_format = "@"
        w.freeze_panes = "A2"
        w.auto_filter.ref = w.dimensions
    path = C.P("SQARMA_Step_2_STAR3D_inputs.xlsx")
    wb.save(path)
    # reconcile
    wb2 = openpyxl.load_workbook(path, read_only=True)
    ok = True
    for name, f in SHEETS:
        rows = C.read_tsv(C.P(f))
        got = list(wb2[name].iter_rows(values_only=True))
        cols = list(rows[0]) if rows else ["empty"]
        ok &= list(got[0]) == cols and len(got) - 1 == len(rows)
        ok &= all([("" if v is None else v) for v in g] == [(r.get(c) or "")[:32000] for c in cols]
                  for g, r in zip(got[1:], rows))
    return ok


def gates(fs):
    G = collections.OrderedDict()
    rows = {}
    with gzip.open(os.path.join(C.STEP1, "seed_rows.jsonl.gz"), "rt") as f:
        for line in f:
            r = json.loads(line)
            rows[(r["family"], r["seed_row_id"])] = r
    # 1 families
    G["g1_all_families_one_status"] = len(fs) == 4227 and len({r["family"] for r in fs}) == 4227 and all(r["outcome"] for r in fs)
    G["g1_outcome_counts"] = dict(collections.Counter(r["outcome"] for r in fs))
    # 2 pairs
    pairs = C.read_tsv(C.P("pairs_to_run.tsv"))
    cr = {r["candidate_id"]: r for r in C.read_tsv(C.P("chain_review.tsv"))}
    p_ok, p_bad = True, []
    for p in pairs:
        q, t = cr[p["query_candidate"]], cr[p["target_candidate"]]
        rq, rt = rows.get((p["family"], p["query_seed_row"])), rows.get((p["family"], p["target_seed_row"]))
        cond = dict(rows_exist=bool(rq and rt),
                    distinct_strings=bool(rq and rt) and rq["ungapped_sequence"].upper() != rt["ungapped_sequence"].upper(),
                    species_differ=q["species_taxid"] and t["species_taxid"] and q["species_taxid"] != t["species_taxid"],
                    entries_differ=q["pdb_id"] != t["pdb_id"],
                    verified=q["decision"] == OK_DEC and t["decision"] == OK_DEC,
                    technical=q["technical_status"] in TECH and t["technical_status"] in TECH)
        if not all(cond.values()):
            p_ok = False
            p_bad.append((p["pair_id"], [k for k, v in cond.items() if not v]))
    G["g2_primary_jobs_meet_all_rules"] = p_ok
    G["g2_violations"] = p_bad[:10]
    # 3 staged files
    man = C.read_tsv(C.P("inputs_star3d", "file_manifest.tsv")) if os.path.exists(C.P("inputs_star3d", "file_manifest.tsv")) else []
    xw = C.read_tsv(C.P("residue_crosswalk.tsv.gz")) if os.path.exists(C.P("residue_crosswalk.tsv.gz")) else []
    xw_by = collections.defaultdict(list)
    for x in xw:
        xw_by[x["candidate_id"]].append(x)
    f_ok, f_bad = True, []
    for m in man:
        path = os.path.join(C.ROOT, m["file"])
        st = gemmi.read_structure(path)
        chains = [ch.name for ch in st[0]]
        staged_ids = [f"{r.seqid.num}{r.seqid.icode.strip()}" for ch in st[0] for r in ch]
        conv = C.read_tsv(os.path.join(C.P("inputs_star3d"), f"{m['stem']}.conversion.tsv"))
        lo, hi = map(int, m["label_interval"].split("-"))
        xmap = collections.Counter(x["staged_resid"] for x in xw_by[m["candidate_id"]] if x["staged_resid"])
        cond = dict(sha=C.sha256(path) == m["sha256"], chain=chains == [m["staged_chain"]],
                    reader=m["reader_ok"] == "True",
                    each_staged_residue_maps_once=all(xmap[s] == 1 for s in staged_ids) and len(staged_ids) == sum(xmap.values()),
                    inside_interval=all(lo <= int(c["label_seq_id"]) <= hi for c in conv))
        if not all(cond.values()):
            f_ok = False
            f_bad.append((m["stem"], [k for k, v in cond.items() if not v]))
    G["g3_staged_files_valid"] = f_ok
    G["g3_violations"] = f_bad[:10]
    # 4 reference round-trip
    ref = C.read_tsv(C.P("rfam_pair_reference.tsv.gz")) if os.path.exists(C.P("rfam_pair_reference.tsv.gz")) else []
    pmap = {p["pair_id"]: p for p in pairs}
    xi = {(x["candidate_id"], x["staged_resid"]): x for x in xw if x["staged_resid"]}
    rt_ok, unobs_ok = True, True
    for r in ref:
        p = pmap[r["pair_id"]]
        for side, cand in (("query", p["query_candidate"]), ("target", p["target_candidate"])):
            if r[f"{side}_has_coords"] == "yes":
                x = xi.get((cand, r[f"{side}_staged_resid"]))
                rt_ok &= bool(x) and int(x["rfam_column"]) == int(r["rfam_column"])
            if r[f"{side}_has_coords"] == "no":
                unobs_ok &= r[f"{side}_char"] not in "-." and r[f"{side}_ungapped_position"] != ""
            if r[f"{side}_char"] in "-.":
                unobs_ok &= r[f"{side}_has_coords"] == ""
    G["g4_reference_endpoints_round_trip"] = rt_ok
    G["g4_unobserved_distinct_from_gap"] = unobs_ok
    # 5 ranking after eligibility
    pc = C.read_tsv(C.P("pair_candidates.tsv"))
    by = collections.defaultdict(list)
    for r in pc:
        by[(r["family"], r["seq_pair"])].append(r)
    rk_ok = True
    for k, rs in by.items():
        sel = [r for r in rs if r["status"] == "selected"]
        if not sel:
            continue
        sk = json.loads(sel[0]["key"].replace("(", "[").replace(")", "]"))
        for r in rs:
            if r["status"] == "eligible_alternate":
                rk_ok &= json.loads(r["key"].replace("(", "[").replace(")", "]")) <= sk
            if r["status"] == "pruned_cannot_outrank_selected":
                rk_ok &= json.loads(r["upper_bound"].replace("(", "[").replace(")", "]")) < sk
    G["g5_selected_outranks_all_alternatives_and_bounds"] = rk_ok
    # 6 tests
    t = subprocess.run([sys.executable, "-m", "pytest", "-q", os.path.join(C.ROOT, "tests")], capture_output=True, text=True)
    open(C.P("pytest_results.txt"), "w").write(t.stdout + t.stderr)
    G["g6_tests_pass"] = t.returncode == 0
    G["g6_tests_summary"] = (t.stdout.strip().splitlines() or [""])[-1]
    # 8 step 1 unchanged
    rec = C.jload(C.P("input_reconciliation.json"))
    G["g8_step1_unchanged"] = all(C.sha256(os.path.join(C.STEP1, n)) == h for n, h in rec["step1_file_sha256"].items())
    start = min(r["retrieved_utc"] for r in C.read_tsv(C.MANIFEST))
    import datetime as dt
    t0 = dt.datetime.strptime(start, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc).timestamp()
    newer = []
    for d in ("../v3_full Rfam", "../rfam_star3d_pilot_v2"):
        for p in glob.glob(os.path.join(C.ROOT, d, "**", "*"), recursive=True):
            if os.path.isfile(p) and ".venv" not in p and "__pycache__" not in p and os.path.getmtime(p) > t0:
                newer.append(os.path.relpath(p, C.ROOT))
    G["g8_earlier_study_files_unchanged"] = not newer
    G["g8_earlier_files_modified"] = newer[:10]
    return G


def trace(pairs):
    if not pairs:
        return "No primary pair passed; no end-to-end trace possible."
    p = next((x for x in pairs if x["family"] == "RF00522"), pairs[0])
    ref = [r for r in C.read_tsv(C.P("rfam_pair_reference.tsv.gz")) if r["pair_id"] == p["pair_id"]]
    r = next(r for r in ref if r["correspondence"] == "residue_pair" and r["query_has_coords"] == "yes" and r["target_has_coords"] == "yes")
    xw = [x for x in C.read_tsv(C.P("residue_crosswalk.tsv.gz")) if x["candidate_id"] in (p["query_candidate"], p["target_candidate"])
          and x["rfam_column"] == r["rfam_column"]]
    lines = [f"Pair {p['pair_id']} ({p['family']} {p['family_id']}): query {p['query_pdb_chain']} "
             f"[{p['query_species']}] seed row {p['query_seed_row']}; target {p['target_pdb_chain']} "
             f"[{p['target_species']}] seed row {p['target_seed_row']}.",
             f"Rfam column {r['rfam_column']} ({r['rf_class']}): query {r['query_char']} = {r['query_residue']} "
             f"(staged {p['query_stem']} {p['query_chain']}:{r['query_staged_resid']}), target {r['target_char']} = "
             f"{r['target_residue']} (staged {p['target_stem']} {p['target_chain']}:{r['target_staged_resid']})."]
    for x in xw:
        lines.append(f"  crosswalk {x['candidate_id']}: col {x['rfam_column']} -> ungapped {x['ungapped_position']} -> "
                     f"{x['label_asym_id']}/{x['label_seq_id']} -> auth {x['auth_chain']}:{x['auth_seq_num']}{x['auth_ins_code']} "
                     f"({x['comp_id']}) -> {x['staged_file']} {x['staged_chain']}:{x['staged_resid']}")
    lines.append(f"Origin: query {p['query_origin']}; target {p['target_origin']}.")
    return "\n".join(lines)


def summary(fs, G, interim):
    pairs = C.read_tsv(C.P("pairs_to_run.tsv"))
    cr = [r for r in C.read_tsv(C.P("chain_review.tsv")) if r.get("technical_status")]
    oc = collections.Counter(r["outcome"] for r in fs)
    pf = collections.Counter(p["family"] for p in pairs)
    ids = {r["family"]: r["family_id"] for r in fs}
    reasons = collections.Counter(r["reason"] for r in cr if r["decision"] not in ("verified_natural",))
    tech = collections.Counter(r["technical_status"] for r in cr if r["decision"] == "verified_natural")
    fam_txt = ", ".join(f"{f} {ids[f]} ({n})" for f, n in pf.most_common())
    ex = G.get("g7_end_to_end_trace", "").splitlines()
    lines = [
        f"# SQARMA Step 2 — structure pairs prepared for original STAR3D ({'INTERIM' if interim else 'FINAL'})", "",
        f"*Automated pipeline, {C.now()[:10]}; no human review has occurred. **STAR3D has not yet been run in this phase.***", "",
        f"**Scope.** All {len(fs):,} Rfam 15.1 families were searched (explicit seed tags, Step 1 leads, and both Rfam 15.1 "
        f"family-to-PDB mapping files). {sum(1 for r in fs if r['discovery_status'] == 'exact_seed_match_found')} families have "
        f"a PDB chain matching an entire seed row exactly; {oc.get('no_candidate_in_searched_sources', 0):,} had no candidate "
        f"in the searched sources and {oc.get('candidate_found_no_exact_seed_match', 0)} had candidates without an exact seed match.", "",
        f"**Result.** {len(pf)} families yield eligible pairs; **{len(pairs)} jobs are ready** ({fam_txt}). "
        f"Each job pairs two verified-natural RNAs from different species, different PDB entries and different seed-row "
        f"sequences; inputs: {tech.get('ready', 0)} verified chains fully modelled, {tech.get('ready_with_coverage_flags', 0)} "
        f"with coverage flags (recorded per residue).", "",
        "**Why others did not qualify** (chain decisions; one chain can block many pairs):", ""]
    for k, v in reasons.most_common(6):
        lines.append(f"- {k}: {v}")
    lines += ["", f"Family outcomes: " + "; ".join(f"{k} {v}" for k, v in oc.most_common()) + ".", "",
              "**Traceable example.**", ""] + [f"> {x}" for x in ex[:4]] + [
              "", "**Limits.** Species and natural origin are automated decisions from deposited annotations and NCBI "
              "records; tRNA pairs dominate the job list; pairs within a family share structures (flagged); a pair "
              "being eligible says nothing about whether the methods will agree.", "",
              f"Validation: {sum(1 for v in G.values() if v is True)}/{sum(1 for v in G.values() if isinstance(v, bool))} gates pass "
              f"(validation_report.md)."]
    open(C.P("professor_summary.md"), "w").write("\n".join(lines) + "\n")


def main():
    fs = family_status()
    structures()
    interim = any(r["family_complete"] == "no" for r in fs)
    wb_ok = workbook(interim)
    G = gates(fs)
    G["g7_workbook_reconciles_with_tsvs"] = wb_ok
    pairs = C.read_tsv(C.P("pairs_to_run.tsv"))
    G["g7_end_to_end_trace"] = trace(pairs)
    gate_bools = {k: v for k, v in G.items() if isinstance(v, bool)}
    C.jdump(dict(gates=G, passed=all(gate_bools.values()), interim=interim, utc=C.now()), C.P("validation.json"))
    lines = [f"# Step 2 validation report ({'INTERIM' if interim else 'FINAL'})", "", f"Generated {C.now()} (automated).", "",
             "| Gate | Result |", "|---|---|"]
    for k, v in G.items():
        if isinstance(v, bool):
            lines.append(f"| {k} | {'PASS' if v else 'FAIL'} |")
    lines += ["", "## Details", "```", json.dumps({k: v for k, v in G.items() if not isinstance(v, bool)}, indent=1), "```"]
    open(C.P("validation_report.md"), "w").write("\n".join(lines) + "\n")
    summary(fs, G, interim)
    print("\n".join(lines[:30]))


if __name__ == "__main__":
    main()
