"""Step 2 / stage 5: chain review decisions, pair formation, branch-and-bound representative selection, staged
STAR3D inputs, residue crosswalk and the frozen Rfam pair reference. STAR3D is NOT executed.

Resumable: reviews are cached per candidate (cache/review/<candidate>.json, invalidated when its genome evidence
changes); staged files are rewritten only if missing or their source changed. Candidates whose (parent-resolved
sequence, species) evidence key is missing are written to origin_queue_extra.tsv for stage 3 and marked pending.

Outputs: chain_review.tsv, pair_candidates.tsv (every evaluated/pruned combination, with reason),
pairs_to_run.tsv, inputs_star3d/<stem>.pdb + <stem>.conversion.tsv + file_manifest.tsv,
residue_crosswalk.tsv.gz, rfam_pair_reference.tsv.gz, selection_state.json
"""
import collections
import gzip
import hashlib
import itertools
import json
import os
import sys

import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import s04_coordinates as X  # noqa: E402

OKSCREEN = {"exact_unique", "resolved_by_explicit_tag", "resolved_by_row_organism"}
USABLE_TECH = {"ready", "ready_with_coverage_flags"}
TRNA = {"RF00005"}


# ------------------------------------------------------------------ inputs
def load():
    iv = C.read_tsv(C.P("interval_screening.tsv"))
    dep = {r["candidate_id"]: r["deposited_sequence"] for r in C.read_tsv(C.P("candidate_intervals.tsv"))}
    for r in iv:
        r["deposited_sequence"] = dep.get(r["candidate_id"], "")
    meta = {r["pdb_id"]: r for r in C.read_tsv(C.P("entry_metadata.tsv"))}
    extra = C.read_tsv(C.P("wildcard_rescreen.tsv")) if os.path.exists(C.P("wildcard_rescreen.tsv")) else []
    sc = {r["candidate_id"]: r for r in iv}
    for r in extra:              # wildcard chains re-screened on the mmCIF parent sequence
        if r["candidate_id"] in sc or r["screening_status"] not in OKSCREEN:
            continue
        m = meta.get(r["pdb_id"], {})
        r.update(screen_status="ok", screen_species_taxid="", screen_observed_fraction="1",
                 method=m.get("method", ""), resolution=(m.get("resolution") or "").split(";")[0])
        sc[r["candidate_id"]] = r
    rows = X.load_rows()
    tracks = {}
    with gzip.open(os.path.join(C.STEP1, "family_tracks.jsonl.gz"), "rt") as f:
        for line in f:
            t = json.loads(line)
            tracks[t["family"]] = t
    return sc, rows, tracks, meta


# ------------------------------------------------------------------ origin decision
def origin_decision(rv, fam):
    """Returns (decision, reason_code, evidence_text, evidence_key or None if pending)."""
    if rv.get("mmcif_status") != "ok":
        return "provisional_origin", "MMCIF_CHAIN_UNRESOLVED", rv.get("note", ""), None
    if not rv["species_taxid"]:
        return "provisional_origin", "SPECIES_UNRESOLVED", rv["species_basis"], None
    seq = rv["deposited_parent_sequence"]
    if not set(seq) <= C.STD:
        return "provisional_origin", "UNKNOWN_RESIDUES_IN_CHAIN", "".join(sorted(set(seq) - C.STD)), None
    key = (rv["deposited_parent_sha256"], rv["species_taxid"])
    path = C.CP("origin", f"{key[0][:20]}_{key[1]}.json")
    ev = C.jload(path)
    pos = C.jload(C.CP("origin_exact_local", f"{key[0][:20]}_{key[1]}.json"))
    if pos and pos.get("exact_local_resolves") and (ev is None or ev["outcome"] != "whole_chain_found"):
        # found in an independent record (exact local search, re-confirmed on the Mac) - union rule
        h = pos["whole_chain_hits_confirmed"][0]
        ev = dict(outcome="whole_chain_found", routes_tried=(ev or {}).get("routes_tried", []) + pos["routes_tried"],
                  interval_hits=(ev or {}).get("interval_hits", []),
                  whole_chain_hits=[dict(h, independent="yes", exact_confirmed="yes",
                                         source="exact_local:" + h["source"], independence_basis=h.get("independence_basis", ""))],
                  ncbi_outcome=(ev or {}).get("outcome", "pending"))
    if ev is None:
        closed = C.jload(C.P("STEP2_CLOSED.json"))
        if closed:   # R7: Step 2 closed; unresolved search is provisional, never 'engineered'
            note = f"Exact local search: {pos['outcome']}" if pos else "no completed search"
            return "provisional_origin", "SEARCH_INCOMPLETE", note + "; NCBI search not completed (Step 2 closed, R7)", key
        return "pending", "GENOME_EVIDENCE_PENDING", "", key
    doc = rv["doc_documented_engineering"] in (True, "True")
    why = (f"pdbx_mutation={rv['doc_pdbx_mutation']}" if rv["doc_mutation_documented"] in (True, "True") else "") + \
          (f" words={rv['doc_engineering_words']}" if rv["doc_engineering_words"] else "")
    t7 = "T7" in (rv["doc_engineering_words"] or "")
    good = [h for h in ev["whole_chain_hits"] if h.get("independent") == "yes" and h.get("exact_confirmed") == "yes"]
    evtxt = "; ".join(f"{h['record']}:{h['start']}-{h['end']}{h['strand']} [{h.get('source')}; {h.get('independence_basis')}]"
                      for h in good[:3])
    if ev["outcome"] == "whole_chain_found":
        if doc:
            return "provisional_origin", "CONFLICTING_EVIDENCE", f"genome: {evtxt} | documentation: {why}", key
        return "verified_natural", "WHOLE_CHAIN_INDEPENDENT_RECORD", evtxt, key
    if doc:
        return "excluded_engineered", "T7_5PRIME_LEADER" if t7 else "DOCUMENTED_CONSTRUCT_CHANGE", why.strip(), key
    if ev["outcome"] == "interval_found_flank_differs":
        loc = "; ".join(f"{h['record']}:{h['start']}-{h['end']}{h['strand']}" for h in ev["interval_hits"][:3])
        if fam in TRNA and seq.endswith("CCA"):
            return "provisional_origin", "TRNA_3PRIME_CCA_POSSIBLY_NOT_ENCODED", f"interval at {loc}", key
        both = {x["route"] for x in ev["routes_tried"]} >= {"blast_core_nt"}
        natural_source = (rv.get("doc_src_method") or "").lower() == "nat"
        if natural_source:     # R4: RNA isolated from the organism cannot carry an engineered flank
            return "provisional_origin", "FLANK_DIFFERS_FROM_REFERENCE_NATURAL_SOURCE", f"interval at {loc}", key
        if both:
            return "excluded_engineered", "NON_NATIVE_FLANK_AT_LOCUS", f"interval found at {loc}; whole chain in no searched record", key
        return "provisional_origin", "FLANK_DIFFERS_FROM_REFERENCE", f"interval at {loc}", key
    if ev["outcome"] == "search_failed":
        if C.jload(C.P("STEP2_CLOSED.json")):
            return "provisional_origin", "SEARCH_INCOMPLETE", "search failed; Step 2 closed (R7)", key
        return "pending", "SEARCH_FAILED_RETRY", ev.get("note", ""), key
    return "provisional_origin", "GENOME_NOT_FOUND", "no confirmed whole-chain match in searched records", key


def review(c, rows, tax, ph):
    """Cached per candidate; recomputed when the evidence file appears/changes."""
    path = C.CP("review", c["candidate_id"].replace("/", "_") + ".json")
    rv = C.jload(path)
    if rv is None:
        try:
            rv = X.interval_review(c, rows[c["family"]], tax, ph)
        except Exception as e:     # download/parse failure: technical, not biological; retried next run
            return dict(candidate_id=c["candidate_id"], mmcif_status="retrieval_failed", note=repr(e)[:200],
                        technical_status="technical_unusable_or_unresolved", decision="pending",
                        reason="MMCIF_RETRIEVAL_FAILED")
        C.jdump(rv, path)
    # row-organism resolution must be consistent with the mmCIF species
    if c["screening_status"] == "resolved_by_row_organism" and rv.get("species_taxid") != c.get("screen_species_taxid"):
        rv = dict(rv, technical_status="technical_unusable_or_unresolved",
                  technical_problems=(rv.get("technical_problems", "") + ";row_organism_resolution_not_confirmed").strip(";"))
    d, why, ev, key = origin_decision(rv, c["family"])
    rv = dict(rv, decision=d, reason=why, origin_evidence=ev, evidence_key="|".join(key) if key else "")
    return rv


# ------------------------------------------------------------------ ranking keys
def res_val(x):
    try:
        return float(str(x).split(";")[0])
    except ValueError:
        return None


def combo_key(a, b, fa, fb):
    """Numeric key (bigger = better) and stable id. fa/fb are dicts with obs, comp, res."""
    ra, rb = res_val(fa["res"]), res_val(fb["res"])
    nores = (ra is None) + (rb is None)
    worst = max([r for r in (ra, rb) if r is not None], default=99.0)
    return (min(fa["obs"], fb["obs"]), min(fa["comp"], fb["comp"]), -nores, -worst), "|".join(sorted((a, b)))


def stem_for(c):
    ch = "".join(f"{x}{'u' if x.isupper() else 'l' if x.islower() else 'd'}" for x in c["auth_chain"])
    return f"{c['family'].lower()}_{c['pdb_id'].lower()}_{ch}_{c['seqres_start']}_{c['seqres_end']}"


def walk_pair(fam, combos, reviews, eligible, facts):
    """Branch-and-bound over combinations [(upper_bound, stable_id, a, b)]: every combination whose upper bound is not
    below the best VERIFIED key is evaluated; eligibility is decided before ranking. Returns best, key, pending, recs."""
    combos = sorted(combos, key=lambda x: (tuple(-v for v in x[0]), x[1]))
    best, best_key, pend, recs = None, None, False, []
    for ub, sid, a, b in combos:
        if best_key is not None and tuple(ub) < tuple(best_key[0]):
            recs.append(dict(family=fam, combo=sid, status="pruned_cannot_outrank_selected", upper_bound=ub))
            continue
        ea, eb = eligible(a), eligible(b)
        ra, rb = reviews[a], reviews[b]
        if "pending" in (ra["decision"], rb["decision"]):
            pend = True
            recs.append(dict(family=fam, combo=sid, status="pending_evidence", upper_bound=ub))
            continue
        if not (ea and eb):
            why = [f"{x}:{reviews[x]['decision']}/{reviews[x].get('technical_status')}" for x in (a, b) if not eligible(x)]
            recs.append(dict(family=fam, combo=sid, status="member_not_eligible", reason=";".join(why), upper_bound=ub))
            continue
        if ra["species_taxid"] == rb["species_taxid"]:
            recs.append(dict(family=fam, combo=sid, status="same_species_mmcif", upper_bound=ub))
            continue
        k = combo_key(a, b, facts(a, True), facts(b, True))
        recs.append(dict(family=fam, combo=sid, status="eligible", key=k[0], upper_bound=ub))
        if best_key is None or k[0] > best_key[0] or (k[0] == best_key[0] and k[1] < best_key[1]):
            best, best_key = (a, b), k
    return best, best_key, pend, recs


def star3d_read(path):
    """Residues as STAR3D PDB.java reads them: first model, chain stops at TER, fixed columns, ATOM+HETATM.
    Returns ([((chain, num, icode), resname)], parse_ok)."""
    order, seen, ok, ter = [], None, True, set()
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if line.startswith("TER"):
            if seen:
                ter.add(seen[0])
            continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                key = (line[21], int(line[22:26]), line[26].strip())
                int(line[6:11])
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                ok = False
                continue
            if key[0] in ter:
                continue
            if key != seen:
                order.append((key, line[17:20].strip()))
                seen = key
    return order, ok


# ------------------------------------------------------------------ staging
def stage(c, rv, rows, out_dir):
    pdb, asym = c["pdb_id"], rv["label_asym_id"]
    lo, hi = int(c["seqres_start"]), int(c["seqres_end"])
    stem = stem_for(c)
    path = os.path.join(out_dir, f"{stem}.pdb")
    src = X.mmcif_path(pdb)
    st = gemmi.read_structure(src)
    st.setup_entities()
    while len(st) > 1:
        del st[1]
    st.remove_alternative_conformations()
    auth = rv["auth_chain"]
    alias = auth if len(auth) == 1 else "A"
    out = gemmi.Structure()
    out.name = stem[:4]
    model, newch = gemmi.Model("1"), gemmi.Chain(alias)
    conv = []
    for ch in st[0]:
        for r in ch:
            if r.subchain != asym or r.label_seq is None or not lo <= r.label_seq <= hi:
                continue
            newch.add_residue(r.clone())
            conv.append(dict(staged_index0=len(conv), label_asym_id=asym, label_seq_id=r.label_seq, auth_chain=auth,
                             staged_chain=alias, auth_seq_num=r.seqid.num, auth_ins_code=r.seqid.icode.strip(),
                             comp_id=r.name, record=("HETATM" if r.het_flag == "H" else "ATOM"), n_atoms=len(r)))
    model.add_chain(newch)
    out.add_model(model)
    out.setup_entities()
    opts = gemmi.PdbWriteOptions(minimal=True)
    opts.ter_records = True
    tmp = path + ".tmp"
    out.write_pdb(tmp, opts)
    os.replace(tmp, path)
    C.write_tsv(os.path.join(out_dir, f"{stem}.conversion.tsv"), conv)
    order, ok_cols = star3d_read(path)
    want = [(alias, x["auth_seq_num"], x["auth_ins_code"]) for x in conv]
    checks = dict(fixed_columns_parse=ok_cols, single_chain=len({k[0] for k, _ in order}) == 1 and (not order or order[0][0][0] == alias),
                  residue_order_identity=[(k, n) for k, n in order] == [((alias, x["auth_seq_num"], x["auth_ins_code"]), x["comp_id"]) for x in conv],
                  residue_ids_unique=len(set(want)) == len(want))
    return dict(stem=stem, file=os.path.relpath(path, C.ROOT), sha256=C.sha256(path), staged_chain=alias,
                chain_alias_used="yes" if alias != auth else "no", residues_written=len(conv),
                atoms_written=sum(x["n_atoms"] for x in conv), source_mmcif=os.path.relpath(src, C.ROOT),
                source_mmcif_sha256=C.sha256(src), reader_checks=";".join(f"{k}={'pass' if v else 'FAIL'}" for k, v in checks.items()),
                reader_ok=all(checks.values())), conv


def crosswalk(c, rv, rows, conv, stem):
    row = rows[c["family"]][c["seed_row_id"]]
    d = X.entry_data(c["pdb_id"])
    pos = {p["seq_id"]: p for p in d["chains"][rv["label_asym_id"]]["positions"]}
    staged = {x["label_seq_id"]: x for x in conv}
    lo = int(c["seqres_start"])
    out, k = [], 0
    for col, ch in enumerate(row["original_gapped_sequence"], 1):
        if ch in "-.":
            continue
        k += 1
        lab = lo + k - 1
        p = pos[lab]
        s = staged.get(lab)
        miss = [a for a in X.BACKBONE if a not in p["atoms"]] if p["observed"] else list(X.BACKBONE)
        out.append(dict(candidate_id=c["candidate_id"], family=c["family"], seed_row_id=c["seed_row_id"], rfam_column=col,
                        row_char=ch, ungapped_position=k, pdb_id=c["pdb_id"], entity_id=rv["entity_id"],
                        label_asym_id=rv["label_asym_id"], label_seq_id=lab, auth_chain=rv["auth_chain"],
                        auth_seq_num=p["auth_num"], auth_ins_code=p["icode"], comp_id=p["mon"], parent=p["parent"],
                        parent_source=p.get("parent_source", ""), modified="yes" if p["mon"] not in C.STD else "no",
                        model=1, observed="yes" if p["observed"] else "no",
                        missing_backbone_atoms=",".join(miss) if p["observed"] else "unobserved",
                        staged_file=f"inputs_star3d/{stem}.pdb" if s else "",
                        staged_chain=s["staged_chain"] if s else "",
                        staged_resid=f"{s['auth_seq_num']}{s['auth_ins_code']}" if s else "",
                        transformation=("chain alias " + rv["auth_chain"] + "->" + s["staged_chain"]) if s and s["staged_chain"] != rv["auth_chain"] else ""))
    return out


# ------------------------------------------------------------------ main
def main():
    tax = C.Taxonomy()
    ph = {str(x) for x in C.policy()["species"]["placeholders_unresolved"]}
    sc, rows, tracks, meta = load()
    partner = C.jload(C.P("possible_partner.json"), {})
    yields = {r["family"]: r for r in C.read_tsv(C.P("yield_by_family.tsv"))}
    fams_q = sorted({c["family"] for c in sc.values() if c.get("screen_status") == "ok"
                     and c["screening_status"] in OKSCREEN and (partner.get(c["candidate_id"]) or "wildcard" in c.get("row_basis", ""))})
    reviews, pend_keys, combos_out, selected, fam_state = {}, {}, [], [], {}

    def rv_of(cid):
        if cid not in reviews:
            reviews[cid] = review(sc[cid], rows, tax, ph)
            r = reviews[cid]
            if r["decision"] == "pending" and r.get("evidence_key"):
                pend_keys[cid] = r
        return reviews[cid]

    def facts(cid, reviewed):
        c = sc[cid]
        if reviewed:
            r = reviews[cid]
            return dict(obs=float(r.get("observed_fraction") or 0), comp=float(r.get("backbone_complete_fraction") or 0),
                        res=r.get("resolution", ""), sp=r.get("species_taxid", ""))
        ub = float(c.get("screen_observed_fraction") or 0)
        return dict(obs=ub, comp=ub, res=c.get("resolution", ""), sp=c.get("screen_species_taxid", ""))

    def eligible(cid):
        r = rv_of(cid)
        return r["decision"] == "verified_natural" and r.get("technical_status") in USABLE_TECH

    for fam in fams_q:
        cands = [c for c in sc.values() if c["family"] == fam and c.get("screen_status") == "ok"
                 and c["screening_status"] in OKSCREEN]
        groups = collections.defaultdict(list)
        for c in cands:
            groups[c["interval_sequence"]].append(c["candidate_id"])
        exhaustive = False
        for attempt in range(2):
            fam_combos, fam_sel, pending_pairs = [], [], 0
            for sa, sb in itertools.combinations(sorted(groups), 2):
                combos = []
                for a in groups[sa]:
                    for b in groups[sb]:
                        if sc[a]["pdb_id"] == sc[b]["pdb_id"]:
                            continue
                        fa, fb = facts(a, a in reviews), facts(b, b in reviews)
                        if fa["sp"] and fb["sp"] and fa["sp"] == fb["sp"]:
                            continue          # same species at species rank (screening; asserted below)
                        ub = (1.0, 1.0, 0, 0) if exhaustive else combo_key(a, b, fa, fb)[0]
                        combos.append((ub, "|".join(sorted((a, b))), a, b))
                if not combos:
                    continue
                best, best_key, pend, recs = walk_pair(fam, combos, reviews, eligible, facts)
                if best and not pend:
                    for r in recs:
                        if r["combo"] == best_key[1]:
                            r["status"] = "selected"
                        elif r["status"] == "eligible":
                            r["status"] = "eligible_alternate"
                    fam_sel.append(dict(family=fam, seq_a=sa, seq_b=sb, members=best, key=best_key))
                elif pend:
                    pending_pairs += 1
                fam_combos += [dict(r, seq_pair=f"{hashlib.sha256(sa.encode()).hexdigest()[:8]}|{hashlib.sha256(sb.encode()).hexdigest()[:8]}")
                               for r in recs]
            # assert the screening facts used for pruning against the mmCIF for every reviewed structure
            bad = [cid for cid in groups_ids(groups) if cid in reviews and reviews[cid].get("mmcif_status") == "ok" and (
                abs(float(reviews[cid]["observed_fraction"]) - float(sc[cid].get("screen_observed_fraction") or 0)) > 1e-3
                or (sc[cid].get("screen_species_taxid") and reviews[cid]["species_taxid"] != sc[cid]["screen_species_taxid"]))]
            if bad and not exhaustive:
                exhaustive = True
                print(f"  {fam}: screening facts differ from mmCIF for {len(bad)} structures -> exhaustive re-walk", flush=True)
                continue
            break
        combos_out += fam_combos
        selected += fam_sel
        fam_state[fam] = dict(pairs_selected=len(fam_sel), seq_pairs_pending=pending_pairs, exhaustive=exhaustive,
                              screening_mismatch_structures=len(bad))
        print(f"{fam}: groups={len(groups)} selected={len(fam_sel)} pending={pending_pairs} reviewed="
              f"{sum(1 for c in cands if c['candidate_id'] in reviews)}/{len(cands)}", flush=True)

    # pending evidence -> stage 3 queue
    extra = []
    for cid, r in pend_keys.items():
        c = sc[cid]
        extra.append(dict(candidate_id=cid, family=c["family"], deposited_seq_sha256=r["deposited_parent_sha256"],
                          deposited_sequence=r["deposited_parent_sequence"], interval_sequence=c["interval_sequence"],
                          species_taxid=r["species_taxid"], source_taxids=r.get("species_source_taxa", "")))
    C.write_tsv(C.P("origin_queue_extra.tsv"), extra,
                ["candidate_id", "family", "deposited_seq_sha256", "deposited_sequence", "interval_sequence",
                 "species_taxid", "source_taxids"])

    # chain review table (every screened candidate; unreviewed ones get an explicit status)
    cr = []
    for cid, c in sorted(sc.items()):
        r = reviews.get(cid)
        if r is None:
            if c["screening_status"] not in OKSCREEN:
                st = "not_primary_" + c["screening_status"]
            elif c["family"] not in fams_q or not partner.get(cid):
                st = "not_evaluated_no_possible_partner"
            else:
                st = "not_evaluated_cannot_outrank_selected"
            cr.append(dict(candidate_id=cid, family=c["family"], pdb_id=c["pdb_id"], auth_chain=c["auth_chain"],
                           seed_row_id=c.get("seed_row_id", ""), matching_rows=c.get("matching_rows", ""),
                           screening_status=c["screening_status"], seqres_interval=f"{c['seqres_start']}-{c['seqres_end']}",
                           screen_species=c.get("screen_species_name", ""), decision=st, reason="", technical_status=""))
            continue
        row = {k: v for k, v in r.items() if k != "deposited_parent_sequence"}
        row.update(candidate_id=cid, family=c["family"], pdb_id=c["pdb_id"], seed_row_id=c.get("seed_row_id", ""),
                   matching_rows=c.get("matching_rows", ""), screening_status=c["screening_status"],
                   seqres_interval=f"{c['seqres_start']}-{c['seqres_end']}", screen_species=c.get("screen_species_name", ""),
                   decision_mode="automated", decision_utc=C.now())
        row.setdefault("auth_chain", c["auth_chain"])
        cr.append(row)
    C.write_tsv(C.P("chain_review.tsv"), cr, sorted({k for r in cr for k in r},
                                                       key=lambda k: (k not in CR_ORDER, CR_ORDER.index(k) if k in CR_ORDER else 0, k)))
    C.write_tsv(C.P("pair_candidates.tsv"), combos_out,
                ["family", "seq_pair", "combo", "status", "reason", "key", "upper_bound"])

    # staging + crosswalk + reference for selected pairs
    out_dir = C.P("inputs_star3d")
    os.makedirs(out_dir, exist_ok=True)
    staged, xw_all, ref_all, pairs = {}, [], [], []
    use = collections.Counter(m for s in selected for m in s["members"])
    for s in sorted(selected, key=lambda s: (s["family"], s["key"][1])):
        a, b = sorted(s["members"], key=lambda x: stem_for(sc[x]))
        for m in (a, b):
            if m not in staged:
                info, conv = stage(sc[m], reviews[m], rows, out_dir)
                staged[m] = (info, conv)
                xw_all += crosswalk(sc[m], reviews[m], rows, conv, info["stem"])
        qa, ta = staged[a][0], staged[b][0]
        pid = f"{s['family']}__{qa['stem']}__vs__{ta['stem']}"
        qrow, trow = rows[s["family"]][sc[a]["seed_row_id"]], rows[s["family"]][sc[b]["seed_row_id"]]
        rf = tracks[s["family"]]["gc"].get("RF", "")
        qx = {x["rfam_column"]: x for x in xw_all if x["candidate_id"] == a}
        tx = {x["rfam_column"]: x for x in xw_all if x["candidate_id"] == b}
        both, same = 0, 0
        for col in range(1, len(qrow["original_gapped_sequence"]) + 1):
            q, t = qx.get(col), tx.get(col)
            if q and t:
                both += 1
                same += q["parent"] == t["parent"]
            ref_all.append(dict(pair_id=pid, rfam_column=col,
                                rf_class=("unknown" if not rf else "insert" if rf[col - 1] in ".-~" else "consensus"),
                                rf_char=rf[col - 1] if rf else "", query_char=qrow["original_gapped_sequence"][col - 1],
                                target_char=trow["original_gapped_sequence"][col - 1],
                                query_ungapped_position=q["ungapped_position"] if q else "",
                                target_ungapped_position=t["ungapped_position"] if t else "",
                                query_residue=f"{q['auth_chain']}:{q['auth_seq_num']}{q['auth_ins_code']}" if q else "",
                                target_residue=f"{t['auth_chain']}:{t['auth_seq_num']}{t['auth_ins_code']}" if t else "",
                                query_staged_resid=q["staged_resid"] if q else "", target_staged_resid=t["staged_resid"] if t else "",
                                query_has_coords=("yes" if q and q["staged_resid"] else "no") if q else "",
                                target_has_coords=("yes" if t and t["staged_resid"] else "no") if t else "",
                                correspondence=("residue_pair" if q and t else "both_gap" if not q and not t else
                                                "query_only" if q else "target_only")))
        ra, rb = reviews[a], reviews[b]
        flags = []
        for m, r in ((a, ra), (b, rb)):
            if r["technical_status"] == "ready_with_coverage_flags":
                flags.append(f"{'query' if m == a else 'target'}_coverage:{r['observed']}/{r['interval_length']} observed,"
                             f" {r['backbone_complete']} full backbone")
            if use[m] > 1:
                flags.append(f"{'query' if m == a else 'target'}_reused_in_{use[m]}_jobs")
            if r.get("modified_residues"):
                flags.append(f"{'query' if m == a else 'target'}_modified:{r['modified_residues']}")
        pairs.append(dict(pair_id=pid, family=s["family"], family_id=yields[s["family"]]["family_id"],
                          query_candidate=a, target_candidate=b,
                          query_seed_row=sc[a]["seed_row_id"], target_seed_row=sc[b]["seed_row_id"],
                          query_species=f"{ra['species_name']} ({ra['species_taxid']})",
                          target_species=f"{rb['species_name']} ({rb['species_taxid']})",
                          query_pdb_chain=f"{sc[a]['pdb_id']}_{ra['auth_chain']}", target_pdb_chain=f"{sc[b]['pdb_id']}_{rb['auth_chain']}",
                          query_label_interval=f"{ra['label_asym_id']}:{sc[a]['seqres_start']}-{sc[a]['seqres_end']}",
                          target_label_interval=f"{rb['label_asym_id']}:{sc[b]['seqres_start']}-{sc[b]['seqres_end']}",
                          query_stem=qa["stem"], query_chain=qa["staged_chain"], query_file=qa["file"], query_sha256=qa["sha256"],
                          target_stem=ta["stem"], target_chain=ta["staged_chain"], target_file=ta["file"], target_sha256=ta["sha256"],
                          query_technical=ra["technical_status"], target_technical=rb["technical_status"],
                          query_origin=f"{ra['decision']} ({ra['origin_evidence'][:120]})",
                          target_origin=f"{rb['decision']} ({rb['origin_evidence'][:120]})",
                          ranking_key=json.dumps(s["key"][0]),
                          pair_identity=f"{same}/{both} columns where both rows have a residue ({(same / both * 100 if both else 0):.1f}%)",
                          flags=";".join(flags),
                          step3_preprocess=f"java -cp STAR3D.jar Preprocess {qa['stem']} {qa['staged_chain']} && "
                                           f"java -cp STAR3D.jar Preprocess {ta['stem']} {ta['staged_chain']}",
                          step3_align=f"java -jar STAR3D.jar -o out/{pid}__forward.aln -p {qa['stem']} {qa['staged_chain']} "
                                      f"{ta['stem']} {ta['staged_chain']}",
                          star3d_run="no"))
    C.write_tsv(C.P("pairs_to_run.tsv"), pairs, list(pairs[0]) if pairs else ["pair_id"])
    C.write_tsv(C.P("residue_crosswalk.tsv.gz"), xw_all, list(xw_all[0]) if xw_all else ["candidate_id"])
    C.write_tsv(C.P("rfam_pair_reference.tsv.gz"), ref_all, list(ref_all[0]) if ref_all else ["pair_id"])
    man = [dict(candidate_id=m, **{k: v for k, v in info.items()}, pdb_id=sc[m]["pdb_id"],
                auth_chain=reviews[m]["auth_chain"], label_asym_id=reviews[m]["label_asym_id"],
                label_interval=f"{sc[m]['seqres_start']}-{sc[m]['seqres_end']}",
                rcsb_url=f"https://files.rcsb.org/download/{sc[m]['pdb_id']}.cif.gz", model=1,
                altloc_rule="first conformer (gemmi remove_alternative_conformations)",
                transformations=("chain alias " + reviews[m]["auth_chain"] + "->" + info["staged_chain"]) if info["chain_alias_used"] == "yes" else "none")
           for m, (info, conv) in sorted(staged.items())]
    C.write_tsv(C.P("inputs_star3d", "file_manifest.tsv"), man, list(man[0]) if man else ["candidate_id"])
    C.jdump(dict(families=fam_state, pending_evidence_keys=len({r["evidence_key"] for r in pend_keys.values()}),
                 pairs_selected=len(pairs), staged_files=len(staged), utc=C.now()), C.P("selection_state.json"))
    print(f"selected pairs {len(pairs)} | staged {len(staged)} | pending evidence keys "
          f"{len({r['evidence_key'] for r in pend_keys.values()})} | reviewed {len(reviews)}")


def groups_ids(groups):
    return [x for v in groups.values() for x in v]


CR_ORDER = ["candidate_id", "family", "pdb_id", "auth_chain", "label_asym_id", "entity_id", "seed_row_id", "matching_rows",
            "screening_status", "seqres_interval", "mmcif_interval_confirmed", "species_taxid", "species_name",
            "species_basis", "decision", "reason", "origin_evidence", "technical_status", "technical_problems",
            "observed_fraction", "backbone_complete_fraction", "interval_length", "cropped_5prime", "cropped_3prime"]

if __name__ == "__main__":
    main()
