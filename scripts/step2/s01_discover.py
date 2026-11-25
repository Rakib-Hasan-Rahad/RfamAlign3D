"""Step 2 / stage 1: reconcile Step 1 inputs, collect candidate evidence for ALL families, and screen every candidate
chain against every seed row of its family using the deposited polymer sequences (pdb_seqres).

Evidence sources: explicit #=GR row tags (Step 1), Step 1 leads, and the Rfam 15.1 family-to-PDB mapping (union of
Rfam.pdb.gz and database_files/pdb_full_region.txt.gz). Screening is exact substring matching of the whole ungapped
row; the result is confirmed later against the mmCIF (stage 3). Nothing here is an eligibility decision.

Outputs (outputs/02_structure_preparation/):
  input_reconciliation.json, discovery_evidence.tsv, candidate_chains.tsv, candidate_intervals.tsv,
  family_discovery.tsv
"""
import collections
import gzip
import itertools
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

EXPECTED = dict(families=4227, seed_rows=109390, links=383, linked_rows=172, link_entries=306, link_chains=381,
                leads=195, families_ge2_linked_seqs=31, possible_unordered_sequence_pairs=440)


def load_step1():
    fams = C.read_tsv(os.path.join(C.STEP1, "families.tsv"))
    links = C.read_tsv(os.path.join(C.STEP1, "row_pdb_links.tsv"))
    leads = C.read_tsv(os.path.join(C.STEP1, "pdb_leads.tsv"))
    rows = collections.defaultdict(dict)
    with gzip.open(os.path.join(C.STEP1, "seed_rows.jsonl.gz"), "rt") as f:
        for line in f:
            r = json.loads(line)
            rows[r["family"]][r["seed_row_id"]] = r
    return fams, links, leads, rows


def reconcile(fams, links, leads, rows):
    by = collections.defaultdict(set)
    for x in links:
        by[x["family"]].add(x["sequence_group_hash"])
    got = dict(families=len(fams), seed_rows=sum(len(v) for v in rows.values()), links=len(links),
               linked_rows=len({(x["family"], x["seed_row_id"]) for x in links}),
               link_entries=len({x["pdb_id"] for x in links}),
               link_chains=len({(x["pdb_id"], x["chain_id"]) for x in links}), leads=len(leads),
               families_ge2_linked_seqs=sum(1 for v in by.values() if len(v) >= 2),
               possible_unordered_sequence_pairs=sum(len(v) * (len(v) - 1) // 2 for v in by.values()))
    prov = C.jload(os.path.join(C.STEP1, "provenance.json"), {})
    return dict(expected=EXPECTED, observed=got, all_match=got == EXPECTED,
                step1_provenance=prov,
                step1_file_sha256={n: C.sha256(os.path.join(C.STEP1, n)) for n in
                                   ("families.tsv", "row_pdb_links.tsv", "pdb_leads.tsv", "seed_rows.jsonl.gz",
                                    "family_tracks.jsonl.gz", "provenance.json")})


def load_seqres(path):
    chains = {}
    cur = None
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                p = line[1:].split(None, 3)
                cur = p[0]
                chains[cur] = dict(mol=p[1].replace("mol:", ""), length=int(p[2].replace("length:", "")),
                                   name=p[3].strip() if len(p) > 3 else "", seq=[])
            else:
                chains[cur]["seq"].append(line.strip())
    for v in chains.values():
        v["seq"] = "".join(v["seq"]).upper()
    by_entry = collections.defaultdict(list)
    for k in chains:
        pdb, ch = k.split("_", 1)
        by_entry[pdb.upper()].append(ch)
    return chains, by_entry


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def wildcard_find(hay, row):
    """Occurrences where unknown letters (non-ACGU) on EITHER side are allowed to match anything. Diagnostic only."""
    n, out = len(row), []
    for i in range(len(hay) - n + 1):
        ok = True
        for a, b in zip(row, hay[i:i + n]):
            if a != b and a in C.STD and b in C.STD:
                ok = False
                break
        if ok:
            out.append(i)
    return out


def match_chain(seq, fam_rows, tag_rows):
    """Exact screening of one deposited sequence against every row of one family.
    Returns (intervals, diagnostics). Each interval: dict(start,end (1-based), rows=[ids], status, basis)."""
    occ = []
    by_seq = collections.defaultdict(list)
    for rid, r in fam_rows.items():
        by_seq[r["ungapped_sequence"].upper()].append(rid)
    for s, rids in by_seq.items():
        if not set(s) <= C.STD or len(s) > len(seq):
            continue
        for i in find_all(seq, s):
            occ.append((i + 1, i + len(s), s, sorted(rids)))
    intervals = []
    occ.sort()
    clusters = []
    for o in occ:
        if clusters and o[0] <= clusters[-1]["end"]:
            c = clusters[-1]
            c["end"] = max(c["end"], o[1])
            c["occ"].append(o)
        else:
            clusters.append(dict(start=o[0], end=o[1], occ=[o]))
    seq_count = collections.Counter(o[2] for o in occ)
    for c in clusters:
        rids = sorted({r for o in c["occ"] for r in o[3]})
        strings = {o[2] for o in c["occ"]}
        repeated = any(seq_count[s] > 1 for s in strings)
        if len(c["occ"]) == 1 and len(rids) == 1 and not repeated:
            st, basis, chosen = "exact_unique", "single_row_single_occurrence", rids[0]
        elif repeated:
            st, basis, chosen = "ambiguous_multiple_intervals", "row string occurs more than once in chain", ""
        else:
            tagged = [r for r in rids if r in tag_rows]
            if len(tagged) == 1 and len(strings) == 1:
                st, basis, chosen = "resolved_by_explicit_tag", "#=GR tag names one of the matching rows", tagged[0]
            else:
                st, basis, chosen = "ambiguous_rows", f"{len(rids)} row IDs match this interval", ""
        o = next((o for o in c["occ"] if chosen in o[3]), c["occ"][0]) if chosen else None
        intervals.append(dict(start=o[0] if o else c["start"], end=o[1] if o else c["end"], rows=rids,
                              status=st, basis=basis, chosen_row=chosen,
                              n_distinct_row_strings=len(strings)))
    diag = ""
    if not intervals:
        best = []
        for s, rids in by_seq.items():
            if len(s) <= len(seq) and (not set(s) <= C.STD or not set(seq) <= C.STD):
                w = wildcard_find(seq, s)
                if w:
                    best.append((rids[0], w[0] + 1, "".join(sorted(set(s) - C.STD)), "".join(sorted(set(seq) - C.STD))))
        if best:
            diag = "wildcard_only:" + ";".join(f"{r}@{p} row_unknown={a or '-'} deposited_unknown={b or '-'}"
                                               for r, p, a, b in best[:5])
    return intervals, diag


def main():
    os.makedirs(C.OUT, exist_ok=True)
    fams, links, leads, rows = load_step1()
    rec = reconcile(fams, links, leads, rows)
    C.jdump(rec, C.P("input_reconciliation.json"))
    print("Step 1 reconciliation:", "OK" if rec["all_match"] else "DIFFERENCE", rec["observed"])

    for key, url, path in [
        ("rfam_pdb_gz", "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.pdb.gz", "inputs/rfam/15.1/Rfam.pdb.gz"),
        ("rfam_pdb_full_region", "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/pdb_full_region.txt.gz",
         "inputs/rfam/15.1/database_files/pdb_full_region.txt.gz"),
        ("rfam_pdb_full_region_sql", "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/pdb_full_region.sql",
         "inputs/rfam/15.1/database_files/pdb_full_region.sql"),
        ("pdb_seqres", "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz", "inputs/pdb/pdb_seqres.txt.gz"),
        ("ncbi_taxdump", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", "inputs/ncbi_taxonomy/taxdump.tar.gz")]:
        done = os.path.exists(C.MANIFEST) and any(r["key"] == key for r in C.read_tsv(C.MANIFEST))
        if not done:
            C.log_source(key, url, os.path.join(C.ROOT, path), "downloaded_2026-11-13",
                         "header Last-Modified recorded in inputs/*/*.headers where available")

    seqres, by_entry = load_seqres(os.path.join(C.ROOT, "inputs", "pdb", "pdb_seqres.txt.gz"))

    # ---------------- evidence records
    ev = []
    for x in links:
        ev.append(dict(family=x["family"], source="explicit_seed_row_tag", pdb_id=x["pdb_id"], chain_hint=x["chain_id"],
                       seed_row_id=x["seed_row_id"], detail=x["raw_tag"], searched="yes"))
    for x in leads:
        ctx = x["original_evidence"]
        neg = ""
        if x["evidence_source"] == "#=GF CC":
            i = ctx.upper().find(x["candidate_pdb_id"])
            tail = re.split(r"[\[.]", ctx[i:], maxsplit=1)[0] if i >= 0 else ""
            if re.search(r"not included", tail, re.I):
                neg = "negative_statement:not included"
        ev.append(dict(family=x["family"], source="step1_lead" + ("_GR_chainless" if x["evidence_source"] == "#=GR" else "_GF_CC"),
                       pdb_id=x["candidate_pdb_id"], chain_hint=x["candidate_chain"], seed_row_id=x["seed_row_id"],
                       detail=(neg + " | " if neg else "") + ctx[:200], searched="no" if neg else "yes"))
    for fname, label in [("inputs/rfam/15.1/Rfam.pdb.gz", "rfam_pdb_gz"),
                         ("inputs/rfam/15.1/database_files/pdb_full_region.txt.gz", "pdb_full_region")]:
        with gzip.open(os.path.join(C.ROOT, fname), "rt") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                ev.append(dict(family=p[0], source=f"rfam_mapping:{label}", pdb_id=p[1].upper(), chain_hint=p[2],
                               seed_row_id="", searched="yes",
                               detail=f"pdb_start={p[3]} pdb_end={p[4]} bits={p[5]} evalue={p[6]} cm={p[7]}-{p[8]}"
                                      + (f" is_significant={p[10]}" if len(p) > 10 else "")))

    # ---------------- candidate chains
    cand = collections.defaultdict(lambda: dict(sources=set(), ev=[]))
    entry_missing = collections.defaultdict(set)
    tag_rows = collections.defaultdict(set)          # (family, PDB, chain or '*') -> rows named by tags
    for i, e in enumerate(ev):
        e["evidence_id"] = f"E{i + 1:06d}"
        if e["searched"] != "yes":
            continue
        pdb = e["pdb_id"]
        if e["source"] in ("explicit_seed_row_tag",) or e["source"].startswith("rfam_mapping"):
            chains = [e["chain_hint"]]
        else:
            chains = [c for c in by_entry.get(pdb, []) if seqres[f"{pdb.lower()}_{c}"]["mol"] == "na"]
            if pdb not in by_entry:
                entry_missing[e["family"]].add(pdb)
        if e["seed_row_id"]:
            tag_rows[(e["family"], pdb, e["chain_hint"] if e["source"] == "explicit_seed_row_tag" else "*")].add(e["seed_row_id"])
        for ch in chains:
            k = (e["family"], pdb, ch)
            cand[k]["sources"].add(e["source"].split(":")[0])
            cand[k]["ev"].append(e["evidence_id"])
    C.write_tsv(C.P("discovery_evidence.tsv"), ev,
                ["evidence_id", "family", "source", "pdb_id", "chain_hint", "seed_row_id", "searched", "detail"])

    chain_rows, interval_rows = [], []
    for (fam, pdb, ch), v in sorted(cand.items()):
        key = f"{pdb.lower()}_{ch}"
        s = seqres.get(key)
        base = dict(family=fam, pdb_id=pdb, auth_chain=ch, sources=";".join(sorted(v["sources"])),
                    evidence_ids=";".join(v["ev"]))
        if s is None:
            st = "entry_absent_from_seqres" if pdb not in by_entry else "chain_absent_from_seqres"
            chain_rows.append(dict(base, seqres_status=st, screening="search_incomplete"))
            continue
        if s["mol"] != "na":
            chain_rows.append(dict(base, seqres_status="not_nucleic_acid", seqres_length=s["length"],
                                   screening="no_exact_seed_match"))
            continue
        trows = tag_rows.get((fam, pdb, ch), set()) | tag_rows.get((fam, pdb, "*"), set())
        ivs, diag = match_chain(s["seq"], rows[fam], trows)
        dep_hash = C.hashlib.sha256(s["seq"].encode()).hexdigest()
        chain_rows.append(dict(base, seqres_status="present", seqres_length=s["length"], seqres_name=s["name"],
                               deposited_seq_sha256=dep_hash, deposited_has_non_ACGU="".join(sorted(set(s["seq"]) - C.STD)),
                               n_intervals=len(ivs), screening=("exact_seed_match" if any(
                                   i["status"] in ("exact_unique", "resolved_by_explicit_tag") for i in ivs)
                                   else "ambiguous_seed_match" if ivs else
                                   "wildcard_only_needs_mmcif" if diag else "no_exact_seed_match"),
                               diagnostic=diag, explicit_tag_rows=";".join(sorted(trows))))
        for iv in ivs:
            cid = f"{fam}__{pdb}_{ch}__{iv['start']}-{iv['end']}"
            interval_rows.append(dict(candidate_id=cid, family=fam, pdb_id=pdb, auth_chain=ch,
                                      seqres_start=iv["start"], seqres_end=iv["end"], screening_status=iv["status"],
                                      row_basis=iv["basis"], seed_row_id=iv["chosen_row"], matching_rows=";".join(iv["rows"]),
                                      n_matching_rows=len(iv["rows"]), explicit_tag_rows=";".join(sorted(trows)),
                                      interval_sequence=s["seq"][iv["start"] - 1:iv["end"]],
                                      deposited_sequence=s["seq"], deposited_length=s["length"],
                                      deposited_seq_sha256=dep_hash, sources=base["sources"]))
    C.write_tsv(C.P("candidate_chains.tsv"), chain_rows,
                ["family", "pdb_id", "auth_chain", "sources", "evidence_ids", "seqres_status", "seqres_length",
                 "seqres_name", "deposited_seq_sha256", "deposited_has_non_ACGU", "n_intervals", "screening",
                 "diagnostic", "explicit_tag_rows"])
    C.write_tsv(C.P("candidate_intervals.tsv"), interval_rows)

    # ---------------- family discovery status (all families)
    bych = collections.defaultdict(list)
    for r in chain_rows:
        bych[r["family"]].append(r)
    fam_rows = []
    for f in fams:
        a = f["family"]
        cs = bych.get(a, [])
        n_exact = sum(1 for r in cs if r["screening"] == "exact_seed_match")
        incomplete = [r for r in cs if r["screening"] == "search_incomplete"] or entry_missing.get(a)
        wild = sum(1 for r in cs if r["screening"] == "wildcard_only_needs_mmcif")
        if n_exact:
            st = "exact_seed_match_found"
        elif not cs and not entry_missing.get(a):
            st = "no_candidate_in_searched_sources"
        elif incomplete or wild:
            st = "search_incomplete"
        else:
            st = "candidate_found_no_exact_seed_match"
        fam_rows.append(dict(family=a, family_id=f["family_id"], discovery_status=st,
                             evidence_records=sum(1 for e in ev if e["family"] == a),
                             candidate_entries=len({r["pdb_id"] for r in cs}), candidate_chains=len(cs),
                             chains_exact_seed_match=n_exact,
                             chains_ambiguous=sum(1 for r in cs if r["screening"] == "ambiguous_seed_match"),
                             chains_wildcard_only=wild,
                             chains_search_incomplete=sum(1 for r in cs if r["screening"] == "search_incomplete"),
                             entries_missing_from_seqres=";".join(sorted(
                                 {r["pdb_id"] for r in cs if r["seqres_status"] == "entry_absent_from_seqres"}
                                 | entry_missing.get(a, set())))))
    C.write_tsv(C.P("family_discovery.tsv"), fam_rows)
    print("evidence", len(ev), "candidate chains", len(chain_rows), "intervals", len(interval_rows))
    print(collections.Counter(r["screening"] for r in chain_rows))
    print(collections.Counter(r["screening_status"] for r in interval_rows))
    print(collections.Counter(r["discovery_status"] for r in fam_rows))


if __name__ == "__main__":
    main()
