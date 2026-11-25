"""Step 2 / stage 2: entry metadata (RCSB GraphQL, derived from the mmCIF; screening/ranking only), obsolete status,
screening-level species (pinned NCBI taxdump), row-organism resolution of ambiguous intervals, and the early
yield check (yield_by_family.tsv). Nothing here is a final eligibility decision.

Outputs: entry_metadata.tsv, interval_screening.tsv, yield_by_family.tsv, cache/rcsb/graphql_*.json,
         cache/rcsb/holdings_removed.json, cache/ncbi/row_accession_taxa.json
"""
import collections
import json
import os
import re
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

GQL = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} exptl{method}
 rcsb_entry_info{resolution_combined deposited_model_count}
 rcsb_accession_info{deposit_date initial_release_date revision_date major_revision minor_revision}
 polymer_entities{ rcsb_id
  entity_poly{type pdbx_seq_one_letter_code_can}
  rcsb_polymer_entity{pdbx_description pdbx_mutation details src_method}
  rcsb_entity_source_organism{ncbi_taxonomy_id ncbi_scientific_name source_type provenance_source}
  rcsb_polymer_entity_container_identifiers{entity_id asym_ids auth_asym_ids}
  polymer_entity_instances{ rcsb_polymer_entity_instance_container_identifiers{asym_id auth_asym_id}
   rcsb_polymer_instance_feature{type feature_positions{beg_seq_id end_seq_id}} } } } }"""


def graphql(ids):
    out = {}
    todo = sorted(ids)
    for i in range(0, len(todo), 40):
        batch = todo[i:i + 40]
        path = C.CP("rcsb", "graphql", f"{batch[0]}_{batch[-1]}_{len(batch)}.json")
        if not os.path.exists(path):
            body = json.dumps({"query": GQL, "variables": {"ids": batch}}).encode()
            data = C.http("https://data.rcsb.org/graphql", data=body, headers={"Content-Type": "application/json"},
                          timeout=300, min_interval=0.5)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "wb").write(data)
        for e in (json.load(open(path)).get("data") or {}).get("entries") or []:
            if e:
                out[e["rcsb_id"]] = e
        print(f"  graphql {min(i + 40, len(todo))}/{len(todo)}", end="\r")
    print()
    return out


def holdings_removed(ids):
    path = C.CP("rcsb", "holdings_removed.json")
    d = C.jload(path, {})
    for pid in sorted(ids):
        if pid in d:
            continue
        try:
            d[pid] = json.loads(C.http(f"https://data.rcsb.org/rest/v1/holdings/removed/{pid}"))
        except Exception as e:  # 404 = not removed
            d[pid] = {"error": str(e)}
        C.jdump(d, path)
    return d


def row_taxa(accs):
    """Organism taxid of seed-row source accessions (EMBL/GenBank) via NCBI esummary; cached."""
    path = C.CP("ncbi", "row_accession_taxa.json")
    d = C.jload(path, {})
    todo = sorted(a for a in accs if a not in d)
    for i in range(0, len(todo), 150):
        b = todo[i:i + 150]
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id="
               + urllib.parse.quote(",".join(b)))
        r = json.loads(C.http(url, timeout=180))
        res = r.get("result", {})
        found = {}
        for uid in res.get("uids", []):
            x = res[uid]
            found[x.get("accessionversion")] = dict(taxid=str(x.get("taxid", "")), title=x.get("title", ""))
        for a in b:
            d[a] = found.get(a, {"taxid": "", "title": "", "note": "not returned by esummary"})
        C.jdump(d, path)
    return d


def unobs(features, kind):
    s = set()
    for f in features or []:
        if f["type"] == kind:
            for p in f.get("feature_positions") or []:
                if p.get("beg_seq_id") is not None:
                    s.update(range(p["beg_seq_id"], (p.get("end_seq_id") or p["beg_seq_id"]) + 1))
    return s


def main():
    tax = C.Taxonomy()
    pol = C.policy()
    placeholders = {str(x) for x in pol["species"]["placeholders_unresolved"]}
    chains = C.read_tsv(C.P("candidate_chains.tsv"))
    ivs = C.read_tsv(C.P("candidate_intervals.tsv"))
    missing = {r["pdb_id"] for r in chains if r["seqres_status"] == "entry_absent_from_seqres"}
    removed = holdings_removed(missing)
    need = {r["pdb_id"] for r in ivs} | {r["pdb_id"] for r in chains if r["screening"] == "wildcard_only_needs_mmcif"}
    meta = graphql(need)

    ent_rows = []
    for pid in sorted(need | missing):
        e = meta.get(pid)
        if e is None:
            rm = removed.get(pid, {}).get("rcsb_repository_holdings_removed")
            ent_rows.append(dict(pdb_id=pid, status="obsolete" if rm else "metadata_unavailable",
                                 replaced_by=";".join(rm.get("id_codes_replaced_by") or []) if rm else "",
                                 remove_date=(rm or {}).get("remove_date", ""), note=(rm or {}).get("details", "")))
            continue
        acc = e.get("rcsb_accession_info") or {}
        ent_rows.append(dict(pdb_id=pid, status="current", title=(e.get("struct") or {}).get("title", ""),
                             method=";".join(x["method"] for x in e.get("exptl") or []),
                             resolution=";".join(str(x) for x in (e.get("rcsb_entry_info") or {}).get("resolution_combined") or []),
                             models=(e.get("rcsb_entry_info") or {}).get("deposited_model_count", ""),
                             deposit_date=acc.get("deposit_date", "")[:10],
                             initial_release_date=acc.get("initial_release_date", "")[:10],
                             revision_date=acc.get("revision_date", "")[:10],
                             revision=f"{acc.get('major_revision')}.{acc.get('minor_revision')}"))
    C.write_tsv(C.P("entry_metadata.tsv"), ent_rows,
                ["pdb_id", "status", "replaced_by", "remove_date", "note", "title", "method", "resolution", "models",
                 "deposit_date", "initial_release_date", "revision_date", "revision"])

    # ---------------- per interval: entity, species, coverage (screening level)
    rows_by_fam = collections.defaultdict(dict)
    import gzip
    with gzip.open(os.path.join(C.STEP1, "seed_rows.jsonl.gz"), "rt") as f:
        for line in f:
            r = json.loads(line)
            rows_by_fam[r["family"]][r["seed_row_id"]] = r
    accs = set()
    for iv in ivs:
        if iv["screening_status"] == "ambiguous_rows":
            for rid in iv["matching_rows"].split(";"):
                if not rid.startswith("URS"):
                    accs.add(rows_by_fam[iv["family"]][rid]["source_accession"])
    rtax = row_taxa(accs)

    out = []
    for iv in ivs:
        e = meta.get(iv["pdb_id"])
        r = dict(iv)
        r.pop("deposited_sequence", None)
        if e is None:
            r.update(screen_status="entry_metadata_unavailable")
            out.append(r)
            continue
        ent = inst = None
        for pe in e.get("polymer_entities") or []:
            for pi in pe.get("polymer_entity_instances") or []:
                ids = pi["rcsb_polymer_entity_instance_container_identifiers"]
                if ids["auth_asym_id"] == iv["auth_chain"]:
                    ent, inst = pe, pi
        if ent is None:
            r.update(screen_status="auth_chain_not_found_in_metadata")
            out.append(r)
            continue
        can = re.sub(r"\s", "", (ent.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can") or "").upper()
        orgs = ent.get("rcsb_entity_source_organism") or []
        sp = {}
        for o in orgs:
            t = str(o.get("ncbi_taxonomy_id") or "")
            if t and t not in placeholders:
                s = tax.species(t)
                sp[t] = s
        species = {v[0] for v in sp.values() if v[0]}
        pe_info = ent.get("rcsb_polymer_entity") or {}
        lo, hi = int(iv["seqres_start"]), int(iv["seqres_end"])
        un_res = unobs(inst.get("rcsb_polymer_instance_feature"), "UNOBSERVED_RESIDUE_XYZ")
        un_atom = unobs(inst.get("rcsb_polymer_instance_feature"), "UNOBSERVED_ATOM_XYZ")
        n = hi - lo + 1
        obs = sum(1 for i in range(lo, hi + 1) if i not in un_res)
        full = sum(1 for i in range(lo, hi + 1) if i not in un_res and i not in un_atom)
        st_sp = (species.pop() if len(species) == 1 else "")
        sp_note = ("" if st_sp else "multiple_source_species" if len({v[0] for v in sp.values() if v[0]}) > 1
                   else "no_informative_taxon")
        if not st_sp:
            text = f"{pe_info.get('pdbx_description') or ''} ; {(e.get('struct') or {}).get('title', '')}"
            b = {tax.species(t)[0] for t in tax.binomials(text)} - {None}
            if len(b) == 1:
                st_sp, sp_note = b.pop(), "description_binomial"
        r.update(screen_status="ok", entity_id=ent["rcsb_polymer_entity_container_identifiers"]["entity_id"],
                 label_asym_id=inst["rcsb_polymer_entity_instance_container_identifiers"]["asym_id"],
                 seqres_equals_entity_can=can == iv["deposited_sequence"].upper(),
                 source_taxids=";".join(str(o.get("ncbi_taxonomy_id") or "") for o in orgs),
                 source_names=";".join(o.get("ncbi_scientific_name") or "" for o in orgs),
                 source_types=";".join(o.get("source_type") or "" for o in orgs),
                 src_method=pe_info.get("src_method") or "", rna_description=pe_info.get("pdbx_description") or "",
                 pdbx_mutation=pe_info.get("pdbx_mutation") or "", entity_details=pe_info.get("details") or "",
                 screen_species_taxid=st_sp, screen_species_name=tax.name.get(st_sp, "") if st_sp else "",
                 screen_species_note=sp_note or "rcsb_source_organism",
                 superkingdom=tax.superkingdom(st_sp) if st_sp else "",
                 screen_observed_fraction=round(obs / n, 4), screen_complete_fraction=round(full / n, 4),
                 method=";".join(x["method"] for x in e.get("exptl") or []),
                 resolution=min((e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None], key=lambda x: x or 99))
        # row identity by organism for ambiguous intervals
        if iv["screening_status"] == "ambiguous_rows" and st_sp:
            cand = []
            for rid in iv["matching_rows"].split(";"):
                m = re.match(r"URS[0-9A-F]+_(\d+)/", rid)
                t = m.group(1) if m else rtax.get(rows_by_fam[iv["family"]][rid]["source_accession"], {}).get("taxid", "")
                if t and tax.species(t)[0] == st_sp:
                    cand.append(rid)
            r["rows_same_species"] = ";".join(cand)
            strings = {rows_by_fam[iv["family"]][x]["ungapped_sequence"].upper() for x in cand}
            if len(cand) == 1:
                r.update(screening_status="resolved_by_row_organism", seed_row_id=cand[0],
                         row_basis="only matching row whose source organism is the chain's species")
        out.append(r)
    C.write_tsv(C.P("interval_screening.tsv"), out)

    # ---------------- early yield check
    ok = {"exact_unique", "resolved_by_explicit_tag", "resolved_by_row_organism"}
    fam_disc = {r["family"]: r for r in C.read_tsv(C.P("family_discovery.tsv"))}
    by = collections.defaultdict(list)
    for r in out:
        by[r["family"]].append(r)
    yrows = []
    partner = {}
    for fam in sorted(fam_disc):
        L = by.get(fam, [])
        good = [r for r in L if r["screening_status"] in ok and r.get("screen_status") == "ok"]
        seqs = collections.defaultdict(list)
        for r in good:
            seqs[r["interval_sequence"]].append(r)
        known_diff, unresolved_only = 0, 0
        for a, b in __import__("itertools").combinations(sorted(seqs), 2):
            kd = any(x["pdb_id"] != y["pdb_id"] and x["screen_species_taxid"] and y["screen_species_taxid"]
                     and x["screen_species_taxid"] != y["screen_species_taxid"] for x in seqs[a] for y in seqs[b])
            un = any(x["pdb_id"] != y["pdb_id"] and (not x["screen_species_taxid"] or not y["screen_species_taxid"])
                     for x in seqs[a] for y in seqs[b])
            known_diff += kd
            unresolved_only += (not kd and un)
        for r in good:
            partner[r["candidate_id"]] = any(
                s != r["interval_sequence"] and any(
                    y["pdb_id"] != r["pdb_id"] and not (y["screen_species_taxid"] and r["screen_species_taxid"]
                                                        and y["screen_species_taxid"] == r["screen_species_taxid"])
                    for y in ys) for s, ys in seqs.items())
        yrows.append(dict(family=fam, family_id=fam_disc[fam]["family_id"], discovery_status=fam_disc[fam]["discovery_status"],
                          candidate_entries=fam_disc[fam]["candidate_entries"], candidate_chains=fam_disc[fam]["candidate_chains"],
                          exact_unambiguous_intervals=len(good),
                          ambiguous_intervals=sum(1 for r in L if r["screening_status"] not in ok),
                          distinct_matched_sequences=len(seqs),
                          seq_pairs_known_different_species=known_diff,
                          seq_pairs_species_unresolved=unresolved_only,
                          different_sequence_pair_possible="yes" if (known_diff or unresolved_only) else "no",
                          intervals_with_possible_partner=sum(1 for r in good if partner[r["candidate_id"]])))
    C.write_tsv(C.P("yield_by_family.tsv"), yrows)
    C.jdump(partner, C.P("possible_partner.json"))
    y = [r for r in yrows if r["different_sequence_pair_possible"] == "yes"]
    print("entries", len(ent_rows), collections.Counter(r["status"] for r in ent_rows))
    print("interval screening", collections.Counter(r["screening_status"] for r in out))
    print("families with exact unambiguous match", sum(1 for r in yrows if int(r["exact_unambiguous_intervals"])),
          "| different-sequence pair possible", len(y),
          "| seq pairs known diff species", sum(r["seq_pairs_known_different_species"] for r in yrows),
          "| unresolved", sum(r["seq_pairs_species_unresolved"] for r in yrows),
          "| intervals w/ partner", sum(r["intervals_with_possible_partner"] for r in yrows))


if __name__ == "__main__":
    main()
