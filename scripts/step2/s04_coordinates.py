"""Step 2 / stage 4: coordinate review from the experimental mmCIF (source record), once per PDB entry.

entry_data(PDB) downloads cache/mmcif/<PDB>.cif.gz once (validated: gzip, parses, entry ID), and caches a parsed
summary cache/mmcif_parsed/<PDB>.json: provenance (method, resolution, revision history), every polymer entity with
its RNA-source categories (_entity_src_nat / _entity_src_gen gene source / _pdbx_entity_src_syn; host kept apart),
pdbx_mutation and description, and for every nucleic-acid instance the full deposited polymer (pdbx_poly_seq_scheme:
seq_id, mon_id, parent, author number, insertion code, observed) with model-1 / first-altloc atom names.

interval_review(...) gives the per-candidate-interval facts used for chain review: mmCIF sequence confirmation of the
seed-row match, species from mmCIF, construct documentation, coverage and technical status.

main(): (1) prefetch entries most likely needed by stage 5 (per sequence-species group, all structures tied at the
best screening coverage), in parallel with stage 3; (2) re-screen wildcard-only chains with parent-resolved
sequences (origin_queue_extra.tsv); (3) the 10-candidate validation sample (validation_sample.tsv).
"""
import collections
import concurrent.futures as cf
import gzip
import json
import os
import random
import re
import shutil
import sys
import threading

import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

BACKBONE = ("P", "O5'", "C5'", "C4'", "C3'", "O3'")
DNA_PARENT = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U"}
LOCKS = collections.defaultdict(threading.Lock)


def s(v):
    return gemmi.cif.as_string(v) if v not in (None, "?", ".") else ""


def table(block, prefix, cols):
    t = block.find_mmcif_category(prefix)
    if len(t.tags) == 0:
        return []
    tags = [x[len(prefix):] for x in t.tags]
    idx = {c: tags.index(c) for c in cols if c in tags}
    return [{c: (s(r[idx[c]]) if c in idx else "") for c in cols} for r in t]


def mmcif_path(pdb):
    return C.CP("mmcif", f"{pdb}.cif.gz")


def download(pdb):
    path = mmcif_path(pdb)
    if os.path.exists(path):
        return path
    for d in ("../v3_full Rfam/rfam_phase2_3_verification_15_1/cache/mmcif", "../rfam_star3d_pilot_v2/inputs/structures/mmcif"):
        src = os.path.normpath(os.path.join(C.ROOT, d, f"{pdb}.cif.gz"))
        if os.path.exists(src):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = os.path.join(os.path.dirname(path), "tmp_" + os.path.basename(path))
            shutil.copy2(src, tmp)
            if gemmi.cif.read(tmp).sole_block().name.upper() == pdb:
                os.replace(tmp, path)
                C.log_source(f"mmcif:{pdb}", f"https://files.rcsb.org/download/{pdb}.cif.gz", path,
                             "copied_from_earlier_SQARMA_cache", os.path.relpath(src, C.ROOT))
                return path
    url = f"https://files.rcsb.org/download/{pdb}.cif.gz"

    def val(p):
        C.gz_ok(p)
        if gemmi.cif.read(p).sole_block().name.upper() != pdb:
            raise ValueError("entry id mismatch")
    return C.fetch_cached(url, path, f"mmcif:{pdb}", val)


def parse_entry(pdb):
    path = download(pdb)
    block = gemmi.cif.read(path).sole_block()
    prov = dict(pdb_id=pdb, rcsb_page=f"https://www.rcsb.org/structure/{pdb}",
                download_url=f"https://files.rcsb.org/download/{pdb}.cif.gz", mmcif=os.path.relpath(path, C.ROOT),
                mmcif_sha256=C.sha256(path), mmcif_bytes=os.path.getsize(path),
                title=s(block.find_value("_struct.title")),
                method=";".join(r["method"] for r in table(block, "_exptl.", ["method"])))
    res = [s(block.find_value(k)) for k in ("_refine.ls_d_res_high", "_em_3d_reconstruction.resolution",
                                            "_reflns.d_resolution_high")]
    res = [float(x) for x in res if re.match(r"^[0-9.]+$", x or "")]
    prov["resolution"] = min(res) if res else ""
    rev = table(block, "_pdbx_audit_revision_history.", ["ordinal", "major_revision", "minor_revision", "revision_date"])
    if rev:
        last = max(rev, key=lambda r: int(r["ordinal"] or 0))
        prov["revision"] = f"{last['major_revision']}.{last['minor_revision']}"
        prov["revision_date"] = last["revision_date"]
    prov["deposit_date"] = s(block.find_value("_pdbx_database_status.recvd_initial_deposition_date"))
    ents = {}
    for r in table(block, "_entity.", ["id", "type", "src_method", "pdbx_description", "pdbx_mutation", "details",
                                       "pdbx_fragment"]):
        ents[r["id"]] = dict(r, sources=[], hosts=[])
    for r in table(block, "_entity_poly.", ["entity_id", "type", "pdbx_seq_one_letter_code_can"]):
        if r["entity_id"] in ents:
            ents[r["entity_id"]]["poly_type"] = r["type"]
            ents[r["entity_id"]]["can"] = re.sub(r"\s", "", r["pdbx_seq_one_letter_code_can"])
    for r in table(block, "_entity_src_nat.", ["entity_id", "pdbx_ncbi_taxonomy_id", "pdbx_organism_scientific", "strain"]):
        ents.get(r["entity_id"], {}).setdefault("sources", []).append(
            dict(category="entity_src_nat", taxid=r["pdbx_ncbi_taxonomy_id"], name=r["pdbx_organism_scientific"], strain=r["strain"]))
    for r in table(block, "_entity_src_gen.", ["entity_id", "pdbx_gene_src_ncbi_taxonomy_id", "pdbx_gene_src_scientific_name",
                                               "pdbx_host_org_ncbi_taxonomy_id", "pdbx_host_org_scientific_name",
                                               "pdbx_gene_src_gene"]):
        e = ents.get(r["entity_id"], {})
        e.setdefault("sources", []).append(dict(category="entity_src_gen", taxid=r["pdbx_gene_src_ncbi_taxonomy_id"],
                                                name=r["pdbx_gene_src_scientific_name"], strain=r["pdbx_gene_src_gene"]))
        e.setdefault("hosts", []).append(dict(taxid=r["pdbx_host_org_ncbi_taxonomy_id"], name=r["pdbx_host_org_scientific_name"]))
    for r in table(block, "_pdbx_entity_src_syn.", ["entity_id", "ncbi_taxonomy_id", "organism_scientific", "details"]):
        ents.get(r["entity_id"], {}).setdefault("sources", []).append(
            dict(category="pdbx_entity_src_syn", taxid=r["ncbi_taxonomy_id"], name=r["organism_scientific"], strain=r["details"]))
    parents = {r["id"]: (r["mon_nstd_parent_comp_id"], r["type"]) for r in
               table(block, "_chem_comp.", ["id", "mon_nstd_parent_comp_id", "type"])}
    scheme = collections.defaultdict(list)
    for r in table(block, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "pdb_seq_num",
                                                     "pdb_ins_code", "pdb_strand_id", "pdb_mon_id", "hetero"]):
        scheme[r["asym_id"]].append(r)
    st = gemmi.read_structure(path)
    st.setup_entities()
    n_models = len(st)
    st.remove_alternative_conformations()
    atoms = {}
    for ch in st[0]:
        for rr in ch:
            if rr.label_seq is not None:
                atoms.setdefault((rr.subchain, rr.label_seq), []).append(
                    dict(name=rr.name, num=rr.seqid.num, icode=rr.seqid.icode.strip(), het=rr.het_flag,
                         atoms=sorted({a.name for a in rr})))
    chains = {}
    for asym, rows in scheme.items():
        eid = rows[0]["entity_id"]
        pt = ents.get(eid, {}).get("poly_type", "")
        if "ribonucleotide" not in pt and "nucleotide" not in pt:
            continue
        seen = set()
        pos = []
        for r in rows:
            sid = int(r["seq_id"])
            if sid in seen:     # microheterogeneity: keep the first listed monomer
                continue
            seen.add(sid)
            mon = r["mon_id"]
            par = mon if mon in C.STD else DNA_PARENT.get(mon) or (parents.get(mon, ("", ""))[0] or "")
            par = par if par in C.STD or par == "T" else (DNA_PARENT.get(par, "X") if par else "X")
            a = atoms.get((asym, sid), [])
            pos.append(dict(seq_id=sid, mon=mon, parent=par, auth_num=r["pdb_seq_num"], icode=r["pdb_ins_code"],
                            observed=bool(r["pdb_mon_id"]) and bool(a), het=a[0]["het"] if a else "",
                            atoms=a[0]["atoms"] if a else [], model_num=a[0]["num"] if a else None,
                            model_icode=a[0]["icode"] if a else "",
                            parent_source="standard" if mon in C.STD else ("chem_comp" if par != "X" else "")))
        # parents for modified residues without a chem_comp parent: wwPDB canonical code of the entity, position by
        # position (only when its length equals the polymer length), else gemmi's tabulated one-letter code.
        can = ents.get(eid, {}).get("can", "")
        for i, p in enumerate(pos):
            if p["parent"] != "X":
                continue
            c = can[i].upper() if len(can) == len(pos) else ""
            if c in C.STD:
                p["parent"], p["parent_source"] = c, "entity_poly_can"
                continue
            ri = gemmi.find_tabulated_residue(p["mon"])
            if ri is not None and ri.one_letter_code.upper() in C.STD:
                p["parent"], p["parent_source"] = ri.one_letter_code.upper(), "gemmi_table"
        chains[asym] = dict(entity_id=eid, auth_chain=rows[0]["pdb_strand_id"], positions=pos,
                            parent_seq="".join(p["parent"] for p in pos))
    return dict(provenance=dict(prov, models=n_models), entities=ents, chains=chains, parsed_utc=C.now())


def ccd_parent(comp):
    """mon_nstd_parent_comp_id from the wwPDB Chemical Component Dictionary (RCSB), cached per component."""
    path = C.CP("ccd", f"{comp}.json")
    d = C.jload(path)
    if d is None:
        try:
            block = gemmi.cif.read_string(C.http(f"https://files.rcsb.org/ligands/view/{comp}.cif").decode()).sole_block()
            par = s(block.find_value("_chem_comp.mon_nstd_parent_comp_id"))
            d = dict(comp=comp, parent=par, type=s(block.find_value("_chem_comp.type")), retrieved=C.now())
        except Exception as e:
            return None          # not cached: retried later
        C.jdump(d, path)
    par = (d.get("parent") or "").split(",")[0].strip()
    return par if par in C.STD else None


def fix_parents(d):
    """Resolve remaining X parents of modified residues through the CCD (unknown 'N'/'UNK' stay unknown)."""
    changed = False
    for ch in d["chains"].values():
        for p in ch["positions"]:
            if p["parent"] == "X" and p["mon"] not in ("N", "UNK", "X"):
                par = ccd_parent(p["mon"])
                if par:
                    p["parent"], p["parent_source"] = par, "wwPDB_CCD"
                    changed = True
        if changed:
            ch["parent_seq"] = "".join(p["parent"] for p in ch["positions"])
    return changed


def entry_data(pdb):
    pdb = pdb.upper()
    path = C.CP("mmcif_parsed", f"{pdb}.json.gz")
    with LOCKS[pdb]:
        if os.path.exists(path):
            d = json.load(gzip.open(path, "rt"))
            if not d.get("ccd_checked") and fix_parents(d) is not None:
                d["ccd_checked"] = True
                tmp = os.path.join(os.path.dirname(path), "tmp_" + os.path.basename(path))
                with gzip.open(tmp, "wt") as f:
                    json.dump(d, f)
                os.replace(tmp, path)
            return d
        d = parse_entry(pdb)
        fix_parents(d)
        d["ccd_checked"] = True
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = os.path.join(os.path.dirname(path), "tmp_" + os.path.basename(path))
        with gzip.open(tmp, "wt") as f:
            json.dump(d, f)
        os.replace(tmp, path)
        return d


# ------------------------------------------------------------------ species and documentation (mmCIF)
ENG = re.compile(r"\b(mutant|mutation|mutated|engineered|chimer\w*|fusion|fused|linker|circular(ly)? permut\w*|"
                 r"tagged|tag|T7 leader|aptamer[- ]?tag|stabiliz\w+ (loop|mutation)|replaced|substitut\w+)\b", re.I)


def species_from_mmcif(ent, tax, placeholders, title):
    srcs = []
    for x in ent.get("sources", []):
        if x.get("taxid"):
            srcs.append(x)
        elif x.get("name") and tax.sci2id.get(x["name"].strip().lower()):
            # organism given by scientific name only: exact NCBI scientific-name match in the pinned snapshot
            srcs.append(dict(x, taxid=tax.sci2id[x["name"].strip().lower()], category=x["category"] + "_by_name"))
    informative = [x for x in srcs if x["taxid"].split(",")[0].strip() not in placeholders
                   and not tax.is_placeholder(x["taxid"].split(",")[0].strip())]
    sp = {}
    for x in informative:
        for t in x["taxid"].split(","):
            t = t.strip()
            if t:
                sp[t] = tax.species(t)
    species = {v[0] for v in sp.values() if v[0]}
    if len(species) == 1:
        t = species.pop()
        return t, "mmcif_" + "+".join(sorted({x["category"] for x in informative})), sorted(sp)
    if len(species) > 1:
        return "", "multiple_source_species", sorted(sp)
    b = {tax.species(t)[0] for t in tax.binomials(f"{ent.get('pdbx_description', '')} ; {title}")} - {None}
    if len(b) == 1:
        return b.pop(), "description_binomial", []
    return "", ("placeholder_taxon_only" if srcs else "no_source_annotation"), sorted(x["taxid"] for x in srcs)


def documentation(ent):
    mut = ent.get("pdbx_mutation", "")
    mut_real = bool(mut) and not re.match(r"^(no|none|n/?a|wild[- ]?type|wt)$", mut.strip(), re.I)
    text = " ; ".join(x for x in (ent.get("pdbx_description", ""), ent.get("details", ""), ent.get("pdbx_fragment", "")) if x)
    words = sorted({m.group(0).lower() for m in ENG.finditer(text)})
    return dict(pdbx_mutation=mut, mutation_documented=mut_real, engineering_words=";".join(words),
                documented_engineering=mut_real or bool(words), description=ent.get("pdbx_description", ""),
                details=ent.get("details", ""), src_method=ent.get("src_method", ""))


# ------------------------------------------------------------------ interval review
def find_asym(d, auth):
    hits = [a for a, c in d["chains"].items() if c["auth_chain"] == auth]
    return hits


def interval_review(cand, fam_rows, tax, placeholders):
    """cand: interval_screening row (+ seed_row_id). Returns a dict of mmCIF-verified facts."""
    d = entry_data(cand["pdb_id"])
    out = dict(candidate_id=cand["candidate_id"], mmcif_sha256=d["provenance"]["mmcif_sha256"])
    asyms = find_asym(d, cand["auth_chain"])
    if len(asyms) != 1:
        return dict(out, mmcif_status="auth_chain_not_unique_or_absent", note=f"asyms={asyms}")
    asym = asyms[0]
    ch = d["chains"][asym]
    ent = d["entities"].get(ch["entity_id"], {})
    seq = ch["parent_seq"]
    lo, hi = int(cand["seqres_start"]), int(cand["seqres_end"])
    row = fam_rows[cand["seed_row_id"]]
    ug = row["ungapped_sequence"].upper()
    occ = [m.start() + 1 for m in re.finditer(f"(?={re.escape(ug)})", seq)]
    others = [rid for rid, r in fam_rows.items() if rid != cand["seed_row_id"]
              and set(r["ungapped_sequence"].upper()) <= C.STD and r["ungapped_sequence"].upper() in seq[lo - 1:hi]]
    pos = ch["positions"][lo - 1:hi]
    seq_ok = seq[lo - 1:hi] == ug and occ == [lo] and [p["seq_id"] for p in pos] == list(range(lo, hi + 1))
    n = hi - lo + 1
    obs = [p for p in pos if p["observed"]]
    full = [p for p in obs if all(a in p["atoms"] for a in BACKBONE)]
    nobb = [p for p in obs if not any(a in p["atoms"] for a in BACKBONE)]
    # longest internal unobserved stretch
    longest, cur = 0, 0
    first = next((i for i, p in enumerate(pos) if p["observed"]), None)
    last = max((i for i, p in enumerate(pos) if p["observed"]), default=None)
    for p in pos[first:last + 1] if first is not None else []:
        cur = 0 if p["observed"] else cur + 1
        longest = max(longest, cur)
    mods = [f"{p['mon']}{p['auth_num']}{p['icode'] or ''}" for p in pos if p["mon"] not in C.STD]
    problems = []
    if not obs:
        problems.append("no_observed_residue_in_interval")
    if nobb:
        problems.append(f"{len(nobb)}_observed_residues_without_any_backbone_atom")
    if any(len(p["mon"]) > 3 for p in obs):
        problems.append("residue_name_longer_than_3")
    if any(p["model_num"] is not None and not -999 <= p["model_num"] <= 9999 for p in obs):
        problems.append("author_number_outside_-999..9999")
    keys = [(p["model_num"], p["model_icode"]) for p in obs]
    if len(set(keys)) != len(keys):
        problems.append("author_residue_id_collision")
    if not seq_ok:
        problems.append("mmcif_sequence_does_not_confirm_unique_interval")
    if problems:
        tech = "technical_unusable_or_unresolved"
    elif len(full) == n:
        tech = "ready"
    else:
        tech = "ready_with_coverage_flags"
    sp, sp_basis, sp_taxa = species_from_mmcif(ent, tax, placeholders, d["provenance"]["title"])
    nums = [p["model_num"] for p in obs]
    return dict(out, mmcif_status="ok", label_asym_id=asym, entity_id=ch["entity_id"], auth_chain=ch["auth_chain"],
                mmcif_interval_confirmed=seq_ok, mmcif_row_occurrences=";".join(map(str, occ)),
                other_rows_inside_interval=";".join(others[:10]),
                deposited_parent_sequence=seq, deposited_length=len(seq),
                deposited_parent_sha256=C.hashlib.sha256(seq.encode()).hexdigest(),
                deposited_equals_seqres=seq == cand.get("deposited_sequence", seq),
                cropped_5prime=lo - 1, cropped_3prime=len(seq) - hi, interval_length=n,
                observed=len(obs), observed_fraction=round(len(obs) / n, 4),
                backbone_complete=len(full), backbone_complete_fraction=round(len(full) / n, 4),
                residues_missing_backbone_atoms=";".join(f"{p['mon']}{p['auth_num']}{p['icode'] or ''}"
                                                         for p in obs if p not in full)[:400],
                unobserved_residues=";".join(f"{p['mon']}{p['auth_num']}{p['icode'] or ''}" for p in pos if not p["observed"])[:400],
                longest_internal_missing=longest, modified_residues=";".join(mods),
                insertion_codes=";".join(sorted({f"{p['model_num']}{p['model_icode']}" for p in obs if p["model_icode"]})),
                author_numbers_monotonic=all(a < b for a, b in zip(nums, nums[1:])) if nums else "",
                author_number_range=f"{min(nums)}..{max(nums)}" if nums else "",
                technical_status=tech, technical_problems=";".join(problems),
                species_taxid=sp, species_name=tax.name.get(sp, "") if sp else "", species_basis=sp_basis,
                species_source_taxa=";".join(sp_taxa),
                rna_sources=json.dumps(ent.get("sources", [])), expression_hosts=json.dumps(ent.get("hosts", [])),
                method=d["provenance"]["method"], resolution=d["provenance"]["resolution"],
                **{f"doc_{k}": v for k, v in documentation(ent).items()})


# ------------------------------------------------------------------ main: prefetch, wildcard re-screen, sample
def load_rows():
    rows = collections.defaultdict(dict)
    with gzip.open(os.path.join(C.STEP1, "seed_rows.jsonl.gz"), "rt") as f:
        for line in f:
            r = json.loads(line)
            rows[r["family"]][r["seed_row_id"]] = r
    return rows


def main():
    import s01_discover as D
    iv = C.read_tsv(C.P("interval_screening.tsv"))
    partner = C.jload(C.P("possible_partner.json"), {})
    q = [r for r in iv if partner.get(r["candidate_id"])]
    groups = collections.defaultdict(list)
    for r in q:
        groups[(r["family"], r["interval_sequence"], r["screen_species_taxid"])].append(r)
    first = set()
    for g in groups.values():
        best = max(float(r["screen_observed_fraction"] or 0) for r in g)
        first |= {r["pdb_id"] for r in g if float(r["screen_observed_fraction"] or 0) == best}
    chains = C.read_tsv(C.P("candidate_chains.tsv"))
    yields = {r["family"]: r for r in C.read_tsv(C.P("yield_by_family.tsv"))}
    wild = [r for r in chains if r["screening"] == "wildcard_only_needs_mmcif" and "deposited_unknown=-" not in r["diagnostic"]]
    order = sorted(first) + sorted({r["pdb_id"] for r in wild} - first)
    cached = [p for p in order if os.path.exists(C.CP("mmcif_parsed", f"{p}.json.gz"))]
    print(f"prefetch: {len(order)} entries ({len(first)} tied-best in {len(groups)} groups, {len(wild)} wildcard chains);"
          f" already parsed {len(cached)}", flush=True)
    import time
    t0, done = time.time(), [0]
    lock = threading.Lock()

    def job(p):
        try:
            entry_data(p)
            st = "ok"
        except Exception as e:
            st = f"FAILED {e!r}"[:200]
        with lock:
            done[0] += 1
            if done[0] % 10 == 0 or st != "ok":
                print(f"[{done[0]}/{len(order)} {(time.time() - t0) / 60:.1f} min] {p} {st}", flush=True)
        return p, st
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        status = dict(ex.map(job, order))
    C.jdump(status, C.P("prefetch_status.json"))

    # wildcard chains: re-screen with parent-resolved deposited sequence
    rows = load_rows()
    extra = []
    tax = C.Taxonomy()
    for r in wild:
        if status.get(r["pdb_id"]) != "ok":
            continue
        d = entry_data(r["pdb_id"])
        asyms = find_asym(d, r["auth_chain"])
        if len(asyms) != 1:
            continue
        seq = d["chains"][asyms[0]]["parent_seq"]
        trows = set(r["explicit_tag_rows"].split(";")) - {""}
        ivs, diag = D.match_chain(seq, rows[r["family"]], trows)
        for x in ivs:
            cid = f"{r['family']}__{r['pdb_id']}_{r['auth_chain']}__{x['start']}-{x['end']}"
            extra.append(dict(candidate_id=cid, family=r["family"], pdb_id=r["pdb_id"], auth_chain=r["auth_chain"],
                              seqres_start=x["start"], seqres_end=x["end"], screening_status=x["status"],
                              row_basis=x["basis"] + " (mmCIF parent-resolved sequence)", seed_row_id=x["chosen_row"],
                              matching_rows=";".join(x["rows"]), interval_sequence=seq[x["start"] - 1:x["end"]],
                              deposited_sequence=seq, deposited_seq_sha256=C.hashlib.sha256(seq.encode()).hexdigest(),
                              sources=r["sources"], origin="wildcard_rescreen"))
    C.write_tsv(C.P("wildcard_rescreen.tsv"), extra)
    print("wildcard re-screen: new intervals", len(extra), collections.Counter(x["screening_status"] for x in extra))

    # validation sample: mapping-only candidates
    rnd = random.Random(20261113)
    maponly = [r for r in iv if r["sources"] == "rfam_mapping"]
    fams = sorted({r["family"] for r in maponly})
    rnd.shuffle(fams)
    sample = []
    for f in fams:
        rs = sorted((r for r in maponly if r["family"] == f), key=lambda r: r["candidate_id"])
        sample.append(rnd.choice(rs))
        if len(sample) == 8:
            break
    amb = sorted((r for r in maponly if r["screening_status"].startswith("ambiguous")), key=lambda r: r["candidate_id"])
    sample.append(rnd.choice(amb))
    rej = sorted((r for r in chains if r["sources"] == "rfam_mapping" and r["screening"] == "no_exact_seed_match"
                  and int(r["seqres_length"] or 0) < 400), key=lambda r: (r["family"], r["pdb_id"]))
    rj = rnd.choice(rej)
    vs = []
    for r in sample:
        d = entry_data(r["pdb_id"])
        a = find_asym(d, r["auth_chain"])
        seq = d["chains"][a[0]]["parent_seq"] if len(a) == 1 else ""
        ids = r["matching_rows"].split(";")
        occ = {rid: [m.start() + 1 for m in re.finditer(f"(?={rows[r['family']][rid]['ungapped_sequence'].upper()})", seq)]
               for rid in ids}
        lo, hi = int(r["seqres_start"]), int(r["seqres_end"])
        if r["screening_status"].startswith("ambiguous"):
            # reproduced = the same row IDs occur in the mmCIF polymer, overlapping the screened cluster
            reproduced = all(any(lo <= o <= hi or o <= lo <= o + len(rows[r['family']][x]['ungapped_sequence']) - 1
                                 for o in occ[x]) for x in ids) and len(ids) > 1
        else:
            reproduced = occ[ids[0]] == [lo]
        vs.append(dict(kind="ambiguous" if r["screening_status"].startswith("ambiguous") else "mapping_only_match",
                       candidate_id=r["candidate_id"], family=r["family"], pdb_id=r["pdb_id"], auth_chain=r["auth_chain"],
                       screening_status=r["screening_status"], matching_rows=r["matching_rows"],
                       mmcif_occurrences=json.dumps(occ), seqres_equals_mmcif_parent=seq == r.get("deposited_sequence", seq) or "",
                       reproduced=reproduced))
    d = entry_data(rj["pdb_id"])
    a = find_asym(d, rj["auth_chain"])
    seq = d["chains"][a[0]]["parent_seq"] if len(a) == 1 else ""
    hit = [rid for rid, x in rows[rj["family"]].items() if set(x["ungapped_sequence"].upper()) <= C.STD
           and x["ungapped_sequence"].upper() in seq]
    vs.append(dict(kind="rejected_no_exact_match", candidate_id=f"{rj['family']}__{rj['pdb_id']}_{rj['auth_chain']}",
                   family=rj["family"], pdb_id=rj["pdb_id"], auth_chain=rj["auth_chain"], screening_status=rj["screening"],
                   matching_rows=";".join(hit), mmcif_occurrences="", reproduced=not hit))
    C.write_tsv(C.P("validation_sample.tsv"), vs)
    print("validation sample:", [(v["candidate_id"], v["reproduced"]) for v in vs])


if __name__ == "__main__":
    main()
