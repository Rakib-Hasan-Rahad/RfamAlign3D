"""Step 2 / stage 3: independent natural-sequence evidence, ONCE per (full deposited sequence, species).

Routes (policy origin + revision R1):
  A. Bacteria/Archaea: NCBI Datasets genome assembly for the entity's own taxon (strain) if one exists, otherwise
     the species reference genome; downloaded once per species and searched locally on both strands.
  B. Everything else, and A-misses: NCBI BLAST URL API (blastn, core_nt, txid<species>[Organism:exp]), one batched
     job per species; each full-length identical hit is confirmed by exact string search of the efetched record
     region, and the record must pass the independence rule (esummary).
For every key the family interval(s) are searched in the same records, so a flank difference at a locus can be
recognised. Evidence is cached per key in cache/origin/<sha>_<species>.json; reruns skip finished keys.

Outputs: sequence_evidence.tsv (one row per key); queue progress printed.
"""
import collections
import concurrent.futures as cf
import gzip
import io
import json
import os
import re
import sys
import socket
import threading
import time
import urllib.parse
import xml.etree.ElementTree as ET
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
BLAST = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
LOCK = threading.Lock()
TAX = [None]
COMP = str.maketrans("ACGUTN", "UGCAAN")


def rc(s):
    return s.translate(COMP)[::-1]


def to_rna(s):
    return s.upper().replace("T", "U")


def key_path(seq_sha, sp):
    return C.CP("origin", f"{seq_sha[:20]}_{sp}.json")


# ------------------------------------------------------------------ queue
def build_queue():
    iv = C.read_tsv(C.P("interval_screening.tsv"))
    partner = C.jload(C.P("possible_partner.json"), {})
    extra = C.read_tsv(C.P("origin_queue_extra.tsv")) if os.path.exists(C.P("origin_queue_extra.tsv")) else []
    q = {}
    dep = {r["deposited_seq_sha256"]: r for r in C.read_tsv(C.P("candidate_intervals.tsv"))}
    for r in iv + extra:
        if not partner.get(r["candidate_id"]) and r not in extra:
            continue
        sp = r.get("species_taxid") or r.get("screen_species_taxid")
        if not sp or TAX[0].is_placeholder(sp):
            continue
        seq = r.get("deposited_sequence") or dep[r["deposited_seq_sha256"]]["deposited_sequence"]
        k = (r["deposited_seq_sha256"], sp)
        e = q.setdefault(k, dict(seq_sha=k[0], species=sp, seq=seq, intervals=set(), strains=set(), families=set(),
                                 candidates=set()))
        e["intervals"].add(r["interval_sequence"])
        e["families"].add(r["family"])
        e["candidates"].add(r["candidate_id"])
        for t in (r.get("source_taxids") or "").split(";"):
            if t.strip():
                e["strains"].add(t.strip())
    return q


# ------------------------------------------------------------------ route A: genome assemblies
def assemblies_for(taxa, species, limit=8):
    """R6: strain taxa first (reference, then complete), then the species reference, then species complete
    genomes and chromosomes (stable accession order); at most `limit` distinct assemblies."""
    out = []
    for t in list(taxa) + [species]:
        for flt in ("filters.reference_only=true", "filters.assembly_level=complete_genome",
                    "filters.assembly_level=chromosome"):
            url = f"{DS}genome/taxon/{t}/dataset_report?{flt}&filters.assembly_version=current&page_size={limit}"
            try:
                d = json.loads(C.http(url, timeout=120))
            except Exception:
                continue
            reps = sorted(d.get("reports") or [], key=lambda r: (not r["accession"].startswith("GCF_"), r["accession"]))
            for r in reps:
                if r["accession"] not in [x[0] for x in out]:
                    out.append((r["accession"], r.get("organism", {}).get("tax_id"), f"taxon {t} {flt}"))
            if len(out) >= limit:
                return out[:limit]
    return out[:limit]


def assembly_for(taxid):
    for flt in ("filters.reference_only=true", "filters.assembly_level=complete_genome",
                "filters.assembly_level=chromosome"):
        url = f"{DS}genome/taxon/{taxid}/dataset_report?{flt}&filters.assembly_version=current&page_size=5"
        try:
            d = json.loads(C.http(url, timeout=120))
        except Exception:
            continue
        reps = d.get("reports") or []
        reps.sort(key=lambda r: (not r["accession"].startswith("GCF_"), r["accession"]))
        if reps:
            return reps[0]["accession"], reps[0].get("organism", {}).get("tax_id"), flt
    return None, None, None


def genome_records(acc):
    path = C.CP("genomes", f"{acc}.fna.gz")
    if not os.path.exists(path):
        data = C.http(f"{DS}genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA", timeout=900)
        z = zipfile.ZipFile(io.BytesIO(data))
        fna = [n for n in z.namelist() if n.endswith(".fna")]
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with gzip.open(path + ".part", "wb") as out:
            for n in fna:
                out.write(z.read(n))
        os.replace(path + ".part", path)
        C.log_source(f"genome:{acc}", f"{DS}genome/accession/{acc}/download", path, "downloaded")
    recs, cur, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    recs[cur] = to_rna("".join(buf))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        recs[cur] = to_rna("".join(buf))
    return recs, C.sha256(path)


def search_records(recs, s):
    hits = []
    for name, seq in recs.items():
        for strand, q in (("+", s), ("-", rc(s))):
            i = seq.find(q)
            while i != -1 and len(hits) < 20:
                hits.append(dict(record=name, start=i + 1, end=i + len(q), strand=strand))
                i = seq.find(q, i + 1)
    return hits


# ------------------------------------------------------------------ route B: BLAST + efetch confirmation
def blast(species, queries):
    """queries: {qid: seq}. Returns parsed hits {qid: [(acc, start, end, strand)]} (full-length identical only)."""
    tag = C.hashlib.sha256(json.dumps(sorted(queries.items())).encode()).hexdigest()[:16]
    path = C.CP("blast", f"{species}_{tag}.xml")
    if not os.path.exists(path):
        fasta = "".join(f">{k}\n{v.replace('U', 'T')}\n" for k, v in sorted(queries.items()))
        short = min(len(v) for v in queries.values()) < 30
        body = urllib.parse.urlencode(dict(CMD="Put", PROGRAM="blastn", DATABASE="core_nt", QUERY=fasta,
                                           ENTREZ_QUERY=f"txid{species}[Organism:exp]", WORD_SIZE=7 if short else 11,
                                           EXPECT=1000 if short else 10, HITLIST_SIZE=50, FORMAT_TYPE="XML",
                                           FILTER="F", MEGABLAST="off")).encode()
        ridf = path[:-4] + ".rid.json"           # resume: reuse a submitted job (NCBI keeps RIDs ~24 h)
        rj = C.jload(ridf)
        if rj and time.time() - rj["t"] < 20 * 3600:
            rid = rj["rid"]
        else:
            r = C.http(BLAST, data=body, timeout=180, min_interval=11).decode()
            rid = re.search(r"RID = (\S+)", r).group(1)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            C.jdump({"rid": rid, "t": time.time(), "species": species, "n_queries": len(queries)}, ridf)
        t0 = time.time()
        while True:
            time.sleep(45)
            s = C.http(f"{BLAST}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}", timeout=120, min_interval=11).decode()
            if "Status=READY" in s:
                break
            if "Status=FAILED" in s or "Status=UNKNOWN" in s:
                if os.path.exists(ridf):
                    os.remove(ridf)       # expired/failed job: resubmit next run
                raise RuntimeError(f"BLAST {rid} {re.search(r'Status=(\\w+)', s)}")
            if time.time() - t0 > 5400:   # still WAITING: keep the RID (queue position) for the next run / poll_blast.py
                raise RuntimeError(f"BLAST {rid} still waiting after 90 min; RID kept for resume")
        x = C.http(f"{BLAST}?CMD=Get&FORMAT_TYPE=XML&RID={rid}", timeout=600, min_interval=11)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        open(tmp, "wb").write(x)
        ET.parse(tmp)                     # must be complete XML before it is cached
        os.replace(tmp, path)
        C.log_source(f"blast:{species}", f"{BLAST}?RID={rid}", path, "downloaded", f"{len(queries)} queries")
    hits = collections.defaultdict(list)
    root = ET.parse(path).getroot()
    for it in root.iter("Iteration"):
        qid = it.findtext("Iteration_query-def").split()[0]
        qlen = int(it.findtext("Iteration_query-len"))
        for h in it.iter("Hit"):
            acc = h.findtext("Hit_accession")
            for hsp in h.iter("Hsp"):
                if (int(hsp.findtext("Hsp_identity")) == qlen and int(hsp.findtext("Hsp_align-len")) == qlen
                        and int(hsp.findtext("Hsp_query-from")) == 1 and int(hsp.findtext("Hsp_query-to")) == qlen):
                    a, b = int(hsp.findtext("Hsp_hit-from")), int(hsp.findtext("Hsp_hit-to"))
                    hits[qid].append((acc, min(a, b), max(a, b), "+" if a <= b else "-"))
    return hits, path


def esummary(accs):
    path = C.CP("ncbi", "esummary_nuccore.json")
    with LOCK:
        d = C.jload(path, {})
    todo = [a for a in accs if a not in d]
    for i in range(0, len(todo), 100):
        b = todo[i:i + 100]
        r = json.loads(C.http(f"{EU}esummary.fcgi?db=nuccore&retmode=json&id={','.join(b)}", timeout=180))
        res = r.get("result", {})
        got = {}
        for uid in res.get("uids", []):
            x = res[uid]
            got[x.get("accessionversion")] = got[x.get("caption")] = dict(
                accver=x.get("accessionversion"), title=x.get("title"), slen=x.get("slen"), taxid=str(x.get("taxid")),
                sourcedb=x.get("sourcedb"), createdate=x.get("createdate"), biomol=x.get("biomol"))
        with LOCK:
            d = C.jload(path, {})
            for a in b:
                d[a] = got.get(a) or got.get(a.split(".")[0]) or {"missing": True}
            C.jdump(d, path)
    with LOCK:
        d = C.jload(path, {})
    return {a: d.get(a, {"missing": True}) for a in accs}


def efetch_region(acc, start, end, pad=60):
    path = C.CP("ncbi", "regions", f"{acc}_{start}_{end}_{pad}.fasta")
    if not os.path.exists(path):
        s0 = max(1, start - pad)
        t = C.http(f"{EU}efetch.fcgi?db=nuccore&id={acc}&seq_start={s0}&seq_stop={end + pad}&rettype=fasta&retmode=text",
                   timeout=120).decode()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not t.startswith(">"):
            raise RuntimeError(f"efetch returned no FASTA for {acc}")
        open(path + ".tmp", "w").write(t)
        os.replace(path + ".tmp", path)
    lines = open(path).read().splitlines()
    return to_rna("".join(l for l in lines if not l.startswith(">"))), lines[0] if lines else ""


# ------------------------------------------------------------------ evaluate one key
def independent(meta, chain_len, species, tax, deposit_date):
    if meta.get("missing"):
        return False, "esummary_missing"
    if (meta.get("sourcedb") or "").lower() == "pdb":
        return False, "pdb_derived_record"
    t = (meta.get("title") or "").lower()
    if "synthetic construct" in t or re.match(r"^chain ", t):
        return False, "synthetic_or_structure_record"
    if tax.species(meta.get("taxid"))[0] != species:
        return False, f"record_species_{meta.get('taxid')}"
    long_enough = int(meta.get("slen") or 0) >= chain_len + 50
    older = bool(meta.get("createdate") and deposit_date and meta["createdate"][:10].replace("/", "-") < deposit_date)
    if long_enough or older:
        return True, "long_record" if long_enough else "record_predates_structure"
    return False, "short_record_not_older_than_structure"


def assembly_phase(e, tax, res):
    """Route A for Bacteria/Archaea (R6): up to 8 assemblies, searched locally on both strands."""
    for acc, otax, flt in assemblies_for(sorted(e["strains"] - {e["species"]}), e["species"]):
        recs, fsha = genome_records(acc)
        h = search_records(recs, e["seq"])
        res["routes_tried"].append(dict(route="assembly", assembly=acc, assembly_taxid=otax, filter=flt,
                                        fasta_sha256=fsha, whole_chain_hits=len(h)))
        if h:
            res["whole_chain_hits"] += [dict(x, source=f"assembly:{acc}", independent="yes",
                                             exact_confirmed="yes", independence_basis="genome_assembly") for x in h]
            return True
        for sq in sorted(e["intervals"]):
            res["interval_hits"] += [dict(x, interval=sq, source=f"assembly:{acc}", exact_confirmed="yes")
                                     for x in search_records(recs, sq)]
    return False


def blast_phase(e, tax, res, hits, bpath, prefix, queries, deposit_date):
    res["routes_tried"].append(dict(route="blast_core_nt", species=e["species"], xml=os.path.relpath(bpath, C.ROOT),
                                    whole_chain_candidates=len(hits.get(prefix + "chain", []))))
    seq = e["seq"]
    for qid, target in queries.items():
        hl = hits.get(qid, [])
        if not hl:
            continue
        metas = esummary(sorted({h[0] for h in hl}))
        is_chain = qid == prefix + "chain"
        for acc, a, b, strand in hl[:25]:
            m = metas.get(acc, {})
            ok, why = independent(m, len(seq), e["species"], tax, deposit_date)
            if not ok and not is_chain:
                continue
            region, _ = efetch_region(m.get("accver") or acc, a, b)
            exact = (target in region) or (rc(target) in region)
            rec = dict(record=m.get("accver") or acc, start=a, end=b, strand=strand, source="blast_core_nt",
                       exact_confirmed="yes" if exact else "no", independent="yes" if ok else "no",
                       independence_basis=why, title=(m.get("title") or "")[:120], biomol=m.get("biomol", ""))
            if is_chain:
                res["whole_chain_hits"].append(rec)
                if ok and exact:
                    break
            elif ok and exact:
                res["interval_hits"].append(dict(rec, interval=target))
                break


def finalize(res, path):
    good = [h for h in res["whole_chain_hits"] if h.get("independent") == "yes" and h.get("exact_confirmed") == "yes"]
    ivh = [h for h in res["interval_hits"] if h.get("exact_confirmed") == "yes"]
    res["outcome"] = "whole_chain_found" if good else ("interval_found_flank_differs" if ivh else "not_found")
    C.jdump(res, path)
    return res


def new_res(e, tax):
    return dict(seq_sha256=e["seq_sha"], species_taxid=e["species"], species_name=tax.name.get(e["species"], ""),
                deposited_sequence=e["seq"], deposited_length=len(e["seq"]), families=sorted(e["families"]),
                candidates=sorted(e["candidates"]), strains=sorted(e["strains"]), routes_tried=[],
                whole_chain_hits=[], interval_hits=[], evaluated_utc=C.now())


def run_species(sp, es, tax, deposit_dates, report):
    """All keys of one species: assembly route first (prokaryotes), then ONE batched BLAST for the remainder."""
    pending = []
    for e in es:
        res = new_res(e, tax)
        if not set(e["seq"]) <= C.STD:
            res.update(outcome="unknown_residues_in_chain", note="".join(sorted(set(e["seq"]) - C.STD)))
            C.jdump(res, key_path(e["seq_sha"], sp))
            report(e, res)
            continue
        if tax.superkingdom(sp) in ("Bacteria", "Archaea"):
            try:
                if assembly_phase(e, tax, res):
                    report(e, finalize(res, key_path(e["seq_sha"], sp)))
                    continue
            except Exception as ex:
                res["routes_tried"].append(dict(route="assembly", error=repr(ex)[:200]))
        pending.append((e, res))
    if not pending:
        return
    queries, per = {}, {}
    for e, res in pending:
        pre = e["seq_sha"][:10] + "_"
        qs = {pre + "chain": e["seq"]}
        for i, s in enumerate(sorted(e["intervals"])):
            if s != e["seq"]:
                qs[f"{pre}interval{i}"] = s
        per[e["seq_sha"]] = (pre, qs)
        queries.update(qs)
    try:
        hits, bpath = blast(sp, queries)
    except Exception as ex:
        for e, res in pending:
            res.update(outcome="search_failed", note=repr(ex)[:300])
            report(e, res)          # not cached -> retried next run
        return
    for e, res in pending:
        pre, qs = per[e["seq_sha"]]
        dd = min(deposit_dates.get(c.split("__")[1].split("_")[0], "9999") for c in e["candidates"])
        blast_phase(e, tax, res, hits, bpath, pre, qs, dd)
        report(e, finalize(res, key_path(e["seq_sha"], sp)))


def main():
    tax = C.Taxonomy()
    TAX[0] = tax
    q = build_queue()
    dd = {r["pdb_id"]: r.get("deposit_date", "") for r in C.read_tsv(C.P("entry_metadata.tsv"))}
    todo = [e for k, e in sorted(q.items()) if not os.path.exists(key_path(*k))]
    print(f"sequence-species keys: {len(q)}; cached {len(q) - len(todo)}; to evaluate {len(todo)}", flush=True)
    by_sp = collections.defaultdict(list)
    for e in todo:
        by_sp[e["species"]].append(e)
    t0, done = time.time(), [0]

    last = [time.time()]
    socket.setdefaulttimeout(300)          # backstop for any socket without its own timeout

    def watchdog():                        # a silent network stall must not block the run forever
        while True:
            time.sleep(60)
            if time.time() - last[0] > 60 * 60:
                print(f"WATCHDOG: no completed check for 60 min at {C.now()}; exiting so the driver restarts "
                      f"from cache", flush=True)
                os._exit(3)
    threading.Thread(target=watchdog, daemon=True).start()

    def report(e, r):
        with LOCK:
            last[0] = time.time()
            done[0] += 1
            print(f"[{done[0]}/{len(todo)} {(time.time() - t0) / 60:.1f} min] {tax.name.get(e['species'], e['species'])} "
                  f"{e['seq_sha'][:8]} -> {r.get('outcome')}", flush=True)

    # R6 phase 1: local assembly search for every prokaryotic key (fast), before any BLAST queue wait
    pro = [e for e in todo if tax.superkingdom(e["species"]) in ("Bacteria", "Archaea")]
    print(f"phase 1: assembly search for {len(pro)} prokaryotic keys", flush=True)

    def asm_only(e):
        res = new_res(e, tax)
        try:
            if set(e["seq"]) <= C.STD and assembly_phase(e, tax, res):
                report(e, finalize(res, key_path(e["seq_sha"], e["species"])))
        except Exception as ex:
            print(f"  assembly route error {e['seq_sha'][:8]}: {ex!r}"[:200], flush=True)
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        list(ex.map(asm_only, pro))
    rest = [e for e in todo if not os.path.exists(key_path(e["seq_sha"], e["species"]))]
    by_sp = collections.defaultdict(list)
    for e in rest:
        by_sp[e["species"]].append(e)
    print(f"phase 2: {len(rest)} keys in {len(by_sp)} species need BLAST", flush=True)
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(run_species, sp, by_sp[sp], tax, dd, report)
                for sp in sorted(by_sp, key=lambda s: -len(by_sp[s]))]
        for f in futs:
            f.result()
    rows = []
    for k in sorted(q):
        r = C.jload(key_path(*k)) or dict(seq_sha256=k[0], species_taxid=k[1], outcome="search_failed")
        rows.append(dict(seq_sha256=r["seq_sha256"], species_taxid=r["species_taxid"],
                         species_name=tax.name.get(k[1], ""), deposited_length=len(q[k]["seq"]),
                         outcome=r["outcome"], families=";".join(sorted(q[k]["families"])),
                         n_candidates=len(q[k]["candidates"]),
                         whole_chain_evidence=";".join(f"{h['record']}:{h['start']}-{h['end']}{h['strand']}({h.get('source')})"
                                                       for h in r.get("whole_chain_hits", [])
                                                       if h.get("independent") == "yes" and h.get("exact_confirmed", "yes") == "yes")[:500],
                         interval_evidence=";".join(f"{h['record']}:{h['start']}-{h['end']}{h['strand']}"
                                                    for h in r.get("interval_hits", [])[:5]),
                         routes=";".join(x["route"] + ":" + str(x.get("assembly") or x.get("species")) for x in r.get("routes_tried", [])),
                         note=r.get("note", "")))
    C.write_tsv(C.P("sequence_evidence.tsv"), rows)
    print(collections.Counter(r["outcome"] for r in rows))


if __name__ == "__main__":
    main()
