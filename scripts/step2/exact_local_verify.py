#!/usr/bin/env python3
"""Independent exact-match verification of pending Step 2 natural-sequence checks (stdlib only).

  python3 exact_local_verify.py pending_checks.tsv [--workers 4]

Per species, once: download versioned NCBI Datasets records and search the FULL deposited sequence (U->T) exactly
on both strands, plus every family interval (for flank diagnostics only). No BLAST, no remote search.
  Bacteria/Archaea : all current complete-genome/chromosome assemblies of the species (if none: any level), <= 60
  Eukaryota        : the species reference assembly (genome) + its RefSeq RNA transcripts (mature-transcript records)
  Viruses          : all NCBI Virus genome records of the species
Evidence -> evidence/<check_id>.json in the Mac's cache/origin schema, with source=exact_local.
Not found = unresolved (never 'engineered'); the Mac applies the frozen decision rules on import.
Resumable: finished checks and downloaded archives are kept.
"""
import concurrent.futures as cf
import csv
import datetime as dt
import gzip
import io
import json
import os
import sys
import threading
import time
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.abspath(__file__))
DS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
UA = {"User-Agent": "research-pipeline/1.0"}
LOCK = threading.Lock()
COMP = str.maketrans("ACGTN", "TGCAN")


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def get(url, timeout=1800, tries=4):
    err = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception as e:
            err = e
            time.sleep(10 * (i + 1))
    raise err


def jdump(o, p):
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(o, f, indent=1, sort_keys=True)
    os.replace(tmp, p)


def log(msg):
    with LOCK:
        print(f"{now()} {msg}", flush=True)
        with open(os.path.join(ROOT, "exact_local_verify.log"), "a") as f:
            f.write(f"{now()} {msg}\n")


def reports(taxon, extra):
    out, token = [], ""
    while True:
        url = f"{DS}genome/taxon/{taxon}/dataset_report?filters.assembly_version=current&page_size=500{extra}" + \
              (f"&page_token={token}" if token else "")
        d = json.loads(get(url, 300))
        out += d.get("reports") or []
        token = d.get("next_page_token")
        if not token or len(out) >= 2000:
            return out


def download(kind, ident, what):
    """Archive cached under refs/<kind>_<ident>_<what>.zip (kept for provenance)."""
    os.makedirs(os.path.join(ROOT, "refs"), exist_ok=True)
    p = os.path.join(ROOT, "refs", f"{kind}_{ident}_{what}.zip")
    if not os.path.exists(p):
        if kind == "virus":
            url = f"{DS}virus/taxon/{ident}/genome/download?include_sequence=GENOME"
        else:
            url = f"{DS}genome/accession/{ident}/download?include_annotation_type={what}"
        data = get(url)
        with open(p + ".part", "wb") as f:
            f.write(data)
        zipfile.ZipFile(p + ".part").testzip()
        os.replace(p + ".part", p)
    return p


def fasta_records(zpath):
    z = zipfile.ZipFile(zpath)
    for n in z.namelist():
        if not (n.endswith(".fna") or n.endswith(".fa") or n.endswith(".fasta")):
            continue
        name, buf = None, []
        with z.open(n) as fh:
            for line in io.TextIOWrapper(fh, encoding="ascii", errors="replace"):
                if line.startswith(">"):
                    if name:
                        yield n, name, "".join(buf).upper()
                    name, buf = line[1:].strip(), []
                else:
                    buf.append(line.strip())
        if name:
            yield n, name, "".join(buf).upper()


def search(zpath, label, queries):
    """queries: {qid: DNA seq}. Returns {qid: [hit dicts]} (both strands, exact)."""
    rc = {q: s.translate(COMP)[::-1] for q, s in queries.items()}
    hits = {q: [] for q in queries}
    for fname, hdr, seq in fasta_records(zpath):
        acc = hdr.split()[0]
        for q, s in queries.items():
            for strand, t in (("+", s), ("-", rc[q])):
                i = seq.find(t)
                while i != -1 and len(hits[q]) < 20:
                    hits[q].append(dict(record=acc, start=i + 1, end=i + len(t), strand=strand, source=label,
                                        title=hdr[:150], file=os.path.basename(fname)))
                    i = seq.find(t, i + 1)
    return hits


def references(sp, kingdom):
    """[(label, zip path)] downloaded once per species; failures are recorded, not fatal."""
    refs, notes = [], []
    try:
        if kingdom in ("Bacteria", "Archaea"):
            rs = [r for r in reports(sp, "") if (r.get("assembly_info") or {}).get("assembly_level") in ("Complete Genome", "Chromosome")]
            if not rs:
                rs = reports(sp, "")
            rs = sorted(rs, key=lambda r: (not r["accession"].startswith("GCF_"), r["accession"]))[:60]
            for i in range(0, len(rs), 10):
                ids = ",".join(r["accession"] for r in rs[i:i + 10])
                refs.append((f"assemblies:{ids}", download("genome", ids, "GENOME_FASTA")))
        elif kingdom in ("Riboviria", "Varidnaviria", "Viruses", "Monodnaviria", "Duplodnaviria"):
            refs.append((f"ncbi_virus_taxon:{sp}", download("virus", sp, "GENOME")))
        else:
            rs = reports(sp, "&filters.reference_only=true") or reports(sp, "&filters.assembly_level=chromosome")
            for r in rs[:1]:
                refs.append((f"assembly:{r['accession']}", download("genome", r["accession"], "GENOME_FASTA")))
                try:
                    refs.append((f"refseq_rna:{r['accession']}", download("genome", r["accession"], "RNA_FASTA")))
                except Exception as e:
                    notes.append(f"RNA_FASTA unavailable for {r['accession']}: {e!r}"[:200])
            if not rs:
                notes.append("no reference/chromosome assembly for species")
    except Exception as e:
        notes.append(f"reference retrieval failed: {e!r}"[:300])
    return refs, notes


def run_species(sp, checks):
    out_dir = os.path.join(ROOT, "evidence")
    todo = [c for c in checks if not os.path.exists(os.path.join(out_dir, c["check_id"] + ".json"))]
    if not todo:
        return
    log(f"species {sp} {checks[0]['species_name']}: {len(todo)} checks; fetching references")
    refs, notes = references(sp, checks[0]["kingdom"])
    queries = {}
    for c in todo:
        queries[c["check_id"] + "|chain"] = c["deposited_sequence"].replace("U", "T")
        for i, iv in enumerate(sorted(set(c["intervals"].split(";")) - {""})):
            if iv != c["deposited_sequence"]:
                queries[c["check_id"] + f"|interval{i}"] = iv.replace("U", "T")
    allhits = {q: [] for q in queries}
    for label, zp in refs:
        h = search(zp, label, queries)
        for q in queries:
            allhits[q] += h[q]
    for c in todo:
        chain = [dict(x, independent="yes", exact_confirmed="yes",
                      independence_basis="genome_assembly" if x["source"].startswith(("assembl", "ncbi_virus")) else "refseq_transcript")
                 for x in allhits[c["check_id"] + "|chain"]]
        iv = [dict(x, interval=queries[q].replace("T", "U"), exact_confirmed="yes")
              for q in queries if q.startswith(c["check_id"] + "|interval") for x in allhits[q]]
        outcome = "whole_chain_found" if chain else ("interval_found_flank_differs" if iv else
                                                     ("not_found" if refs else "search_failed"))
        rec = dict(check_id=c["check_id"], seq_sha256=c["seq_sha256"], species_taxid=sp, species_name=c["species_name"],
                   deposited_sequence=c["deposited_sequence"], deposited_length=len(c["deposited_sequence"]),
                   families=c["families"].split(";"), outcome=outcome, whole_chain_hits=chain, interval_hits=iv,
                   routes_tried=[dict(route="exact_local", reference=l, archive=os.path.basename(z),
                                      archive_bytes=os.path.getsize(z)) for l, z in refs],
                   notes=notes, method="exact substring search, both strands, full deposited sequence (U->T)",
                   host=os.uname().nodename, evaluated_utc=now(), source="exact_local")
        if outcome != "search_failed":
            jdump(rec, os.path.join(out_dir, c["check_id"] + ".json"))
        log(f"  {c['check_id']} {c['species_name']} -> {outcome} ({len(chain)} whole-chain hits)")


def main():
    pend = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    os.makedirs(os.path.join(ROOT, "evidence"), exist_ok=True)
    by = {}
    for c in pend:
        by.setdefault(c["species_taxid"], []).append(c)
    log(f"start: {len(pend)} checks in {len(by)} species, {workers} workers")
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_species, sp, cs) for sp, cs in sorted(by.items(), key=lambda x: -len(x[1]))]
        for f in futs:
            try:
                f.result()
            except Exception as e:
                log(f"species worker error: {e!r}")
    done = len(os.listdir(os.path.join(ROOT, "evidence")))
    log(f"finished: {done}/{len(pend)} checks have exact-local-search evidence records")


if __name__ == "__main__":
    main()
