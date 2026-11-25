"""R8 fallback: rerun still-queued species BLAST jobs in megablast mode (same database, taxon filter and queries).

  python scripts/step2/megablast_fallback.py [--hours 6]

Why: blastn (word size 11/7) jobs for a few species stayed WAITING in NCBI's queue for hours. The decision only needs
full-length 100%-identical whole-chain hits; any such hit contains megablast's 28-nt seed whenever the query is at
least 28 nt, so megablast finds every hit the decision can use. Every hit is still confirmed by exact string search of
the efetched record (unchanged rules in s03). Short interval queries (<28 nt) only label flank differences and may be
missed; this is recorded.
Writes cache/blast/<species>_<tag>.megablast.xml; s03_sequence_evidence.py uses it when the blastn XML is absent.
"""
import argparse
import collections
import json
import os
import re
import sys
import time
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import s03_sequence_evidence as S  # noqa: E402
import xml.etree.ElementTree as ET  # noqa: E402

B = S.BLAST


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--hours", type=float, default=6); a = ap.parse_args()
    tax = C.Taxonomy(); S.TAX[0] = tax
    q = S.build_queue()
    todo = [e for k, e in sorted(q.items()) if not os.path.exists(S.key_path(*k))]
    by_sp = collections.defaultdict(list)
    for e in todo:
        by_sp[e["species"]].append(e)
    jobs = {}
    for sp, es in by_sp.items():
        queries = {}
        for e in es:
            if not set(e["seq"]) <= C.STD:
                continue
            pre = e["seq_sha"][:10] + "_"
            queries[pre + "chain"] = e["seq"]
            for i, s in enumerate(sorted(e["intervals"])):
                if s != e["seq"]:
                    queries[f"{pre}interval{i}"] = s
        if not queries:
            continue
        tag = C.hashlib.sha256(json.dumps(sorted(queries.items())).encode()).hexdigest()[:16]
        path = C.CP("blast", f"{sp}_{tag}.xml")
        mpath = path[:-4] + ".megablast.xml"
        if os.path.exists(path) or os.path.exists(mpath):
            continue
        fasta = "".join(f">{k}\n{v.replace('U', 'T')}\n" for k, v in sorted(queries.items()))
        body = urllib.parse.urlencode(dict(CMD="Put", PROGRAM="blastn", MEGABLAST="on", WORD_SIZE=28, DATABASE="core_nt",
                                           QUERY=fasta, ENTREZ_QUERY=f"txid{sp}[Organism:exp]", EXPECT=10,
                                           HITLIST_SIZE=50, FORMAT_TYPE="XML", FILTER="F")).encode()
        r = C.http(B, data=body, timeout=180, min_interval=11).decode()
        rid = re.search(r"RID = (\S+)", r).group(1)
        short = sorted(k for k, v in queries.items() if len(v) < 28)
        jobs[mpath] = rid
        print(f"{time.strftime('%H:%M:%S')} {tax.name.get(sp, sp)} ({sp}): {len(queries)} queries -> megablast RID {rid}"
              + (f"; queries <28 nt: {short}" if short else ""), flush=True)
    t0 = time.time()
    while jobs and time.time() - t0 < a.hours * 3600:
        time.sleep(60)
        for mpath, rid in list(jobs.items()):
            s = C.http(f"{B}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}", timeout=120, min_interval=11).decode()
            st = (re.search(r"Status=(\w+)", s) or [None, "?"])[1]
            if st == "READY":
                x = C.http(f"{B}?CMD=Get&FORMAT_TYPE=XML&RID={rid}", timeout=600, min_interval=11)
                open(mpath + ".tmp", "wb").write(x); ET.parse(mpath + ".tmp"); os.replace(mpath + ".tmp", mpath)
                C.log_source(f"blast:{os.path.basename(mpath)}", f"{B}?RID={rid}", mpath, "downloaded", "megablast fallback (R8)")
                print(f"{time.strftime('%H:%M:%S')} {rid} READY -> {os.path.basename(mpath)}", flush=True)
                del jobs[mpath]
            elif st in ("FAILED", "UNKNOWN"):
                print(f"{time.strftime('%H:%M:%S')} {rid} {st}", flush=True); del jobs[mpath]
        if jobs:
            print(f"{time.strftime('%H:%M:%S')} waiting: {len(jobs)} ({(time.time() - t0) / 60:.0f} min)", flush=True)
    print("remaining:", sorted(jobs.values()) or "none")


if __name__ == "__main__":
    main()
