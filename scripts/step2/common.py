"""Shared helpers for Step 2: paths, cached HTTP with a sources manifest, TSV/JSON I/O, taxonomy."""
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import ssl
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request

import certifi
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "outputs", "02_structure_preparation")
STEP1 = os.path.join(ROOT, "outputs", "01_inventory")
CACHE = os.path.join(ROOT, "cache")
P = lambda *a: os.path.join(OUT, *a)  # noqa: E731
CP = lambda *a: os.path.join(CACHE, *a)  # noqa: E731
SSL = ssl.create_default_context(cafile=certifi.where())
UA = "research-pipeline/1.0 (academic RNA structure study)"
MANIFEST = P("sources_manifest.tsv")
STD = set("ACGU")
csv.field_size_limit(1 << 30)


def policy():
    return yaml.safe_load(open(P("step2_policy.yaml")))


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_tsv(path, rows, cols=None):
    rows = list(rows)
    if cols is None:
        cols = list(rows[0]) if rows else ["empty"]
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})
    return path


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def jdump(obj, path):
    """Atomic JSON write (temp file + rename): an interruption never leaves a truncated cache entry."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)


def jload(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path) as f:
        return json.load(f)


def log_source(key, url, path, status, note=""):
    new = not os.path.exists(MANIFEST)
    with open(MANIFEST, "a") as f:
        if new:
            f.write("key\turl\tlocal_path\tstatus\tbytes\tsha256\tretrieved_utc\tnote\n")
        b = os.path.getsize(path) if path and os.path.exists(path) else ""
        s = sha256(path) if path and os.path.exists(path) else ""
        f.write(f"{key}\t{url}\t{os.path.relpath(path, ROOT) if path else ''}\t{status}\t{b}\t{s}\t{now()}\t{note}\n")


_last = {"t": 0.0}


def http(url, data=None, headers=None, timeout=120, tries=4, min_interval=0.34):
    """GET/POST with retries and a polite minimum interval. Returns bytes; raises on final failure."""
    err = None
    for i in range(tries):
        wait = min_interval - (time.time() - _last["t"])
        if wait > 0:
            time.sleep(wait)
        _last["t"] = time.time()
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": UA, **(headers or {})})
            with urllib.request.urlopen(req, context=SSL, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            err = e
            if e.code in (400, 404, 410):
                raise
        except Exception as e:  # network
            err = e
        time.sleep(2 * (i + 1) ** 2)
    raise err


def fetch_cached(url, path, key, validate=None):
    """Download once into path (temp file + validate + rename). Never overwrites an existing valid file."""
    if os.path.exists(path):
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = http(url, timeout=600)
    tmp = os.path.join(os.path.dirname(path), "tmp_" + os.path.basename(path))   # keeps the extension
    with open(tmp, "wb") as f:
        f.write(data)
    if validate:
        validate(tmp)
    os.replace(tmp, path)
    log_source(key, url, path, "downloaded")
    return path


def gz_ok(p):
    with gzip.open(p) as f:
        while f.read(1 << 20):
            pass


# ------------------------------------------------------------------ taxonomy (pinned NCBI taxdump)
class Taxonomy:
    def __init__(self, path=os.path.join(ROOT, "inputs", "ncbi_taxonomy", "taxdump.tar.gz")):
        cache = CP("ncbi", "taxonomy_index.json.gz")
        if os.path.exists(cache) and jload(cache + ".src", {}).get("sha256") == sha256(path):
            d = json.load(gzip.open(cache, "rt"))
        else:
            d = dict(parent={}, rank={}, name={}, merged={}, deleted=[], sci2id={})
            with tarfile.open(path) as tf:
                for line in tf.extractfile("nodes.dmp").read().decode().splitlines():
                    f = [x.strip() for x in line.split("|")]
                    d["parent"][f[0]] = f[1]
                    d["rank"][f[0]] = f[2]
                for line in tf.extractfile("names.dmp").read().decode().splitlines():
                    f = [x.strip() for x in line.split("|")]
                    if f[3] == "scientific name":
                        d["name"][f[0]] = f[1]
                        d["sci2id"].setdefault(f[1].lower(), f[0])
                for line in tf.extractfile("merged.dmp").read().decode().splitlines():
                    f = [x.strip() for x in line.split("|")]
                    d["merged"][f[0]] = f[1]
                d["deleted"] = [line.split("|")[0].strip() for line in
                                tf.extractfile("delnodes.dmp").read().decode().splitlines()]
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            with gzip.open(cache, "wt") as f:
                json.dump(d, f)
            jdump({"sha256": sha256(path)}, cache + ".src")
        self.parent, self.rank, self.name = d["parent"], d["rank"], d["name"]
        self.merged, self.deleted, self.sci2id = d["merged"], set(d["deleted"]), d["sci2id"]

    def current(self, taxid):
        t = str(taxid).strip()
        if t in self.merged:
            return self.merged[t], "merged"
        if t in self.deleted:
            return None, "deleted"
        if t not in self.parent:
            return None, "unknown"
        return t, "current"

    def lineage(self, taxid):
        out, t = [], str(taxid)
        while t in self.parent and t not in [x for x in out]:
            out.append(t)
            if self.parent[t] == t:
                break
            t = self.parent[t]
        return out

    def species(self, taxid):
        """Return (species_taxid, species_name, status). status: species|merged->species|no_species_rank|deleted|unknown"""
        cur, st = self.current(taxid)
        if cur is None:
            return None, None, st
        for t in self.lineage(cur):
            if self.rank.get(t) == "species":
                return t, self.name.get(t), ("merged_then_species" if st == "merged" else "species")
        return None, None, "no_species_rank"

    def is_placeholder(self, taxid):
        """R5: taxa that do not identify an organism: lineage under 'unclassified sequences' (12908) or
        'environmental samples' (48479), or names with metagenome/unidentified/uncultured/synthetic construct."""
        cur, _ = self.current(taxid)
        if not cur:
            return True
        lin = self.lineage(cur)
        if "12908" in lin or "48479" in lin or "32630" in lin:
            return True
        n = (self.name.get(cur) or "").lower()
        return any(w in n for w in ("metagenome", "unidentified", "uncultured", "synthetic construct"))

    def superkingdom(self, taxid):
        cur, _ = self.current(taxid)
        if not cur:
            return ""
        for t in self.lineage(cur):
            if self.rank.get(t) in ("superkingdom", "domain", "acellular root", "realm"):
                return self.name.get(t, "")
        lin = [self.name.get(t, "") for t in self.lineage(cur)]
        return "Viruses" if "Viruses" in lin else ""

    def binomials(self, text):
        """Exact NCBI scientific names (species rank or below) of the form 'Genus species[ ...]' found in text."""
        hits = set()
        words = re.findall(r"[A-Za-z][A-Za-z0-9\-\.]*", text or "")
        for i in range(len(words) - 1):
            for n in (4, 3, 2):
                cand = " ".join(words[i:i + n]).lower()
                tid = self.sci2id.get(cand)
                if tid and words[i][0].isupper():
                    sp = self.species(tid)[0]
                    if sp:
                        hits.add(tid)
                        break
        return sorted(hits)
