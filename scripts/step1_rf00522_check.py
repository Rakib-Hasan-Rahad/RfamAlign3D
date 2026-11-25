"""Independent check of RF00522 against the raw Rfam.seed.gz (does NOT import step1_inventory.py).

Slices the raw RF00522 record by plain text search, counts rows/links with its own simple logic, saves the raw
record for inspection, and compares with (a) the previous run's counts and (b) the production outputs.
"""
import csv
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "01_inventory")
EXPECTED_PREVIOUS = dict(seed_rows=43, linked_rows=8, distinct_linked_sequences=7, linked_pdb_chains=16)

text = gzip.open(os.path.join(ROOT, "inputs", "rfam", "15.1", "Rfam.seed.gz")).read().decode("latin-1")
records = text.split("# STOCKHOLM 1.0")
rec = [r for r in records if "\n#=GF AC   RF00522\n" in r]
assert len(rec) == 1, f"RF00522 records found: {len(rec)}"
rec = "# STOCKHOLM 1.0" + rec[0].split("\n//")[0] + "\n//\n"
open(os.path.join(OUT, "RF00522_raw_record.sto"), "w").write(rec)

seqs, tags = {}, []
for line in rec.splitlines():
    if not line or line.startswith("# STOCKHOLM") or line.startswith("//"):
        continue
    if line.startswith("#=GR"):
        _, rid, tag, _ann = line.split(None, 3)
        if tag.endswith("_SS") and tag.count("_") == 2 and len(tag.split("_")[0]) == 4:
            tags.append((rid, tag.split("_")[0].upper(), tag.split("_")[1]))
    elif not line.startswith("#"):
        rid, s = line.split()
        seqs[rid] = seqs.get(rid, "") + s

linked_rows = {r for r, _, _ in tags}
indep = dict(
    seed_rows=len(seqs), linked_rows=len(linked_rows),
    distinct_linked_sequences=len({seqs[r].replace("-", "").replace(".", "").upper() for r in linked_rows}),
    linked_pdb_chains=len({(p, c) for _, p, c in tags}),
    all_link_rows_are_sequence_rows=linked_rows <= set(seqs),
    alignment_width=sorted({len(s) for s in seqs.values()}),
    links=sorted(f"{r}|{p}|{c}" for r, p, c in tags))

fam = next(r for r in csv.DictReader(open(os.path.join(OUT, "families.tsv")), delimiter="\t") if r["family"] == "RF00522")
plinks = sorted(f"{r['seed_row_id']}|{r['pdb_id']}|{r['chain_id']}"
                for r in csv.DictReader(open(os.path.join(OUT, "row_pdb_links.tsv")), delimiter="\t")
                if r["family"] == "RF00522")
production = dict(seed_rows=int(fam["total_seed_rows"]), linked_rows=int(fam["linked_rows"]),
                  distinct_linked_sequences=int(fam["distinct_linked_sequences"]),
                  linked_pdb_chains=int(fam["distinct_pdb_chains"]))

res = dict(raw_record_sha256=hashlib.sha256(rec.encode("latin-1")).hexdigest(), independent=indep,
           production=production, expected_previous=EXPECTED_PREVIOUS,
           independent_equals_previous={k: indep[k] == v for k, v in EXPECTED_PREVIOUS.items()},
           independent_equals_production={k: indep[k] == v for k, v in production.items()},
           link_list_identical=indep["links"] == plinks)
json.dump(res, open(os.path.join(OUT, "rf00522_independent_check.json"), "w"), indent=2)
print(json.dumps({k: v for k, v in res.items() if k != "independent"}, indent=1))
ok = all(res["independent_equals_previous"].values()) and all(res["independent_equals_production"].values()) \
    and res["link_list_identical"] and indep["all_link_rows_are_sequence_rows"]
print("RF00522 independent check:", "PASS" if ok else "DIFFERENCE - investigate")
sys.exit(0 if ok else 1)
