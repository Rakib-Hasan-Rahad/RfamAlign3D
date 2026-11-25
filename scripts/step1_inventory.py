"""Step 1: master inventory of explicit Rfam seed-row <-> PDB-chain links (Rfam 15.1, every family).

Reads the pinned Rfam.seed.gz, preserves every seed row (gapped + ungapped) and every column track, extracts
explicit row-level PDB-chain tags, keeps unresolved PDB mentions as leads, validates, and exports a workbook
backed by lossless TSV / JSONL.gz files in outputs/01_inventory/.

Definitions (documented in outputs/01_inventory/README.md):
  explicit link   #=GR <row_id> <PDB>_<chain>_SS <annotation>  with tag matching LINK_RE; row_id must be a
                  sequence row of the same family. PDB ID uppercased; chain preserved exactly.
  PDB lead        (a) a #=GR tag that starts like a PDB ID but does not match LINK_RE (e.g. no chain);
                  (b) a PDB-ID-like token in #=GF CC free text of a family whose CC mentions "PDB".
  ungapped        gapped row with only '-' and '.' removed (case and all other symbols kept).
  sequence hash   sha256 of ungapped.upper() encoded as ASCII. No T->U or other change; T and non-ACGU
                  symbols are flagged per row.
Exit status 0 only if every validation check passes.
"""
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import sys

import openpyxl
import yaml
from openpyxl.styles import Font

PARSER_VERSION = "step1-inventory 1.0.0"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "01_inventory")
LINK_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})_([A-Za-z0-9]{1,4})_SS$")
PDBLIKE_TAG_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:$|[_\-.])")
FREE_PDB_RE = re.compile(r"(?<![0-9A-Z])([1-9][A-Za-z0-9]{3})(?:_([A-Za-z0-9]{1,4}))?(?![A-Za-z0-9])")
UNIT_RE = re.compile(r"^[0-9]+[a-z]+$")
RANGE_RE = re.compile(r"^(.*)/(\d+)-(\d+)$")
GAPS = "-."
EXCEL_MAX = 32767
COMPANION = "SEE_COMPANION_FILE"


# ------------------------------------------------------------------ parsing
def decode_lines(raw_bytes, anomalies):
    for n, b in enumerate(raw_bytes.split(b"\n"), 1):
        b = b.rstrip(b"\r")
        try:
            yield n, b.decode("utf-8")
        except UnicodeDecodeError:
            t = b.decode("latin-1")
            kind = t[:7] if t.startswith("#=") else "non-markup"
            anomalies.append(dict(family="", line=n, kind="latin1_fallback_decode",
                                  detail=f"{kind}: not valid UTF-8; decoded as Latin-1"))
            yield n, t


def new_record(start_line):
    return dict(start_line=start_line, gf=collections.OrderedDict(), gs=[], seq=collections.OrderedDict(),
                gr=collections.OrderedDict(), gc=collections.OrderedDict(), problems=[])


def parse_stockholm(lines):
    """Yield one record per '# STOCKHOLM' ... '//' block. Fragments are joined by key (wrapped alignments);
    a key repeated inside one blank-line-separated block is a duplicate, not a wrap."""
    rec = None
    block_seen = set()
    for n, line in lines:
        if line.startswith("# STOCKHOLM"):
            if rec is not None:
                rec["problems"].append(("missing_terminator", f"new header at line {n}"))
                yield rec
            rec = new_record(n)
            block_seen = set()
            continue
        if rec is None:
            if line.strip():
                yield dict(orphan_line=(n, line))
            continue
        if line.startswith("//"):
            rec["end_line"] = n
            yield rec
            rec = None
            continue
        if not line.strip():
            block_seen = set()
            continue
        if line[0].isspace():
            rec["problems"].append(("leading_whitespace_line", f"line {n}"))
            continue

        def frag(store, key, text, kind):
            if key in block_seen:
                rec["problems"].append(("duplicate_key_in_block", f"{kind} {key!r} line {n}"))
            block_seen.add(key)
            store[key] = store.get(key, "") + text

        if line.startswith("#=GF"):
            p = line.split(None, 2)
            tag = p[1] if len(p) > 1 else ""
            rec["gf"].setdefault(tag, []).append(p[2] if len(p) > 2 else "")
        elif line.startswith("#=GS"):
            p = line.split(None, 3)
            rec["gs"].append(dict(row_id=p[1], tag=p[2] if len(p) > 2 else "", text=p[3] if len(p) > 3 else "",
                                  line=n))
        elif line.startswith("#=GR"):
            p = line.split(None, 3)
            if len(p) < 4:
                rec["problems"].append(("malformed_GR", f"line {n}"))
                continue
            frag(rec["gr"], (p[1], p[2]), p[3].rstrip(), "GR")
            rec.setdefault("gr_line", {}).setdefault((p[1], p[2]), n)
        elif line.startswith("#=GC"):
            p = line.split(None, 2)
            if len(p) < 3:
                rec["problems"].append(("malformed_GC", f"line {n}"))
                continue
            frag(rec["gc"], p[1], p[2].rstrip(), "GC")
        elif line.startswith("#"):
            rec["problems"].append(("unrecognized_markup", f"line {n}: {line[:40]}"))
        else:
            p = line.split()
            if len(p) != 2:
                rec["problems"].append(("malformed_sequence_line", f"line {n}: {len(p)} fields"))
                continue
            frag(rec["seq"], p[0], p[1], "sequence")
    if rec is not None:
        rec["problems"].append(("missing_terminator", "end of file"))
        yield rec


# ------------------------------------------------------------------ helpers
def ungap(s):
    return s.translate({ord(c): None for c in GAPS})


def seq_hash(ungapped):
    return hashlib.sha256(ungapped.upper().encode("ascii")).hexdigest()


def split_range(row_id):
    m = RANGE_RE.match(row_id)
    if not m:
        return row_id, "", "", ""
    a, b = int(m.group(2)), int(m.group(3))
    return m.group(1), a, b, ("+" if a <= b else "-")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gf1(rec, tag):
    v = rec["gf"].get(tag)
    return v[0].strip() if v else ""


def write_tsv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ------------------------------------------------------------------ build
def build(seed_path):
    anomalies = []
    raw = gzip.open(seed_path).read()
    families, seed_rows, tracks, links, leads, census = [], [], [], {}, [], collections.Counter()
    alphabet = collections.Counter()
    for rec in parse_stockholm(decode_lines(raw, anomalies)):
        if "orphan_line" in rec:
            n, l = rec["orphan_line"]
            anomalies.append(dict(family="", line=n, kind="line_outside_record", detail=l[:80]))
            continue
        acc, fid = gf1(rec, "AC"), gf1(rec, "ID")
        for k, d in rec["problems"]:
            anomalies.append(dict(family=acc, line="", kind=k, detail=d))
        for tag, vals in rec["gf"].items():
            census[("GF", tag)] += len(vals)
        for g in rec["gs"]:
            census[("GS", g["tag"])] += 1
        for (_, tag) in rec["gr"]:
            census[("GR", "<PDB>_<chain>_SS" if LINK_RE.match(tag) else tag)] += 1
        for tag in rec["gc"]:
            census[("GC", tag)] += 1

        widths = {len(s) for s in rec["seq"].values()}
        width = widths.pop() if len(widths) == 1 else None
        if width is None:
            anomalies.append(dict(family=acc, line="", kind="unequal_row_lengths", detail=str(sorted(widths))))
        for tag, s in rec["gc"].items():
            if len(s) != width:
                anomalies.append(dict(family=acc, line="", kind="GC_width_mismatch", detail=f"{tag}: {len(s)} vs {width}"))
        for (rid, tag), s in rec["gr"].items():
            if rid not in rec["seq"]:
                anomalies.append(dict(family=acc, line=rec["gr_line"][(rid, tag)], kind="orphan_GR",
                                      detail=f"{rid} {tag}"))
            if len(s) != width:
                anomalies.append(dict(family=acc, line="", kind="GR_width_mismatch",
                                      detail=f"{rid} {tag}: {len(s)} vs {width}"))
        for g in rec["gs"]:
            if g["row_id"] not in rec["seq"]:
                anomalies.append(dict(family=acc, line=g["line"], kind="orphan_GS", detail=g["row_id"]))

        # rows
        row_info = {}
        for rid, gapped in rec["seq"].items():
            ug = ungap(gapped)
            alphabet.update(gapped)
            src, a, b, strand = split_range(rid)
            nonacgu = sorted(set(ug.upper()) - set("ACGU"))
            info = dict(family=acc, seed_row_id=rid, source_accession=src, source_start=a, source_end=b,
                        source_strand=strand, original_gapped_sequence=gapped, ungapped_sequence=ug,
                        ungapped_length=len(ug), sequence_group_hash=seq_hash(ug),
                        contains_T="T" in ug.upper(), non_ACGU_symbols="".join(nonacgu),
                        has_lowercase=any(c.islower() for c in ug))
            row_info[rid] = info
            seed_rows.append(dict(info, row_annotations=[dict(tag=t, annotation=s) for (r, t), s in rec["gr"].items()
                                                         if r == rid],
                                  row_gs=[dict(tag=g["tag"], text=g["text"]) for g in rec["gs"] if g["row_id"] == rid]))
        tracks.append(dict(family=acc, family_id=fid, alignment_width=width, n_rows=len(rec["seq"]),
                           gc=dict(rec["gc"]), gf={k: v for k, v in rec["gf"].items()},
                           record_lines=[rec["start_line"], rec.get("end_line")]))

        # explicit links + GR leads
        fam_links = []
        for (rid, tag), ann in rec["gr"].items():
            m = LINK_RE.match(tag)
            if m:
                if rid not in row_info:
                    leads.append(dict(family=acc, seed_row_id=rid, candidate_pdb_id=m.group(1).upper(),
                                      candidate_chain=m.group(2), raw_tag=tag, evidence_source="#=GR",
                                      original_evidence=f"#=GR {rid} {tag}",
                                      reason="link tag names a row id that is not a sequence row of this family"))
                    continue
                key = (acc, rid, m.group(1).upper(), m.group(2))
                if key in links:
                    links[key]["raw_tag"] += ";" + tag
                    anomalies.append(dict(family=acc, line="", kind="duplicate_link_key_after_normalization",
                                          detail=f"{rid} {tag}"))
                    continue
                ri = row_info[rid]
                links[key] = dict(family=acc, family_id=fid, seed_row_id=rid,
                                  source_accession=ri["source_accession"], source_start=ri["source_start"],
                                  source_end=ri["source_end"], source_strand=ri["source_strand"],
                                  original_gapped_sequence=ri["original_gapped_sequence"],
                                  ungapped_sequence=ri["ungapped_sequence"], ungapped_length=ri["ungapped_length"],
                                  sequence_group_hash=ri["sequence_group_hash"], pdb_id=m.group(1).upper(),
                                  chain_id=m.group(2), evidence_type="explicit_GR_row_tag", raw_tag=tag,
                                  structure_annotation=ann)
                fam_links.append(key)
            elif PDBLIKE_TAG_RE.match(tag):
                leads.append(dict(family=acc, seed_row_id=rid, candidate_pdb_id=tag[:4].upper(), candidate_chain="",
                                  raw_tag=tag, evidence_source="#=GR",
                                  original_evidence=f"#=GR {rid} {tag} {ann}",
                                  reason="row-level tag looks like a PDB ID but does not match <PDB>_<chain>_SS "
                                         "(chain not stated)"))

        # free-text leads: #=GF CC of families that mention PDB
        cc = " ".join(x.strip() for x in rec["gf"].get("CC", []))
        if re.search(r"\bPDB\b", cc, re.I):
            explicit_pdbs = {k[2] for k in fam_links}
            seen = set()
            for m in FREE_PDB_RE.finditer(cc):
                tok, ch = m.group(1).upper(), m.group(2) or ""
                if not re.search(r"[A-Z]", tok[1:]) or (tok, ch) in seen:
                    continue
                if UNIT_RE.match(m.group(1)):   # e.g. "22nt": a length, not a PDB ID
                    continue
                seen.add((tok, ch))
                ctx = cc[max(0, m.start() - 80): m.end() + 40]
                leads.append(dict(family=acc, seed_row_id="", candidate_pdb_id=tok, candidate_chain=ch,
                                  raw_tag=m.group(0), evidence_source="#=GF CC", original_evidence=ctx,
                                  also_explicitly_linked_in_family="yes" if tok in explicit_pdbs else "no",
                                  reason=("family free text names PDB ID and chain but no seed row" if ch else
                                          "family free text names PDB ID only; no seed row or chain")))

        fl = [links[k] for k in fam_links]
        families.append(dict(
            family=acc, family_id=fid, description=gf1(rec, "DE"), rfam_type=gf1(rec, "TP"),
            alignment_width=width, total_seed_rows=len(rec["seq"]),
            distinct_row_sequences=len({r["sequence_group_hash"] for r in row_info.values()}),
            linked_rows=len({x["seed_row_id"] for x in fl}),
            distinct_linked_sequences=len({x["sequence_group_hash"] for x in fl}),
            distinct_pdb_entries=len({x["pdb_id"] for x in fl}),
            distinct_pdb_chains=len({(x["pdb_id"], x["chain_id"]) for x in fl}),
            unresolved_leads=sum(1 for x in leads if x["family"] == acc),
            link_status="explicit link(s) found" if fl else "no explicit link found in this seed record"))
    for x in leads:
        x.setdefault("also_explicitly_linked_in_family", "")
    return dict(families=families, seed_rows=seed_rows, tracks=tracks, links=list(links.values()), leads=leads,
                census=census, alphabet=alphabet, anomalies=anomalies, raw_len=len(raw),
                raw_sha=hashlib.sha256(raw).hexdigest())


# ------------------------------------------------------------------ export
FAM_COLS = ["family", "family_id", "description", "rfam_type", "alignment_width", "total_seed_rows",
            "distinct_row_sequences", "linked_rows", "distinct_linked_sequences", "distinct_pdb_entries",
            "distinct_pdb_chains", "unresolved_leads", "link_status"]
LINK_COLS = ["family", "family_id", "seed_row_id", "source_accession", "source_start", "source_end", "source_strand",
             "pdb_id", "chain_id", "evidence_type", "raw_tag", "sequence_group_hash", "ungapped_length",
             "ungapped_sequence", "original_gapped_sequence", "structure_annotation"]
LEAD_COLS = ["family", "seed_row_id", "candidate_pdb_id", "candidate_chain", "raw_tag", "evidence_source",
             "also_explicitly_linked_in_family", "reason", "original_evidence"]
TEXT_COLS = {"family", "family_id", "seed_row_id", "source_accession", "pdb_id", "chain_id", "raw_tag",
             "candidate_pdb_id", "candidate_chain", "sequence_group_hash", "ungapped_sequence",
             "original_gapped_sequence", "structure_annotation", "description", "rfam_type", "original_evidence"}


def excel_value(v, key):
    if isinstance(v, str) and len(v) > EXCEL_MAX:
        return f"{COMPANION}|{key}"
    return v


def write_sheet(ws, rows, cols, keyfn):
    ws.append(cols)
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in rows:
        vals = [excel_value(r.get(c, ""), keyfn(r, c)) for c in cols]
        ws.append(vals)
        for c, cell in zip(cols, ws[ws.max_row]):
            if c in TEXT_COLS or isinstance(cell.value, str):
                cell.data_type = "s"          # never a formula, never a number
                cell.number_format = "@"
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions


def export(d, prov):
    os.makedirs(OUT, exist_ok=True)
    fams = sorted(d["families"], key=lambda r: r["family"])
    links = sorted(d["links"], key=lambda r: (r["family"], r["seed_row_id"], r["pdb_id"], r["chain_id"]))
    leads = sorted(d["leads"], key=lambda r: (r["family"], r["evidence_source"], r["seed_row_id"],
                                              r["candidate_pdb_id"], r["candidate_chain"]))
    write_tsv(os.path.join(OUT, "families.tsv"), fams, FAM_COLS)
    write_tsv(os.path.join(OUT, "row_pdb_links.tsv"), links, LINK_COLS)
    write_tsv(os.path.join(OUT, "pdb_leads.tsv"), leads, LEAD_COLS)
    write_tsv(os.path.join(OUT, "anomalies.tsv"), d["anomalies"], ["family", "line", "kind", "detail"])
    write_tsv(os.path.join(OUT, "feature_census.tsv"),
              [dict(markup=k[0], tag=k[1], lines=v) for k, v in sorted(d["census"].items())],
              ["markup", "tag", "lines"])
    with gzip.open(os.path.join(OUT, "seed_rows.jsonl.gz"), "wt", encoding="utf-8") as f:
        for r in d["seed_rows"]:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    with gzip.open(os.path.join(OUT, "family_tracks.jsonl.gz"), "wt", encoding="utf-8") as f:
        for r in d["tracks"]:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "About"
    for k, v in prov.items():
        ws.append([k, str(v)])
    ws.append([])
    for line in ABOUT:
        ws.append([line])
    write_sheet(wb.create_sheet("Families"), fams, FAM_COLS, lambda r, c: f"family_tracks.jsonl.gz|{r['family']}|{c}")
    write_sheet(wb.create_sheet("Row_PDB_links"), links, LINK_COLS,
                lambda r, c: f"seed_rows.jsonl.gz|{r['family']}|{r['seed_row_id']}|{c}")
    write_sheet(wb.create_sheet("PDB_leads"), leads, LEAD_COLS, lambda r, c: f"pdb_leads.tsv|{r['family']}|{c}")
    wb.save(os.path.join(OUT, "rfam15.1_step1_inventory.xlsx"))
    return fams, links, leads


ABOUT = [
    "Step 1 documents what every Rfam 15.1 seed record says. It does not judge structure suitability; no family is",
    "excluded and no pair is selected. Explicit link = #=GR <seed row> <PDB>_<chain>_SS. PDB_leads = PDB mentions",
    "without an unambiguous row and chain. Full sequences for every row: seed_rows.jsonl.gz; column tracks:",
    "family_tracks.jsonl.gz. A cell showing SEE_COMPANION_FILE|<key> holds a value longer than Excel's 32,767 limit.",
]


# ------------------------------------------------------------------ validation
def validate(d, fams, links, leads, readme_path, prev_dir):
    V = collections.OrderedDict()
    accs = [f["family"] for f in fams]
    m = re.search(r"Rfam\s+15\.1\s+contains\s+(\d+)\s+families", open(readme_path, encoding="latin-1").read())
    release_n = int(m.group(1)) if m else None
    V["families_parsed"] = len(accs)
    V["unique_accessions"] = len(set(accs))
    V["release_readme_family_count"] = release_n
    V["coverage_each_family_once"] = len(accs) == len(set(accs)) and all(a for a in accs)
    V["coverage_matches_release_count"] = release_n == len(accs)

    kinds = collections.Counter(a["kind"] for a in d["anomalies"])
    V["anomaly_counts"] = dict(kinds)
    integrity_kinds = {"unequal_row_lengths", "GC_width_mismatch", "GR_width_mismatch", "duplicate_key_in_block",
                       "orphan_GR", "orphan_GS", "malformed_sequence_line", "malformed_GR", "malformed_GC",
                       "missing_terminator", "line_outside_record", "leading_whitespace_line"}
    V["alignment_integrity_no_structural_anomalies"] = not any(k in integrity_kinds for k in kinds)
    V["all_rows_have_family_width"] = all(t["alignment_width"] is not None for t in d["tracks"])
    V["every_family_has_RF_and_SS_cons"] = all("RF" in t["gc"] and "SS_cons" in t["gc"] for t in d["tracks"])
    fb = [a for a in d["anomalies"] if a["kind"] == "latin1_fallback_decode"]
    V["decode_fallback_lines"] = len(fb)
    V["decode_fallbacks_only_in_GF_free_text"] = all(a["detail"].startswith("#=GF") for a in fb)
    V["gapped_alphabet"] = dict(sorted(d["alphabet"].items()))

    rows = {(r["family"], r["seed_row_id"]): r for r in d["seed_rows"]}
    V["seed_rows_total"] = len(d["seed_rows"])
    V["seed_row_keys_unique"] = len(rows) == len(d["seed_rows"])
    V["links_total"] = len(links)
    V["links_resolve_to_exact_seed_row"] = all((x["family"], x["seed_row_id"]) in rows for x in links)
    V["link_sequences_equal_seed_row"] = all(
        rows[(x["family"], x["seed_row_id"])]["original_gapped_sequence"] == x["original_gapped_sequence"]
        for x in links)
    V["link_keys_unique"] = len({(x["family"], x["seed_row_id"], x["pdb_id"], x["chain_id"]) for x in links}) == len(links)
    V["ungapped_equals_gapped_minus_gaps"] = all(ungap(r["original_gapped_sequence"]) == r["ungapped_sequence"]
                                                for r in d["seed_rows"])
    V["rows_with_T"] = sum(r["contains_T"] for r in d["seed_rows"])
    V["rows_with_non_ACGU"] = sum(bool(r["non_ACGU_symbols"]) for r in d["seed_rows"])
    V["rows_with_lowercase_residues"] = sum(r["has_lowercase"] for r in d["seed_rows"])
    rows_per_link = collections.Counter((x["family"], x["seed_row_id"]) for x in links)
    V["rows_linked_to_multiple_chains"] = sum(1 for v in rows_per_link.values() if v > 1)
    fam_per_chain = collections.defaultdict(set)
    for x in links:
        fam_per_chain[(x["pdb_id"], x["chain_id"])].add(x["family"])
    V["chains_linked_in_multiple_families"] = sorted(f"{p}_{c}:{','.join(sorted(f))}"
                                                     for (p, c), f in fam_per_chain.items() if len(f) > 1)
    V["families_with_links"] = sum(1 for f in fams if f["linked_rows"])
    V["families_without_links"] = sum(1 for f in fams if not f["linked_rows"])
    V["leads_total"] = len(leads)
    V["leads_by_source"] = dict(collections.Counter(x["evidence_source"] for x in leads))

    # summary counts agree with detail
    by = collections.defaultdict(list)
    for x in links:
        by[x["family"]].append(x)
    lead_n = collections.Counter(x["family"] for x in leads)
    ok = True
    for f in fams:
        L = by.get(f["family"], [])
        ok &= f["linked_rows"] == len({x["seed_row_id"] for x in L})
        ok &= f["distinct_linked_sequences"] == len({x["sequence_group_hash"] for x in L})
        ok &= f["distinct_pdb_entries"] == len({x["pdb_id"] for x in L})
        ok &= f["distinct_pdb_chains"] == len({(x["pdb_id"], x["chain_id"]) for x in L})
        ok &= f["unresolved_leads"] == lead_n.get(f["family"], 0)
    V["family_summary_agrees_with_details"] = ok

    # export integrity: reopen everything
    V.update(export_roundtrip(d, fams, links, leads))

    # informational comparison with the previous (Nov 2026) inventory, if present
    prev = os.path.join(prev_dir, "outputs", "structure_links.csv")
    if os.path.exists(prev):
        with open(prev, newline="", encoding="utf-8") as f:
            pr = list(csv.DictReader(f))
        cols = pr[0].keys() if pr else []
        fk = next((c for c in cols if c.lower() in ("family", "accession", "rfam_acc", "family_acc", "rfam_accession")), None)
        rk = next((c for c in cols if "row" in c.lower() and "id" in c.lower()), None)
        pk = next((c for c in cols if c.lower() in ("pdb_id", "pdb")), None)
        ck = next((c for c in cols if c.lower() in ("chain", "chain_id", "pdb_chain", "auth_chain")), None)
        if all((fk, rk, pk, ck)):
            old = {(r[fk], r[rk], r[pk].upper(), r[ck]) for r in pr}
            new = {(x["family"], x["seed_row_id"], x["pdb_id"], x["chain_id"]) for x in links}
            V["previous_inventory_comparison"] = dict(previous_links=len(old), current_links=len(new),
                                                      only_previous=sorted(map("|".join, old - new))[:20],
                                                      only_current=sorted(map("|".join, new - old))[:20])
        else:
            V["previous_inventory_comparison"] = dict(note="columns not recognized", columns=list(cols))
    return V


def export_roundtrip(d, fams, links, leads):
    R = {}
    t_f, t_l, t_d = (read_tsv(os.path.join(OUT, n)) for n in ("families.tsv", "row_pdb_links.tsv", "pdb_leads.tsv"))
    R["tsv_roundtrip_families"] = t_f == [{c: str(r.get(c, "")) for c in FAM_COLS} for r in fams]
    R["tsv_roundtrip_links"] = t_l == [{c: str(r.get(c, "")) for c in LINK_COLS} for r in links]
    R["tsv_roundtrip_leads"] = t_d == [{c: str(r.get(c, "")) for c in LEAD_COLS} for r in leads]
    with gzip.open(os.path.join(OUT, "seed_rows.jsonl.gz"), "rt", encoding="utf-8") as f:
        back = [json.loads(x) for x in f]
    R["jsonl_seed_rows_roundtrip"] = back == json.loads(json.dumps(d["seed_rows"]))
    with gzip.open(os.path.join(OUT, "family_tracks.jsonl.gz"), "rt", encoding="utf-8") as f:
        backt = [json.loads(x) for x in f]
    R["jsonl_family_tracks_roundtrip"] = backt == json.loads(json.dumps(d["tracks"]))

    wb = openpyxl.load_workbook(os.path.join(OUT, "rfam15.1_step1_inventory.xlsx"), read_only=False)
    comp = {(r["family"], r["seed_row_id"]): r for r in back}
    text_ok, val_ok, n_comp = True, True, 0
    for name, rows, cols in (("Families", fams, FAM_COLS), ("Row_PDB_links", links, LINK_COLS),
                             ("PDB_leads", leads, LEAD_COLS)):
        ws = wb[name]
        got = list(ws.iter_rows(values_only=False))
        val_ok &= [c.value for c in got[0]] == cols and len(got) - 1 == len(rows)
        for r, cells in zip(rows, got[1:]):
            for c, cell in zip(cols, cells):
                want = r.get(c, "")
                v = "" if cell.value is None else cell.value
                if isinstance(v, str) and v.startswith(COMPANION + "|"):
                    n_comp += 1
                    if name == "Row_PDB_links":
                        val_ok &= comp[(r["family"], r["seed_row_id"])].get(c, want) == want
                    continue
                if c in TEXT_COLS and want != "":
                    text_ok &= cell.data_type == "s"
                val_ok &= v == want
    R["xlsx_identifiers_are_text"] = text_ok
    R["xlsx_values_equal_records"] = val_ok
    R["xlsx_cells_redirected_to_companion"] = n_comp
    return R


def main():
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    src = next(s for s in cfg["sources"] if s["key"] == "rfam_seed")
    seed = os.path.join(ROOT, src["path"])
    file_sha = sha256_file(seed)
    if file_sha != src["sha256"]:
        sys.exit(f"Rfam.seed.gz sha256 {file_sha} != pinned {src['sha256']}")
    run_utc = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    d = build(seed)
    prov = collections.OrderedDict(
        step="Step 1 - master inventory of explicit seed-row to PDB-chain links",
        rfam_release=cfg["reference"]["release"], source_url=src["url"], local_path=src["path"],
        seed_gz_sha256=file_sha, seed_gz_bytes=os.path.getsize(seed), decompressed_bytes=d["raw_len"],
        decompressed_sha256=d["raw_sha"], run_utc=run_utc, parser_version=PARSER_VERSION,
        parser_script_sha256=sha256_file(os.path.abspath(__file__)),
        link_rule=LINK_RE.pattern, sequence_hash_rule="sha256(ascii(uppercase(ungapped))); no T->U",
        ungap_rule="remove '-' and '.' only", python=sys.version.split()[0], openpyxl=openpyxl.__version__)
    fams, links, leads = export(d, prov)
    V = validate(d, fams, links, leads, os.path.join(ROOT, "inputs", "rfam", "15.1", "README"),
                 os.path.normpath(os.path.join(ROOT, "..", "v3_full Rfam", "rfam_seed_inventory_15_1")))
    gates = [k for k, v in V.items() if isinstance(v, bool)]
    failed = [k for k in gates if not V[k]]
    json.dump(dict(provenance=prov, checks=V, failed_checks=failed), open(os.path.join(OUT, "validation.json"), "w"),
              indent=2, default=str)
    json.dump(prov, open(os.path.join(OUT, "provenance.json"), "w"), indent=2)
    print(json.dumps({k: V[k] for k in V if k not in ("gapped_alphabet",)}, indent=1, default=str)[:4000])
    print(f"\n{len(gates) - len(failed)}/{len(gates)} gate checks passed" + (f"; FAILED {failed}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
