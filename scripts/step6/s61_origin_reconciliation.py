"""Step 6.1: every origin check ever cached, its R7 state, its R8 state and what it decides for the cohort.

  python scripts/step6/s61_origin_reconciliation.py [<repo_root> <out_dir>]
Read-only on the project. Writes origin_check_reconciliation.tsv and origin_summary.json.
"""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import collections
import csv
import glob
import json
import os
import sys

ROOT, OUT = _ARGS[0] if len(_ARGS) > 0 else _REPO, _ARGS[1] if len(_ARGS) > 1 else _os.path.join(_REPO, "outputs", "06_audit", "step2_r8")
os.makedirs(OUT, exist_ok=True)
S2 = os.path.join(ROOT, "outputs", "02_structure_preparation")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


cur = {(r["seq_sha256"][:20] + "_" + r["species_taxid"]): r for r in tsv(f"{S2}/snapshots/step2_r8/sequence_evidence.tsv")}
cache = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob(f"{ROOT}/cache/origin/*.json")}
loc = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob(f"{ROOT}/cache/origin_exact_local/*.json")}
closed = {u["check_id"]: u for u in json.load(open(f"{S2}/STEP2_CLOSED.json"))["unresolved_checks"]}
CR = {r["candidate_id"]: r for r in tsv(f"{S2}/snapshots/step2_r8/chain_review.tsv")}
CF = {r["candidate_id"]: r for r in tsv(f"{S2}/snapshots/step2_final/chain_review.tsv")}
jobs = tsv(f"{ROOT}/outputs/03_star3d/step2_r8/step3_jobs.tsv")
used = {j[f"{s}_candidate"] for j in jobs for s in "ab"}
mega = {os.path.basename(p).split("_")[0] for p in glob.glob(f"{ROOT}/cache/blast/*.megablast.xml")}

by_key = collections.defaultdict(list)
for c, r in CR.items():
    if r["evidence_key"]:
        sha, sp = r["evidence_key"].split("|")
        by_key[sha[:20] + "_" + sp].append(c)

rows = []
for k in sorted(set(cache) | set(cur) | set(closed)):
    d = cache.get(k, {})
    l = loc.get(k, {})
    r8_eval = d.get("evaluated_utc", "") >= "2026-11-25"
    routes = sorted({r.get("route", "") for r in d.get("routes_tried", [])})
    sp = k.split("_")[1]
    route_r8 = ("megablast" if (r8_eval and sp in mega and "blast_core_nt" in routes) else
                ("blastn" if "blast_core_nt" in routes else ("assembly" if "assembly" in routes else "")))
    if k in closed:
        r7 = "SEARCH_INCOMPLETE (R7 open)"
    elif l.get("exact_local_resolves"):
        r7 = "whole_chain_found (exact local genome search)"
    else:
        r7 = d.get("outcome", "")
    chains = by_key.get(k, [])
    dec_r8 = collections.Counter(CR[c]["decision"] + ":" + CR[c]["reason"] for c in chains)
    dec_r7 = collections.Counter(CF[c]["decision"] + ":" + CF[c]["reason"] for c in chains)
    in_cohort = sorted(c for c in chains if c in used)
    good = [h for h in d.get("whole_chain_hits", []) if h.get("independent") == "yes" and h.get("exact_confirmed") == "yes"]
    lgood = l.get("whole_chain_hits_confirmed", []) if l.get("exact_local_resolves") else []
    basis = sorted({h.get("independence_basis", "") for h in good} | {"genome_assembly(local)" for _ in lgood})
    recs = [h.get("record") for h in good][:3] + [h.get("record") for h in lgood][:2]
    if k not in cur:
        disp = "removed_from_queue: chain has no possible partner or key superseded; cached result kept; chain not used"
    elif any(v.startswith("verified_natural") for v in dec_r8):
        disp = "exact_independent_sequence_support"
    elif any("engineered" in v for v in dec_r8):
        disp = "documented_engineering_or_non_native_flank"
    elif cur[k]["outcome"] == "unknown_residues_in_chain":
        disp = "not_searchable_unknown_residues"
    elif cur[k]["outcome"] in ("not_found", "interval_found_flank_differs"):
        disp = "no_whole_chain_match_in_searched_references (origin unresolved; not used)"
    elif cur[k]["outcome"] == "whole_chain_found":
        disp = "exact_independent_sequence_support (no chain decision uses it)"
    else:
        disp = cur[k]["outcome"]
    if not chains:
        disp += " [no chain currently references this key]"
    rows.append(dict(
        check_id=k, species_taxid=sp, species=(d or l).get("species_name", cur.get(k, {}).get("species_name", "")),
        families=";".join(d.get("families", [])), deposited_length=d.get("deposited_length", ""),
        query_sha256=d.get("seq_sha256", ""), in_current_queue=k in cur, open_at_R7=k in closed,
        state_at_R7=r7, search_completed=bool(d) and d.get("outcome") is not None,
        reevaluated_R8=r8_eval, route=route_r8, ncbi_outcome_now=cur.get(k, {}).get("outcome", d.get("outcome", "")),
        local_exact_outcome=l.get("outcome", ""), supporting_records=";".join(x for x in recs if x),
        independence_basis=";".join(basis), n_chains=len(chains), chains_in_cohort=len(in_cohort),
        decisions_R7=json.dumps(dict(dec_r7)), decisions_R8=json.dumps(dict(dec_r8)), disposition=disp))

with open(f"{OUT}/origin_check_reconciliation.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# evidence strength for the 138 cohort chains
strength = collections.Counter()
weak = []
for c in sorted(used):
    ev = CR[c]["origin_evidence"]
    b = ("genome_assembly" if ("genome_assembly" in ev or "assembly:" in ev) else
         "long_record(>=chain+50nt)" if "long_record" in ev else
         "older_record_only" if "record_predates_structure" in ev else "other")
    strength[b] += 1
    if b == "older_record_only":
        weak.append((c, ev))
summ = dict(
    keys_total=len(rows), keys_in_current_queue=len(cur), open_at_R7=len(closed),
    open_at_R7_now=dict(collections.Counter(cur[k]["outcome"] for k in closed)),
    current_outcomes=dict(collections.Counter(r["outcome"] for r in cur.values())),
    dispositions=dict(collections.Counter(r["disposition"].split(" [")[0] for r in rows)),
    reevaluated_R8=sum(r["reevaluated_R8"] for r in rows),
    cohort_chains=len(used), cohort_chain_evidence=dict(strength), cohort_chains_older_record_only=weak,
    keys_used_by_cohort=sum(1 for r in rows if r["chains_in_cohort"]))
json.dump(summ, open(f"{OUT}/origin_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != "cohort_chains_older_record_only"}, indent=1))
print("older-record-only cohort chains:", len(weak)); [print(" ", c, e[:90]) for c, e in weak]
