"""Read-only audit of the frozen step2_interim_440 results (no STAR3D rerun; Step 2 untouched).
Writes outputs/04_comparison/step2_interim_440/audit/ and prints every figure with its source."""
import collections
import csv
import glob
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s41_compare as K  # noqa: E402

C = sys.argv[1] if len(sys.argv) > 1 else "step2_interim_440"
R = K.ROOT
B = os.path.join(R, "outputs", "04_comparison", C)
M = os.path.join(B, "motif")
RET = os.path.join(R, "outputs", "03_star3d", C, "returned")
SNAP = os.path.join(R, "outputs", "02_structure_preparation", "snapshots", C)
A = os.path.join(B, "audit")
os.makedirs(A, exist_ok=True)
out = {}

res = K.tsv(os.path.join(B, "residue_comparisons.tsv.gz"))
edges = K.tsv(os.path.join(B, "alignment_edges.tsv.gz"))
oc = K.tsv(os.path.join(B, "run_outcomes.tsv"))
obs = K.tsv(os.path.join(M, "disagreement_observations.tsv"))
ev = K.tsv(os.path.join(M, "overlap_evidence.tsv.gz"))
pairs = {p["pair_id"]: p for p in K.tsv(os.path.join(SNAP, "pairs_to_run.tsv"))}

# ---------------------------------------------------------------- Q1 denominators
a = [r for r in res if r["side"] == "A"]
cat = collections.Counter(r["category"] for r in a)
both = [r for r in a if r["category"] in ("same_partner", "different_partner")]
ass_same = [r for r in a if r["category"] == "same_partner" and r["coords"] == "yes" and r["rfam_partner_coords"] == "yes"]
ass_diff = [r for r in a if r["scored_different"] == "yes"]
# assessable for a different-partner row needs A coords + Rfam partner coords + STAR3D partner coords
diff_not = [r for r in a if r["category"] == "different_partner" and r["scored_different"] != "yes"]
em = collections.Counter(e["method"] for e in edges)
shared = sum(1 for e in edges if e["method"] == "rfam" and e["in_both"] == "yes")
rg = K.tsv(os.path.join(B, "regions.tsv"))
out["Q1"] = dict(
    source="residue_comparisons.tsv.gz, side=A (411 compared pairs)",
    a_side_nucleotides=len(a), categories=dict(cat),
    mapped_by_both=len(both), same_partner=cat["same_partner"], different_partner=cat["different_partner"],
    jointly_assessable=len(ass_same) + len(ass_diff), assessable_same=len(ass_same), assessable_different=len(ass_diff),
    different_not_assessable=len(diff_not), same_not_assessable=cat["same_partner"] - len(ass_same),
    disagreement_rate=f"{len(ass_diff)}/{len(ass_same) + len(ass_diff)}",
    edges=dict(source="alignment_edges.tsv.gz", rows=len(edges), rfam_rows=em["rfam"], star3d_rows=em["star3d"],
               shared_pairs=shared, union=em["rfam"] + em["star3d"] - shared,
               definition="65,919 = rows = Rfam edges + STAR3D edges (SUM; a shared residue pair appears once per method)"),
    regions=dict(source="regions.tsv", total=len(rg), by_side=dict(collections.Counter(r["side"] for r in rg))))

# ---------------------------------------------------------------- Q2 four-cell table
evo = collections.defaultdict(list)
for e in ev:
    evo[e["observation_id"]].append(e)
cells = collections.Counter()
loop_only = collections.Counter()
for o in obs:
    es = evo[o["observation_id"]]
    rf = any(e["rfam_motif"] == "hit" for e in es)
    at = any(e["bgsu_group"] == "hit" for e in es)
    cells[("Rfam" if rf else "-") + "/" + ("Atlas" if at else "-")] += 1
    if o["classification"] == "structural_loop_only":
        st = {e["role"]: e["bgsu_group"] for e in es}
        loop_only["all 3 endpoints in Atlas-covered chains (searched, no group)" if all(v == "assessed_no_group" for v in st.values())
                  else "no endpoint in an Atlas-covered chain (Atlas did not cover)" if all(v == "not_covered" for v in st.values())
                  else "mixed coverage across endpoints"] += 1
out["Q2"] = dict(source="motif/overlap_evidence.tsv.gz grouped by observation_id; hit = ANY of the 3 endpoints",
                 rfam_only=cells["Rfam/-"], atlas_only=cells["-/Atlas"], both=cells["Rfam/Atlas"], neither=cells["-/-"],
                 total=sum(cells.values()), any_rfam=cells["Rfam/-"] + cells["Rfam/Atlas"],
                 any_atlas=cells["-/Atlas"] + cells["Rfam/Atlas"], loop_only_breakdown=dict(loop_only),
                 classification=dict(collections.Counter(o["classification"] for o in obs)))

# ---------------------------------------------------------------- Q3 baseline, recomputed per source, query-only
qd = {o["observation_id"]: next(e for e in evo[o["observation_id"]] if e["role"] == "A_residue") for o in obs}
# agreements: rebuild query-endpoint evidence with the same lookup code (s42.endpoint not importable cheaply) ->
# use motif_baseline.tsv rows, and recompute disagreements directly for the exact rule.
bl = K.tsv(os.path.join(M, "motif_baseline.tsv"))
pooled = {(b["source"], b["endpoints"]): b for b in bl if b["scope"] == "pooled"}
d_rfam_q = sum(1 for e in qd.values() if e["rfam_motif"] == "hit")
d_atlas_cov = [e for e in qd.values() if e["bgsu_group"] in ("hit", "assessed_no_group")]
d_atlas_q = sum(1 for e in d_atlas_cov if e["bgsu_group"] == "hit")
out["Q3"] = dict(
    source="motif/motif_baseline.tsv (pooled) + overlap_evidence.tsv.gz A_residue rows",
    rfam_query_only=dict(disagreement=pooled[("rfam_motif", "query_only")]["disagreement_hits"],
                         agreement=pooled[("rfam_motif", "query_only")]["agreement_hits"],
                         recomputed_disagreement=f"{d_rfam_q}/{len(qd)}"),
    bgsu_loop_query_only=dict(disagreement=pooled[("bgsu_loop", "query_only")]["disagreement_hits"],
                              agreement=pooled[("bgsu_loop", "query_only")]["agreement_hits"]),
    atlas_query_only_baseline=dict(disagreement=pooled[("bgsu_atlas_group", "query_only")]["disagreement_hits"],
                                   agreement=pooled[("bgsu_atlas_group", "query_only")]["agreement_hits"],
                                   meaning="member = Atlas group at the query endpoint; assessable = query endpoint in an "
                                           "Atlas-covered chain (the earlier mixed Rfam-or-Atlas row was removed)"),
    atlas_query_only_disagreement=f"{d_atlas_q}/{len(d_atlas_cov)}",
    note_39_vs_81="39 = observations with an Atlas group at ANY of 3 endpoints, no coverage filter; 81 = observations whose "
                  "QUERY endpoint lies in an Atlas-covered chain. Different endpoint rules -> 39/81 is not a valid rate.")

# ---------------------------------------------------------------- Q4 per family
fam = collections.defaultdict(lambda: collections.Counter())
for p in pairs.values():
    fam[p["family"]]["pairs"] += 1
for o in oc:
    f = o["family"]
    fam[f]["valid"] += o["step3_status"] == "valid_alignment"
    fam[f]["no_alignment"] += o["step3_status"] == "no_alignment"
for r in ass_same + ass_diff:
    f = pairs[r["pair_id"]]["family"]
    fam[f]["assessable"] += 1
    fam[f]["different"] += r["scored_different"] == "yes"
for o in obs:
    es = evo[o["observation_id"]]
    fam[o["family"]]["rfam_hit"] += any(e["rfam_motif"] == "hit" for e in es)
    fam[o["family"]]["atlas_hit"] += any(e["bgsu_group"] == "hit" for e in es)
cols = ["pairs", "valid", "no_alignment", "assessable", "different", "rfam_hit", "atlas_hit"]
rows = []
for f in sorted(fam):
    rows.append(dict(family=f, family_id=next(p["family_id"] for p in pairs.values() if p["family"] == f), **{c: fam[f][c] for c in cols}))
sub = {c: sum(r[c] for r in rows if r["family"] not in ("RF00005", "RF00001")) for c in cols}
rows.append(dict(family="SUBTOTAL excl. RF00005+RF00001", family_id="", **sub))
rows.append(dict(family="TOTAL", family_id="", **{c: sum(r[c] for r in rows[:-1]) for c in cols}))
K.wtsv(os.path.join(A, "family_table.tsv"), rows)
na = [o for o in oc if o["step3_status"] == "no_alignment"]
struct = collections.Counter()
for o in na:
    p = pairs[o["pair_id"]]
    struct[p["query_pdb_chain"]] += 1
    struct[p["target_pdb_chain"]] += 1
out["Q4"] = dict(source="snapshot pairs_to_run.tsv + run_outcomes.tsv + residue_comparisons + overlap_evidence",
                 table="audit/family_table.tsv", no_alignment_by_family=dict(collections.Counter(o["family"] for o in na)),
                 no_alignment_structures_top=struct.most_common(6), no_alignment_distinct_structures=len(struct))

# ---------------------------------------------------------------- Q6 input evidence
pre = [json.load(open(p)) for p in glob.glob(os.path.join(RET, "preprocess", "*", "preprocess.json"))]
dl = [p for p in glob.glob(os.path.join(RET, "preprocess", "*", "preprocess.stdout")) if "Downloading" in open(p).read()]
dl += [p for p in glob.glob(os.path.join(RET, "runs", "*", "align.stdout")) if "Downloading" in open(p).read()]
manifest = {m["stem"].lower(): m for m in K.tsv(os.path.join(SNAP, "inputs_star3d", "file_manifest.tsv"))}
conv = {}
for m in manifest.values():
    rows_ = K.tsv(os.path.join(SNAP, "inputs_star3d", f"{m['stem']}.conversion.tsv"))
    conv[m["stem"].lower()] = [(r["staged_chain"], r["auth_seq_num"] + r["auth_ins_code"]) for r in rows_]
mca_ok, mca_bad = 0, []
for p in pre:
    sid = p["key"].split("__")[0]
    mca = glob.glob(os.path.join(RET, "preprocess", p["key"], "STAR3D_struct_info", "*.mca"))[0]
    lines = open(mca).read().split("Residue conformations")[1].split("---\n", 1)[1].split("\n\n")[0]
    got = []
    for l in lines.splitlines():
        m = re.match(r"^'?([A-Za-z0-9])'?(-?\d+)(?:\.([A-Za-z]))?\s*:", l.strip())
        if m:
            got.append((m.group(1), m.group(2) + (m.group(3) or "")))
    want = conv[sid]
    if got == want:
        mca_ok += 1
    else:
        mca_bad.append((p["key"], len(got), len(want), got[:2], want[:2]))
ids = set()
for r in glob.glob(os.path.join(RET, "runs", "*", "result.json")):
    j = json.load(open(r))
    ids |= {x.split(":")[0] for x in j["identity"]["direction"].split(" -> ")}
ps = K.tsv(os.path.join(B, "pair_scores.tsv"))
out["Q6"] = dict(source="returned/preprocess/*/preprocess.json, *.stdout, STAR3D_struct_info/*.mca, runs/*/result.json; "
                        "snapshot inputs_star3d/*.conversion.tsv",
                 preprocess_records=len(pre), preprocess_input_sha_equals_frozen=sum(
                     p["identity"]["input_sha256"] == manifest[p["key"].split("__")[0]]["sha256"] for p in pre),
                 download_messages_found=len(dl),
                 mca_residue_list_equals_staged_file=f"{mca_ok}/{len(pre)}", mca_mismatches=mca_bad[:3],
                 star3d_ids_all_lowercase=all(i == i.lower() for i in ids), distinct_ids=len(ids),
                 jobs_inputs_changed=sum(bool(json.load(open(r)).get("inputs_changed")) for r in glob.glob(os.path.join(RET, "runs", "*", "result.json"))),
                 rf_track=dict(collections.Counter(p["rf_available"] for p in ps)),
                 rf_track_families=sorted({p["family"] for p in ps if p["rf_available"] == "yes"}))
# ---------------------------------------------------------------- Q7 diagnostics: offsets, patterns, hairpins
dg = K.tsv(os.path.join(B, "diagnostics", "disagreement_diagnostics.tsv.gz"))
bres = {(r["pair_id"], r["unit_id"]): r for r in res if r["side"] == "B"}
# offset recomputed independently from B ungapped positions (crosswalk), never from author numbers
off_ok = sum(int(d["offset"]) == int(bres[(d["pair_id"], next(r["star3d_partner"] for r in ass_diff if r["pair_id"] == d["pair_id"] and r["unit_id"] == d["unit_id"]))]["ungapped_position"])
             - int(bres[(d["pair_id"], next(r["rfam_partner"] for r in ass_diff if r["pair_id"] == d["pair_id"] and r["unit_id"] == d["unit_id"]))]["ungapped_position"])
             for d in dg[:300])
lead = [d for d in dg if d["unit_id"] == "3RG5|1|A|U|47|||C"]
lr = [r for r in ass_diff if r["unit_id"] == "3RG5|1|A|U|47|||C"]
lead_info = {}
if lr:
    r = lr[0]
    lead_info = dict(rfam_partner=r["rfam_partner"], rfam_partner_pos=bres[(r["pair_id"], r["rfam_partner"])]["ungapped_position"],
                     star3d_partner=r["star3d_partner"], star3d_partner_pos=bres[(r["pair_id"], r["star3d_partner"])]["ungapped_position"],
                     offset=lead[0]["offset"] if lead else None,
                     author_number_difference=0, note="47B and 47C share author number 47; ungapped positions differ by 1")
pat = collections.Counter(d["pattern"] for d in dg)
gp = collections.defaultdict(collections.Counter)
for d in dg:
    if d["pattern"] == "global_misregistration":
        gp[d["pair_id"]][int(d["offset"])] += 1
modal = sum(c.most_common(1)[0][1] for c in gp.values())
tot_g = sum(sum(c.values()) for c in gp.values())
ass_by_ctx = collections.Counter(r["ss_context"] for r in ass_same + ass_diff)
dif_by_ctx = collections.Counter(r["ss_context"] for r in ass_diff)
hp_n, hp_d = dif_by_ctx["hairpin_loop"], ass_by_ctx["hairpin_loop"]
oth_n, oth_d = len(ass_diff) - hp_n, len(ass_same) + len(ass_diff) - hp_d
out["Q7"] = dict(source="diagnostics/disagreement_diagnostics.tsv.gz + residue_comparisons.tsv.gz",
                 offset_definition="ungapped position (Step 2 crosswalk, insertion codes preserved) of STAR3D's B partner minus that of Rfam's B partner",
                 offset_recomputed_matches=f"{off_ok}/{min(300, len(dg))} (first 300 rows)", lead_example=lead_info,
                 patterns={k: f"{v}/{len(dg)}" for k, v in pat.items()},
                 patterns_mutually_exclusive="yes: first matching rule assigns exactly one pattern per observation",
                 pairwide_pairs=len(gp), pairwide_modal_offset_share=f"{modal}/{tot_g}",
                 hairpin=dict(hairpin=f"{hp_n}/{hp_d}", outside_hairpin=f"{oth_n}/{oth_d}",
                              rule="disagreements / jointly assessable A-side observations whose Rfam column is in a consensus hairpin loop (SS_cons)"),
                 by_context={k: f"{dif_by_ctx[k]}/{ass_by_ctx[k]}" for k in sorted(ass_by_ctx)})
hl = json.load(open(os.path.join(R, "outputs", "05_results", C, "headline.json"))) if os.path.exists(os.path.join(R, "outputs", "05_results", C, "headline.json")) else {}
out["Q8_website_reconciliation"] = dict(
    pairs=(len(oc), hl.get("pairs")), valid=(sum(o["step3_status"] == "valid_alignment" for o in oc), hl.get("valid")),
    no_stacks=(len(na), hl.get("no_stacks")), assessable=(len(ass_same) + len(ass_diff), hl.get("assessable")),
    different=(len(ass_diff), hl.get("different")),
    all_match=bool(hl) and len(oc) == hl.get("pairs") and len(ass_diff) == hl.get("different") and len(ass_same) + len(ass_diff) == hl.get("assessable"))
json.dump(out, open(os.path.join(A, "audit.json"), "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
for r in rows:
    print("\t".join(str(r[k]) for k in ["family", "family_id"] + cols))
