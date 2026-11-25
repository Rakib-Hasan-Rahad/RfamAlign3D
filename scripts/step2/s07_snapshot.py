"""Step 2 snapshots and their comparison with a later export. Reads outputs only; never touches running stages.

  python scripts/step2/s07_snapshot.py snapshot   # freeze the CURRENT on-disk export into snapshots/interim_<UTC>/
  python scripts/step2/s07_snapshot.py compare    # compare current pairs_to_run.tsv with the latest snapshot

Pair identity for comparison = family + unordered pair of the two seed-row ungapped sequences (seq_pair_key);
'replaced' = same key, different representative structures. Snapshot files are made read-only.
"""
import collections
import glob
import gzip
import hashlib
import json
import os
import shutil
import stat
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

COPY = ["SQARMA_Step_2_STAR3D_inputs.xlsx", "pairs_to_run.tsv", "rfam_pair_reference.tsv.gz", "residue_crosswalk.tsv.gz",
        "chain_review.tsv", "family_status.tsv", "yield_by_family.tsv", "pair_candidates.tsv", "selection_state.json",
        "validation.json", "validation_report.md", "professor_summary.md", "step2_policy.yaml", "structures.tsv",
        "sequence_evidence.tsv", "input_reconciliation.json"]


def seed_seqs():
    d = {}
    with gzip.open(os.path.join(C.STEP1, "seed_rows.jsonl.gz"), "rt") as f:
        for line in f:
            r = json.loads(line)
            d[(r["family"], r["seed_row_id"])] = r["ungapped_sequence"].upper()
    return d


def key_of(p, seqs):
    a = hashlib.sha256(seqs[(p["family"], p["query_seed_row"])].encode()).hexdigest()[:12]
    b = hashlib.sha256(seqs[(p["family"], p["target_seed_row"])].encode()).hexdigest()[:12]
    return f"{p['family']}|{'|'.join(sorted((a, b)))}"


def members(p):
    return "|".join(sorted((p["query_candidate"], p["target_candidate"])))


def snapshot(name=None):
    src_time = C.jload(C.P("validation.json"))["utc"]
    ts = C.now().replace("-", "").replace(":", "")
    out = C.P("snapshots", name or f"interim_{ts}")
    os.makedirs(out)
    for f in COPY:
        if os.path.exists(C.P(f)):
            shutil.copy2(C.P(f), os.path.join(out, f))
    pairs = C.read_tsv(C.P("pairs_to_run.tsv"))
    fam = {r["family"]: r for r in C.read_tsv(C.P("family_status.tsv"))}
    seqs = seed_seqs()
    groups = collections.defaultdict(list)
    for p in pairs:
        g = "frozen_family_complete" if fam[p["family"]]["family_complete"] == "yes" else "eligible_family_incomplete"
        groups[g].append(dict(p, snapshot_group=g, seq_pair_key=key_of(p, seqs)))
    cols = list(pairs[0]) + ["snapshot_group", "seq_pair_key"] if pairs else ["pair_id"]
    C.write_tsv(os.path.join(out, "pairs_frozen.tsv"), groups["frozen_family_complete"], cols)
    C.write_tsv(os.path.join(out, "pairs_selection_may_change.tsv"), groups["eligible_family_incomplete"], cols)
    # Rfam references split by group
    gof = {p["pair_id"]: g for g, ps in groups.items() for p in ps}
    ref = C.read_tsv(C.P("rfam_pair_reference.tsv.gz"))
    for g in groups:
        C.write_tsv(os.path.join(out, f"rfam_pair_reference__{g}.tsv.gz"), [r for r in ref if gof[r["pair_id"]] == g],
                    list(ref[0]))
    # staged inputs: copy, re-hash, group membership
    sdir = os.path.join(out, "inputs_star3d")
    os.makedirs(sdir)
    man = {m["stem"]: m for m in C.read_tsv(C.P("inputs_star3d", "file_manifest.tsv"))}
    shutil.copy2(C.P("inputs_star3d", "file_manifest.tsv"), os.path.join(sdir, "file_manifest.tsv"))
    used = collections.defaultdict(set)
    for g, ps in groups.items():
        for p in ps:
            used[p["query_stem"]].add(g)
            used[p["target_stem"]].add(g)
    hashes, bad = [], []
    for stem, m in sorted(man.items()):
        for ext in (".pdb", ".conversion.tsv"):
            shutil.copy2(C.P("inputs_star3d", stem + ext), os.path.join(sdir, stem + ext))
        h = C.sha256(os.path.join(sdir, stem + ".pdb"))
        if h != m["sha256"]:
            bad.append(stem)
        hashes.append(dict(stem=stem, groups=";".join(sorted(used.get(stem, {"unused"}))), file=f"inputs_star3d/{stem}.pdb",
                           sha256=h, manifest_sha256=m["sha256"], matches_manifest=h == m["sha256"],
                           pdb_id=m["pdb_id"], auth_chain=m["auth_chain"], staged_chain=m["staged_chain"],
                           label_interval=m["label_interval"], source_mmcif_sha256=m["source_mmcif_sha256"]))
    C.write_tsv(os.path.join(out, "input_file_hashes.tsv"), hashes)
    counts = {g: dict(pairs=len(ps), staged_files=len({s for p in ps for s in (p["query_stem"], p["target_stem"])}),
                      families=sorted({p["family"] for p in ps})) for g, ps in groups.items()}
    v = C.jload(C.P("validation.json"))
    fs = C.read_tsv(C.P("family_status.tsv"))
    info = dict(snapshot_utc=C.now(), source_export_utc=src_time, pairs_total=len(pairs), groups=counts,
                staged_files_total=len(man), staged_hash_mismatches=bad,
                validation_passed=v["passed"],
                gates={k: x for k, x in v["gates"].items() if isinstance(x, bool)},
                families_complete=sum(1 for r in fs if r["family_complete"] == "yes"),
                families_incomplete=[r["family"] for r in fs if r["family_complete"] == "no"],
                star3d_run=False)
    C.jdump(info, os.path.join(out, "snapshot_info.json"))
    md = [f"# Step 2 INTERIM snapshot {info['snapshot_utc']}", "",
          f"Frozen copy of the Step 2 export validated at {src_time} (all {len(info['gates'])} acceptance gates: "
          f"{'PASS' if v['passed'] else 'NOT ALL PASS'}). Evidence checks were still running; STAR3D has not been run.", "",
          "| Group | Pairs | Staged files | Families |", "|---|---:|---:|---|"]
    for g, c in counts.items():
        md.append(f"| {g} | {c['pairs']} | {c['staged_files']} | {', '.join(c['families'])} |")
    md += ["", "- **frozen_family_complete**: the family has an outcome for every sequence pair; these selections are final.",
           "- **eligible_family_incomplete**: both members are verified and usable now, and every combination that could "
           "outrank the selection was already evaluated; but the family still has pending sequence pairs, so the final "
           "export re-derives it and reports any change against this snapshot.",
           "- Chains with no confirmed whole-chain record (not_found) are provisional_origin and are not in any pair.", "",
           "Files: `pairs_frozen.tsv`, `pairs_selection_may_change.tsv`, `rfam_pair_reference__<group>.tsv.gz`, "
           "`residue_crosswalk.tsv.gz`, `inputs_star3d/` + `input_file_hashes.tsv`, workbook, `validation_report.md`, "
           "`snapshot_info.json`, `snapshot_manifest.tsv` (SHA-256 of every file here)."]
    open(os.path.join(out, "SNAPSHOT.md"), "w").write("\n".join(md) + "\n")
    allf = sorted(p for p in glob.glob(os.path.join(out, "**", "*"), recursive=True) if os.path.isfile(p))
    C.write_tsv(os.path.join(out, "snapshot_manifest.tsv"),
                [dict(file=os.path.relpath(p, out), sha256=C.sha256(p), bytes=os.path.getsize(p)) for p in allf])
    for p in glob.glob(os.path.join(out, "**", "*"), recursive=True) + [out]:
        mode = os.stat(p).st_mode
        os.chmod(p, mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    print(json.dumps(info, indent=1))
    print("snapshot:", os.path.relpath(out, C.ROOT))


def compare():
    snaps = sorted(glob.glob(C.P("snapshots", "*")), key=os.path.getmtime)
    if not snaps:
        print("no snapshot")
        return
    snap = snaps[-1]
    seqs = seed_seqs()
    old = {}
    for f in ("pairs_frozen.tsv", "pairs_selection_may_change.tsv"):
        for p in C.read_tsv(os.path.join(snap, f)):
            old[p["seq_pair_key"]] = p
    new = {key_of(p, seqs): p for p in C.read_tsv(C.P("pairs_to_run.tsv"))}
    rows = []
    for k in sorted(set(old) | set(new)):
        o, n = old.get(k), new.get(k)
        st = ("added" if not o else "removed" if not n else
              "replaced" if members(o) != members(n) else "unchanged")
        rows.append(dict(seq_pair_key=k, family=(n or o)["family"], change=st,
                         snapshot_group=o["snapshot_group"] if o else "", snapshot_pair_id=o["pair_id"] if o else "",
                         final_pair_id=n["pair_id"] if n else "",
                         snapshot_members=members(o) if o else "", final_members=members(n) if n else "",
                         files_changed=("" if not (o and n) or st != "unchanged" else
                                        "yes" if (o["query_sha256"], o["target_sha256"]) != (n["query_sha256"], n["target_sha256"]) else "no")))
    C.write_tsv(C.P("snapshot_comparison.tsv"), rows)
    c = collections.Counter(r["change"] for r in rows)
    res = dict(snapshot=os.path.relpath(snap, C.ROOT), counts=dict(c),
               frozen_pairs_changed=[r["seq_pair_key"] for r in rows if r["snapshot_group"] == "frozen_family_complete"
                                     and r["change"] != "unchanged"],
               unchanged_but_file_hash_changed=[r["seq_pair_key"] for r in rows if r["files_changed"] == "yes"])
    C.jdump(res, C.P("snapshot_comparison.json"))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "snapshot":
        snapshot(sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        compare()
