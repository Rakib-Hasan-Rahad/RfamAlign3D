"""One-command check that the published results are intact and that the headline numbers follow from the data.

  python scripts/verify.py [--cohort step2_final]

Needs only the Python standard library and the files in this repository (release archives are checked when
they have been unpacked; otherwise those checks are reported as SKIPPED, never as passed).

1. Integrity: every file listed in each sha256 manifest has the recorded checksum
   - Step 2 snapshot (structures, pairs, crosswalk, staged STAR3D inputs)   snapshot_manifest.tsv
   - STAR3D returns (raw alignments, superpositions, preprocessing)          RETURN.sha256   [release archive]
   - Verification archive (exact local search evidence, STAR3D outputs)      COLLECT.sha256  [release archive]
2. Recount: the headline counts are recomputed independently from residue_comparisons.tsv.gz and
   pair_scores.tsv and compared with the published overall_metric.json and docs/results.json.
Exit status 0 only if every executed check passes.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = []


def report(name, ok, detail=""):
    RESULTS.append((name, ok))
    tag = {True: "PASS", False: "FAIL", None: "SKIPPED"}[ok]
    print(f"[{tag:7s}] {name}{(' - ' + detail) if detail else ''}")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check_manifest(name, base, entries):
    if not os.path.isdir(base):
        return report(name, None, f"{os.path.relpath(base, ROOT)} not present (unpack the release archive)")
    missing = [p for p, _ in entries if not os.path.exists(os.path.join(base, p))]
    if len(missing) == len(entries):
        return report(name, None, "files not present (unpack the release archive)")
    bad = [p for p, h in entries if os.path.exists(os.path.join(base, p)) and sha(os.path.join(base, p)) != h]
    present = len(entries) - len(missing)
    if bad:
        return report(name, False, f"{len(bad)} of {present} present files do not match, e.g. {bad[:3]}")
    if missing:
        return report(name, None, f"{present} present files match; {len(missing)} not present "
                                  "(unpack the release archive to check them)")
    report(name, True, f"{present}/{len(entries)} files match")


def sha_list(path):
    out = []
    for line in open(path):
        if line.strip():
            h, p = line.rstrip("\n").split(None, 1)
            out.append((p.lstrip("*").strip(), h))
    return out


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="step2_r8")
    c = ap.parse_args().cohort
    snap = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", c)
    man = [(r["file"], r["sha256"]) for r in tsv(os.path.join(snap, "snapshot_manifest.tsv"))]
    check_manifest("Step 2 snapshot checksums", snap, man)
    ret = os.path.join(ROOT, "outputs", "03_star3d", c, "returned")
    check_manifest("STAR3D return checksums", ret, sha_list(os.path.join(ret, "RETURN.sha256")))
    pa = os.path.join(ROOT, "outputs", "verification_archive")
    check_manifest("Verification archive checksums", pa, sha_list(os.path.join(pa, "COLLECT.sha256")))

    # ---- independent recount from the residue-level table (RNA A view; assessability as in s45)
    base = os.path.join(ROOT, "outputs", "04_comparison", c)
    res = tsv(os.path.join(base, "residue_comparisons.tsv.gz"))
    bco = {(r["pair_id"], r["unit_id"]): r["coords"] == "yes" for r in res if r["side"] == "B"}
    n = collections.Counter()
    for r in res:
        if r["side"] != "A" or r["category"] == "neither":
            continue
        a = r["coords"] == "yes"
        rp = bco.get((r["pair_id"], r["rfam_partner"]), False)
        sp = bco.get((r["pair_id"], r["star3d_partner"]), False)
        need = {"same_partner": rp, "different_partner": rp and sp, "rfam_only": rp, "star3d_only": sp}[r["category"]]
        n[r["category"] if a and need else "unassessable"] += 1
    A, D, R, S = n["same_partner"], n["different_partner"], n["rfam_only"], n["star3d_only"]
    pub = json.load(open(os.path.join(base, "overall", "overall_metric.json")))["overall"]
    report("Recount: outcome counts equal published overall_metric.json",
           (A, D, R, S) == (pub["same_partner"], pub["different_partner"], pub["rfam_only"], pub["star3d_only"]),
           f"same {A}, different {D}, Rfam only {R}, STAR3D only {S}; difference rate {(D + R + S)}/{(A + D + R + S)}")
    ps = tsv(os.path.join(base, "pair_scores.tsv"))
    sh, ro, so = (sum(int(p[k]) for p in ps) for k in ("shared_edges", "rfam_only_edges", "star3d_only_edges"))
    exp_path = os.path.join(ROOT, "docs", "results.json")
    if os.path.exists(exp_path):
        exp = json.load(open(exp_path))["headline"]
        got = dict(pairs_aligned=len(ps), same_partner=A, different_partner=D, rfam_only=R, star3d_only=S,
                   shared_edges=sh, rfam_only_edges=ro, star3d_only_edges=so)
        diff = {k: (got[k], exp[k]) for k in got if k in exp and got[k] != exp[k]}
        report("Recount: equals the numbers reported in docs/results.json", not diff, str(diff) if diff else
               f"{len(got)} headline values match")
    else:
        report("Recount: docs/results.json", None, "file not present")
    ran = [ok for _, ok in RESULTS if ok is not None]
    print(f"\n{sum(ran)}/{len(ran)} executed checks passed; {sum(ok is None for _, ok in RESULTS)} skipped")
    sys.exit(0 if all(ran) else 1)


if __name__ == "__main__":
    main()
