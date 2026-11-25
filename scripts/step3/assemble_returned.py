"""Assemble a cohort's STAR3D return folder from reused runs of an earlier cohort plus newly run pairs.

  python scripts/step3/assemble_returned.py <new_cohort> <old_cohort> <new_runs_dir>

<new_runs_dir> holds the runner's runs/ and preprocess/ for the newly run pairs (copied back from the run host).
For every job in outputs/03_star3d/<new_cohort>/step3_jobs.tsv:
  reuse  if the reconciliation says reuse AND the old result's identity (direction, both input sha256, STAR3D jar,
         parameters, command) equals the new job's; the Java path is not compared (native vs container runs)
  new    the run from <new_runs_dir>; its identity must equal the job's as well
Writes outputs/03_star3d/<new_cohort>/returned/: runs/, preprocess/, step3_jobs.tsv, run_provenance.tsv, ledger.json,
RETURN.sha256. Fails if any job has no matching result.
"""
import ast
import collections
import csv
import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KEYS = ("direction", "a_sha256", "b_sha256", "jar", "parameters", "cmd")
# Only outputs are kept (same set as every earlier return); the runner's private STAR3D working copy (jar, bundled
# binaries, libraries, staged PDB copies) is never kept or redistributed.
RUN_KEEP = lambda rel: rel in ("result.json", "align.stdout", "align.stderr") or rel.startswith("out/")  # noqa: E731
PRE_KEEP = lambda rel: rel in ("preprocess.stdout", "preprocess.stderr", "preprocess.json") or rel.startswith("STAR3D_struct_info/")  # noqa: E731


def copy_kept(src, dst, keep):
    for d, _, fs in os.walk(src):
        for fn in fs:
            rel = os.path.relpath(os.path.join(d, fn), src)
            if keep(rel):
                os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
                shutil.copy2(os.path.join(d, fn), os.path.join(dst, rel))


def tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def ident_of(rec):
    i = rec.get("identity")
    i = ast.literal_eval(i) if isinstance(i, str) else (i or {})
    return {k: i.get(k) for k in KEYS}


def main():
    new, old, nd = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    base = os.path.join(ROOT, "outputs", "03_star3d", new)
    oret = os.path.join(ROOT, "outputs", "03_star3d", old, "returned")
    out = os.path.join(base, "returned")
    if os.path.exists(out):
        sys.exit(f"{out} exists; remove it first")
    jobs = tsv(os.path.join(base, "step3_jobs.tsv"))
    jar = json.load(open(os.path.join(oret, "tools_local.json")))["star3d_jar_sha256"]
    rec_dir = os.path.join(ROOT, "outputs", "reconciliation", f"{old}__to__{new}")
    reuse = {l.strip() for l in open(os.path.join(rec_dir, "reuse_list.txt")) if l.strip()}
    oprov = {r["pair_id"]: r for r in tsv(os.path.join(oret, "run_provenance.tsv"))}
    os.makedirs(os.path.join(out, "runs")); os.makedirs(os.path.join(out, "preprocess"))
    prov, problems, status = [], [], collections.Counter()
    for j in jobs:
        pid = j["pair_id"]
        want = dict(direction=f"{j['a_star3d_id']}:{j['a_chain']} -> {j['b_star3d_id']}:{j['b_chain']}",
                    a_sha256=j["a_sha256"], b_sha256=j["b_sha256"], jar=jar, parameters=j["parameters"], cmd=j["align_cmd"])
        src = os.path.join(oret, "runs", pid) if pid in reuse else os.path.join(nd, "runs", pid)
        pre = os.path.join(oret, "preprocess") if pid in reuse else os.path.join(nd, "preprocess")
        rp = os.path.join(src, "result.json")
        if not os.path.exists(rp):
            problems.append(f"{pid}: no result"); continue
        rec = json.load(open(rp))
        if ident_of(rec) != want:
            diff = [k for k in KEYS if ident_of(rec)[k] != want[k]]
            problems.append(f"{pid}: identity differs in {diff}"); continue
        copy_kept(src, os.path.join(out, "runs", pid), RUN_KEEP)
        for side in ("a", "b"):
            p = (rec.get("preprocess") or {}).get(side)
            if p and not os.path.exists(os.path.join(out, "preprocess", p)):
                copy_kept(os.path.join(pre, p), os.path.join(out, "preprocess", p), PRE_KEEP)
        status[rec["status"]] += 1
        if pid in reuse:
            o = oprov.get(pid, {})
            prov.append(dict(pair_id=pid, source=f"reused from {old}", execution_environment=o.get("execution_environment", ""),
                             status=rec["status"], identity_matches_final_job=True))
        else:
            prov.append(dict(pair_id=pid, source=f"new run for {new}", execution_environment="batch-native",
                             status=rec["status"], identity_matches_final_job=True))
    if problems:
        sys.exit("FAILED:\n" + "\n".join(problems[:20]))
    shutil.copy(os.path.join(base, "step3_jobs.tsv"), out)
    for f in ("tools_local.json",) + tuple(x for x in os.listdir(oret) if x.startswith("tools_") and x != "tools_local.json"):
        shutil.copy(os.path.join(oret, f), out)
    with open(os.path.join(out, "run_provenance.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prov[0]), delimiter="\t"); w.writeheader(); w.writerows(prov)
    src_n = collections.Counter("reused" if p["source"].startswith("reused") else "new" for p in prov)
    json.dump(dict(total=len(prov), by_status=dict(status), by_source=dict(src_n)), open(os.path.join(out, "ledger.json"), "w"), indent=1)
    lines = []
    for d, _, fs in os.walk(out):
        for fn in sorted(fs):
            if fn == "RETURN.sha256":
                continue
            p = os.path.join(d, fn)
            h = hashlib.sha256(open(p, "rb").read()).hexdigest()
            lines.append(f"{h}  ./{os.path.relpath(p, out)}\n")
    open(os.path.join(out, "RETURN.sha256"), "w").writelines(sorted(lines, key=lambda x: x.split("  ", 1)[1]))
    print(json.dumps(dict(total=len(prov), by_status=dict(status), by_source=dict(src_n), files=len(lines)), indent=1))


if __name__ == "__main__":
    main()
