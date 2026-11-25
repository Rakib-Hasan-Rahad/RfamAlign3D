#!/usr/bin/env python3
"""Step 3 runner (Linux batch host; Python >= 3.9, standard library only). Original STAR3D v1.2, one forward run/pair.

  python3 step3_runner.py setup                       # verify bundle, unpack STAR3D once, record tool hashes
  python3 step3_runner.py run --cpus 100-115 [--only PAIR ...] [--limit N]
  python3 step3_runner.py ledger                      # print progress

Layout under the production root (this directory):
  inputs/<id>.pdb                    staged Step 2 inputs (hash-verified against step3_jobs.tsv)
  software/STAR3D_source/            read-only template (unpacked from the pinned tarball)
  preprocess/<id>__<chain>__<sha12>/ cached Preprocess output (+ preprocess.json); one per input+chain+tool hashes
  runs/<pair_id>/                    private working directory per job; result.json is the job record
  ledger.json                        counts by status (rewritten after every job)
Resumable: a job with result.json whose execution identity matches is never re-run. Nothing is retried
automatically and parameters are never changed per job.
"""
import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import threading
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
JAVA = os.environ.get("STAR3D_JAVA", "/usr/lib/jvm/jre-1.8.0-openjdk/bin/java")
DOCKER = os.environ.get("STAR3D_DOCKER_IMAGE", "")   # local fallback: run java inside this linux/amd64 image
ENV_LABEL = os.environ.get("STAR3D_ENV_LABEL", "batch-native")
TAR_SHA = "c9da44059d9cfc40828e798ebd4f6194c55415e17aaae7d90ac0f849ec2ed946"
TEMPLATE = os.path.join(ROOT, "software", "STAR3D_source")
LOCK = threading.Lock()
PP_LOCKS = {}


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jdump(o, p):
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(o, f, indent=1, sort_keys=True)
    os.replace(tmp, p)


def jload(p):
    try:
        with open(p) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def jobs():
    with open(os.path.join(ROOT, "step3_jobs.tsv"), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def tools():
    return jload(os.path.join(ROOT, "software", "tools.json"))


# ------------------------------------------------------------------ setup
def setup(_):
    bad = []
    for line in open(os.path.join(ROOT, "MANIFEST.sha256")):
        h, p = line.rstrip("\n").split("  ", 1)
        if sha(os.path.join(ROOT, p)) != h:
            bad.append(p)
    if bad:
        sys.exit(f"MANIFEST mismatch: {bad}")
    tar = os.path.join(ROOT, "software", "STAR3D_v1.2.tar.gz")
    if sha(tar) != TAR_SHA:
        sys.exit("tarball sha256 mismatch")
    if not os.path.isdir(TEMPLATE):
        with tarfile.open(tar) as tf:
            tf.extractall(os.path.join(ROOT, "software"))
    for t in ("MC-Annotate", "RemovePseudoknots"):
        os.chmod(os.path.join(TEMPLATE, "tools", t), 0o755)
    jv = subprocess.run((["docker", "run", "--rm", "--platform", "linux/amd64", DOCKER] if DOCKER else []) + [JAVA, "-version"],
                        capture_output=True, text=True).stderr.strip().splitlines()
    info = dict(tarball_sha256=TAR_SHA, star3d_jar_sha256=sha(os.path.join(TEMPLATE, "STAR3D.jar")),
                mc_annotate_sha256=sha(os.path.join(TEMPLATE, "tools", "MC-Annotate")),
                remove_pseudoknots_sha256=sha(os.path.join(TEMPLATE, "tools", "RemovePseudoknots")),
                java=JAVA, java_version=jv, host=os.uname().nodename, execution_environment=ENV_LABEL, docker_image=DOCKER, kernel=os.uname().release,
                python=sys.version.split()[0], setup_utc=now())
    jdump(info, os.path.join(ROOT, "software", "tools.json"))
    n = 0
    for j in jobs():
        for s in ("a", "b"):
            p = os.path.join(ROOT, j[f"{s}_file"])
            if sha(p) != j[f"{s}_sha256"]:
                sys.exit(f"input hash mismatch {p}")
            n += 1
    print(json.dumps(info, indent=1))
    print(f"bundle verified; {n} job inputs hash-checked")


# ------------------------------------------------------------------ working directory from the template
def workdir(path):
    """Fresh STAR3D home: jar/lib/tools hard-linked from the read-only template; private PDB/ struct_info/ out/."""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    shutil.copytree(TEMPLATE, path, dirs_exist_ok=True, copy_function=os.link,
                    ignore=shutil.ignore_patterns("src", "class", "PDB", "STAR3D_struct_info", "out"))
    for d in ("PDB", "STAR3D_struct_info", "out"):
        os.makedirs(os.path.join(path, d), exist_ok=True)


def run_cmd(cmd, cwd, cpus, log):
    t0 = time.time()
    if DOCKER:
        # STAR3D derives its home from the JAR URL; host paths with spaces become %20 and break it, so mount at /w
        full = ["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{cwd}:/w", "-w", "/w", DOCKER] + cmd
    else:
        full = ["taskset", "-c", cpus, "/usr/bin/time", "-v"] + cmd
    with open(log + ".stdout", "w") as o, open(log + ".stderr", "w") as e:
        rc = subprocess.call(full, cwd=cwd, stdout=o, stderr=e)
    err = open(log + ".stderr").read()
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", err)
    return rc, round(time.time() - t0, 2), int(m.group(1)) if m else None


# ------------------------------------------------------------------ preprocessing cache
def preprocess(sid, chain, in_sha, cpus):
    tl = tools()
    key = f"{sid}__{chain}__{in_sha[:12]}"
    d = os.path.join(ROOT, "preprocess", key)
    with LOCK:
        lk = PP_LOCKS.setdefault(key, threading.Lock())
    with lk:
        rec = jload(os.path.join(d, "preprocess.json"))
        ident = dict(input_sha256=in_sha, chain=chain, jar=tl["star3d_jar_sha256"], mca=tl["mc_annotate_sha256"],
                     rpk=tl["remove_pseudoknots_sha256"], java=tl["java"])
        if rec and rec.get("identity") == ident:
            return rec, d
        workdir(d)
        pdb = os.path.join(d, "PDB", f"{sid}.pdb")
        shutil.copy2(os.path.join(ROOT, "inputs", f"{sid}.pdb"), pdb)
        rc, el, rss = run_cmd([JAVA, "-cp", "STAR3D.jar", "Preprocess", sid, chain], d, cpus, os.path.join(d, "preprocess"))
        files = [f"{sid}.mca", f"{sid}_{chain}.ct", f"{sid}_{chain}.npk.ct"]
        si = os.path.join(d, "STAR3D_struct_info")
        produced = {f: sha(os.path.join(si, f)) for f in files if os.path.exists(os.path.join(si, f))
                    and os.path.getsize(os.path.join(si, f)) > 0}
        out = open(os.path.join(d, "preprocess.stdout")).read()
        problems = []
        if rc != 0:
            problems.append(f"exit {rc}")
        if sha(pdb) != in_sha:
            problems.append("input PDB changed (download fallback?)")
        if "Downloading" in out:
            problems.append("Preprocess attempted a download")
        missing = [f for f in files if f not in produced]
        if missing:
            problems.append("missing/empty " + ",".join(missing))
        rec = dict(identity=ident, key=key, status="ok" if not problems else "failed", problems=problems,
                   intermediates=produced, elapsed_s=el, max_rss_kb=rss, utc=now())
        jdump(rec, os.path.join(d, "preprocess.json"))
        return rec, d


# ------------------------------------------------------------------ one job
def parse_aln(path):
    txt = open(path).read()
    decl = re.search(r"^#Aligned nucleotide:\s*(\d+)", txt, re.M)
    rmsd = re.search(r"^#Alignment RMSD:\s*([0-9.]+)", txt, re.M)
    maps = [l.strip() for l in txt.splitlines() if "<->" in l and not l.startswith("#")]
    a = [m.split("<->")[0] for m in maps]
    b = [m.split("<->")[1] for m in maps]
    return dict(declared=int(decl.group(1)) if decl else None, rmsd=float(rmsd.group(1)) if rmsd else None,
                parsed=len(maps), one_to_one=len(set(a)) == len(a) and len(set(b)) == len(b),
                header_ok=txt.startswith("#STAR3d alignment for"))


def run_job(j, cpus):
    pid = j["pair_id"]
    d = os.path.join(ROOT, "runs", pid)
    tl = tools()
    ident = dict(direction=f"{j['a_star3d_id']}:{j['a_chain']} -> {j['b_star3d_id']}:{j['b_chain']}",
                 a_sha256=j["a_sha256"], b_sha256=j["b_sha256"], jar=tl["star3d_jar_sha256"], java=JAVA,
                 parameters=j["parameters"], cmd=j["align_cmd"])
    prev = jload(os.path.join(d, "result.json"))
    if prev and prev.get("identity") == ident and prev.get("status") not in (None, "running"):
        return prev
    pre = {}
    for s in ("a", "b"):
        rec, pdir = preprocess(j[f"{s}_star3d_id"], j[f"{s}_chain"], j[f"{s}_sha256"], cpus)
        pre[s] = (rec, pdir)
    if any(r["status"] != "ok" for r, _ in pre.values()):
        res = dict(identity=ident, pair_id=pid, status="preprocess_failed", utc=now(),
                   preprocess={s: dict(key=r["key"], status=r["status"], problems=r["problems"]) for s, (r, _) in pre.items()})
        os.makedirs(d, exist_ok=True)
        jdump(res, os.path.join(d, "result.json"))
        return res
    workdir(d)
    jdump(dict(identity=ident, pair_id=pid, status="running", utc=now()), os.path.join(d, "result.json"))
    before = {}
    for s, (rec, pdir) in pre.items():
        sid = j[f"{s}_star3d_id"]
        shutil.copy2(os.path.join(ROOT, "inputs", f"{sid}.pdb"), os.path.join(d, "PDB", f"{sid}.pdb"))
        for f, h in rec["intermediates"].items():
            shutil.copy2(os.path.join(pdir, "STAR3D_struct_info", f), os.path.join(d, "STAR3D_struct_info", f))
            before[f] = h
    aln = os.path.join("out", f"{pid}__forward.aln")
    cmd = [JAVA, "-jar", "STAR3D.jar", "-o", aln, "-p", j["a_star3d_id"], j["a_chain"], j["b_star3d_id"], j["b_chain"]]
    rc, el, rss = run_cmd(cmd, d, cpus, os.path.join(d, "align"))
    stdout = open(os.path.join(d, "align.stdout")).read()
    stderr = open(os.path.join(d, "align.stderr")).read()
    mutated = [f for f, h in before.items() if sha(os.path.join(d, "STAR3D_struct_info", f)) != h]
    inputs_changed = [s for s in ("a", "b") if sha(os.path.join(d, "PDB", f"{j[s + '_star3d_id']}.pdb")) != j[f"{s}_sha256"]]
    ap = os.path.join(d, aln)
    outs = {}
    for p in (ap, ap + ".pdb"):
        if os.path.exists(p):
            outs[os.path.relpath(p, d)] = dict(sha256=sha(p), bytes=os.path.getsize(p))
    parsed = parse_aln(ap) if os.path.exists(ap) and os.path.getsize(ap) > 0 else None
    if rc != 0 and not parsed:
        status = "runtime_failure"
    elif "No similar stacks" in stdout + stderr and (not parsed or parsed["parsed"] == 0):
        status = "no_alignment"
    elif parsed and parsed["header_ok"] and parsed["declared"] == parsed["parsed"] and parsed["one_to_one"] \
            and parsed["parsed"] > 0 and not inputs_changed:
        status = "valid_alignment"
    elif parsed:
        status = "invalid_output_quarantined"
    else:
        status = "runtime_failure"
    sup = "available" if (aln + ".pdb") in outs and outs[aln + ".pdb"]["bytes"] > 0 else "unavailable"
    res = dict(identity=ident, pair_id=pid, family=j["family"], status=status, exit_code=rc, elapsed_s=el,
               execution_environment=ENV_LABEL, docker_image=DOCKER,
               max_rss_kb=rss, cpus=cpus, host=os.uname().nodename, utc=now(), outputs=outs, parsed=parsed,
               superposition=sup, intermediates_mutated=mutated, inputs_changed=inputs_changed,
               preprocess={s: r["key"] for s, (r, _) in pre.items()},
               stdout_tail=stdout[-400:], stderr_tail=stderr[-400:] if rc else "")
    jdump(res, os.path.join(d, "result.json"))
    return res


# ------------------------------------------------------------------ batch + ledger
def ledger():
    js = jobs()
    c = {}
    for j in js:
        r = jload(os.path.join(ROOT, "runs", j["pair_id"], "result.json"))
        st = r["status"] if r else "waiting"
        c[st] = c.get(st, 0) + 1
    led = dict(utc=now(), total=len(js), by_status=c)
    jdump(led, os.path.join(ROOT, "ledger.json"))
    return led


def run(a):
    js = jobs()
    if a.only:
        js = [j for j in js if j["pair_id"] in set(a.only) or j["family"] in set(a.only)]
    # long jobs first (sum of staged residues), stable by pair id
    js.sort(key=lambda j: (-(int(j["a_residues"]) + int(j["b_residues"])), j["pair_id"]))
    if a.limit:
        js = js[:a.limit]
    n_cpu = len(expand(a.cpus))
    workers = a.workers or n_cpu
    print(f"{now()} running {len(js)} jobs, {workers} workers inside CPU set {a.cpus}", flush=True)

    def one(j):
        try:
            r = run_job(j, a.cpus)
        except Exception as e:  # infrastructure error: recorded, never silently retried
            d = os.path.join(ROOT, "runs", j["pair_id"])
            os.makedirs(d, exist_ok=True)
            r = dict(pair_id=j["pair_id"], status="runtime_failure", error=repr(e)[:400], utc=now())
            jdump(r, os.path.join(d, "result.json"))
        with LOCK:
            led = ledger()
            print(f"{now()} {j['pair_id']} -> {r['status']} ({r.get('elapsed_s', '')} s) | {led['by_status']}", flush=True)
        return r
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one, js))
    print(json.dumps(ledger(), indent=1))


def expand(spec):
    out = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("setup")
    r = sub.add_parser("run")
    r.add_argument("--cpus", required=True, help="CPU set, e.g. 100-115 (bounds total CPU use)")
    r.add_argument("--workers", type=int, default=0)
    r.add_argument("--only", nargs="*")
    r.add_argument("--limit", type=int, default=0)
    sub.add_parser("ledger")
    a = ap.parse_args()
    {"setup": setup, "run": run, "ledger": lambda _: print(json.dumps(ledger(), indent=1))}[a.cmd](a)
