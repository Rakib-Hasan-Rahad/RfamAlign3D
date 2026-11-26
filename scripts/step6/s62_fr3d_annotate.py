"""Step 6.2: independent base-pair annotation of every cohort chain with FR3D (pinned commit 288f98c, as in config v1).

  <project>/.venv/bin/python scripts/step6/s62_fr3d_annotate.py [<repo_root> <out_dir> <scratch_dir>]
Annotates the deposited mmCIF (first model, FR3D defaults), restricted to the cohort chains of each entry.
Writes <audit_dir>/fr3d/<PDB>_basepair.txt and fr3d/fr3d_manifest.tsv (sha256 of input and output).
"""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import collections
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import urllib.request
import sys

ROOT, OUT = _ARGS[0] if len(_ARGS) > 0 else _REPO, _ARGS[1] if len(_ARGS) > 1 else _os.path.join(_REPO, "outputs", "06_audit", "step2_r8")
TMP = _ARGS[2] if len(_ARGS) > 2 else _os.path.join(OUT, '_tmp_fr3d')
jobs = list(csv.DictReader(open(f"{ROOT}/outputs/03_star3d/step2_r8/step3_jobs.tsv"), delimiter="\t"))
chains = collections.defaultdict(set)
for j in jobs:
    for s in "ab":
        pdb, ch = j[f"{s}_pdb_chain"].split("_", 1)
        chains[pdb].add(ch)
os.makedirs(f"{OUT}/fr3d", exist_ok=True)
os.makedirs(f"{TMP}/in", exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()  # noqa: E731
mpath = f"{OUT}/fr3d/fr3d_manifest.tsv"
prev = {r["pdb"]: r for r in csv.DictReader(open(mpath), delimiter="\t")} if os.path.exists(mpath) else {}
man = []
for pdb in sorted(chains):
    dst = f"{OUT}/fr3d/{pdb}_basepair.txt"
    src = f"{ROOT}/cache/mmcif/{pdb}.cif.gz"
    if os.path.exists(dst) and pdb in prev and not os.path.exists(src):
        # published annotation reused; the mmCIF it was made from is identified by the recorded sha256
        if prev[pdb]["output_sha256"] != sha(dst):
            sys.exit(f"{pdb}: FR3D output does not match fr3d_manifest.tsv")
        man.append(prev[pdb]); continue
    if not os.path.exists(src):
        os.makedirs(os.path.dirname(src), exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.cif.gz", src)
        want = prev.get(pdb, {}).get("mmcif_gz_sha256")
        if want and sha(src) != want:
            print(f"WARNING {pdb}: downloaded mmCIF differs from the one annotated (sha256); RCSB file was updated")
    if not os.path.exists(dst):
        cif = f"{TMP}/in/{pdb}.cif"
        with gzip.open(src, "rb") as fi, open(cif, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        r = subprocess.run([sys.executable, "-m", "fr3d.classifiers.NA_pairwise_interactions", "-i", f"{TMP}/in",
                            "-o", f"{OUT}/fr3d", "-c", "basepair", "--chain", ",".join(sorted(chains[pdb])),
                            f"{pdb}.cif"], capture_output=True, text=True)
        os.remove(cif)
        if r.returncode or not os.path.exists(dst):
            print("FAILED", pdb, r.stderr[-400:]); continue
    man.append(dict(pdb=pdb, chains=",".join(sorted(chains[pdb])), mmcif_gz_sha256=sha(src), output_sha256=sha(dst),
                    n_lines=sum(1 for _ in open(dst))))
    print(pdb, man[-1]["chains"], man[-1]["n_lines"], flush=True)
with open(f"{OUT}/fr3d/fr3d_manifest.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(man[0]), delimiter="\t"); w.writeheader(); w.writerows(man)
print("entries", len(man), "of", len(chains))
