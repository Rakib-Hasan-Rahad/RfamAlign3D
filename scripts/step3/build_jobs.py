"""Step 3 / local: build the production job manifest and transfer bundle from a FROZEN Step 2 snapshot.

  python scripts/step3/build_jobs.py <cohort>

Reads only outputs/02_structure_preparation/snapshots/<cohort>/ (immutable). Writes outputs/03_star3d/<cohort>/:
  step3_jobs.tsv            one row per approved pair (identity, inputs, execution settings; result fields empty)
  bundle/inputs/<id>.pdb    staged inputs renamed to lower-case STAR3D IDs (STAR3D lower-cases IDs; see README)
  bundle/software/STAR3D_v1.2.tar.gz, bundle/step3_runner.py, bundle/step3_jobs.tsv, bundle/MANIFEST.sha256
"""
import csv
import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TARBALL_SHA = "c9da44059d9cfc40828e798ebd4f6194c55415e17aaae7d90ac0f849ec2ed946"
JAVA = "/usr/lib/jvm/jre-1.8.0-openjdk/bin/java"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    cohort = sys.argv[1]
    snap = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    out = os.path.join(ROOT, "outputs", "03_star3d", cohort)
    bdir = os.path.join(out, "bundle")
    os.makedirs(os.path.join(bdir, "inputs"), exist_ok=True)
    os.makedirs(os.path.join(bdir, "software"), exist_ok=True)
    pairs = tsv(os.path.join(snap, "pairs_frozen.tsv")) + tsv(os.path.join(snap, "pairs_selection_may_change.tsv"))
    man = {m["stem"]: m for m in tsv(os.path.join(snap, "inputs_star3d", "file_manifest.tsv"))}
    hashes = {h["stem"]: h for h in tsv(os.path.join(snap, "input_file_hashes.tsv"))}
    info = json.load(open(os.path.join(snap, "snapshot_info.json")))
    manifest_sha = sha(os.path.join(snap, "pairs_to_run.tsv"))
    tar = os.path.join(ROOT, "inputs", "software", "STAR3D_v1.2.tar.gz")
    assert sha(tar) == TARBALL_SHA
    shutil.copy2(tar, os.path.join(bdir, "software", "STAR3D_v1.2.tar.gz"))

    # unique inputs: lower-case STAR3D id must stay unique (case-safe chain markers guarantee it; asserted)
    ids = {}
    for stem, m in man.items():
        sid = stem.lower()
        assert sid not in ids.values(), f"lower-case collision {sid}"
        ids[stem] = sid
        src = os.path.join(snap, "inputs_star3d", f"{stem}.pdb")
        assert sha(src) == m["sha256"] == hashes[stem]["sha256"], stem
        dst = os.path.join(bdir, "inputs", f"{sid}.pdb")
        shutil.copy2(src, dst)
        assert sha(dst) == m["sha256"]

    jobs = []
    for p in sorted(pairs, key=lambda p: p["pair_id"]):
        q, t = man[p["query_stem"]], man[p["target_stem"]]
        assert q["sha256"] == p["query_sha256"] and t["sha256"] == p["target_sha256"], p["pair_id"]
        assert q["candidate_id"] == p["query_candidate"] and t["candidate_id"] == p["target_candidate"], "direction"
        assert p["query_stem"] < p["target_stem"], "frozen direction is lexicographic query<target"
        qi, ti = ids[p["query_stem"]], ids[p["target_stem"]]
        jobs.append(dict(
            cohort=cohort, step2_manifest_sha256=manifest_sha, step2_snapshot_utc=info["snapshot_utc"],
            snapshot_group=p["snapshot_group"], pair_id=p["pair_id"], family=p["family"], family_id=p["family_id"],
            direction="A->B (query=A, target=B; frozen in Step 2)",
            a_candidate=p["query_candidate"], b_candidate=p["target_candidate"],
            a_seed_row=p["query_seed_row"], b_seed_row=p["target_seed_row"],
            a_species=p["query_species"], b_species=p["target_species"],
            a_pdb_chain=p["query_pdb_chain"], b_pdb_chain=p["target_pdb_chain"],
            a_label_interval=p["query_label_interval"], b_label_interval=p["target_label_interval"],
            a_staged_stem=p["query_stem"], b_staged_stem=p["target_stem"],
            a_star3d_id=qi, b_star3d_id=ti, a_chain=p["query_chain"], b_chain=p["target_chain"],
            a_file=f"inputs/{qi}.pdb", b_file=f"inputs/{ti}.pdb", a_sha256=q["sha256"], b_sha256=t["sha256"],
            a_residues=q["residues_written"], b_residues=t["residues_written"],
            a_rcsb=f"https://www.rcsb.org/structure/{p['query_pdb_chain'].split('_')[0]}",
            b_rcsb=f"https://www.rcsb.org/structure/{p['target_pdb_chain'].split('_')[0]}",
            crosswalk_keys=f"{p['query_candidate']};{p['target_candidate']}", rfam_reference_key=p["pair_id"],
            star3d_package_sha256=TARBALL_SHA, java=JAVA,
            preprocess_cmd=f"{JAVA} -cp STAR3D.jar Preprocess <id> <chain>",
            align_cmd=f"{JAVA} -jar STAR3D.jar -o out/{p['pair_id']}__forward.aln -p {qi} {p['query_chain']} {ti} {p['target_chain']}",
            parameters="defaults (-r 4.0 -s 3 -g -5 -e -2 -m 3 -i 0); -t not passed (thread_num=1)",
            status="", outputs="", elapsed_s="", validation=""))
    cols = list(jobs[0])
    for path in (os.path.join(out, "step3_jobs.tsv"), os.path.join(bdir, "step3_jobs.tsv")):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(jobs)
    shutil.copy2(os.path.join(ROOT, "scripts", "step3", "step3_runner.py"), os.path.join(bdir, "step3_runner.py"))
    lines = []
    for dp, _, fs in os.walk(bdir):
        for fn in sorted(fs):
            if fn == "MANIFEST.sha256":
                continue
            p = os.path.join(dp, fn)
            lines.append(f"{sha(p)}  {os.path.relpath(p, bdir)}")
    open(os.path.join(bdir, "MANIFEST.sha256"), "w").write("\n".join(sorted(lines, key=lambda l: l[66:])) + "\n")
    print(f"{len(jobs)} jobs, {len(ids)} inputs -> {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
