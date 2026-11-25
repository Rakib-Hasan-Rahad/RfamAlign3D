"""Reconcile a newer frozen Step 2 snapshot against an older one (both read-only).

  python scripts/reconcile.py step2_interim_440 step2_interim_v2

Pair identity = family + unordered pair of seed-row ungapped sequences. For every pair in either cohort:
  change            added / removed / replaced (same sequence pair, different structures) / unchanged
  star3d_reuse      reuse only if unchanged AND frozen direction, both staged-file hashes, STAR3D IDs/chains, package
                    and parameters are identical; otherwise run (new) or archive (superseded)
  reference_changed the pair's Rfam reference rows or crosswalk rows differ -> Step 4 must be recomputed (no STAR3D rerun)
Outputs outputs/reconciliation/<old>__to__<new>/: pair_reconciliation.tsv, summary.json, reuse_list.txt
"""
import collections
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "step4"))
import s41_compare as K  # noqa: E402

ROOT = K.ROOT
SN = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots")


def load(c):
    d = os.path.join(SN, c)
    pairs = {}
    for f in ("pairs_frozen.tsv", "pairs_selection_may_change.tsv"):
        for p in K.tsv(os.path.join(d, f)):
            pairs[p["seq_pair_key"]] = p
    ref = collections.defaultdict(list)
    for r in K.tsv(os.path.join(d, "rfam_pair_reference.tsv.gz")):
        ref[r["pair_id"]].append("\t".join(r.values()))
    xw = collections.defaultdict(list)
    for x in K.tsv(os.path.join(d, "residue_crosswalk.tsv.gz")):
        xw[x["candidate_id"]].append("\t".join(v for k, v in x.items() if k != "staged_file"))
    h = lambda xs: hashlib.sha256("\n".join(xs).encode()).hexdigest()  # noqa: E731
    return pairs, {k: h(v) for k, v in ref.items()}, {k: h(v) for k, v in xw.items()}


def ident(p):
    return (p["query_candidate"], p["target_candidate"], p["query_sha256"], p["target_sha256"], p["query_stem"].lower(),
            p["query_chain"], p["target_stem"].lower(), p["target_chain"])


def main():
    old, new = sys.argv[1], sys.argv[2]
    po, ro, xo = load(old)
    pn, rn, xn = load(new)
    out = os.path.join(ROOT, "outputs", "reconciliation", f"{old}__to__{new}")
    os.makedirs(out, exist_ok=True)
    rows = []
    for k in sorted(set(po) | set(pn)):
        o, n = po.get(k), pn.get(k)
        if o and n:
            same_members = {o["query_candidate"], o["target_candidate"]} == {n["query_candidate"], n["target_candidate"]}
            change = "unchanged" if same_members else "replaced"
        else:
            change = "added" if n else "removed"
        reuse = bool(o and n and change == "unchanged" and ident(o) == ident(n) and o["pair_id"] == n["pair_id"])
        refch = bool(o and n and (ro.get(o["pair_id"]) != rn.get(n["pair_id"]) or
                                  xo.get(o["query_candidate"]) != xn.get(n["query_candidate"]) or
                                  xo.get(o["target_candidate"]) != xn.get(n["target_candidate"])))
        rows.append(dict(seq_pair_key=k, family=(n or o)["family"], change=change,
                         old_pair_id=o["pair_id"] if o else "", new_pair_id=n["pair_id"] if n else "",
                         star3d_action=("reuse" if reuse else "run_new" if n else "archive_superseded"),
                         reference_or_crosswalk_changed="yes" if refch else ("" if not (o and n) else "no"),
                         new_group=n["snapshot_group"] if n else ""))
        if change == "replaced" and o:
            rows.append(dict(seq_pair_key=k, family=o["family"], change="replaced(old)", old_pair_id=o["pair_id"], new_pair_id="",
                             star3d_action="archive_superseded", reference_or_crosswalk_changed="", new_group=""))
    K.wtsv(os.path.join(out, "pair_reconciliation.tsv"), rows)
    open(os.path.join(out, "reuse_list.txt"), "w").write("\n".join(r["new_pair_id"] for r in rows if r["star3d_action"] == "reuse") + "\n")
    s = dict(old=old, new=new, old_pairs=len(po), new_pairs=len(pn),
             change=dict(collections.Counter(r["change"] for r in rows if r["change"] != "replaced(old)")),
             star3d=dict(collections.Counter(r["star3d_action"] for r in rows)),
             reference_or_crosswalk_changed=sum(r["reference_or_crosswalk_changed"] == "yes" for r in rows),
             added_by_family=dict(collections.Counter(r["family"] for r in rows if r["change"] == "added")))
    json.dump(s, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
