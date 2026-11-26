"""Step 6.4: a row-level correction candidate for the Rfam 15.1 RF00001 seed (mycobacterial 5S rRNA, helix IV/loop E),
built from STAR3D's correspondences and tested with evidence that does NOT come from STAR3D.

  python scripts/step6/s64_correction_case.py [<repo_root> <out_dir>]

1. Proposal: for each residue of the two mycobacterial seed rows (X75601.1/91-203 M. tuberculosis,
   X55251.1/3-116 M. smegmatis), the seed column of its STAR3D counterpart is collected over every cohort partner;
   where a strict majority of partners agree on one column that differs from the seed column, the residue is moved
   there. Residues that STAR3D leaves unmapped in most partners are placed in the first free column between their
   neighbours. Two variants: core (residues whose majority is >=2/3 of partners) and extended (strict majority).
2. Integrity: width unchanged, ungapped sequences identical, residue order preserved, only the two rows change.
3. Tests independent of STAR3D:
   a. consensus secondary structure: canonical/G-U pairs of the row at SS_cons base-paired columns, before/after;
   b. 3D column consistency: every FR3D base pair of each mycobacterial chain is placed at its column pair; support =
      number of OTHER 5S seed rows (non-mycobacterial cohort structures, their own unchanged seed columns) whose FR3D
      annotation has the same LW class at the same column pair. Before vs after, per interaction.
Writes in <audit_dir>/correction_case/: RF00001.seed.original.sto, RF00001.seed.proposed_<variant>.sto, *.diff,
patch_check.json, interaction_column_support.tsv, ss_cons_check.json.
"""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import collections
import csv
import difflib
import glob
import gzip
import json
import os
import sys

ROOT, OUT = _ARGS[0] if len(_ARGS) > 0 else _REPO, _ARGS[1] if len(_ARGS) > 1 else _os.path.join(_REPO, "outputs", "06_audit", "step2_r8")
D = os.path.join(OUT, "correction_case")
os.makedirs(D, exist_ok=True)
ROWS = {"X75601.1/91-203": "Mycobacterium tuberculosis", "X55251.1/3-116": "Mycolicibacterium smegmatis"}

# ---------------------------------------------------------------- seed (RF00001) from the frozen Rfam 15.1 input
lines, cur, keep = [], [], False
with gzip.open(f"{ROOT}/inputs/rfam/15.1/Rfam.seed.gz", "rt", encoding="latin-1") as f:
    for line in f:
        cur.append(line)
        if line.startswith("#=GF AC") and line.split()[2] == "RF00001":
            keep = True
        if line.startswith("//"):
            if keep:
                lines = cur
                break
            cur = []
open(f"{D}/RF00001.seed.original.sto", "w", encoding="latin-1").writelines(lines)
seqs, gc = collections.OrderedDict(), collections.defaultdict(str)
for l in lines:
    if l.startswith("#=GC"):
        _, k, v = l.split(None, 2); gc[k] += v.strip()
    elif l.strip() and not l.startswith(("#", "//")):
        n, s = l.split(); seqs[n] = seqs.get(n, "") + s
W = len(gc["SS_cons"])


def colmap(s):          # residue position (1-based) -> column (1-based)
    m, p = {}, 0
    for i, c in enumerate(s):
        if c not in ".-":
            p += 1; m[p] = i + 1
    return m


def sspairs(ss):
    st, out = collections.defaultdict(list), {}
    op, cl = "<([{", ">)]}"
    for i, c in enumerate(ss, 1):
        if c in op:
            st[c].append(i)
        elif c in cl:
            j = st[op[cl.index(c)]].pop(); out[i] = j; out[j] = i
    return out


SSP = sspairs(gc["SS_cons"])


def key(u):
    f = u.split("|")
    return (f[0], f[1], f[2], f[4], f[7] if len(f) > 7 else "")


# ---------------------------------------------------------------- STAR3D votes from every cohort pair of the rows
votes = {r: collections.defaultdict(collections.Counter) for r in ROWS}
unit_pos = {}            # chain -> {unit key: (seed row, position, seed column, nt)}
partners = collections.defaultdict(set)
for f in sorted(glob.glob(f"{ROOT}/outputs/05_results/step2_r8/app/data/pairs/RF00001__*/pair.json")):
    p = json.load(open(f))
    if not p.get("res"):
        continue
    for s, o in (("A", "B"), ("B", "A")):
        ch, row = p[s.lower()]["chain"], p[s.lower()]["seed"]
        unit_pos.setdefault(ch, {key(x["u"]): (row, x["p"], x["col"], x["n"]) for x in p["res"][s]})
        if row in ROWS:
            partners[row].add(p[o.lower()]["chain"])
            oth = {x["u"]: x for x in p["res"][o]}
            for x in p["res"][s]:
                votes[row][x["p"]][oth[x["sp"]]["col"] if x["sp"] else "-"] += 1


def propose(row, frac):
    cm = colmap(seqs[row])
    new = dict(cm)
    moved = {}
    for pos, c in votes[row].items():
        n = sum(c.values()); m, k = c.most_common(1)[0]
        if k / n > 0.5 and k / n >= frac and m != cm[pos]:
            new[pos] = m; moved[pos] = (cm[pos], m, f"{k}/{n}")
    # residues left unmapped by STAR3D: first free column after the previous residue
    ps = sorted(new)
    for i, pos in enumerate(ps):
        if new[pos] == "-":
            lo = new[ps[i - 1]] + 1
            hi = next(new[q] for q in ps[i + 1:] if new[q] != "-")
            new[pos] = lo if lo < hi else cm[pos]
            moved[pos] = (cm[pos], new[pos], "unmapped by STAR3D in most partners; placed in the first free column")
    # monotonicity / collision check; anything violating is reverted to the seed placement
    while True:
        ps = sorted(new)
        bad = [q for a_, q in zip(ps, ps[1:]) if new[q] <= new[a_]]
        if not bad:
            break
        q = bad[0]
        new[q] = cm[q]; moved.pop(q, None)
    moved = {q: v for q, v in moved.items() if new[q] != cm[q]}       # placements that did not change are not moves
    ungapped = seqs[row].replace("-", "").replace(".", "")
    assert len(set(new.values())) == len(new), "column collision"
    s = ["-"] * W
    for pos, c in new.items():
        s[c - 1] = ungapped[pos - 1]
    return "".join(s), moved


def ss_check(s):
    pr = can = 0
    for i, j in SSP.items():
        if i < j and s[i - 1] not in ".-" and s[j - 1] not in ".-":
            pr += 1; can += (s[i - 1] + s[j - 1]).upper() in {"AU", "UA", "GC", "CG", "GU", "UG"}
    return dict(both_columns_occupied=pr, canonical_or_GU=can)


# ---------------------------------------------------------------- FR3D column-pair reference from other rows
fr3d = collections.defaultdict(list)
for f in glob.glob(f"{OUT}/fr3d/*_basepair.txt"):
    for line in open(f):
        x = line.rstrip("\n").split("\t")
        if len(x) >= 3 and not x[1].startswith("n"):
            fr3d[(key(x[0])[0], key(x[0])[2])].append((key(x[0]), key(x[2]), x[1][0] + x[1][1:3].upper()))
ref = collections.defaultdict(set)      # (c1, c2, class) -> set of seed rows
for ch, um in unit_pos.items():
    row = next(iter(um.values()))[0]
    if row in ROWS:
        continue
    for k1, k2, c in fr3d[tuple(ch.split("_", 1))]:
        if k1 in um and k2 in um:
            ref[(um[k1][2], um[k2][2], c)].add(row)
n_ref_rows = len({next(iter(um.values()))[0] for um in unit_pos.values()} - set(ROWS))

variants = {}
diffs = {}
for vname, frac in (("core", 2 / 3), ("extended", 0.5 + 1e-9)):
    newseqs = dict(seqs)
    info = {}
    for row in ROWS:
        s, moved = propose(row, frac)
        newseqs[row] = s
        info[row] = moved
    # integrity
    chk = dict(width_original=W, width_proposed={r: len(newseqs[r]) for r in ROWS},
               ungapped_identical={r: newseqs[r].replace("-", "") == seqs[r].replace("-", "").replace(".", "") for r in ROWS},
               rows_changed=sorted(r for r in seqs if newseqs[r] != seqs[r]),
               rows_total=len(seqs), residues_moved={r: len(info[r]) for r in ROWS},
               moves={r: {str(k): v for k, v in sorted(info[r].items())} for r in ROWS})
    chk["ss_cons"] = {r: dict(original=ss_check(seqs[r]), proposed=ss_check(newseqs[r])) for r in ROWS}
    # per interaction column support
    rows_out = []
    for ch, um in unit_pos.items():
        row = next(iter(um.values()))[0]
        if row not in ROWS:
            continue
        cmn = colmap(newseqs[row])
        for k1, k2, c in fr3d[tuple(ch.split("_", 1))]:
            if k1 in um and k2 in um and (um[k1][1], um[k2][1]) and um[k1][1] < um[k2][1]:
                p1, p2 = um[k1][1], um[k2][1]
                if not (p1 in info[row] or p2 in info[row]):
                    continue
                b = len(ref[(um[k1][2], um[k2][2], c)]); a_ = len(ref[(cmn[p1], cmn[p2], c)])
                rows_out.append(dict(variant=vname, chain=ch, seed_row=row, u1=f"{um[k1][3]}{p1}", u2=f"{um[k2][3]}{p2}",
                                     fr3d_class=c, cols_original=f"{um[k1][2]}-{um[k2][2]}", cols_proposed=f"{cmn[p1]}-{cmn[p2]}",
                                     rows_supporting_original=b, rows_supporting_proposed=a_,
                                     type="canonical_WC" if c == "cWW" and um[k1][3] + um[k2][3] in ("AU", "UA", "GC", "CG")
                                     else "GU_wobble" if c == "cWW" and um[k1][3] + um[k2][3] in ("GU", "UG") else "noncanonical"))
    agg = collections.defaultdict(collections.Counter)
    for r in rows_out:
        t = agg[(r["chain"], r["type"])]
        t["n"] += 1; t["supp_orig"] += r["rows_supporting_original"] > 0; t["supp_prop"] += r["rows_supporting_proposed"] > 0
        t["better"] += r["rows_supporting_proposed"] > r["rows_supporting_original"]
        t["worse"] += r["rows_supporting_proposed"] < r["rows_supporting_original"]
    chk["column_support_summary"] = {f"{k[0]}|{k[1]}": dict(v) for k, v in sorted(agg.items())}
    chk["reference_rows"] = n_ref_rows
    variants[vname] = (newseqs, chk, rows_out)
    # proposed Stockholm: replace only the two sequence lines, keep everything else byte-identical
    out = []
    for l in lines:
        nm = l.split()[0] if l.strip() and not l.startswith(("#", "//")) else None
        if nm in ROWS:
            pad = l[:len(l) - len(l.lstrip()[len(nm):].lstrip()) - 1] if False else None
            seqstart = l.index(seqs[nm][:10], len(nm))
            out.append(l[:seqstart] + newseqs[nm] + "\n")
        else:
            out.append(l)
    open(f"{D}/RF00001.seed.proposed_{vname}.sto", "w", encoding="latin-1").writelines(out)
    diffs[vname] = "".join(difflib.unified_diff(lines, out, "RF00001.seed.original.sto", f"RF00001.seed.proposed_{vname}.sto", n=0))
    open(f"{D}/RF00001.seed.proposed_{vname}.diff", "w").write(diffs[vname])

json.dump({v: x[1] for v, x in variants.items()}, open(f"{D}/patch_check.json", "w"), indent=1)
allrows = [r for v in variants.values() for r in v[2]]
with open(f"{D}/interaction_column_support.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(allrows[0]), delimiter="\t"); w.writeheader(); w.writerows(allrows)
for v, (ns, chk, _) in variants.items():
    print("==", v, "rows changed", chk["rows_changed"], "moved", chk["residues_moved"], "width", chk["width_proposed"],
          "ungapped identical", chk["ungapped_identical"])
    print("  ss_cons", json.dumps(chk["ss_cons"]))
    for k, s in chk["column_support_summary"].items():
        print("  ", k, s)
    a, b = colmap(seqs["X75601.1/91-203"])[70] - 2, colmap(seqs["X75601.1/91-203"])[112]
    print("  SS_cons ", gc["SS_cons"][a:b])
    for r in ROWS:
        print("  orig    ", seqs[r][a:b], r); print("  proposed", ns[r][a:b])
print("partners", {r: len(v) for r, v in partners.items()}, "reference rows", n_ref_rows)
