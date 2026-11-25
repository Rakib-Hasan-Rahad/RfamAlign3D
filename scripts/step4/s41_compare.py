"""Step 4 / part 1: compare Rfam seed correspondences with original STAR3D mappings (no motif data is read here).

  python scripts/step4/s41_compare.py <cohort>

Inputs (read-only): the frozen Step 2 snapshot (pairs, residue_crosswalk, rfam_pair_reference), the returned Step 3
results, and Step 1 family tracks (#=GC RF, SS_cons). Outputs in outputs/04_comparison/<cohort>/:
  run_outcomes.tsv         one row per pair: Step 3 status + Step 4 validation (valid alignments only are scored)
  alignment_edges.tsv.gz   every Rfam and STAR3D edge, raw and assessable, endpoint columns and RF classes
  pair_scores.tsv          direction-independent scores (numerator / denominator), RF strata, endpoint breakdown
  side_scores.tsv          A-side and B-side residue category counts (two views of one comparison; never summed)
  residue_comparisons.tsv.gz  one row per nucleotide per side
  regions.tsv, region_members.tsv.gz   side-specific disagreement regions
  freeze_before_motifs.json            sha256 of every table above (motif-independent)

Definitions (frozen; see README):
  Rfam edge   = (A residue, B residue) for every original seed column holding a nucleotide in both rows.
  STAR3D edge = one raw '<->' mapping line resolved through the crosswalk (staged chain + residue ID -> residue).
  coordinate-available residue = observed in the staged file (crosswalk observed=yes); defined from Step 2 only,
  never from either method's output. Assessable edge = both endpoints coordinate-available.
  Stable residue ID = BGSU-style unit ID PDB|1|chain|comp|number[|||ins].
"""
import collections
import csv
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
csv.field_size_limit(1 << 30)


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def wtsv(p, rows, cols=None):
    rows = list(rows)
    cols = cols or (list(rows[0]) if rows else ["empty"])
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def frac(n, d):
    return f"{n}/{d}" if d else "NA"


# ------------------------------------------------------------------ consensus structure
OPEN, CLOSE = "<([{", ">)]}"


def parse_ss(ss):
    """Pair-matching parse of a WUSS SS_cons. Returns per column (1-based) dict(char, partner, pk, context) and an
    error string. Bracket types <>()[]{} are nested structure; letter pairs Aa..Zz are pseudoknots."""
    n = len(ss)
    partner, pk, err = [None] * (n + 1), [False] * (n + 1), ""
    stacks = collections.defaultdict(list)
    for i, c in enumerate(ss, 1):
        if c in OPEN:
            stacks[c].append(i)
        elif c in CLOSE:
            o = OPEN[CLOSE.index(c)]
            if not stacks[o]:
                err = f"unmatched {c} at {i}"
                continue
            j = stacks[o].pop()
            partner[i], partner[j] = j, i
        elif c.isalpha() and c.isupper():
            stacks["pk" + c].append(i)
        elif c.isalpha() and c.islower():
            k = "pk" + c.upper()
            if not stacks[k]:
                err = f"unmatched pseudoknot {c} at {i}"
                continue
            j = stacks[k].pop()
            partner[i], partner[j] = j, i
            pk[i] = pk[j] = True
    if any(v for v in stacks.values()):
        err = err or "unclosed pairs"
    ctx = [""] * (n + 1)
    nested = [p if (p and not pk[i]) else None for i, p in enumerate(partner)]
    for i in range(1, n + 1):
        if nested[i]:
            ctx[i] = "paired"
            continue
        if pk[i]:
            ctx[i] = "pseudoknot_paired"
            continue
        # enclosing nested pair: nearest j<i opening with partner>i
        enc = None
        for j in range(i - 1, 0, -1):
            if nested[j] and nested[j] > i and j < nested[j]:
                enc = j
                break
        if enc is None:
            ctx[i] = "exterior"
            continue
        k, branches, left, right = enc + 1, 0, 0, 0
        end = nested[enc]
        seen_branch = False
        while k < end:
            if nested[k] and nested[k] > k:
                branches += 1
                seen_branch = True
                k = nested[k] + 1
                continue
            if not seen_branch:
                left += 1
            else:
                right += 1
            k += 1
        ctx[i] = ("hairpin_loop" if branches == 0 else "multiloop" if branches > 1 else
                  "bulge" if (left == 0 or right == 0) else "internal_loop")
    out = {i: dict(char=ss[i - 1], partner=partner[i], pk=pk[i], context=ctx[i] if not err else "unknown")
           for i in range(1, n + 1)}
    return out, err


def rf_class(ch):
    if ch is None:
        return "unknown"
    return "insert" if ch in ".-~" else "consensus"


# ------------------------------------------------------------------ helpers
def unit_id(x):
    u = f"{x['pdb_id']}|1|{x['auth_chain']}|{x['comp_id']}|{x['auth_seq_num']}"
    return u + (f"|||{x['auth_ins_code']}" if x["auth_ins_code"] else "")


def endpoint_group(a, b):
    s = {a, b}
    if "unknown" in s:
        return "unknown_involved"
    if s == {"consensus"}:
        return "both_consensus"
    if s == {"insert"}:
        return "both_insert"
    return "mixed_consensus_insert"


def compare_pair(pair, xa, xb, aln_lines, tracks):
    """Core comparison for one pair. xa/xb: crosswalk rows (one per nucleotide) of A and B, in ungapped order.
    aln_lines: raw mapping lines. Returns dict of edges, scores, side rows, problems."""
    fam = pair["family"]
    rf = tracks.get("RF")
    ss_info, ss_err = parse_ss(tracks["SS_cons"]) if tracks.get("SS_cons") else ({}, "absent")
    by_col_a = {int(x["rfam_column"]): x for x in xa}
    by_col_b = {int(x["rfam_column"]): x for x in xb}
    stA = {(x["staged_chain"], x["staged_resid"]): x for x in xa if x["staged_resid"]}
    stB = {(x["staged_chain"], x["staged_resid"]): x for x in xb if x["staged_resid"]}
    problems = []
    cov = lambda x: x["observed"] == "yes" and bool(x["staged_resid"])  # noqa: E731
    colclass = lambda c: rf_class(rf[c - 1]) if rf else "unknown"  # noqa: E731
    R = {}
    for c in sorted(set(by_col_a) & set(by_col_b)):
        a, b = by_col_a[c], by_col_b[c]
        R[(unit_id(a), unit_id(b))] = dict(method="rfam", a=a, b=b, rfam_column=c, raw_line="")
    S = {}
    for line in aln_lines:
        l, r = line.split("<->")
        ca, ra = l.split(":", 1)
        cb, rb = r.split(":", 1)
        a, b = stA.get((ca, ra)), stB.get((cb, rb))
        if a is None or b is None:
            problems.append(f"unresolved mapping {line}")
            continue
        S[(unit_id(a), unit_id(b))] = dict(method="star3d", a=a, b=b, rfam_column="", raw_line=line)
    if len(S) != len(aln_lines) - sum(1 for p in problems if p.startswith("unresolved")):
        problems.append("duplicate STAR3D edges after resolution")
    edges = []
    for kind, D, other in (("rfam", R, S), ("star3d", S, R)):
        for k, e in D.items():
            a, b = e["a"], e["b"]
            ass = cov(a) and cov(b)
            ra_, rb_ = colclass(int(a["rfam_column"])), colclass(int(b["rfam_column"]))
            edges.append(dict(pair_id=pair["pair_id"], method=kind, a_unit=k[0], b_unit=k[1],
                              a_column=a["rfam_column"], b_column=b["rfam_column"], rfam_column=e["rfam_column"],
                              raw_line=e["raw_line"], a_coords="yes" if cov(a) else "no", b_coords="yes" if cov(b) else "no",
                              assessable="yes" if ass else "no",
                              exclusion_reason="" if ass else ";".join(x for x, ok in (("A endpoint unobserved", cov(a)),
                                                                                        ("B endpoint unobserved", cov(b))) if not ok),
                              a_rf=ra_, b_rf=rb_, endpoint_rf_group=endpoint_group(ra_, rb_),
                              in_both="yes" if k in other else "no"))
    # ---------------- direction-independent scores
    Ra = {k for k in R if cov(R[k]["a"]) and cov(R[k]["b"])}
    Sa = {k for k in S if cov(S[k]["a"]) and cov(S[k]["b"])}
    sh = Ra & Sa
    sc = dict(pair_id=pair["pair_id"], family=fam,
              raw_rfam_edges=len(R), raw_star3d_edges=len(S), raw_shared=len(set(R) & set(S)),
              rfam_edges=len(Ra), star3d_edges=len(Sa), shared_edges=len(sh),
              rfam_recovery=frac(len(sh), len(Ra)), star3d_agreement=frac(len(sh), len(Sa)),
              rfam_only_edges=len(Ra - Sa), star3d_only_edges=len(Sa - Ra),
              rf_available="yes" if rf else "no", ss_cons_parse=ss_err or "ok")
    for stratum in ("consensus", "insert", "unknown"):
        es = {k for k in Ra if colclass(R[k]["rfam_column"]) == stratum}
        sc[f"rfam_{stratum}_reproduced"] = frac(len(es & Sa), len(es))
    for g in ("both_consensus", "mixed_consensus_insert", "both_insert", "unknown_involved"):
        es = {k for k in Sa if endpoint_group(colclass(int(S[k]["a"]["rfam_column"])),
                                              colclass(int(S[k]["b"]["rfam_column"]))) == g}
        sc[f"star3d_{g}_in_rfam"] = frac(len(es & Ra), len(es))
    rc = {k for k in Ra if colclass(int(R[k]["a"]["rfam_column"])) == "consensus" == colclass(int(R[k]["b"]["rfam_column"]))}
    scs = {k for k in Sa if colclass(int(S[k]["a"]["rfam_column"])) == "consensus" == colclass(int(S[k]["b"]["rfam_column"]))}
    sc["consensus_only_rfam_recovery"] = frac(len(rc & scs), len(rc))
    sc["consensus_only_star3d_agreement"] = frac(len(rc & scs), len(scs))
    # ---------------- residue tables: A view and B view of the same mapping sets
    side_rows, residue_rows = [], []
    for side, X, Y, idx in (("A", xa, xb, 0), ("B", xb, xa, 1)):
        rmap = {k[idx]: k[1 - idx] for k in R}
        smap = {k[idx]: k[1 - idx] for k in S}
        yby = {unit_id(y): y for y in Y}
        cnt = collections.Counter()
        for x in X:
            u = unit_id(x)
            rp, sp = rmap.get(u), smap.get(u)
            cat = ("same_partner" if rp and sp and rp == sp else "different_partner" if rp and sp else
                   "rfam_only" if rp else "star3d_only" if sp else "neither")
            scored = cat == "different_partner" and cov(x) and cov(yby[rp]) and cov(yby[sp])
            cnt[cat] += 1
            cnt["observed"] += cov(x)
            if rp and sp:
                cnt["both_mapped"] += 1
            if scored:
                cnt["scored_different"] += 1
            col = int(x["rfam_column"])
            si = ss_info.get(col, {})
            k = (u, rp) if side == "A" else (rp, u)
            ks = (u, sp) if side == "A" else (sp, u)
            residue_rows.append(dict(pair_id=pair["pair_id"], side=side, unit_id=u, ungapped_position=x["ungapped_position"],
                                     nucleotide=x["parent"], comp_id=x["comp_id"], rfam_column=col,
                                     rf_class=colclass(col), ss_char=si.get("char", ""), ss_context=si.get("context", "unknown"),
                                     coords="yes" if cov(x) else "no", category=cat,
                                     rfam_partner=rp or "", rfam_partner_column=yby[rp]["rfam_column"] if rp else "",
                                     rfam_partner_coords=("yes" if cov(yby[rp]) else "no") if rp else "",
                                     star3d_partner=sp or "", star3d_partner_column=yby[sp]["rfam_column"] if sp else "",
                                     raw_line=S[ks]["raw_line"] if sp else "",
                                     scored_different="yes" if scored else "no"))
        side_rows.append(dict(pair_id=pair["pair_id"], family=fam, side=side,
                              pdb_chain=pair["query_pdb_chain"] if side == "A" else pair["target_pdb_chain"],
                              nucleotides=len(X), observed=cnt["observed"], unobserved=len(X) - cnt["observed"],
                              same_partner=cnt["same_partner"], different_partner=cnt["different_partner"],
                              rfam_only=cnt["rfam_only"], star3d_only=cnt["star3d_only"], neither=cnt["neither"],
                              both_mapped=cnt["both_mapped"], scored_different=cnt["scored_different"],
                              different_of_both_mapped=frac(cnt["different_partner"], cnt["both_mapped"])))
    return dict(edges=edges, scores=sc, sides=side_rows, residues=residue_rows, problems=problems)


def regions(residue_rows):
    """Maximal runs of consecutive ungapped positions (per pair and side) carrying a scored different_partner."""
    out, mem = [], []
    by = collections.defaultdict(list)
    for r in residue_rows:
        by[(r["pair_id"], r["side"])].append(r)
    for (pid, side), rs in sorted(by.items()):
        rs.sort(key=lambda r: int(r["ungapped_position"]))
        run, n = [], 0

        def flush():
            nonlocal n
            if not run:
                return
            n += 1
            rid = f"{pid}:{side}:region_{n:03d}"
            out.append(dict(region_id=rid, pair_id=pid, side=side, length=len(run),
                            start_position=run[0]["ungapped_position"], end_position=run[-1]["ungapped_position"],
                            start_unit=run[0]["unit_id"], end_unit=run[-1]["unit_id"],
                            rf_composition=";".join(f"{k}:{v}" for k, v in sorted(collections.Counter(x["rf_class"] for x in run).items())),
                            ss_composition=";".join(f"{k}:{v}" for k, v in sorted(collections.Counter(x["ss_context"] for x in run).items()))))
            for x in run:
                mem.append(dict(region_id=rid, pair_id=pid, side=side, unit_id=x["unit_id"],
                                ungapped_position=x["ungapped_position"], rfam_partner=x["rfam_partner"],
                                star3d_partner=x["star3d_partner"], rfam_column=x["rfam_column"], raw_line=x["raw_line"]))
            run.clear()
        prev = None
        for r in rs:
            pos = int(r["ungapped_position"])
            if r["scored_different"] == "yes" and (not run or pos == prev + 1):
                run.append(r)
            else:
                flush()
                if r["scored_different"] == "yes":
                    run.append(r)
            prev = pos
        flush()
    return out, mem


def main():
    cohort = sys.argv[1]
    snap = os.path.join(ROOT, "outputs", "02_structure_preparation", "snapshots", cohort)
    ret = os.path.join(ROOT, "outputs", "03_star3d", cohort, "returned")
    out = os.path.join(ROOT, "outputs", "04_comparison", cohort)
    os.makedirs(out, exist_ok=True)
    pairs = {p["pair_id"]: p for p in tsv(os.path.join(snap, "pairs_to_run.tsv"))}
    jobs = {j["pair_id"]: j for j in tsv(os.path.join(ROOT, "outputs", "03_star3d", cohort, "step3_jobs.tsv"))}
    xw = collections.defaultdict(list)
    for x in tsv(os.path.join(snap, "residue_crosswalk.tsv.gz")):
        xw[x["candidate_id"]].append(x)
    for v in xw.values():
        v.sort(key=lambda x: int(x["ungapped_position"]))
    tracks = {}
    with gzip.open(os.path.join(ROOT, "outputs", "01_inventory", "family_tracks.jsonl.gz"), "rt") as f:
        for line in f:
            t = json.loads(line)
            if t["family"] in {p["family"] for p in pairs.values()}:
                tracks[t["family"]] = t["gc"]
    outcomes, E, SC, SD, RES = [], [], [], [], []
    for pid in sorted(pairs):
        p, j = pairs[pid], jobs[pid]
        res = json.load(open(os.path.join(ret, "runs", pid, "result.json")))
        o = dict(pair_id=pid, family=p["family"], step3_status=res["status"], declared=(res.get("parsed") or {}).get("declared", ""),
                 rmsd=(res.get("parsed") or {}).get("rmsd", ""), superposition=res.get("superposition", ""),
                 aln_sha256="", step4_status="", problems="")
        if res["status"] != "valid_alignment":
            o["step4_status"] = "not_scored_" + res["status"]
            outcomes.append(o)
            continue
        alnp = os.path.join(ret, "runs", pid, "out", f"{pid}__forward.aln")
        o["aln_sha256"] = sha(alnp)
        if o["aln_sha256"] != res["outputs"][f"out/{pid}__forward.aln"]["sha256"]:
            o.update(step4_status="quarantined_hash_mismatch")
            outcomes.append(o)
            continue
        lines = [l.strip() for l in open(alnp) if "<->" in l and not l.startswith("#")]
        c = compare_pair(p, xw[p["query_candidate"]], xw[p["target_candidate"]], lines, tracks[p["family"]])
        if c["problems"]:
            o.update(step4_status="quarantined_unresolved_mapping", problems=";".join(c["problems"])[:300])
            outcomes.append(o)
            continue
        o["step4_status"] = "compared"
        outcomes.append(o)
        E += c["edges"]
        SC.append(c["scores"])
        SD += c["sides"]
        RES += c["residues"]
    RG, RM = regions(RES)
    files = {
        "run_outcomes.tsv": (outcomes, None), "alignment_edges.tsv.gz": (E, None), "pair_scores.tsv": (SC, None),
        "side_scores.tsv": (SD, None), "residue_comparisons.tsv.gz": (RES, None),
        "regions.tsv": (RG, ["region_id", "pair_id", "side", "length", "start_position", "end_position", "start_unit",
                             "end_unit", "rf_composition", "ss_composition"]),
        "region_members.tsv.gz": (RM, None)}
    for n, (rows, cols) in files.items():
        wtsv(os.path.join(out, n), rows, cols)
    fr = {n: sha(os.path.join(out, n)) for n in files}
    json.dump(dict(cohort=cohort, motif_data_read=False, sha256=fr), open(os.path.join(out, "freeze_before_motifs.json"), "w"), indent=1)
    oc = collections.Counter(o["step4_status"] for o in outcomes)
    print(f"pairs {len(outcomes)} {dict(oc)} | edges {len(E)} | regions {len(RG)}")


if __name__ == "__main__":
    main()
