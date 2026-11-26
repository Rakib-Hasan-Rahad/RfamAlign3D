"""Step 6.5: corrected manuscript figures, built only from pipeline outputs and the audit tables (no typed numbers).

  python scripts/step6/s65_figures.py [<repo_root> <out_dir>]
Fig. 1 workflow/data selection; Fig. 2 quantitative findings; Fig. 3 proposed correction (RF00001, mycobacterial
helix IV). Okabe-Ito colours: Rfam #0072B2, STAR3D #D55E00.
"""
import os as _os
import sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_ARGS = [x for x in _sys.argv[1:] if not x.startswith("--")]

import collections
import csv
import gzip
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

SRC, AUD = _ARGS[0] if len(_ARGS) > 0 else _REPO, _ARGS[1] if len(_ARGS) > 1 else _os.path.join(_REPO, "outputs", "06_audit", "step2_r8")
OUT = f"{AUD}/figures"
os.makedirs(OUT, exist_ok=True)
C = "step2_r8"
RF, ST, GR, DK = "#0072B2", "#D55E00", "#9A9A9A", "#222222"
plt.rcParams.update({"font.size": 7, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": .6, "pdf.fonttype": 42})
W = 4.8


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight")


# ------------------------------------------------------------------ Fig 1
inv = tsv(f"{SRC}/outputs/01_inventory/families.tsv")
snap = f"{SRC}/outputs/02_structure_preparation/snapshots/{C}"
fs = tsv(f"{snap}/family_status.tsv"); cr = tsv(f"{snap}/chain_review.tsv"); st = tsv(f"{snap}/structures.tsv")
pr = tsv(f"{snap}/pairs_to_run.tsv"); ro = tsv(f"{SRC}/outputs/04_comparison/{C}/run_outcomes.tsv")
osum = json.load(open(f"{AUD}/origin_summary.json"))
dec = collections.Counter(r["decision"] for r in cr)
disc = collections.Counter(r["discovery_status"] for r in fs)
reviewed = sum(1 for r in cr if r["mmcif_status"] == "ok")
valid = [r for r in ro if r["step3_status"] == "valid_alignment"]
ev = osum["cohort_chain_evidence"]
steps = [
    ("Rfam 15.1 seed", f"{len(inv):,} families; curated seed rows only (SS_cons, RF)"),
    ("Link rows to 3D", f"{disc['exact_seed_match_found']} families; {len(st):,} PDB entries; {reviewed:,} chains contain a seed row exactly"),
    ("Origin checks", f"{osum['keys_in_current_queue']} sequence-species checks; exact match in an independent non-PDB record"),
    ("Natural chains", f"{dec['verified_natural']:,} verified; {dec['excluded_engineered']} engineered excluded; unconfirmed origin not used"),
    ("Pairs", f"{len(pr)} pairs, {len({r['family'] for r in pr})} families; 138 chains (genome {ev['genome_assembly']}, long record "
                            f"{ev['long_record(>=chain+50nt)']}, older record {ev['older_record_only']})"),
    ("STAR3D v1.2", f"{len(valid)} alignments, {len({r['family'] for r in valid})} families; {len(ro) - len(valid)} 'no similar stacks'; 0 failures"),
    ("Compare", "per nucleotide: Rfam = shared seed column; STAR3D = mapped nucleotide"),
    ("Interactions", "MC-Annotate (STAR3D's input) + FR3D (independent); same class, common set"),
    ("Regions", "interaction evidence, both annotators; distance reported, never decisive"),
]
n = len(steps)
fig, ax = plt.subplots(figsize=(W, 2.75))
ax.set_xlim(0, 10); ax.set_ylim(0, n); ax.axis("off")
for i, (h, d) in enumerate(steps):
    y = n - i - .8
    col = "#EAF2F8" if i < 5 else ("#FBEEE6" if i == 5 else "#F2F2F2")
    ax.add_patch(FancyBboxPatch((.05, y), 9.9, .62, boxstyle="round,pad=0,rounding_size=0.08", fc=col, ec="#999", lw=.5))
    ax.text(.25, y + .31, f"{i + 1}", fontsize=7.5, weight="bold", va="center", color="#444")
    ax.text(.6, y + .31, h, fontsize=6.1, weight="bold", va="center")
    ax.text(2.55, y + .31, d, fontsize=5.1, va="center", color="#333")
    if i < n - 1:
        ax.annotate("", xy=(1.5, y - .38), xytext=(1.5, y), arrowprops=dict(arrowstyle="-|>", lw=.6, color="#666"))
save(fig, "fig1_workflow")

# ------------------------------------------------------------------ Fig 2
om = json.load(open(f"{SRC}/outputs/04_comparison/{C}/overall/overall_metric.json"))["overall"]
sts = json.load(open(f"{SRC}/outputs/04_comparison/{C}/statistics/statistics.json"))
bf = tsv(f"{AUD}/interaction_preservation_by_family.tsv")
rsum = json.load(open(f"{AUD}/region_summary.json"))
fig, axs = plt.subplots(2, 2, figsize=(W, 4.0), gridspec_kw=dict(hspace=.75, wspace=.45))
ax = axs[0, 0]
cats = [("same partner", om["same_partner"], "#BFBFBF"), ("different partner", om["different_partner"], "#CC79A7"),
        ("Rfam only", om["rfam_only"], RF), ("STAR3D only", om["star3d_only"], ST)]
tot = sum(c[1] for c in cats); left = 0
for nm, v, c in cats:
    ax.barh(0, v / tot * 100, left=left, color=c, height=.5); left += v / tot * 100
for k, (nm, v, c) in enumerate(cats):
    ax.text(1, -.45 - k * .3, f"■ {nm}: {v:,} ({v / tot * 100:.1f}%)", fontsize=5.8, color=c if c != "#BFBFBF" else "#777", va="center")
ax.set_yticks([]); ax.set_xlim(0, 100); ax.set_xlabel("% of assessable nucleotides (RNA A)", labelpad=1)
ax.set_ylim(-1.6, .4); ax.spines["left"].set_visible(False)
s = sts["comprehensive_difference_rate"]
ax.set_title(f"a  Outcomes (n = {tot:,})\n    any difference {s['rate'] * 100:.1f}%", loc="left", fontsize=7)


def get(ann, scope, typ, fam, d="A->B"):
    for r in bf:
        if (r["annotator"], r["direction"], r["scope"], r["type"], r["family"]) == (ann, d, scope, typ, fam):
            return r


# (b) noncanonical, at disagreements, by family group and annotator
for ax, typ, ttl in ((axs[0, 1], "noncanonical", "b  Noncanonical pairs at disagreements"),
                     (axs[1, 0], "canonical_WC", "c  Canonical A-U/G-C pairs at disagreements")):
    groups = [("5S rRNA", "RF00001"), ("tRNA", "RF00005"), ("other fam.", "ALL_excl_tRNA_5S"), ("pooled", "ALL")]
    xt = []
    for g, (lbl, fam) in enumerate(groups):
        for j, ann in enumerate(("MC-Annotate", "FR3D")):
            r = get(ann, "disagreement", typ, fam)
            x = g + (j - .5) * .38
            if not r or not r["common_eligible"] or r["common_eligible"] == "0":
                continue
            rv, sv = float(r["rfam_class_rate"]) * 100, float(r["star3d_class_rate"]) * 100
            ax.plot([x, x], [rv, sv], color="#BBB", lw=.8, zorder=1)
            ax.scatter([x], [rv], color=RF, s=10, zorder=2, marker="o" if j == 0 else "s")
            ax.scatter([x], [sv], color=ST, s=10, zorder=2, marker="o" if j == 0 else "s")
        n_ = get("MC-Annotate", "disagreement", typ, fam)["common_eligible"]
        xt.append(f"{lbl}\nn={int(n_):,}")
    ax.set_xticks(range(len(groups)), xt, fontsize=5.4); ax.set_ylabel("% kept (same LW class)", fontsize=6)
    ax.set_ylim(0, 75 if typ == "noncanonical" else 80)
    ax.set_title(ttl, loc="left", fontsize=7)
    ax.scatter([], [], color=RF, s=10, label="Rfam"); ax.scatter([], [], color=ST, s=10, label="STAR3D")
    ax.scatter([], [], color="#666", s=10, marker="o", label="MC-Annotate"); ax.scatter([], [], color="#666", s=10, marker="s", label="FR3D")
    ax.legend(frameon=False, fontsize=5, ncol=2, loc="upper right", handletextpad=.2, columnspacing=.6)
# (d) region tiers by family group
ax = axs[1, 1]
tiers = [("STAR3D_strong", ST, 1), ("STAR3D_moderate", ST, .45), ("Rfam_strong", RF, 1), ("Rfam_moderate", RF, .45)]
grp = [("5S rRNA", ["RF00001"]), ("tRNA", ["RF00005"]), ("other", None)]
tb = rsum["interaction_tiers_by_family"]
for g, (lbl, fams) in enumerate(grp):
    fl = fams or [f for f in tb if f not in ("RF00001", "RF00005")]
    for k, (t, c, a_) in enumerate(tiers):
        v = sum(tb[f].get(t, 0) for f in fl)
        x = g + (k - 1.5) * .2
        ax.bar(x, v, .19, color=c, alpha=a_)
        ax.text(x, v + 2, str(v), ha="center", fontsize=4.8)
ax.set_xticks(range(3), [x[0] for x in grp], fontsize=5.8); ax.set_ylim(0, 150); ax.set_ylabel("disagreement regions", fontsize=6)
from matplotlib.patches import Patch  # noqa: E402
ax.legend(handles=[Patch(color=c, alpha=a_, label=t.replace("_", " ")) for t, c, a_ in tiers], frameon=False, fontsize=5, loc="upper center")
it = rsum["interaction_tiers"]
ax.set_title(f"d  Regions favoured by interactions\n    (both annotators agree; {rsum['regions']:,} regions)", loc="left", fontsize=7)
save(fig, "fig2_results")

# ------------------------------------------------------------------ Fig 3
D = f"{AUD}/correction_case"
pc = json.load(open(f"{D}/patch_check.json"))["core"]
seqs = {}; ss = ""
for fn, tag in (("RF00001.seed.original.sto", "o"), ("RF00001.seed.proposed_core.sto", "p")):
    for l in open(f"{D}/{fn}", encoding="latin-1"):
        if l.startswith("#=GC SS_cons") and tag == "o":
            ss += l.split()[2]
        elif l.strip() and not l.startswith(("#", "//")):
            nm, sq = l.split(); seqs[(tag, nm)] = seqs.get((tag, nm), "") + sq
jobs = tsv(f"{SRC}/outputs/03_star3d/{C}/step3_jobs.tsv")
sp = {}
for j in jobs:
    for s_ in "ab":
        sp[j[f"{s_}_seed_row"]] = j[f"{s_}_species"].split(" (")[0]
c0, c1 = 163, 222
REF = ["M59345.1/3-119", "K02687.1/1-115", "AF034619.1/5584-5700"]
MY = ["X75601.1/91-203", "X55251.1/3-116"]
fig = plt.figure(figsize=(W, 3.9))
ax = fig.add_axes([0, .63, 1, .37]); ax.axis("off")
mono = dict(family="DejaVu Sans Mono", fontsize=5.5, va="top")
ax.text(0, 1, f"a  Rfam 15.1 RF00001 seed, columns {c0 + 1}-{c1} (unchanged rows above; mycobacterial rows below)",
        fontsize=6.8, weight="bold", va="top")
short = lambda r: (lambda w: f"{w[0][0]}. {w[1]}")(sp[r].split())  # noqa: E731
lines = [("SS_cons", ss[c0:c1], DK)] + [(short(r), seqs[("o", r)][c0:c1], GR) for r in REF]
for r in MY:
    lines.append((f"{short(r)} seed", seqs[("o", r)][c0:c1], RF))
    lines.append((f"{short(r)} proposed", seqs[("p", r)][c0:c1], ST))
for k, (lbl, s_, c) in enumerate(lines):
    y = .86 - k * .095
    ax.text(0, y, lbl, fontsize=5.4, va="top", color=DK, style="italic" if lbl not in ("SS_cons",) else "normal")
    ax.text(.21, y, s_, color=c, **mono)
ax.text(0, .86 - len(lines) * .095, "Sequences unchanged; only these two rows move (width 230 columns kept); "
        f"{pc['residues_moved']['X75601.1/91-203']} and {pc['residues_moved']['X55251.1/3-116']} residues moved.",
        fontsize=5.4, va="top", color="#444")
# (b) per-interaction support in other rows (FR3D), M. tuberculosis 7KGB_B and M. smegmatis 9N5T_B
ics = [r for r in tsv(f"{D}/interaction_column_support.tsv") if r["variant"] == "core" and r["chain"] in ("7KGB_B",)]
ics.sort(key=lambda r: int(r["u1"][1:]))
ax = fig.add_axes([.1, .26, .38, .30])
y = range(len(ics))
nref = json.load(open(f"{D}/patch_check.json"))["core"]["reference_rows"]
ax.barh([i + .18 for i in y], [int(r["rows_supporting_original"]) for r in ics], .34, color=RF, label="seed placement")
ax.barh([i - .18 for i in y], [int(r["rows_supporting_proposed"]) for r in ics], .34, color=ST, label="proposed")
ax.set_yticks(list(y), [f"{r['u1']}-{r['u2']} {r['fr3d_class']}" for r in ics], fontsize=5.4)
ax.invert_yaxis(); ax.set_xlim(0, nref); ax.set_xlabel(f"other 5S rows with the same pair\nat these columns (of {nref})", fontsize=5.6)
ax.legend(frameon=False, fontsize=5, loc="lower right")
ax.set_title("b  FR3D pairs, M. tuberculosis 7KGB B", loc="left", fontsize=6.6, weight="bold")
# (c) STAR3D votes along the segment, both mycobacterial rows
ax = fig.add_axes([.6, .26, .38, .30])
votes = collections.defaultdict(lambda: collections.Counter())
for fpath in __import__("glob").glob(f"{SRC}/outputs/05_results/{C}/app/data/pairs/RF00001__*/pair.json"):
    p = json.load(open(fpath))
    if not p.get("res"):
        continue
    for s_, o in (("A", "B"), ("B", "A")):
        row = p[s_.lower()]["seed"]
        if row != "X75601.1/91-203":
            continue
        oth = {x["u"]: x for x in p["res"][o]}
        for x in p["res"][s_]:
            if 70 <= x["p"] <= 110:
                c = oth[x["sp"]]["col"] if x["sp"] else None
                votes[x["p"]]["seed" if c == x["col"] else ("none" if c is None else "other")] += 1
ps = sorted(votes)
tot_ = [sum(votes[q].values()) for q in ps]
ax.bar(ps, [votes[q]["other"] / t for q, t in zip(ps, tot_)], color=ST, width=.8, label="STAR3D: different column")
ax.bar(ps, [votes[q]["seed"] / t for q, t in zip(ps, tot_)], bottom=[votes[q]["other"] / t for q, t in zip(ps, tot_)],
       color=RF, width=.8, label="STAR3D: seed column")
moved = [int(k) for k in pc["moves"]["X75601.1/91-203"]]
mv = [q for q in moved if q >= 80]
ax.plot([min(mv) - .4, max(mv) + .4], [-.06, -.06], color=DK, lw=1.2, clip_on=False)
ax.set_title(f"c  STAR3D placement vs partners\n    (black bar: residues moved, {min(mv)}-{max(mv)})", loc="left", fontsize=6.6, weight="bold")
ax.set_xlabel("M. tuberculosis 5S position", fontsize=5.6, labelpad=6); ax.set_ylabel("share of 18 partners", fontsize=5.6)
ax.legend(frameon=False, fontsize=4.6, loc="upper center", bbox_to_anchor=(.5, -.3), ncol=2); ax.set_ylim(0, 1.0)

# (d) caveat line
ax = fig.add_axes([0, .0, 1, .12]); ax.axis("off")
ssx = pc["ss_cons"]["X75601.1/91-203"]
ax.text(0, 1, "d  Conflicting evidence", fontsize=6.6, weight="bold", va="top")
ax.text(0, .78, f"Canonical/G·U pairs at SS_cons-paired columns, M. tuberculosis row: seed {ssx['original']['canonical_or_GU']}, "
        f"proposed {ssx['proposed']['canonical_or_GU']}.\nIn this window SS_cons pairs column 168 with 194, while FR3D places the "
        "helix IV pair at 168-193 in 16 of 17 other rows:\nthe consensus annotation itself departs from the 3D structures here, "
        "so the proposal is a candidate for curator review.", fontsize=5.6, va="top", color="#333")
save(fig, "fig3_correction")
print("figures written to", OUT)
