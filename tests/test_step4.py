"""Step 4 acceptance tests (synthetic data + the RF00522 production fixture when present)."""
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts", "step4"))
import s41_compare as S  # noqa: E402


def xw(pdb, chain, gapped, observed=None, start_col=1):
    """Crosswalk rows for an aligned row string; observed = set of ungapped positions with coordinates."""
    rows, k = [], 0
    for col, c in enumerate(gapped, start_col):
        if c in "-.":
            continue
        k += 1
        obs = observed is None or k in observed
        rows.append(dict(rfam_column=str(col), ungapped_position=str(k), pdb_id=pdb, auth_chain=chain, comp_id=c,
                         auth_seq_num=str(k), auth_ins_code="", parent=c, observed="yes" if obs else "no",
                         staged_chain=chain, staged_resid=str(k) if obs else ""))
    return rows


PAIR = dict(pair_id="P", family="RFX", query_pdb_chain="1AAA_A", target_pdb_chain="2BBB_B")


def run(a_row, b_row, lines, rf="xxxx.xxx", ss="<<..>>..", a_obs=None, b_obs=None):
    xa, xb = xw("1AAA", "A", a_row, a_obs), xw("2BBB", "B", b_row, b_obs)
    tr = {"RF": rf, "SS_cons": ss} if rf is not None else {"SS_cons": ss}
    return S.compare_pair(PAIR, xa, xb, lines, tr), xa, xb


def test_transposing_labels_preserves_scores_and_swaps_sides():
    a, b = "ACGU.ACG", "AC-UGAC-"
    lines = ["A:1<->B:1", "A:2<->B:3", "A:3<->B:2", "A:5<->B:5"]
    c, _, _ = run(a, b, lines)
    rev = [f"{r.replace('B:', 'A:')}<->{l.replace('A:', 'B:')}" for l, r in (x.split("<->") for x in lines)]
    xa, xb = xw("2BBB", "B", b), xw("1AAA", "A", a)
    t = S.compare_pair(dict(PAIR, query_pdb_chain="2BBB_B", target_pdb_chain="1AAA_A"), xa, xb,
                       [f"B:{x.split('<->')[1].split(':')[1]}<->A:{x.split('<->')[0].split(':')[1]}" for x in lines],
                       {"RF": "xxxx.xxx", "SS_cons": "<<..>>.."})
    for k in ("shared_edges", "rfam_edges", "star3d_edges", "rfam_recovery", "star3d_agreement"):
        assert c["scores"][k] == t["scores"][k]
    sa = {s["side"]: s for s in c["sides"]}
    st = {s["side"]: s for s in t["sides"]}
    for k in ("same_partner", "different_partner", "rfam_only", "star3d_only", "neither"):
        assert sa["A"][k] == st["B"][k] and sa["B"][k] == st["A"][k]


def test_edge_vs_residue_disagreement_not_double_counted():
    c, _, _ = run("ACG", "ACG", ["A:1<->B:1", "A:2<->B:3"], rf="xxx", ss="...")
    assert c["scores"]["rfam_only_edges"] == 2 and c["scores"]["star3d_only_edges"] == 1   # A2-B2, A3-B3 vs A2-B3
    a = next(s for s in c["sides"] if s["side"] == "A")
    assert (a["same_partner"], a["different_partner"], a["rfam_only"]) == (1, 1, 1)        # one residue disagrees


def test_mixed_consensus_insert_edge_and_missing_rf():
    c, _, _ = run("AC.", "A.C", ["A:2<->B:2"], rf="xx.", ss="...")
    e = [x for x in c["edges"] if x["method"] == "star3d"][0]
    assert e["endpoint_rf_group"] == "mixed_consensus_insert"          # A2 in col 2 (consensus), B2 in col 3 (insert)
    c2, _, _ = run("AC.", "A.C", ["A:2<->B:2"], rf=None, ss="...")
    assert [x for x in c2["edges"] if x["method"] == "star3d"][0]["endpoint_rf_group"] == "unknown_involved"
    assert c2["scores"]["rfam_unknown_reproduced"] != "NA" or c2["scores"]["rfam_edges"] == 0


def test_missing_coordinates_keep_rfam_edge_but_no_false_disagreement():
    c, _, _ = run("ACG", "ACG", ["A:1<->B:1", "A:2<->B:3"], rf="xxx", ss="...", b_obs={1, 3})
    raw = [e for e in c["edges"] if e["method"] == "rfam"]
    assert len(raw) == 3 and sum(e["assessable"] == "no" for e in raw) == 1
    a = {r["unit_id"].split("|")[-1]: r for r in c["residues"] if r["side"] == "A"}
    assert a["2"]["category"] == "different_partner" and a["2"]["scored_different"] == "no"   # Rfam partner unobserved
    b = {r["unit_id"].split("|")[-1]: r for r in c["residues"] if r["side"] == "B"}
    assert b["2"]["coords"] == "no" and b["2"]["category"] == "rfam_only"                    # not automatically rfam_only? yes: Rfam gave partner


def test_unobserved_without_rfam_partner_is_neither():
    c, _, _ = run("AC-", "A-G", ["A:1<->B:1"], rf="xxx", ss="...", b_obs={1})
    b = {r["unit_id"].split("|")[-1]: r for r in c["residues"] if r["side"] == "B"}
    assert b["2"]["coords"] == "no" and b["2"]["category"] == "neither"


def test_region_boundaries():
    rows = []
    for pos, d in enumerate(["yes", "yes", "no", "yes", "yes", "yes", "no", "yes"], 1):
        rows.append(dict(pair_id="P", side="A", ungapped_position=str(pos), scored_different=d, unit_id=f"u{pos}",
                         rf_class="consensus", ss_context="paired", rfam_partner="", star3d_partner="",
                         rfam_column=str(pos), raw_line=""))
    rg, mem = S.regions(rows)
    assert [r["length"] for r in rg] == [2, 3, 1]                       # singleton survives; breaks at non-disagreeing


def test_ss_cons_pair_matching_and_pseudoknot():
    info, err = S.parse_ss("<<<..AA..>>>..aa.")
    assert not err and info[1]["partner"] == 12 and info[6]["pk"] and info[6]["partner"] == 16
    assert info[4]["context"] == "hairpin_loop" and info[13]["context"] == "exterior"
    info, err = S.parse_ss("<<..>")
    assert err and info[3]["context"] == "unknown"


def test_rf00522_production_fixture():
    p = os.path.join(ROOT, "outputs", "04_comparison", "step2_interim_440", "side_scores.tsv")
    if not os.path.exists(p):
        return
    sides = {r["side"]: r for r in csv.DictReader(open(p), delimiter="\t") if r["family"] == "RF00522"}
    assert [sides["A"][k] for k in ("same_partner", "different_partner", "rfam_only", "star3d_only", "neither")] == ["16", "15", "1", "1", "0"]
    assert [sides["B"][k] for k in ("same_partner", "different_partner", "rfam_only", "star3d_only", "neither")] == ["16", "14", "2", "2", "0"]
    sc = next(r for r in csv.DictReader(open(p.replace("side_scores", "pair_scores")), delimiter="\t") if r["family"] == "RF00522")
    assert (sc["rfam_edges"], sc["star3d_edges"], sc["shared_edges"]) == ("32", "32", "16")


def test_offset_uses_ungapped_positions_not_author_numbers():
    """3W3S G47B and U47C share author number 47 but are distinct nucleotides (ungapped positions 52 and 53)."""
    import gzip
    p = os.path.join(ROOT, "outputs", "04_comparison", "step2_interim_v2", "diagnostics", "disagreement_diagnostics.tsv.gz")
    if not os.path.exists(p):
        return
    d = [r for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t") if r["unit_id"] == "3RG5|1|A|U|47|||C"]
    assert len(d) == 1 and d[0]["offset"] == "1"
