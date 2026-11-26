"""Step 6 regression tests for the rules that replaced s52/s53 scoring (see docs/DECISIONS.md, audit A1-A6)."""
import gzip
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts", "step6"))
import rules as R  # noqa: E402

TGT = {"2BBB|1|B|G|1": dict(n="G", xyz=True), "2BBB|1|B|C|9": dict(n="C", xyz=True),
       "2BBB|1|B|A|5": dict(n="A", xyz=False)}


def test_class_is_compared_in_mapped_order():
    # MC-Annotate lists the target pair w-t; the stored class cWH must be read as cHW from t's side
    assert R.flip("cWH") == "cHW" and R.flip("cWW") == "cWW" and R.flip("tSH") == "tHS"
    assert R.preservation("cHW", R.flip("cWH"), True) == ("yes", "yes", "yes")
    assert R.preservation("cWH", R.flip("cWH"), True)[1] == "no"


def test_unmapped_endpoint_is_never_preserved_or_eligible():
    # s53 bug: an unmapped Rfam endpoint was counted as a STAR3D gain
    u = dict(xyz=True, rp="2BBB|1|B|G|1", sp="2BBB|1|B|G|1")
    v = dict(xyz=True, rp="", sp="2BBB|1|B|C|9")
    assert R.endpoint_state(u, v, "rp", TGT, ("2BBB", "B"))[0] == "unmapped"
    assert R.endpoint_state(u, v, "sp", TGT, ("2BBB", "B"))[0] == "mapped"


def test_missing_coordinates_and_wrong_instance_are_separate_states():
    u = dict(xyz=True, rp="2BBB|1|B|G|1"); v = dict(xyz=True, rp="2BBB|1|B|A|5")
    assert R.endpoint_state(u, v, "rp", TGT, ("2BBB", "B"))[0] == "target_no_coords"
    assert R.endpoint_state(dict(xyz=False, rp="x"), v, "rp", TGT, ("2BBB", "B"))[0] == "src_no_coords"
    w = dict(xyz=True, rp="2BBB|2|B|C|9")              # model 2 is not the analysed instance
    assert R.endpoint_state(u, w, "rp", TGT, ("2BBB", "B"))[0] == "target_no_coords"


def test_canonicality_needs_identities_not_just_cWW():
    assert R.btype("cWW", "G", "C") == "canonical_WC"
    assert R.btype("cWW", "G", "U") == "GU_wobble"
    assert R.btype("cWW", "A", "A") == "noncanonical"
    assert R.btype("tHS", "A", "G") == "noncanonical"
    assert R.btype("", "A", "G") == "no_LW_class"


def test_region_tier_needs_both_annotators_and_ignores_distance():
    assert R.region_tier(0, 0, 0, 0, 0, 0) == "unassessable"
    assert R.region_tier(3, 3, 0, 0, 0, 0) == "no_interaction_difference"
    assert R.region_tier(3, 3, 1, 0, 0, 0) == "STAR3D_weak_single_annotator"   # one contact, one annotator
    assert R.region_tier(3, 3, 1, 0, 1, 0) == "STAR3D_moderate"
    assert R.region_tier(4, 4, 2, 0, 3, 0) == "STAR3D_strong"
    assert R.region_tier(4, 4, 3, 1, 3, 0) == "STAR3D_moderate"               # the other method also gains
    assert R.region_tier(4, 4, 0, 2, 0, 2) == "Rfam_strong"
    assert R.region_tier(4, 4, 1, 0, 0, 1) == "mixed_annotators_disagree"


def test_published_step6_tables_obey_the_rules():
    p = os.path.join(ROOT, "outputs", "06_audit", "step2_r8", "interactions_long.tsv.gz")
    if not os.path.exists(p):
        return
    with gzip.open(p, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            for m in ("rfam", "star3d"):
                if r[f"{m}_state"] != "mapped":
                    assert r[f"{m}_class"] == "" and r[f"{m}_any"] == ""
