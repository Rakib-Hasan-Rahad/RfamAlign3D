"""Parser tests on invented Stockholm data (no real Rfam content)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
import step1_inventory as S  # noqa: E402


def parse(text):
    an = []
    recs = list(S.parse_stockholm(S.decode_lines(text.encode("utf-8"), an)))
    return recs, an


WRAPPED = """# STOCKHOLM 1.0
#=GF AC   RFX0001
#=GF ID   toy
#=GF CC   Structure reported in PDB 9ZZZ and 9YYY_B; length 22nt.

rowA/1-8        ACG-
rowB/10-3       ac.U
#=GR rowA/1-8   9ZZZ_aB_SS ((..
#=GR rowA/1-8   9ZZZ_A_SS  ((..
#=GR rowB/10-3  1QQQ       ....
#=GC SS_cons    ((__
#=GC RF         ACGU

rowA/1-8        UNG.
rowB/10-3       GGUU
#=GR rowA/1-8   9ZZZ_aB_SS ..))
#=GR rowA/1-8   9ZZZ_A_SS  ..))
#=GR rowB/10-3  1QQQ       ....
#=GC SS_cons    __))
#=GC RF         NNNN
//
"""


def test_wrapped_fragments_joined_and_width_kept():
    (rec,), an = parse(WRAPPED)
    assert rec["seq"]["rowA/1-8"] == "ACG-UNG."
    assert rec["gc"]["SS_cons"] == "((____))"
    assert rec["gr"][("rowA/1-8", "9ZZZ_aB_SS")] == "((....))"
    assert not rec["problems"]


def test_ungap_hash_and_range():
    assert S.ungap("ac.U-g") == "acUg"
    assert S.seq_hash("acu") == S.seq_hash("ACU")
    assert S.seq_hash("ACT") != S.seq_hash("ACU")          # no T->U
    assert S.split_range("X.1/10-3") == ("X.1", 10, 3, "-")
    assert S.split_range("noRange") == ("noRange", "", "", "")


def test_links_chain_case_preserved_and_leads(tmp_path):
    p = tmp_path / "seed.gz"
    import gzip
    p.write_bytes(gzip.compress(WRAPPED.encode()))
    d = S.build(str(p))
    keys = sorted((x["seed_row_id"], x["pdb_id"], x["chain_id"]) for x in d["links"])
    assert keys == [("rowA/1-8", "9ZZZ", "A"), ("rowA/1-8", "9ZZZ", "aB")]
    fam = d["families"][0]
    assert (fam["linked_rows"], fam["distinct_linked_sequences"], fam["distinct_pdb_chains"]) == (1, 1, 2)
    leads = {(x["evidence_source"], x["candidate_pdb_id"], x["candidate_chain"]) for x in d["leads"]}
    assert ("#=GR", "1QQQ", "") in leads                     # PDB-like tag without chain
    assert ("#=GF CC", "9YYY", "B") in leads
    assert not any(x["candidate_pdb_id"] == "22NT" for x in d["leads"])   # length, not a PDB ID
    row = next(r for r in d["seed_rows"] if r["seed_row_id"] == "rowB/10-3")
    assert row["original_gapped_sequence"] == "ac.UGGUU" and row["ungapped_sequence"] == "acUGGUU"
    assert row["has_lowercase"]


def test_duplicate_in_block_orphan_and_width_flagged(tmp_path):
    bad = ("# STOCKHOLM 1.0\n#=GF AC   RFX0002\nr1/1-3 ACG\nr1/1-3 ACG\nr2/1-2 AC\n"
           "#=GR ghost/1-3 1ABC_A_SS ...\n#=GC SS_cons ..\n//\n")
    import gzip
    p = tmp_path / "s.gz"
    p.write_bytes(gzip.compress(bad.encode()))
    kinds = {a["kind"] for a in S.build(str(p))["anomalies"]}
    assert {"duplicate_key_in_block", "unequal_row_lengths", "orphan_GR"} <= kinds


def test_latin1_fallback_logged():
    an = []
    lines = list(S.decode_lines(b"#=GF RA   M\xfcller\nok", an))
    assert lines[0][1].endswith("Müller") and an[0]["kind"] == "latin1_fallback_decode"


def test_excel_overlong_value_redirected():
    assert S.excel_value("A" * (S.EXCEL_MAX + 1), "k").startswith(S.COMPANION + "|k")
    assert S.excel_value("ACGU", "k") == "ACGU"
