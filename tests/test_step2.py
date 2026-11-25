"""Step 2 tests on invented data (gate 5 and gate 6 cases)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "step2"))
import s01_discover as D  # noqa: E402
import s03_sequence_evidence as E  # noqa: E402
import s04_coordinates as X  # noqa: E402
import s05_select_prepare as S  # noqa: E402


def rows(**seqs):
    return {k: {"ungapped_sequence": v} for k, v in seqs.items()}


# ---------------------------------------------------------------- ranking happens after eligibility
def test_sequence_in_two_species_not_lost_by_early_choice():
    """Seq B exists in species X (better coverage) and species Y; seq A only in X. Picking B's best structure
    first (X) would block the pair; forming pairs first must select A(X) with B(Y)."""
    reviews = {
        "A_x": dict(decision="verified_natural", technical_status="ready", species_taxid="X"),
        "B_x": dict(decision="verified_natural", technical_status="ready", species_taxid="X"),
        "B_y": dict(decision="verified_natural", technical_status="ready_with_coverage_flags", species_taxid="Y")}
    f = {"A_x": dict(obs=1.0, comp=1.0, res="2.0"), "B_x": dict(obs=1.0, comp=1.0, res="1.5"),
         "B_y": dict(obs=0.9, comp=0.8, res="3.0")}
    facts = lambda c, rv: f[c]  # noqa: E731
    combos = [(S.combo_key(a, b, f[a], f[b])[0], "|".join(sorted((a, b))), a, b) for a, b in (("A_x", "B_x"), ("A_x", "B_y"))]
    best, key, pend, recs = S.walk_pair("RFX", combos, reviews,
                                        lambda c: reviews[c]["decision"] == "verified_natural", facts)
    assert best == ("A_x", "B_y") and not pend
    assert any(r["status"] == "same_species_mmcif" for r in recs)


def test_better_later_alternative_wins_and_pruning_is_sound():
    reviews = {k: dict(decision="verified_natural", technical_status="ready", species_taxid=sp)
               for k, sp in (("a1", "X"), ("b1", "Y"), ("b2", "Z"), ("b3", "W"))}
    f = {"a1": dict(obs=1.0, comp=1.0, res="2.0"), "b1": dict(obs=1.0, comp=0.5, res="2.0"),
         "b2": dict(obs=1.0, comp=1.0, res="2.5"), "b3": dict(obs=0.5, comp=0.5, res="1.0")}
    ub = lambda a, b: (min(f[a]["obs"], f[b]["obs"]),) * 2 + S.combo_key(a, b, f[a], f[b])[0][2:]  # noqa: E731
    combos = [(ub("a1", b), f"a1|{b}", "a1", b) for b in ("b1", "b2", "b3")]
    best, key, pend, recs = S.walk_pair("RFX", combos, reviews, lambda c: True, lambda c, rv: f[c])
    assert best == ("a1", "b2")          # b1 is evaluated first (same upper bound, better resolution) but loses
    st = {r["combo"]: r["status"] for r in recs}
    assert st["a1|b3"] == "pruned_cannot_outrank_selected"


def test_pending_member_blocks_selection():
    reviews = {"a": dict(decision="verified_natural", technical_status="ready", species_taxid="X"),
               "b": dict(decision="pending", technical_status="ready", species_taxid="Y")}
    f = {"a": dict(obs=1, comp=1, res="2"), "b": dict(obs=1, comp=1, res="2")}
    best, key, pend, _ = S.walk_pair("RFX", [((1, 1, 0, -2), "a|b", "a", "b")], reviews,
                                     lambda c: reviews[c]["decision"] == "verified_natural", lambda c, rv: f[c])
    assert pend and best is None


# ---------------------------------------------------------------- seed-row matching
def test_ambiguous_rows_and_tag_resolution():
    fam = rows(**{"r1/1-8": "ACGUACGU", "r2/5-12": "ACGUACGU"})
    ivs, _ = D.match_chain("GGACGUACGUGG", fam, set())
    assert [i["status"] for i in ivs] == ["ambiguous_rows"]
    ivs, _ = D.match_chain("GGACGUACGUGG", fam, {"r2/5-12"})
    assert ivs[0]["status"] == "resolved_by_explicit_tag" and ivs[0]["chosen_row"] == "r2/5-12"


def test_engineered_flank_outside_exact_interval_and_repeat():
    fam = rows(**{"r/1-6": "GCAUGC"})
    ivs, _ = D.match_chain("GGGCAUGCAAAA", fam, set())     # interval exact, extra flanks kept for origin review
    assert ivs[0]["status"] == "exact_unique" and (ivs[0]["start"], ivs[0]["end"]) == (3, 8)
    ivs, _ = D.match_chain("GCAUGCUUGCAUGC", fam, set())
    assert {i["status"] for i in ivs} == {"ambiguous_multiple_intervals"}


def test_unknown_letters_are_never_exact():
    fam = rows(**{"r/1-6": "GCNUGC"})
    ivs, diag = D.match_chain("AAGCAUGCAA", fam, set())
    assert ivs == [] and diag.startswith("wildcard_only")
    fam = rows(**{"r/1-6": "GCAUGC"})
    ivs, diag = D.match_chain("AAGCXUGCAA", fam, set())     # deposited unknown / unresolved modification
    assert ivs == [] and "deposited_unknown=X" in diag


# ---------------------------------------------------------------- origin decisions
def rv(**kw):
    base = dict(mmcif_status="ok", species_taxid="562", species_basis="mmcif_entity_src_nat",
                deposited_parent_sequence="GGGCAUGCAAAA", deposited_parent_sha256="f" * 64,
                doc_documented_engineering=False, doc_mutation_documented=False, doc_pdbx_mutation="",
                doc_engineering_words="")
    base.update(kw)
    return base


def with_evidence(tmp_path, monkeypatch, ev):
    monkeypatch.setattr(S.C, "CACHE", str(tmp_path))
    monkeypatch.setattr(S.C, "CP", lambda *a: os.path.join(str(tmp_path), *a))
    S.C.jdump(ev, os.path.join(str(tmp_path), "origin", f"{'f' * 20}_562.json"))


def test_flank_at_locus_excluded_only_after_broad_search(tmp_path, monkeypatch):
    ev = dict(outcome="interval_found_flank_differs", whole_chain_hits=[], routes_tried=[{"route": "assembly"}],
              interval_hits=[dict(record="NC_1.1", start=10, end=15, strand="+", exact_confirmed="yes")])
    with_evidence(tmp_path, monkeypatch, ev)
    assert S.origin_decision(rv(), "RFX")[1] == "FLANK_DIFFERS_FROM_REFERENCE"
    ev["routes_tried"].append({"route": "blast_core_nt"})
    with_evidence(tmp_path, monkeypatch, ev)
    assert S.origin_decision(rv(), "RFX")[:2] == ("excluded_engineered", "NON_NATIVE_FLANK_AT_LOCUS")
    assert S.origin_decision(rv(deposited_parent_sequence="GCAUGCCCA"), "RF00005")[1] == "TRNA_3PRIME_CCA_POSSIBLY_NOT_ENCODED"


def test_documented_mutation_conflict_and_unknowns(tmp_path, monkeypatch):
    ev = dict(outcome="whole_chain_found", routes_tried=[], interval_hits=[],
              whole_chain_hits=[dict(record="NC_1.1", start=1, end=12, strand="+", independent="yes", exact_confirmed="yes")])
    with_evidence(tmp_path, monkeypatch, ev)
    assert S.origin_decision(rv(), "RFX")[0] == "verified_natural"
    assert S.origin_decision(rv(doc_documented_engineering=True, doc_engineering_words="mutant"), "RFX")[1] == "CONFLICTING_EVIDENCE"
    assert S.origin_decision(rv(deposited_parent_sequence="GGXCA"), "RFX")[1] == "UNKNOWN_RESIDUES_IN_CHAIN"
    assert S.origin_decision(rv(species_taxid=""), "RFX")[1] == "SPECIES_UNRESOLVED"


def test_circular_provenance_rejected():
    class T:
        def species(self, t):
            return ("562", "E. coli", "species") if t in ("562", "83333") else (None, None, "unknown")
    tax = T()
    assert E.independent(dict(sourcedb="pdb", taxid="562", slen=500, title="x"), 30, "562", tax, "2020-01-01") == (False, "pdb_derived_record")
    assert E.independent(dict(sourcedb="insd", taxid="562", slen=30, title="Chain A, RNA"), 30, "562", tax, "2020-01-01")[0] is False
    assert E.independent(dict(sourcedb="insd", taxid="562", slen=32, title="riboswitch", createdate="2021/05/01"),
                         30, "562", tax, "2020-01-01") == (False, "short_record_not_older_than_structure")
    assert E.independent(dict(sourcedb="refseq", taxid="83333", slen=4600000, title="chromosome"), 30, "562", tax, "2020-01-01")[0]


# ---------------------------------------------------------------- coordinates
def fake_entry(positions, can="", sources=None):
    return dict(provenance=dict(mmcif_sha256="0", title="t", method="X-RAY DIFFRACTION", resolution=2.0),
                entities={"1": dict(sources=sources or [dict(category="entity_src_nat", taxid="562", name="E. coli")],
                                    hosts=[], pdbx_description="RNA", can=can)},
                chains={"B": dict(entity_id="1", auth_chain="b", positions=positions,
                                  parent_seq="".join(p["parent"] for p in positions))})


def pos(i, mon, par, observed=True, atoms=X.BACKBONE, num=None, icode=""):
    return dict(seq_id=i, mon=mon, parent=par, auth_num=str(num or i), icode=icode, observed=observed,
                atoms=list(atoms) if observed else [], model_num=(num or i) if observed else None,
                model_icode=icode if observed else "", het="H" if mon not in "ACGU" else "A")


class Tax:
    name = {"562": "Escherichia coli"}

    def species(self, t):
        return ("562", "Escherichia coli", "species") if t == "562" else (None, None, "unknown")

    def binomials(self, text):
        return []

    def is_placeholder(self, t):
        return t == "32630"


def test_missing_coordinates_insertion_code_modified_and_lowercase_chain(monkeypatch):
    p = [pos(1, "G", "G"), pos(2, "PSU", "U", icode="A", num=1), pos(3, "C", "C", observed=False), pos(4, "A", "A")]
    monkeypatch.setattr(X, "entry_data", lambda pdb: fake_entry(p))
    cand = dict(candidate_id="c", pdb_id="1ABC", auth_chain="b", seqres_start="1", seqres_end="4", seed_row_id="r")
    r = X.interval_review(cand, rows(r="GUCA"), Tax(), {"32630"})
    assert r["mmcif_interval_confirmed"] and r["auth_chain"] == "b"
    assert r["technical_status"] == "ready_with_coverage_flags" and r["observed"] == 3
    assert r["unobserved_residues"] == "C3" and r["insertion_codes"] == "1A" and r["modified_residues"] == "PSU1A"
    assert r["species_taxid"] == "562"
    stems = {S.stem_for(dict(family="RF1", pdb_id="1ABC", auth_chain=c, seqres_start=1, seqres_end=4)) for c in ("a", "A", "AB")}
    assert len({s.lower() for s in stems}) == 3          # case-insensitive file systems cannot collide


def test_residue_without_backbone_is_unusable(monkeypatch):
    p = [pos(1, "G", "G"), pos(2, "U", "U", atoms=("N1", "C2"))]
    monkeypatch.setattr(X, "entry_data", lambda pdb: fake_entry(p))
    cand = dict(candidate_id="c", pdb_id="1ABC", auth_chain="b", seqres_start="1", seqres_end="2", seed_row_id="r")
    r = X.interval_review(cand, rows(r="GU"), Tax(), set())
    assert r["technical_status"] == "technical_unusable_or_unresolved"


def test_placeholder_species_unresolved():
    ent = dict(sources=[dict(category="pdbx_entity_src_syn", taxid="32630", name="synthetic construct")], pdbx_description="RNA")
    assert X.species_from_mmcif(ent, Tax(), {"32630"}, "title")[0] == ""


def test_star3d_reader_insertion_code_hetatm_and_ter(tmp_path):
    p = tmp_path / "x.pdb"
    lines = [
        "ATOM      1  P     G A   1       1.000   2.000   3.000  1.00  0.00           P",
        "HETATM    2  P   PSU A   1A      1.000   2.000   3.000  1.00  0.00           P",
        "TER       3        PSU A   1A",
        "ATOM      4  P     C A   2       1.000   2.000   3.000  1.00  0.00           P"]
    p.write_text("\n".join(lines) + "\n")
    order, ok = S.star3d_read(str(p))
    assert ok and order == [(("A", 1, ""), "G"), (("A", 1, "A"), "PSU")]   # chain stops at TER, as in PDB.java


def test_placeholder_taxa_with_species_rank():
    import common as C
    tax = C.Taxonomy()
    assert tax.is_placeholder("408172") and tax.is_placeholder("32630")       # marine metagenome; synthetic construct
    assert not tax.is_placeholder("562") and not tax.is_placeholder("119072")  # E. coli; Caldanaerobacter subsp.
