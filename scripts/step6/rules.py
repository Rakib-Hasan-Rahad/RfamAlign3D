"""Step 6 rules (pure functions, unit-tested in tests/test_step6.py).

Interaction preservation and region support as defined in docs/METHODS.md (section "Step 6"):
  endpoint states   src_no_coords | unmapped | target_no_coords | mapped   (never merged)
  preservation      any_pair | same_class (oriented, compared in mapped endpoint order) | same_class + identity
  type              canonical_WC (cWW AU/UA/GC/CG) | GU_wobble (cWW GU/UG) | noncanonical (all other LW classes,
                    incl. cWW between other bases) | no_LW_class
  region tier       interaction evidence on the common eligible set from two annotators; distance never decides
"""
WC = {"AU", "UA", "GC", "CG"}
GU = {"GU", "UG"}
DIS = {"different_partner", "rfam_only", "star3d_only"}


def key(u):
    """Residue identity without the residue name: PDB|model|chain|number|insertion code."""
    f = u.split("|")
    return (f[0], f[1], f[2], f[4], f[7] if len(f) > 7 else "")


def flip(c):
    """Leontis-Westhof class seen from the other nucleotide: cWH <-> cHW (cWW unchanged)."""
    return c[0] + c[2] + c[1]


def btype(c, i1, i2):
    if not c:
        return "no_LW_class"
    if c == "cWW":
        return "canonical_WC" if i1 + i2 in WC else ("GU_wobble" if i1 + i2 in GU else "noncanonical")
    return "noncanonical"


def endpoint_state(r1, r2, via, target_res, target_chain):
    """r1, r2: source residue records; via: 'rp' (Rfam) or 'sp' (STAR3D); target_res: {unit: record};
    target_chain: (PDB, chain). Returns (state, t1, t2)."""
    if not (r1.get("xyz") and r2.get("xyz")):
        return "src_no_coords", None, None
    t1, t2 = r1.get(via), r2.get(via)
    if not (t1 and t2):
        return "unmapped", t1, t2
    x1, x2 = target_res.get(t1), target_res.get(t2)
    inst = all((t.split("|")[0], t.split("|")[2]) == tuple(target_chain) and t.split("|")[1] == "1" for t in (t1, t2))
    ok = x1 and x2 and x1.get("xyz") and x2.get("xyz") and inst
    return ("mapped" if ok else "target_no_coords"), t1, t2


def preservation(source_class, target_class, ident_ok):
    """target_class: oriented class of t->w in the target annotation (None = no annotated pair).
    Returns (any_pair, same_class, same_class_and_identity) as 'yes'/'no'."""
    anyp = "yes" if target_class is not None else "no"
    cls = "yes" if (target_class and source_class and target_class == source_class) else "no"
    return anyp, cls, ("yes" if cls == "yes" and ident_ok else "no")


def region_tier(elig_mca, elig_fr3d, gain_s_mca, gain_r_mca, gain_s_fr3d, gain_r_fr3d):
    """Gains count source pairs (common eligible set only) whose class is preserved by one method and not the other."""
    if elig_mca == 0 and elig_fr3d == 0:
        return "unassessable"
    net_m, net_f = gain_s_mca - gain_r_mca, gain_s_fr3d - gain_r_fr3d
    if gain_s_mca == gain_r_mca == gain_s_fr3d == gain_r_fr3d == 0:
        return "no_interaction_difference"
    if (net_m > 0 and net_f < 0) or (net_m < 0 and net_f > 0):
        return "mixed_annotators_disagree"
    if net_m == 0 and net_f == 0:
        return "tied"
    side = "STAR3D" if (net_m > 0 or net_f > 0) else "Rfam"
    agree = (net_m > 0 and net_f > 0) if side == "STAR3D" else (net_m < 0 and net_f < 0)
    opp = (gain_r_mca + gain_r_fr3d) if side == "STAR3D" else (gain_s_mca + gain_s_fr3d)
    if agree and min(abs(net_m), abs(net_f)) >= 2 and opp == 0:
        return f"{side}_strong"
    return f"{side}_moderate" if agree else f"{side}_weak_single_annotator"
