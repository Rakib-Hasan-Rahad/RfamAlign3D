"""Single importer for exact-local-search natural-sequence evidence (run from one process only).

  python scripts/step2/import_exact_local.py

1. Pulls exact-local-search evidence/<check_id>.json (read-only copy) into outputs/02_structure_preparation/exact_local_exchange/evidence/.
2. Validates each record against pending_checks.tsv (check id, species, sha256 of the full deposited sequence).
3. Re-confirms every whole-chain hit independently on the Mac: efetch of that exact record region from NCBI
   E-utilities (not BLAST), exact string match of the full deposited sequence on the stated strand, and esummary
   independence rule (not PDB-derived, same species). Unconfirmable hits are kept but not used.
4. Writes cache/origin_exact_local/<sha20>_<species>.json. NCBI evidence (cache/origin) is never overwritten; when both
   exist, both are kept, the check is counted once, and any outcome difference is flagged (import_ledger.tsv).
Resolution rule (frozen): exact local search resolves a check only by a confirmed whole-chain match. 'not_found' and
'interval_found_flank_differs' from a limited reference set leave the check pending while the species-wide NCBI
search is outstanding (absence from limited references is not evidence of engineering).
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import s03_sequence_evidence as E  # noqa: E402

SOCK = os.path.expanduser("~/.ssh/sqarma-control/batch.sock")
REMOTE = "<batch-root>/step2_verify/evidence/"
XDIR = C.P("exact_local_exchange", "evidence")
PDIR = C.CP("origin_exact_local")


def pull():
    os.makedirs(XDIR, exist_ok=True)
    subprocess.run(["rsync", "-a", "-e", f"ssh -S {SOCK} -o BatchMode=yes -o ConnectTimeout=10", f"x:{REMOTE}", XDIR + "/"],
                   check=True)


def confirm(hit, seq, species, tax):
    """Independent re-confirmation of one exact-local-search whole-chain hit via efetch + esummary."""
    acc = hit["record"]
    try:
        region, _ = E.efetch_region(acc, hit["start"], hit["end"], pad=0)
    except Exception as e:
        return False, f"efetch failed: {e!r}"[:120]
    target = seq if hit["strand"] == "+" else E.rc(seq)
    if region != target:
        return False, "efetched region does not equal deposited sequence on stated strand"
    m = E.esummary([acc]).get(acc, {})
    if (m.get("sourcedb") or "").lower() == "pdb":
        return False, "PDB-derived record"
    if tax.species(m.get("taxid"))[0] != species:
        return False, f"record species {m.get('taxid')} != {species}"
    return True, f"efetch {acc}:{hit['start']}-{hit['end']}{hit['strand']} exact; esummary taxid {m.get('taxid')}, sourcedb {m.get('sourcedb')}"


def main():
    tax = C.Taxonomy()
    pull()
    pend = {r["check_id"]: r for r in C.read_tsv(C.P("exact_local_exchange", "pending_checks.tsv"))}
    os.makedirs(PDIR, exist_ok=True)
    led = []
    for fn in sorted(os.listdir(XDIR)):
        if not fn.endswith(".json"):
            continue
        r = C.jload(os.path.join(XDIR, fn))
        cid = r["check_id"]
        p = pend.get(cid)
        problems = []
        if not p:
            problems.append("check id not in pending_checks.tsv")
        elif r["species_taxid"] != p["species_taxid"] or r["seq_sha256"] != p["seq_sha256"]:
            problems.append("species or sequence hash differs from the pending check")
        if C.hashlib.sha256(r["deposited_sequence"].encode()).hexdigest() != r["seq_sha256"]:
            problems.append("sequence does not hash to seq_sha256")
        confirmed = []
        if r["outcome"] == "whole_chain_found" and not problems:
            for h in r["whole_chain_hits"][:5]:
                ok, why = confirm(h, r["deposited_sequence"], r["species_taxid"], tax)
                h["mac_reconfirmation"] = why
                if ok:
                    confirmed.append(h)
                    break
        r["whole_chain_hits_confirmed"] = confirmed
        r["exact_local_resolves"] = bool(confirmed)
        r["import_problems"] = problems
        r["imported_utc"] = C.now()
        C.jdump(r, os.path.join(PDIR, f"{r['seq_sha256'][:20]}_{r['species_taxid']}.json"))
        ncbi = C.jload(C.CP("origin", f"{r['seq_sha256'][:20]}_{r['species_taxid']}.json"))
        conflict = ""
        if ncbi:
            a = ncbi["outcome"] == "whole_chain_found"
            b = bool(confirmed)
            conflict = "" if a == b else f"NCBI={ncbi['outcome']} vs exact_local={r['outcome']} (both kept; found in any independent record wins)"
        led.append(dict(check_id=cid, species=r["species_name"], exact_local_outcome=r["outcome"],
                        whole_chain_hits=len(r["whole_chain_hits"]), mac_reconfirmed=bool(confirmed),
                        resolves_check=bool(confirmed), ncbi_evidence="present" if ncbi else "pending",
                        conflict=conflict, problems=";".join(problems),
                        evidence=(f"{confirmed[0]['record']}:{confirmed[0]['start']}-{confirmed[0]['end']}{confirmed[0]['strand']} "
                                  f"({confirmed[0]['source']})") if confirmed else ""))
    C.write_tsv(C.P("exact_local_exchange", "import_ledger.tsv"), led)
    n = len(led)
    res = sum(r["resolves_check"] for r in led)
    print(f"imported {n} exact-local-search records; {res} resolve their check (whole chain, re-confirmed on Mac); "
          f"{sum(1 for r in led if r['exact_local_outcome'] == 'whole_chain_found' and not r['resolves_check'])} whole-chain hits not re-confirmable; "
          f"conflicts {sum(1 for r in led if r['conflict'])}; problems {sum(1 for r in led if r['problems'])}")


if __name__ == "__main__":
    main()
