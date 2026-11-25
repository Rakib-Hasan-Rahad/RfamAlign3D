"""Resumable Step 2 driver. Safe to rerun after any interruption (network loss, sleep, crash):
every download, BLAST job ID, genome search, mmCIF parse and chain review is cached atomically, and failed network
work is never cached, so a rerun repeats only unfinished work.

  .venv/bin/python scripts/step2/run_all.py            # full run / resume
  .venv/bin/python scripts/step2/run_all.py --from s05  # resume from a stage

Loop: s03 (evidence) and s05 (selection) alternate until no new evidence keys are requested (max 4 rounds).
Progress is written to outputs/02_structure_preparation/checkpoint.json after every stage.
"""
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

PY = sys.executable
S = os.path.dirname(os.path.abspath(__file__))
STAGES = ["s01", "s02", "s03", "s04", "s05", "s06"]
FILES = {"s01": "s01_discover.py", "s02": "s02_metadata_yield.py", "s03": "s03_sequence_evidence.py",
         "s04": "s04_coordinates.py", "s05": "s05_select_prepare.py", "s06": "s06_export_validate.py"}
DONE_MARK = {"s01": "family_discovery.tsv", "s02": "yield_by_family.tsv", "s04": "validation_sample.tsv"}


def run(stage, log):
    t = time.time()
    with open(log, "a") as f:
        f.write(f"\n===== {stage} {C.now()}\n")
        f.flush()
        rc = subprocess.call([PY, "-u", os.path.join(S, FILES[stage])], stdout=f, stderr=subprocess.STDOUT)
    cp = C.jload(C.P("checkpoint.json"), {})
    cp[stage] = dict(rc=rc, finished_utc=C.now(), seconds=round(time.time() - t))
    C.jdump(cp, C.P("checkpoint.json"))
    print(f"{stage}: rc={rc} ({time.time() - t:.0f} s)", flush=True)
    return rc


def run_s03_until_done(log, tries=5):
    """s03 exits 3 when its watchdog sees a stall; every restart resumes from cache."""
    for t in range(tries):
        rc = run("s03", log)
        if rc == 0:
            return rc
        print(f"s03 exited rc={rc}; restart {t + 1}/{tries} from cache", flush=True)
        time.sleep(30)
    return rc


def main():
    start = sys.argv[sys.argv.index("--from") + 1] if "--from" in sys.argv else "s01"
    log = os.path.join(C.ROOT, "logs", "step2", "run_all.log")
    os.makedirs(os.path.dirname(log), exist_ok=True)
    for st in STAGES[STAGES.index(start):STAGES.index("s05")]:
        if st in DONE_MARK and os.path.exists(C.P(DONE_MARK[st])) and start != st:
            print(f"{st}: outputs present, skipped (rerun with --from {st} to force)")
            continue
        run_s03_until_done(log) if st == "s03" else run(st, log)
    for i in range(4):
        n_before = len(C.read_tsv(C.P("origin_queue_extra.tsv"))) if os.path.exists(C.P("origin_queue_extra.tsv")) else -1
        run("s05", log)
        extra = C.read_tsv(C.P("origin_queue_extra.tsv")) if os.path.exists(C.P("origin_queue_extra.tsv")) else []
        def resolved(r):
            k = f"{r['deposited_seq_sha256'][:20]}_{r['species_taxid']}.json"
            pos = C.jload(C.CP("origin_exact_local", k)) or {}
            return os.path.exists(C.CP("origin", k)) or pos.get("exact_local_resolves")
        pending = [r for r in extra if not resolved(r)]
        if not pending:
            break
        print(f"round {i + 1}: {len(pending)} candidates wait for evidence -> s03", flush=True)
        run_s03_until_done(log)
    run("s06", log)
    # report added / removed / replaced pairs against the latest preserved interim snapshot
    with open(log, "a") as f:
        subprocess.call([PY, "-u", os.path.join(S, "s07_snapshot.py"), "compare"], stdout=f, stderr=subprocess.STDOUT)
    print("snapshot comparison -> outputs/02_structure_preparation/snapshot_comparison.{tsv,json}")


if __name__ == "__main__":
    main()
