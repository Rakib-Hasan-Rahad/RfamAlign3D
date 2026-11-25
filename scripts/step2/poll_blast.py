"""Wait on already-submitted NCBI BLAST jobs and save their XML where s03_sequence_evidence.blast() expects it.

  python scripts/step2/poll_blast.py [--hours 20]

For every cache/blast/<species>_<tag>.rid.json without a matching .xml: poll the RID (keeps its queue position;
NCBI keeps RIDs ~24-36 h), and when READY save complete XML to cache/blast/<species>_<tag>.xml. Rerunning
s03_sequence_evidence.py then parses that XML and applies the unchanged confirmation and independence rules.
Added 2026-11-25 (policy R8, execution only): NCBI queue waits exceeded s03's 90-min poll window.
"""
import argparse, glob, json, os, re, sys, time, urllib.request
import ssl
import xml.etree.ElementTree as ET
try:
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
def get(u):
    with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "research-pipeline/1.0"}), timeout=600, context=CTX) as r:
        return r.read()
ap = argparse.ArgumentParser(); ap.add_argument("--hours", type=float, default=20); a = ap.parse_args()
todo = {}
for f in glob.glob(os.path.join(ROOT, "cache", "blast", "*.rid.json")):
    x = f[:-9] + ".xml"
    if not os.path.exists(x):
        todo[x] = json.load(open(f))["rid"]
print(f"{len(todo)} jobs to wait for: {sorted(todo.values())}", flush=True)
t0 = time.time()
while todo and time.time() - t0 < a.hours * 3600:
    for x, rid in list(todo.items()):
        try:
            s = get(f"{B}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}").decode()
        except Exception as e:
            print(f"{time.strftime('%H:%M:%S')} {rid} poll error {e!r}", flush=True); continue
        st = (re.search(r"Status=(\w+)", s) or [None, "?"])[1]
        if st == "READY":
            data = get(f"{B}?CMD=Get&FORMAT_TYPE=XML&RID={rid}")
            open(x + ".tmp", "wb").write(data); ET.parse(x + ".tmp"); os.replace(x + ".tmp", x)
            print(f"{time.strftime('%H:%M:%S')} {rid} READY -> {os.path.basename(x)}", flush=True); del todo[x]
        elif st in ("FAILED", "UNKNOWN"):
            print(f"{time.strftime('%H:%M:%S')} {rid} {st}: needs resubmission by s03", flush=True)
            os.remove(x[:-4] + ".rid.json"); del todo[x]
        time.sleep(12)
    if todo:
        print(f"{time.strftime('%H:%M:%S')} waiting: {len(todo)} ({(time.time() - t0) / 60:.0f} min)", flush=True)
        time.sleep(300)
print("remaining:", sorted(todo.values()) or "none")
sys.exit(1 if todo else 0)
