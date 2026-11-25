"""Place every pinned source from config.yaml into inputs/, verified by sha256.

Order per file: already present and correct -> keep; local mirror copy with the same sha256 -> copy;
otherwise download to a temp file, verify, then move into place. A mismatch is a hard error and an
existing file is never overwritten. Every outcome is appended to logs/sources_manifest.tsv.
"""
import datetime as dt
import hashlib
import os
import shutil
import ssl
import sys
import tempfile
import urllib.request

import certifi
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ROOT, "logs", "sources_manifest.tsv")
COLS = ["key", "url", "local_path", "obtained_by", "bytes", "sha256", "utc"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def log(row):
    new = not os.path.exists(MANIFEST)
    with open(MANIFEST, "a") as f:
        if new:
            f.write("\t".join(COLS) + "\n")
        f.write("\t".join(str(row[c]) for c in COLS) + "\n")


def download(url, dest_dir):
    ctx = ssl.create_default_context(cafile=certifi.where())
    req = urllib.request.Request(url, headers={"User-Agent": "research-pipeline/1.0"})
    fd, tmp = tempfile.mkstemp(dir=dest_dir, suffix=".part")
    with os.fdopen(fd, "wb") as out, urllib.request.urlopen(req, context=ctx, timeout=120) as r:
        shutil.copyfileobj(r, out)
    return tmp


def main():
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    mirror = os.path.normpath(os.path.join(ROOT, cfg.get("local_mirror", "")))
    bad = 0
    for s in cfg["sources"]:
        dest = os.path.join(ROOT, s["path"])
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        how = None
        if os.path.exists(dest):
            if sha256(dest) != s["sha256"]:
                print(f"FAIL {s['key']}: existing {s['path']} has wrong sha256 (left untouched)")
                bad += 1
                continue
            how = "already_present"
        else:
            m = os.path.join(mirror, s.get("mirror", "")) if s.get("mirror") else None
            if m and os.path.exists(m) and sha256(m) == s["sha256"]:
                shutil.copy2(m, dest)
                how = "copied_from_local_mirror"
            else:
                try:
                    tmp = download(s["url"], os.path.dirname(dest))
                except Exception as e:  # network or server failure
                    print(f"FAIL {s['key']}: download error {e}")
                    bad += 1
                    continue
                if sha256(tmp) != s["sha256"]:
                    os.remove(tmp)
                    print(f"FAIL {s['key']}: downloaded file sha256 differs from pinned value")
                    bad += 1
                    continue
                os.replace(tmp, dest)
                how = "downloaded"
        log(dict(key=s["key"], url=s["url"], local_path=s["path"], obtained_by=how,
                 bytes=os.path.getsize(dest), sha256=s["sha256"],
                 utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
        print(f"ok   {s['key']:<18} {how}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
