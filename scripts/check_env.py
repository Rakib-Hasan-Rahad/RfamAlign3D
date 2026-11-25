"""Verify the environment before any data step. Exit 0 only if every check passes.

Checks: Python version and pinned packages (incl. FR3D import), pinned inputs and checksums,
Docker reachable, STAR3D image present, and a STAR3D smoke test inside the container (Java runs,
STAR3D.jar starts, bundled Linux binaries MC-Annotate and RemovePseudoknots execute under emulation).
Writes logs/environment.md.
"""
import hashlib
import importlib.metadata as md
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


def sh(cmd, timeout=600):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout + p.stderr).strip()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    check("python 3.13", sys.version_info[:2] == (3, 13), platform.python_version())

    for line in open(os.path.join(ROOT, "requirements.txt")):
        line = line.split("#")[0].strip()
        if not line:
            continue
        name = line.split("@")[0].split("==")[0].strip()
        want = line.split("==")[1] if "==" in line else None
        try:
            have = md.version(name)
            check(f"pkg {name}", want is None or have == want, have)
        except md.PackageNotFoundError:
            check(f"pkg {name}", False, "not installed")
    try:
        import fr3d  # noqa: F401
        from fr3d.classifiers import NA_pairwise_interactions  # noqa: F401
        check("fr3d import", True)
    except Exception as e:
        check("fr3d import", False, repr(e))

    for s in CFG["sources"]:
        p = os.path.join(ROOT, s["path"])
        check(f"input {s['key']}", os.path.exists(p) and sha256(p) == s["sha256"], s["path"])

    rc, out = sh(["docker", "info", "--format", "{{.Architecture}} {{.ServerVersion}}"], 60)
    check("docker daemon", rc == 0, out.splitlines()[-1] if out else "")
    image = CFG["star3d"]["image"]
    rc, out = sh(["docker", "image", "inspect", image, "--format", "{{.Id}} {{.Architecture}}"], 60)
    check(f"image {image}", rc == 0 and "amd64" in out, out)

    tarball = os.path.join(ROOT, "inputs", "software", "STAR3D_v1.2.tar.gz")
    if rc == 0 and os.path.exists(tarball):
        work = tempfile.mkdtemp(prefix="star3d_smoke_", dir=os.path.join(ROOT, "work"))
        try:
            with tarfile.open(tarball) as tf:
                tf.extractall(work, filter="data")
            src = os.path.join(work, "STAR3D_source")
            for t in ("MC-Annotate", "RemovePseudoknots"):
                os.chmod(os.path.join(src, "tools", t), 0o755)
            base = ["docker", "run", "--rm", "--platform", CFG["star3d"]["platform"],
                    "-v", f"{src}:/star3d", "-w", "/star3d", image, "bash", "-lc"]
            rc, out = sh(base + ["java -version 2>&1 | head -1"])
            check("container java", rc == 0 and "1.8" in out, out)
            rc, out = sh(base + ["java -jar STAR3D.jar 2>&1 | head -3"])
            check("STAR3D.jar starts", "exception" not in out.lower() or "usage" in out.lower(),
                  out.replace("\n", " | ")[:120])
            for t in ("MC-Annotate", "RemovePseudoknots"):
                rc, out = sh(base + [f"./tools/{t} 2>&1 | head -2; echo rc=${{PIPESTATUS[0]}}"])
                ok = "exec format error" not in out and "rc=126" not in out and "rc=127" not in out
                check(f"bundled {t} executes", ok, out.replace("\n", " | ")[:120])
        finally:
            shutil.rmtree(work, ignore_errors=True)

    with open(os.path.join(ROOT, "logs", "environment.md"), "w") as f:
        f.write("# Environment check\n\n")
        _, os_ver = sh(["sw_vers", "-productVersion"], 10)
        f.write(f"- Host: macOS {os_ver}, {platform.machine()}; Python {platform.python_version()}\n")
        _, cv = sh(["colima", "version"], 10)
        f.write(f"- colima: {cv.splitlines()[0] if cv else 'n/a'}\n\n| Check | Result | Detail |\n|---|---|---|\n")
        for n, ok, d in results:
            f.write(f"| {n} | {'PASS' if ok else 'FAIL'} | {d.replace('|', '/')} |\n")
    fails = [n for n, ok, _ in results if not ok]
    print(f"\n{len(results) - len(fails)}/{len(results)} checks passed" + (f"; failed: {fails}" if fails else ""))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
