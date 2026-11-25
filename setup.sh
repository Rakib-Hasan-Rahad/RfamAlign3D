#!/usr/bin/env bash
# One-command environment setup (safe to rerun). Run from anywhere:
#   bash "<path>/ RfamVSTAR3D findings/setup.sh"
set -euo pipefail
cd "$(dirname "$0")"

echo "== 1/4 Python venv (.venv, Python 3.13)"
[ -d .venv ] || python3.13 -m venv .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -r requirements.txt

echo "== 2/4 Docker runtime (colima x86_64 on Apple silicon)"
if ! docker info >/dev/null 2>&1; then
  colima start --arch x86_64 --cpu 4 --memory 6
fi

echo "== 3/4 STAR3D image star3d-runtime:1"
docker image inspect star3d-runtime:1 >/dev/null 2>&1 \
  || docker build --platform linux/amd64 -t star3d-runtime:1 -f env/Dockerfile.star3d env/

echo "== 4/4 Pinned inputs (Rfam 15.1 + STAR3D v1.2, sha256-verified)"
.venv/bin/python scripts/fetch_sources.py

echo "== Environment check"
.venv/bin/python scripts/check_env.py
