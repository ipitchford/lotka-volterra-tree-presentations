#!/usr/bin/env python3
"""Replay every verifier and compare its output byte-for-byte with data/*.json.

Usage (from any directory):
  python3 verify_all.py          # replay into a temporary directory and compare
  python3 verify_all.py --write  # regenerate data/*.json
Exit status is nonzero on any verifier failure or byte mismatch.
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERIFIERS = [
    ('code/verify.py', 'verification.json'),
    ('code/verify_hypergraphs.py', 'direct_hypergraph_verification.json'),
    ('code/lv_core_count.py', 'lv_core_counts.json'),
    ('code/verify_extensions.py', 'extensions_verification.json'),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write', action='store_true', help='regenerate data/*.json')
    args = ap.parse_args()
    data = ROOT / 'data'
    data.mkdir(exist_ok=True)
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        outdir = data if args.write else Path(tmp)
        for script, name in VERIFIERS:
            out = outdir / name
            proc = subprocess.run([sys.executable, str(ROOT / script), '--out', str(out)],
                                  cwd=str(ROOT), capture_output=True, text=True)
            if proc.returncode != 0:
                ok = False
                print(f'FAIL  {script}: exit {proc.returncode}\n{proc.stderr[-2000:]}')
                continue
            status = json.loads(out.read_text())
            status = status.get('status') if isinstance(status, dict) else None
            if status != 'PASS':
                ok = False
            if args.write:
                print(f'WROTE {name}  status={status}')
                continue
            ref = data / name
            same = ref.exists() and ref.read_bytes() == out.read_bytes()
            ok = ok and same
            print(f'{"MATCH" if same else "DIFF "} {name}  status={status}')
    print('verify_all: ' + ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
