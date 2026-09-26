#!/usr/bin/env python3
"""Run a command and report wall time and peak RSS: timed.py LABEL -- CMD..."""
import resource
import subprocess
import sys
import time

label, cmd = sys.argv[1], sys.argv[sys.argv.index("--") + 1 :]
t = time.perf_counter()
rc = subprocess.run(cmd).returncode
wall = time.perf_counter() - t
rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
print(f"[{label}] wall {wall:.2f}s  peak RSS {rss:.0f} MB  exit {rc}", file=sys.stderr)
sys.exit(rc)
