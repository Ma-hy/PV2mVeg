"""Render the figures of the manuscript, each in its own process (output: <output>/figures/<name>.png and .pdf).

usage: python scripts/40_figures.py [names ...] [--tier repo|restricted|all] [-j N]

Without names, every figure available for the active layout is drawn: the repository layout draws the figures of
tier "repo" (pv2mveg.catalog); the project layout draws all 20. Fig. 2 needs numbers/numbers.json
(scripts/31_numbers.py); in the project layout, Figs. 6 and S8 write two summaries that the numbers read.
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))
from pv2mveg.catalog import FIGURES, ORDER  # noqa: E402
from pv2mveg.config import P  # noqa: E402


def render(name):
    script = REPO / 'figures' / FIGURES[name][0]
    t0 = time.time()
    r = subprocess.run([sys.executable, str(script)], cwd=str(REPO), capture_output=True, text=True, encoding='utf-8', errors='replace',
                       env=dict(os.environ, PYTHONIOENCODING='utf-8'))
    P.checks.mkdir(parents=True, exist_ok=True)
    (P.checks / f'{name}.log').write_text(r.stdout + '\n' + r.stderr, encoding='utf-8')
    ok = r.returncode == 0 and (P.figures / f'{name}.png').exists()
    return name, ok, time.time() - t0, r.stderr.strip().split('\n')[-1] if r.returncode else ''


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('names', nargs='*')
    ap.add_argument('--tier', choices=['repo', 'restricted', 'all'])
    ap.add_argument('-j', type=int, default=3)
    a = ap.parse_args()
    tier = a.tier or ('all' if P.layout == 'project' else 'repo')
    names = a.names or [n for n in ORDER if tier == 'all' or FIGURES[n][1] == tier]
    bad = [n for n in names if n not in FIGURES]
    if bad:
        sys.exit(f'unknown figures: {bad}; known: {list(FIGURES)}')
    if P.layout == 'repo' and any(FIGURES[n][1] == 'restricted' for n in names):
        sys.exit('restricted figures need the project layout (README.md, "Data availability")')
    P.mkdirs()
    first = [n for n in names if n in ('Fig6', 'FigS8')]
    rest = [n for n in names if n not in first]
    failed = []
    with ThreadPoolExecutor(max_workers=max(1, a.j)) as ex:
        for group in (first, rest):
            for name, ok, sec, err in ex.map(render, group):
                print(f'{name:7s} {"ok" if ok else "FAILED: " + err}  ({sec:.0f} s)', flush=True)
                if not ok:
                    failed.append(name)
    print('FIGURES_DONE', len(names) - len(failed), 'of', len(names), 'into', P.figures)
    sys.exit(1 if failed else 0)
