"""Analysis step 4 (restricted inputs): the supplementary diagnostics of v1.2, one process per stage.

Runs analysis/experiments/bright_geometry/v12_secondary.py stage by stage, in this order: tables (legacy-format
unit tables, read by the other stages), stats, balance, placebo, subsets, common, crossperiod, temporal, forward.
Each stage runs in a fresh process, as in the archived run (see pv2mveg.analysis). Output:
<output>/analysis/v12/cross_scene/secondary/.

usage: python scripts/13_secondary.py [stage ...]  (project layout; about 10 min)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.analysis import run  # noqa: E402

STAGES = ['tables', 'stats', 'balance', 'placebo', 'subsets', 'common', 'crossperiod', 'temporal', 'forward']

if __name__ == '__main__':
    stages = sys.argv[1:] or STAGES
    bad = [s for s in stages if s not in STAGES]
    if bad:
        sys.exit(f'unknown stages: {bad}; known: {STAGES}')
    for s in stages:
        run('v12_secondary.py', s)
