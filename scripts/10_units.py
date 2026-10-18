"""Analysis step 1 (restricted inputs): units, timing and controls of v1.2.

Runs analysis/experiments/bright_geometry/v12_units.py:
  snapshot   copies the frozen inputs of v1.2 (unit table of the coverage expansion, matched controls on the
             construction-year window, native footprint sums and gap profiles of the phased units, per-unit
             estimates of the expansion) into <output>/analysis/v12/inputs, with a SHA-256 manifest;
  validate   checks that the per-unit primary contrasts of the merged units equal the coverage-expansion
             estimates (units_validation.json).

usage: python scripts/10_units.py        (project layout; about 20 s)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.analysis import run  # noqa: E402

if __name__ == '__main__':
    run('v12_units.py', 'snapshot', 'validate')
