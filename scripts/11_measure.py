"""Analysis step 2 (restricted inputs): native 2 m subset measurements of the units.

Runs analysis/experiments/bright_geometry/v12_measure.py for the three acquisitions: service lanes, centre-distance
bands, gap sides, 0.5 m boundary perturbations and visible bands of the inter-row footprints of the 2,814 phased
units (the original cohort reuses the release measurements), checked unit by unit against the native footprint
counts at 1, 1.5 and 2 m clearance. Output: <output>/analysis/v12/measurements/<tag>/.

usage: python scripts/11_measure.py [tag ...]      (project layout; about 8 min)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.analysis import run  # noqa: E402

if __name__ == '__main__':
    run('v12_measure.py', *sys.argv[1:])
