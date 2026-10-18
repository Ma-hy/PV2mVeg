"""Analysis step 3 (restricted inputs): the main estimates of v1.2.

Runs analysis/experiments/bright_geometry/v12_run.py, all stages in one process as in the archived run:
  estimate   matched cross-scene design at 1.5 (primary), 1 and 2 m clearance: per-unit contrasts of inter-row
             ground and mixed pixel, age groups, station-clustered resampling (cross_scene/)
  buffer     unmatched buffer baseline and its variants (buffer_primary/)
  lanes, spectral, distance, aux (red edge of the 1.2 m imagery), optics (imaging scenarios), annotate, terrain,
  origin     the mechanism tables of Figs. 7-8 and the unit annotations
Output: <output>/analysis/v12/cross_scene/ and buffer_primary/.

usage: python scripts/12_estimate.py [stage ...]   (project layout; about 30 min)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.analysis import run  # noqa: E402

if __name__ == '__main__':
    run('v12_run.py', *sys.argv[1:])
