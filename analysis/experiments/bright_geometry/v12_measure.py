"""v1.2: native 2 m subset measurements of the phased east-west units.

v1.2 adds to the inter-row census the 2,814 east-west fixed-array units whose 30 m pixels share
one hand-drawn construction year but not one Landsat-detected year (coverage expansion, EW_HY).
Their lanes, centre-distance bands, gap sides, 0.5 m boundary perturbations and visible bands are
measured here with the v0.12 code (revision_measure.extra, orientation_diagnostics adapters).

The only change: extra() asserts that its selection equals the release clearance layers, which
cover the original cohort only (value 255 elsewhere). For the new units that assertion is removed
and replaced by an exact per-unit comparison with the native footprint counts of the coverage
expansion (experiments/coverage_expansion/x02_native.py) at 1, 1.5 and 2 m clearance.

Inputs (frozen copies, see v12_units.snapshot): units table and native sums of EW_HY.
Output: analysis/v12/measurements/<tag>/{subsets.parquet (original + new units), new/<cell>.json,
        complete.json}
"""
# Mirrored from experiments/bright_geometry/v12_measure.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '1')
import sys
import json
import time
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_diagnostics as OD
import revision_measure as M
import ecology_revision as R

ROOT = PROJECT  # project data root (pv2mveg configuration; read only)
SRC = ROOT / 'data/eco/releases/fixed_ew_v2_orientation/analysis'
V12 = V12_OUT  # output folder of the v1.2 analysis (pv2mveg configuration)
INP = V12 / 'inputs'
TAGS = ['06', '08', '2026_08']
CLEAR = ['1m', '1p5m', '2m']
_DONE = False


def install():
    """v0.12 adapters, then drop the release-layer assertion of extra()."""
    global _DONE
    if _DONE:
        return
    OD.install()
    old = ("        layer=R.read(R.F.ecology_dir('clearance1p5')/tag/(R.LAYERS[clearance]+'.tif'),r,c,h,w,fill=255)==1\n"
           "        assert np.array_equal(sel,layer&inc),(tag,int(row.cell),clearance,int((sel^(layer&inc)).sum()))\n")
    OD.rewrite(M, 'extra', [(old, '')])
    _DONE = True


def init_worker():
    install()
    M.init('extra', '1p5m')


def units():
    u = pd.read_parquet(INP / 'cells_expansion.parquet')
    return u[u.ew_hy].sort_values('cell').reset_index(drop=True)


def run(tag, workers=8):
    install()
    u = units()
    out = V12 / 'measurements' / tag; new = out / 'new'; new.mkdir(parents=True, exist_ok=True)
    jobs = [(tag, SimpleNamespace(cell=int(r.cell), x=float(r.x), y=float(r.y)), str(new / f'{int(r.cell)}.json')) for r in u.itertuples()]
    todo = [j for j in jobs if not Path(j[2]).exists()]
    print('MEASURE', tag, 'units', len(jobs), 'to measure', len(todo), flush=True)
    t0 = time.time()
    if todo:
        with ProcessPoolExecutor(workers, initializer=init_worker) as pool:
            for i, _ in enumerate(pool.map(M.worker, todo, chunksize=2), 1):
                if i % 200 == 0:
                    print('MEASURE', tag, i, len(todo), f'{time.time() - t0:.0f}s', flush=True)
    rows = []
    for *_, path in jobs:
        rows.extend(json.loads(Path(path).read_text(encoding='utf-8')))
    z = pd.DataFrame(rows)
    # exact check against the native footprint counts of the coverage expansion
    nz = np.load(INP / f'native_ew_hy_{tag}.npz')
    bad = {}
    for j, p in enumerate(CLEAR):
        q = z[z.subset == f'ir_{p}'].set_index('cell').n
        ref = nz['sums'][j][2][q.index.values]
        diff = np.abs(q.values - ref)
        bad[p] = int((diff > 0).sum())
        assert bad[p] == 0, (tag, p, bad[p], int(diff.max()))
    parent = pd.read_parquet(SRC / 'measurements' / tag / 'subsets.parquet')
    assert not set(parent.cell) & set(z.cell)
    both = pd.concat([parent.assign(origin='original'), z.assign(origin='phased')], ignore_index=True)
    both.to_parquet(out / 'subsets.parquet', index=False)
    (out / 'complete.json').write_text(json.dumps(dict(tag=tag, complete=True, original_cells=int(parent.cell.nunique()), phased_cells=int(z.cell.nunique()),
        check='ir_1m/ir_1p5m/ir_2m footprint counts of every phased unit equal the coverage-expansion native sums', mismatches=bad,
        seconds=round(time.time() - t0)), indent=2), encoding='utf-8')
    print('MEASURE_DONE', tag, len(jobs), bad, flush=True)


if __name__ == '__main__':
    for tag in (sys.argv[1:] or TAGS):
        run(tag)
