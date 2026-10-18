"""v1.2 main analyses on the merged east-west units (original cohort + phased units).

Runs the v1.1 estimators unchanged, with the v1.2 loader (v12_units.load_v12) and the output root
analysis/v12/ (same layout as the v1.1 analysis root, so the figure and text code can switch roots).

Stages
  estimate   cross_scene.main at 1.5, 1 and 2 m (primary cross-scene design; the 'original'
             comparison columns refer to the same-image control subset)
  buffer     buffer_primary.main: unmatched buffer baseline and its variants on the primary units
  lanes      service-lane contrasts (cross_scene_lanes rules) from the merged subset measurements
  spectral   cross_scene_mechanisms.spectral (visible bands need the BJ3A1 calibration)
  distance   cross_scene_mechanisms.distance
  aux        auxiliary 1.2 m red-edge indices (e2_spec_analysis with the v1.2 units, construction-year
             ages and hand-drawn station clusters) at 1.5, 1 and 2 m
  optics     imaging-scenario re-aggregation (cross_scene_mechanisms.optics) for the units that have
             forward optical models (original cohort)
  annotate   adds origin, construction year and station label to the unit tables
"""
# Mirrored from experiments/bright_geometry/v12_run.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '1')
import sys
import json
import shutil
import time
import warnings
from functools import partial
from pathlib import Path
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_statistics as N
import v12_units as V

E = N.E; cc = E.cc; es = cc.es; J = E.J
SRC = V.SRC; V12 = V.V12; CS12 = V12 / 'cross_scene'
TAGS = V.TAGS


def route(same='image'):
    N.load = partial(V.load_v12, same=same)
    N.OUT = V12
    CS12.mkdir(parents=True, exist_ok=True)


def estimate():
    route('image')
    import cross_scene as C
    for p in ['1p5m', '1m', '2m']:
        C.OUT = CS12
        C.main(600, p)
        print('ESTIMATE_DONE', p, flush=True)


def buffer():
    route('image')
    import buffer_primary as BP
    BP.OUT = V12 / 'buffer_primary'
    BP.main(600)


def lanes():
    route('image')
    K = ['ndvi', 'nirv', 'msavi']
    GROUPS = [(str(a), a, a) for a in range(15)] + [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('all', 0, 99)]

    def wmean(x, w):
        ok = np.isfinite(x) & (w > 0)
        return float(np.sum(x[ok] * w[ok]) / w[ok].sum()) if ok.any() else np.nan
    rows = []
    for tag in TAGS:
        d = N.load(tag, '1p5m'); tidx, age, ct = d['tidx'], d['age'], d['ct']
        cs = pd.read_csv(CS12 / f'cells_{tag}.csv'); assert np.array_equal(cs.cell.values, tidx)
        sub = pd.read_parquet(V12 / 'measurements' / tag / 'subsets.parquet')
        z = sub[sub.subset == 'lane_1p5m'].drop_duplicates('cell').set_index('cell').reindex(tidx)
        sr, sn, n = z.s_red.values, z.s_nir.values, np.nan_to_num(z.n.values)
        I = cc.index_from_sums(np.nan_to_num(sr), np.nan_to_num(sn), n)
        st = ct.station.to_numpy(); u, inv = np.unique(st, return_inverse=True)
        valid = cs.valid_cross.values & (n >= 50)
        for k in K:
            r, _ = E.A1.per_cell(d, 'all', k)
            x = I[k] - r['ic'] - r['dpre']
            rng = np.random.default_rng(20261010); boots = {g: [] for g, *_ in GROUPS}
            for b in range(2000):
                wt = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u)).astype(float)[inv]
                for g, lo, hi in GROUPS:
                    boots[g].append(wmean(x, wt * valid * ((age >= lo) & (age <= hi))))
            for g, lo, hi in GROUPS:
                m = valid & (age >= lo) & (age <= hi) & np.isfinite(x); nc = len(np.unique(st[m]))
                bb = np.asarray(boots[g]); bb = bb[np.isfinite(bb)]
                l, h = np.percentile(bb, [2.5, 97.5]) if nc >= J.MIN_CL and len(bb) > 50 else (np.nan, np.nan)
                rows.append(dict(tag=tag, index=k, group=g, lane=wmean(x, m.astype(float)), lane_lo=l, lane_hi=h, lane_n=int(m.sum()), lane_clusters=nc))
        print('LANES', tag, int(valid.sum()), flush=True)
    pd.DataFrame(rows).to_csv(CS12 / 'lanes.csv', index=False)


def mech():
    route('image')
    import cross_scene_mechanisms as CM
    CM.OUT = CS12
    return CM


def spectral():
    mech().spectral()


def distance():
    mech().distance()


def optics():
    """cross_scene_mechanisms.optics with the frozen forward-model inputs of the release."""
    CM = mech()
    import inspect
    src = inspect.getsource(CM.optics)
    src = src.replace("pd.read_csv(N.OUT / 'optical_scenario_composition.csv'", "pd.read_csv(SRC / 'optical_scenario_composition.csv'")
    src = src.replace("root = N.folder(clearance) / 'adjacency/ew_only'", "root = SRC / clearance / 'adjacency/ew_only'")
    assert src.count('SRC /') == 2
    env = dict(CM.__dict__); env['SRC'] = SRC
    exec(compile(src, '<v12_optics>', 'exec'), env)
    env['optics']()


def aux():
    """e2_spec_analysis.main with the v1.2 units, construction-year ages and hid clusters."""
    import orientation_diagnostics as OD
    OD.install()
    route('image')
    import e2_spec_analysis as E2
    tidx, T, M, _ = V.merged()
    hid = V.unit_table().hand_hid
    OD.rewrite(E2, 'main', [("    tidx, Td, M = cc.matched('det')\n", "    tidx, Td, M = V12_MATCHED\n"),
                            ("        t['clu'] = np.where(t.station > 0, t.station, 10 ** 9 + blk_all[t.cell.values])\n",
                             "        t['clu'] = V12_HID.reindex(t.cell.values).values.astype(np.int64)\n")],
               {'V12_MATCHED': (tidx, T, M), 'V12_HID': hid})
    dates = ['20250528', '20250628', '20250721', '20250824']
    for p, metres in [('1p5m', 1.5), ('1m', 1.0), ('2m', 2.0)]:
        src = SRC / p / 'adjacency/spectral'
        dest = CS12 / 'aux' if p == '1p5m' else CS12 / 'secondary' / 'aux' / p
        dest.mkdir(parents=True, exist_ok=True)
        for dte in dates:
            for typ in ['prof', 'ctrl']:
                shutil.copy2(src / f'e1_{dte}_{typ}.parquet', dest / f'e1_{dte}_{typ}.parquet')
        shutil.copy2(src / 'e1_veg.csv', dest / 'e1_veg.csv')
        E2.SP = str(dest); E2.CLEARANCE_M = metres; E2.MIN_PIXELS = 50
        root = CS12 if p == '1p5m' else CS12 / p
        elig = {t: set(pd.read_csv(root / f'cells_{t}.csv').query('valid_cross').cell) for t in ['06', '08']}
        E2.ELIGIBLE_CELLS = {dte: elig['06' if dte in dates[:2] else '08'] for dte in dates}
        E2.main()
        if p == '1p5m':
            (CS12 / 'secondary' / 'aux' / '1p5m').mkdir(parents=True, exist_ok=True)
            shutil.copy2(dest / 'e2_indices.csv', CS12 / 'secondary' / 'aux' / '1p5m' / 'e2_indices.csv')
        print('AUX', p, flush=True)


def annotate():
    u = V.unit_table(); tidx, T, M, origin = V.merged()
    org = pd.Series(np.where(origin == 1, 'phased', 'original'), index=tidx)
    for root in [CS12, CS12 / '1m', CS12 / '2m', V12 / 'buffer_primary']:
        for tag in TAGS:
            f = root / f'cells_{tag}.csv'
            if not f.exists():
                continue
            c = pd.read_csv(f)
            c['origin'] = org.reindex(c.cell).values
            c['hand_year'] = u.loc[c.cell, 'hand_year'].values; c['station_label'] = u.loc[c.cell, 'station'].values
            c['det_year'] = u.loc[c.cell, 'det_year'].values
            c.to_csv(f, index=False)
    print('ANNOTATED', flush=True)


def terrain():
    """Terrain and greening of the comparison grounds, August 2025, units aged >=6 years.

    Unit-weighted means of the unit's own slope, the mean slope of its usable matched controls and
    the mean slope of its far buffer ground (>1 km from arrays within 8 km, buffer_primary rules);
    greening = mean Landsat midsummer NDVI 2021-2025 minus 2000-2010 of the same grounds."""
    route('image')
    import buffer_primary as BP
    import buffer_baseline as BB
    cells = pd.read_parquet(os.path.join(cc.TS, 'cells_pool.parquet')); okq = BB.quality(cells)
    Y = np.load(os.path.join(cc.TS, 'cell_y.npz'))['sum_ndvi']; yrs = es.YEARS
    green = np.nanmean(Y[:, (yrs >= 2021) & (yrs <= 2025)], 1) - np.nanmean(Y[:, (yrs >= 2000) & (yrs <= 2010)], 1)
    slope = cells.slope.values
    d = N.load('08', '1p5m'); raw = d['tabs']['raw']
    bp = pd.read_csv(V12 / 'buffer_primary' / 'cells_08.csv'); assert np.array_equal(bp.cell.values, d['tidx'])
    ok = okq & (raw.ogt1000_n.values >= 2000)
    W = BP.matrix(d, cells, ok, 8000, d['native_mode'])
    sel = d['sel']['all']; M = np.clip(d['M'], 0, None)
    use = bp.valid_primary.values & (d['age'] >= 6) & (np.asarray(W.sum(1)).ravel() >= 5)
    def wavg(v):
        num = np.asarray(W @ np.nan_to_num(v)).ravel(); den = np.asarray(W @ np.isfinite(v).astype(float)).ravel()
        return np.where(den > 0, num / np.maximum(den, 1), np.nan)
    def cavg(v):
        x = np.where(sel, v[M], np.nan); return np.nanmean(x, 1)
    out = dict(units=int(use.sum()),
               slope_units=float(np.nanmean(slope[d['tidx']][use])), slope_controls=float(np.nanmean(cavg(slope)[use])), slope_far=float(np.nanmean(wavg(slope)[use])),
               green_units=float(np.nanmean(green[d['tidx']][use])), green_controls=float(np.nanmean(cavg(green)[use])), green_far=float(np.nanmean(wavg(green)[use])))
    (V12 / 'buffer_primary' / 'terrain.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    print('TERRAIN', out, flush=True)


def origin():
    """Primary contrasts by unit origin (original cohort, phased arrays), whole-station resampling.

    Also the v1.1 contrasts of the original cohort timed by the detected year, for the bridge."""
    import buffer_baseline as BB
    rows = []
    groups = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('all', 0, 99)]
    v11 = pd.read_csv(SRC / 'cross_scene' / 'effects.csv', dtype={'tag': str, 'group': str})
    for tag in TAGS:
        c = pd.read_csv(CS12 / f'cells_{tag}.csv')
        st = c.station.values; age = c.age.values
        for org in ['original', 'phased', 'all units']:
            om = np.ones(len(c), bool) if org == 'all units' else (c.origin.values == org)
            for k in ['ndvi', 'nirv', 'msavi']:
                x = c[f'cross_{k}_dadj'].values
                masks = {g: om & (age >= lo) & (age <= hi) for g, lo, hi in groups}
                bs = BB.station_boot(x, st, c.valid_cross.values.astype(float), masks, B=2000)
                for g, lo, hi in groups:
                    m = c.valid_cross.values & masks[g] & np.isfinite(x); nc = len(np.unique(st[m]))
                    bb = bs[g][np.isfinite(bs[g])]
                    l, h = np.percentile(bb, [2.5, 97.5]) if nc >= J.MIN_CL else (np.nan, np.nan)
                    r = dict(tag=tag, origin=org, index=k, group=g, estimate=BB.wmean(x, m.astype(float)), station_lo=l, station_hi=h, n_cells=int(m.sum()), clusters=nc)
                    if org == 'original':
                        z = v11[(v11.tag == tag) & (v11['index'] == k) & (v11.group == g)]
                        r['v11_detected_year'] = float(z.estimate.iloc[0]) if len(z) else np.nan
                        r['v11_units'] = int(z.n_cells.iloc[0]) if len(z) else 0
                    rows.append(r)
    pd.DataFrame(rows).to_csv(CS12 / 'origin.csv', index=False)
    print('ORIGIN done', flush=True)


STAGES = ['estimate', 'buffer', 'lanes', 'spectral', 'distance', 'aux', 'optics', 'annotate', 'terrain', 'origin']

if __name__ == '__main__':
    for s in (sys.argv[1:] or STAGES):
        t0 = time.time(); globals()[s](); print('STAGE_DONE', s, f'{time.time() - t0:.0f}s', flush=True)
