"""v1.2 analysis units and the legacy-format data dictionary.

Units: the 4,730 units of the original 4,759-unit east-west cohort that carry one hand-drawn
construction year, plus the 2,814 phased east-west units of the coverage expansion (EW_HY), whose
30 m pixels share one construction year but not one Landsat-detected year. Both sets follow the
same grid, panel-fraction, east-west and exclusion rules.

Timing: the construction year of the hand-drawn station polygons, assigned by the first author
from year-by-year visual reading of Landsat imagery (the first year at whose end the station is
complete, the previous year incomplete or not started); a unit takes the year shared by at least
20 of its 25 Landsat pixels. Array age = acquisition year - construction year.
Controls: the matching procedure of the mixed-pixel event study (es.match) on the
construction-year window (coverage expansion x03_match); units whose construction and detected
years coincide keep their archived controls.
Clusters: hand-drawn station polygon (hid); the legacy column ct.station holds hid so that every
legacy estimator resamples and averages by hand-drawn station (the label is kept in station_label).
Native footprints: release native sums for the original cohort and coverage-expansion native sums
for the phased units (identical rules, see v12_measure).

snapshot()  copies the coverage-expansion inputs into analysis/v12/inputs with SHA-256 manifest
validate()  per-unit primary contrasts equal the coverage-expansion estimates (parent_hand, ew_hy)
"""
# Mirrored from experiments/bright_geometry/v12_units.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
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
import hashlib
from pathlib import Path
from functools import lru_cache
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_statistics as N

E = N.E; cc = E.cc; es = cc.es; A1 = E.A1
ROOT = PROJECT  # project data root (pv2mveg configuration; read only)
SRC = ROOT / 'data/eco/releases/fixed_ew_v2_orientation/analysis'
V12 = V12_OUT  # output folder of the v1.2 analysis (pv2mveg configuration)
INP = V12 / 'inputs'
CE = ROOT / 'data/eco/coverage_expansion'
TAGS = ['06', '08', '2026_08']
CLEAR = ['1m', '1p5m', '2m']
METRICS = ['ndvi', 'nirv', 'msavi', 'red', 'nir']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def snapshot():
    INP.mkdir(parents=True, exist_ok=True)
    files = {'matched_parent_hand.npz': CE / 'match/matched_parent_hand.npz', 'matched_ew_hy.npz': CE / 'match/matched_ew_hy.npz',
             'cells_expansion.parquet': CE / 'units/cells_expansion.parquet',
             'parent_hand_effects.csv': CE / 'estimates/parent_hand/effects.csv', 'ew_hy_effects.csv': CE / 'estimates/ew_hy/effects.csv'}
    for t in TAGS:
        files[f'native_ew_hy_{t}.npz'] = CE / f'native/ew_hy/{t}/native_sums.npz'
        files[f'profile_ew_{t}.npz'] = CE / f'mount/profile/ew/{t}/profile_sums.npz'
        files[f'parent_hand_cells_{t}.csv'] = CE / f'estimates/parent_hand/cells_{t}.csv'
        files[f'ew_hy_cells_{t}.csv'] = CE / f'estimates/ew_hy/cells_{t}.csv'
    man = {}
    for name, src in files.items():
        dst = INP / name
        if not dst.exists():
            shutil.copy2(src, dst)
        man[name] = dict(source=str(src.relative_to(ROOT)), sha256=sha(dst))
    (INP / 'manifest.json').write_text(json.dumps(man, indent=2), encoding='utf-8')
    print('SNAPSHOT', len(man), flush=True)


@lru_cache(maxsize=1)
def unit_table():
    return pd.read_parquet(INP / 'cells_expansion.parquet').set_index('cell')


@lru_cache(maxsize=1)
def merged():
    a = np.load(INP / 'matched_parent_hand.npz'); b = np.load(INP / 'matched_ew_hy.npz')
    tidx = np.concatenate([a['tidx'], b['tidx']]).astype(np.int64)
    T = np.concatenate([a['T'], b['T']]).astype(int)
    M = np.concatenate([a['M'], b['M']]).astype(np.int64)
    origin = np.r_[np.zeros(len(a['tidx']), int), np.ones(len(b['tidx']), int)]
    assert len(np.unique(tidx)) == len(tidx)
    o = np.argsort(tidx, kind='stable')
    return tidx[o], T[o], M[o], origin[o]


@lru_cache(maxsize=9)
def native(tag, clearance):
    j = CLEAR.index(clearance)
    zr = np.load(SRC.parent / 'native' / tag / 'native_sums.npz'); zn = np.load(INP / f'native_ew_hy_{tag}.npz')
    sums = zr['sums'][j].copy(); votes = zr['scenes'][j].copy()
    tidx, _, _, origin = merged(); new = tidx[origin == 1]
    assert float(np.abs(sums[:, new]).sum()) == 0.0
    sums[:, new] = zn['sums'][j][:, new]; votes[new] = zn['scenes'][j][new]
    return sums, votes


@lru_cache(maxsize=1)
def cells_pool():
    return pd.read_parquet(os.path.join(cc.TS, 'cells_pool.parquet'))


def load_v12(tag, clearance='1p5m', threshold=50, baseline=False, same='primary'):
    """Legacy data dictionary (orientation_statistics.load layout) of the v1.2 units.

    same='primary': the legacy 'same' mode is the primary cross-scene design (as in
    cross_scene_secondary); same='image': it is the subset of controls imaged in the scene of the
    unit's native footprints."""
    assert not baseline
    cells = cells_pool()
    tidx, T, M, origin = merged()
    year = cc.PERIODS[tag]['year']; age = year - T; keep = age >= 0
    tidx, T, M, origin, age = tidx[keep], T[keep], M[keep], origin[keep], age[keep]
    tabs = {v: pd.read_parquet(os.path.join(cc.OUT, f'celltab_{tag}_{v}.parquet')) for v in ('raw', 'fe_red', 'fe_nir')}
    sums, votes = native(tag, clearance)
    for t in tabs.values():
        for f, a in zip(['sr', 'sn', 'n'], sums):
            t['ir_' + f] = a
        t['ir_scene_frac'] = votes.max(1) / np.maximum(votes.sum(1), 1)
        cc.add_index(t, 'ir'); cc.add_index(t, 'ctrl')
    ct = cells.iloc[tidx].reset_index(drop=True)
    u = unit_table().loc[tidx]
    hid = u.hand_hid.values.astype(np.int64); assert (hid > 0).all()
    ct['station_label'] = ct['station'].values; ct['station'] = hid; ct['T_hand'] = T
    blk_all = es.block_id(cells.x.values, cells.y.values)
    R = tabs['raw']
    ctrl_ok = R.ctrl_n.values >= A1.CTRL_MIN_N
    mode = np.where(votes[tidx].sum(1) > 0, votes[tidx].argmax(1), 0)
    allc = ctrl_ok[np.clip(M, 0, None)] & (M >= 0)
    same_img = allc & (R.ctrl_scene.values[np.clip(M, 0, None)] == mode[:, None]) & (mode[:, None] > 0)
    irn = sums[2, tidx]
    scope = cc.scope_filter(cells, tag)[tidx]
    pop = scope & (irn >= threshold) & (allc.sum(1) >= 2)
    pop_img = scope & (irn >= threshold) & (same_img.sum(1) >= 2)
    sel = {'all': allc, 'same': allc if same == 'primary' else same_img, 'same_image': same_img, 'fe_red': allc, 'fe_nir': allc}
    valid = {'all': pop, 'same': pop if same == 'primary' else pop_img, 'same_image': pop_img,
             'fe_red': np.zeros(len(tidx), bool), 'fe_nir': np.zeros(len(tidx), bool)}
    tab = {'all': 'raw', 'same': 'raw', 'same_image': 'raw', 'fe_red': 'fe_red', 'fe_nir': 'fe_nir'}
    return dict(cells=cells, ct=ct, tidx=tidx, T=T, M=M, age=age, tabs=tabs, sel=sel, tab=tab, valid=valid, blk_all=blk_all, bt=blk_all[tidx],
                grp=hid, YL=cc.landsat_Y(), year=year, release='v12', clearance=clearance, parent=np.ones(len(tidx), bool),
                reference_scene=R.ir_scene.values[tidx], native_mode=mode, native_sums=sums[:, tidx], native_scene_votes=votes[tidx],
                origin=origin, design='v1.2: original cohort timed by construction year + phased east-west units; primary cross-scene controls')


def validate():
    """Per-unit primary contrasts equal the coverage-expansion estimates (parent_hand and ew_hy)."""
    rep = {}
    u = unit_table()
    for tag in TAGS:
        d = load_v12(tag)
        r, _ = A1.per_cell(d, 'all', 'ndvi')
        x = r['it'] - r['ic'] - r['dpre']
        v = d['valid']['all'] & np.isfinite(x)
        scope_ce = ~u.loc[d['tidx'], 'excl_2026_08'].values if tag == '2026_08' else np.ones(len(d['tidx']), bool)
        assert np.array_equal(cc.scope_filter(d['cells'], tag)[d['tidx']], scope_ce), tag
        res = {}
        for name, o in [('parent_hand', 0), ('ew_hy', 1)]:
            ref = pd.read_csv(INP / f'{name}_cells_{tag}.csv').set_index('cell')
            m = d['origin'] == o
            ids = d['tidx'][m]
            assert set(ids) == set(ref.index), (tag, name, len(set(ids) ^ set(ref.index)))
            rv = ref.loc[ids, 'valid_cross'].values.astype(bool)
            assert np.array_equal(v[m], rv), (tag, name, int((v[m] ^ rv).sum()))
            diff = np.abs(x[m][rv] - ref.loc[ids, 'cross_ndvi_dadj'].values[rv])
            assert diff.max() < 1e-12, (tag, name, float(diff.max()))
            res[name] = dict(units=int(m.sum()), valid=int(rv.sum()), max_abs_diff=float(diff.max()))
        rep[tag] = res
        print('VALID', tag, res, flush=True)
    (V12 / 'units_validation.json').write_text(json.dumps(rep, indent=2), encoding='utf-8')
    return rep


if __name__ == '__main__':
    what = sys.argv[1:] or ['snapshot', 'validate']
    for w in what:
        globals()[w]()
