"""v1.2: supplementary diagnostics on the merged east-west units (original cohort + phased units).

Copy of cross_scene_secondary.py (v1.1) with three changes: the loader is v12_units.load_v12
(construction-year timing, hand-drawn station clusters, primary controls in the legacy 'same'
mode), the output root is analysis/v12/cross_scene/secondary, and the subset measurements are the
merged ones of v12_measure. Common ground footprints are measured for the original cohort only
(the release clearance layers cover those units); the auxiliary red edge is run by v12_run.aux.

v1.1 text follows.

Supplementary diagnostics re-run under the primary (cross-scene) matched design (v1.1).

v0.7-v0.12 computed the supplementary diagnostics with the original matched controls
imaged in the treated unit's own scene. The primary design of v1.1 (cross_scene.py) uses
every usable original control in any image. This adapter re-runs the same v0.12 code,
with the orientation_diagnostics adapters installed, after one substitution: the control
selection and unit eligibility of the legacy 'same' mode are replaced by those of the
primary design. Every legacy estimator therefore returns primary-design values.

Nothing under the v0.12 analysis tree is overwritten. Outputs use the legacy layout under
analysis/cross_scene/secondary/. In those legacy-format files the column prefix 'same_'
and the flag 'valid_same' denote the primary design, and 'ir_scene' is the modal scene of
the unit's native ground footprints (the original reference scene is kept in
'reference_scene'), so that source-satellite conversions follow the measured ground.

Stages, in order:
  tables       legacy-format unit tables (validated against cross_scene/cells_<tag>.csv)
  stats        revision_statistics.base_stats: groups, regressions, leave-one-station-out,
               station means, pixel-count and control-screen alternatives, mixed bias
  balance      revision_balance: covariate balance and held-out pre-construction checks
  placebo      revision_placebo: leave-one-control-out contrasts
  subsets      boundary perturbations, distance bands, lanes, within-unit contrasts and
               visible-band differences on the primary eligibility (revision_subsets rules)
  common       ≥1/1.5/2 m comparison on units valid at all three clearances, 600 joint
               resamples plus doubled calibration variance, as in cross_scene.py
  aux          auxiliary 1.2 m indices at 1 m and 2 m clearance (1.5 m: cross_scene/aux)
  crossperiod  common ground footprints between acquisitions (new units measured, cached
               v0.12 measurements reused); the footprint scene need not match a reference
  temporal     revision_temporal anchored reconstructions
  forward      revision_forward one-year forward/backward checks
The panchromatic digital-number diagnostics keep same-image controls, because digital
numbers are not calibrated between images; they are not re-run here.
"""
# Mirrored from experiments/bright_geometry/v12_secondary.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
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
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_diagnostics as OD
import orientation_statistics as N
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import v12_units as V

E = N.E; cc = E.cc; es = cc.es; J = E.J
TAGS = ['06', '08', '2026_08']
POL = ['1m', '1p5m', '2m']
SRC = Path(N.OUT)                       # v0.12 analysis root (read-only here)
V12 = V.V12
CS = V12 / 'cross_scene'
SEC = CS / 'secondary'
CS11 = SRC / 'cross_scene'
SEED = 20261010
CAL = np.sqrt(2)
_load = N.load
DESIGN = 'v1.2 primary: original cohort timed by construction year + phased east-west units; ten matched controls in any image, >=2,000 valid >300 m ground pixels, >=2 controls; hand-drawn station clusters'


def load_cross(tag, clearance='1p5m', threshold=50, baseline=False):
    """v1.2 legacy data dictionary whose 'same' mode is the primary (cross-scene) design."""
    return V.load_v12(tag, clearance, threshold, baseline, same='primary')


def route():
    OD.install()
    N.load = load_cross
    N.OUT = SEC
    for _p in POL:
        (SEC / _p / "statistics").mkdir(parents=True, exist_ok=True)
    OD.configure('1p5m')
    assert R.load is load_cross and Path(R.OUT) == SEC and Path(R.folder('1p5m')) == SEC / '1p5m'


def cells(p, tag):
    return pd.read_csv((CS if p == '1p5m' else CS / p) / f'cells_{tag}.csv')


def log(status, **kw):
    path = SEC / 'complete.json'
    state = json.loads(path.read_text()) if path.exists() else dict(design=DESIGN, stages={})
    state['stages'][status] = dict(time=time.strftime('%Y-%m-%d %H:%M:%S'), **kw)
    path.write_text(json.dumps(state, indent=2), encoding='utf-8')
    print('SECONDARY', status, kw, flush=True)


def relabel(path, column, mapping):
    t = pd.read_csv(path)
    if column in t:
        t[column] = t[column].replace(mapping)
    t.to_csv(path, index=False)


# ---------------------------------------------------------------- tables
def tables():
    audit = []
    for p in POL:
        for tag in TAGS:
            t = N.cell_table(tag, p, write=False)
            t['reference_scene'] = t['ir_scene']
            t['ir_scene'] = np.where(t.native_mode > 0, t.native_mode, t.ir_scene)
            c = cells(p, tag)
            assert np.array_equal(t.cell.values, c.cell.values), (p, tag)
            v = t.valid_same.values & np.isfinite(t.same_ndvi_dadj.values)
            assert np.array_equal(v, c.valid_cross.values), (p, tag, int(v.sum()), int(c.valid_cross.sum()))
            err = float(np.max(np.abs(t.same_ndvi_dadj.values[v] - c.cross_ndvi_dadj.values[v])))
            assert err < 1e-9, (p, tag, err)
            root = SEC / p / 'chrono'; root.mkdir(parents=True, exist_ok=True)
            t.to_csv(root / f'a1_{tag}_det_thr50_cells.csv', index=False)
            audit.append(dict(clearance=p, tag=tag, units=int(v.sum()), max_abs_diff_vs_cross_scene=err,
                              native_scene_differs_from_reference=int((v & (t.native_mode.values != t.reference_scene.values)).sum())))
    log('tables', audit=audit)


# ---------------------------------------------------------------- legacy statistics
def stats():
    for p in POL:
        S.base_stats(p)
        out = SEC / p / 'statistics'
        relabel(out / 'period_pairs.csv', 'support', {'same units; own-date pixels and same-scene controls': 'same units; own-date pixels and primary controls'})
        (out / 'base_complete.json').write_text(json.dumps(dict(clearance=p, bootstrap=J.B, design=DESIGN,
            interval='Secondary conditional station bootstrap; controls travel with treatment.')), encoding='utf-8')
    log('stats', clearances=POL)


def balance():
    import revision_balance as RB
    RB.main()
    for p in POL:
        out = SEC / p / 'statistics'
        for name in ['balance.csv', 'holdout.csv']:
            relabel(out / name, 'controls', {'same_scene': 'primary'})
        c = pd.read_csv(out / 'controls.csv').rename(columns={'median_same_controls': 'median_primary_controls', 'median_same_distance': 'median_primary_distance'})
        c.to_csv(out / 'controls.csv', index=False)
    log('balance', clearances=POL)


def placebo():
    import revision_placebo as RP
    RP.main()
    log('placebo')


# ---------------------------------------------------------------- subsets (revision_subsets rules, primary eligibility)
def subsets():
    import k_bg_greenness as K
    for tag in TAGS:
        raw = pd.read_parquet(V12 / 'measurements' / tag / 'subsets.parquet')
        d = load_cross(tag, '1p5m'); t = N.tables(tag, '1p5m').set_index('cell')
        raw = raw[raw.cell.isin(d['tidx'])]
        pos = pd.Series(np.arange(len(d['tidx'])), index=d['tidx'])
        nsel = d['sel']['same'].sum(1)
        B = cc.read_bands(cc.prod_path(tag, '30m_bg'))
        ci, cn = K.cell_vals(B, cc.CTRL_BANDS)
        bgsel = d['sel']['same'] & (cn[np.clip(d['M'], 0, None)] >= 2000) & (d['M'] >= 0)
        ctrl = {k: cc.ctrl_mean(v, d['M'], bgsel)[0] for k, v in ci.items()}
        sensor = {s['index']: s['sensor'] for s in cc.bands_json(tag)['scenes']}
        scope = cc.scope_filter(d['cells'], tag)
        stats_, cells_, changes, gaprows = [], [], [], []
        for name, z in raw.groupby('subset', sort=False):
            z = z.set_index('cell').join(t, rsuffix='_meta').reset_index()
            p = pos.loc[z.cell].values
            ref = d['tabs']['raw'].ir_scene.values[z.cell]
            valid = scope[z.cell] & d['parent'][p] & (z.n.values >= 50) & (nsel[p] >= 2)
            indices = cc.index_from_sums(z.s_red.values, z.s_nir.values, z.n.values)
            isgap = '_gap' in name
            for k in (['ndvi'] if isgap else R.E.METRICS):
                tr = indices[k]; co = z['same_' + k + '_ic'].values
                pre = z['same_' + k + '_d'].values - z['same_' + k + '_dadj'].values
                adj = tr - co - pre
                if isgap:
                    for group, lo, hi in S.GROUPS:
                        keep = valid & (z.age.values >= lo) & (z.age.values <= hi) & np.isfinite(adj)
                        gaprows.append(dict(tag=tag, subset=name, index=k, group=group, n=int(keep.sum()), estimate=float(np.mean(adj[keep])) if keep.any() else np.nan))
                else:
                    stats_.extend(dict(tag=tag, subset=name, index=k, quantity='dadj', **a) for a in S.boot(adj[valid], z.loc[valid]))
                    if name.startswith('boundary_'):
                        keep = valid & z.valid_same.values
                        change = adj - z['same_' + k + '_dadj'].values
                        changes.extend(dict(tag=tag, subset=name, index=k, **a) for a in S.boot(change[keep], z.loc[keep]))
                if k == 'ndvi':
                    cell = z[['cell', 'age', 'grp', 'n']].copy(); cell['subset'] = name; cell['valid'] = valid; cell['ndvi_dadj'] = adj; cell['ndvi'] = tr
                    cells_.append(cell)
            if name.startswith(('ir_', 'lane_')):
                a1 = np.array([sensor.get(int(s)) == 'BJ3A1' for s in ref])
                bgok = valid & (z.n_bg.values >= 50) & (bgsel.sum(1)[p] >= 2) & a1
                ix = K.indices(*[z['bg_' + b].values / np.maximum(z.n_bg.values, 1) for b in K.BANDS])
                for k, tr in ix.items():
                    diff = tr - ctrl[k][p]
                    rows = S.boot(diff[bgok], z.loc[bgok])
                    for a in rows:
                        lo, hi = next((lo, hi) for g, lo, hi in S.GROUPS if g == a['group'])
                        take = bgok & (z.age.values >= lo) & (z.age.values <= hi)
                        a['control'] = float(np.mean(ctrl[k][p][take])); a['treated'] = float(np.mean(tr[take]))
                        a['relative'] = a['estimate'] / a['control'] if k in K.BANDS else np.nan
                    stats_.extend(dict(tag=tag, subset=name, index=k, quantity='direct_spectral_difference', **a) for a in rows)
        out = SEC / 'measurements' / tag; out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(stats_).to_csv(out / 'subset_summary.csv', index=False)
        pd.DataFrame(changes).to_csv(out / 'boundary_paired.csv', index=False)
        pd.DataFrame(gaprows).to_csv(out / 'gap_stratified.csv', index=False)
        allcells = pd.concat(cells_, ignore_index=True); allcells.to_parquet(out / 'subset_cells.parquet', index=False)
        paired = []
        for a, b in [('ir_1p5m', 'lane_1p5m'), ('d1.5_2', 'd2_3'), ('d2_3_side1', 'd2_3_side3')]:
            x = allcells[(allcells.subset == a) & allcells.valid]; y = allcells[(allcells.subset == b) & allcells.valid]
            zz = x.merge(y, on='cell', suffixes=('_a', '_b')); zz['age'] = zz.age_a; zz['grp'] = zz.grp_a
            paired.extend(dict(tag=tag, a=a, b=b, quantity='b_minus_a', **v) for v in S.boot(zz.ndvi_dadj_b - zz.ndvi_dadj_a, zz))
        pd.DataFrame(paired).to_csv(out / 'within_unit.csv', index=False)
        print('SUBSETS', tag, flush=True)
    log('subsets', rule='primary eligibility: parent, scope, >=50 subset footprints, >=2 primary controls; visible bands: BJ3A1 reference scene, >=50 calibrated footprints, >=2 controls with >=2,000 visible pixels')


# ---------------------------------------------------------------- three-threshold common units
def common(B=600):
    rows = []
    groups = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('all', 0, 99)]
    for tag in TAGS:
        t0 = time.time()
        ds = {p: load_cross(tag, p) for p in POL}
        for p in POL:
            assert np.array_equal(ds[p]['tidx'], ds['1p5m']['tidx'])
        d0 = ds['1p5m']; age = d0['age']; ct = d0['ct']
        ground = {p: cc.index_from_sums(*ds[p]['native_sums'])['ndvi'] for p in POL}

        def contrast(p, wp=None):
            r, _ = E.A1.per_cell(ds[p], 'all', 'ndvi', wp)
            return ground[p] - r['ic'] - r['dpre']
        pts = {p: contrast(p) for p in POL}
        valid = {p: ds[p]['valid']['same'] & np.isfinite(pts[p]) for p in POL}
        com = valid['1m'] & valid['1p5m'] & valid['2m']
        st = np.where(ct.station.to_numpy() > 0, ct.station, 10 ** 9 + d0['bt'])
        unique, inv = np.unique(st, return_inverse=True)
        masks = {g: (age >= lo) & (age <= hi) for g, lo, hi in groups}
        boot = {(p, g): [] for p in POL for g in masks}
        rng = np.random.default_rng(SEED)
        for b in range(B):
            wp = es.multiplicities(d0['blk_all'], rng)
            wt = np.bincount(rng.integers(0, len(unique), len(unique)), minlength=len(unique)).astype(float)[inv]
            for p in POL:
                x = contrast(p, wp)
                for g, m in masks.items():
                    w = wt * com * m * np.isfinite(x)
                    boot[p, g].append(float(np.sum(np.where(w > 0, x, 0) * w) / w.sum()) if w.sum() > 0 else np.nan)
        for p in POL:
            for g, m in masks.items():
                ok = com & m
                nc = len(np.unique(st[ok])); est = float(pts[p][ok].mean()) if ok.any() else np.nan
                bb = np.asarray(boot[p, g]); bb = bb[np.isfinite(bb)]
                enough = nc >= J.MIN_CL and len(bb) > 50
                lo, hi = np.percentile(bb, [2.5, 97.5]) if enough else (np.nan, np.nan)
                se = float(bb.std()) if enough else np.nan
                tot = np.sqrt(se ** 2 + (CAL * J.SYS['ndvi']) ** 2)
                rows.append(dict(tag=tag, clearance=p, group=g, n_cells=int(ok.sum()), n_station_clusters=nc, estimate=est, joint_lo=lo, joint_hi=hi, joint_se=se,
                                 joint_plus_calibration_lo=est - 1.96 * tot if enough else np.nan, joint_plus_calibration_hi=est + 1.96 * tot if enough else np.nan))
        print('COMMON', tag, f'{time.time() - t0:.0f}s', flush=True)
    pd.DataFrame(rows).to_csv(SEC / 'threshold_common.csv', index=False)
    log('common', B=B, seed=SEED, rule='units valid at >=1, >=1.5 and >=2 m; joint resamples shared across clearances')


# ---------------------------------------------------------------- common footprints between acquisitions
def crossperiod(workers=6):
    import revision_crossperiod as X
    for a, b in [('06', '08'), ('08', '2026_08')]:
        out = SEC / 'crossperiod' / f'{a}_{b}'; out.mkdir(parents=True, exist_ok=True)
        old = SRC / 'crossperiod' / f'{a}_{b}'; old11 = CS11 / 'secondary' / 'crossperiod' / f'{a}_{b}'
        base = N.tables(a, '1m')
        eligible = {tag: set().union(*[set(N.tables(tag, p).query('valid_same').cell) for p in POL]) for tag in [a, b]}
        tidx_, _, _, origin_ = V.merged(); original = set(tidx_[origin_ == 0])
        t = base[base.cell.isin(eligible[a] & eligible[b] & original)]
        reused = 0
        for c in t.cell:
            f = out / f'{int(c)}.json'
            for src_ in (old11, old):
                if not f.exists() and (src_ / f'{int(c)}.json').exists():
                    shutil.copy2(src_ / f'{int(c)}.json', f); reused += 1
        jobs = [(a, b, SimpleNamespace(cell=int(z.cell), x=z.x, y=z.y), str(out / f'{int(z.cell)}.json')) for z in t.itertuples()]
        todo = [j for j in jobs if not Path(j[3]).exists()]
        print('CROSSPERIOD', a, b, 'units', len(jobs), 'reused', reused, 'to measure', len(todo), flush=True)
        if todo:
            with ProcessPoolExecutor(workers) as pool:
                for i, _ in enumerate(pool.map(OD.cross_worker, todo, chunksize=4), 1):
                    if i % 200 == 0:
                        print('CROSSPERIOD', a, b, i, len(todo), flush=True)
        rows = []
        for *_, path in jobs:
            rows.extend(json.loads(Path(path).read_text()))
        raw = pd.DataFrame(rows); raw.to_parquet(out / 'measurements.parquet', index=False)
        sums, records = [], []
        for p in POL:
            zz = []
            for tag in [a, b]:
                tt = N.tables(tag, p); z = raw[(raw.clearance == p) & (raw.tag == tag)].merge(tt, on='cell')
                z['valid'] = z.valid_same & (z.n >= 50)
                vals = cc.index_from_sums(z.sr.values, z.sn.values, z.n.values)
                for k, v in vals.items():
                    if f'same_{k}_ic' in z:
                        z[k + '_common'] = v - z['same_' + k + '_ic'] - (z['same_' + k + '_d'] - z['same_' + k + '_dadj'])
                zz.append(z)
            z = zz[0][zz[0].valid].merge(zz[1][zz[1].valid], on='cell', suffixes=('_a', '_b')); z['age'] = z.age_a; z['grp'] = z.grp_a
            for k in R.E.METRICS:
                ca = z[k + '_common_a']; cb = z[k + '_common_b']; oa = z['same_' + k + '_dadj_a']; ob = z['same_' + k + '_dadj_b']
                for name, val in [('common_a', ca), ('common_b', cb), ('common_change', cb - ca), ('own_change_same_units', ob - oa), ('support_component', (ob - oa) - (cb - ca))]:
                    sums.extend(dict(a=a, b=b, clearance=p, index=k, quantity=name, **v) for v in S.boot(val, z))
            z['clearance'] = p; records.append(z)
        pd.DataFrame(sums).to_csv(out / 'summary.csv', index=False)
        pd.concat(records).to_parquet(out / 'paired_cells.parquet', index=False)
        (out / 'complete.json').write_text(json.dumps(dict(complete=True, pairs=len(jobs), reused_v012=reused, reference_grid=a, design=DESIGN,
            validity='primary-valid in both acquisitions and >=50 common footprints; original cohort only (release clearance layers); the footprint scene need not match a reference scene')), encoding='utf-8')
    log('crossperiod')


def temporal():
    import revision_temporal as T
    T.temporal()
    proto = SEC / 'temporal' / 'protocol.json'
    p = json.loads(proto.read_text()); p['controls'] = 'Primary design: original matched controls in any image; no rematching.'; p['design'] = DESIGN
    p['conversion'] = 'Frozen outside-array source-satellite-to-OLI coefficients, selected per unit from the scene of its native ground footprints.'
    proto.write_text(json.dumps(p, indent=2), encoding='utf-8')
    log('temporal')


def forward():
    import revision_forward as F
    F.main()
    log('forward')


STAGES = ['tables', 'stats', 'balance', 'placebo', 'subsets', 'common', 'crossperiod', 'temporal', 'forward']

if __name__ == '__main__':
    route()
    for stage in (sys.argv[1:] or STAGES):
        t0 = time.time()
        globals()[stage]()
        print('STAGE_DONE', stage, f'{time.time() - t0:.0f}s', flush=True)
