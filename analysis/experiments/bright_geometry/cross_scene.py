"""Cross-scene matched design: the original ten matched controls, wherever they were imaged.

The primary analysis kept only original controls imaged in the treated unit's own
scene, which removed 1,241-1,718 units per acquisition. Here every original control
with >=2,000 valid 2 m ground pixels more than 300 m from arrays is used, whatever its
scene or satellite (all scenes are calibrated to same-day Sentinel-2 reflectance;
open-ground scene biases -0.001 to +0.003 NDVI). The controls are those of the
mixed-pixel event study, so its held-out, pre-trend and placebo checks (Table 3)
apply unchanged. Because treated and control reflectances may come from two
independently calibrated scenes, the calibration standard error is doubled in
variance (multiplied by sqrt 2) before adding it to the sampling interval.

Population: parent units in the acquisition scope with >=50 native footprints at
>=1.5 m clearance and >=2 usable original controls. Estimand: the primary one.
Uncertainty: 600 joint resamples with the primary seed.
Outputs: analysis/cross_scene/{effects.csv, cells_<tag>.csv, complete.json}
"""
# Mirrored from experiments/bright_geometry/cross_scene.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import json
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_statistics as N

E = N.E; cc = E.cc; es = cc.es; J = E.J
TAGS = ['06', '08', '2026_08']
METRICS = ['ndvi', 'nirv', 'msavi', 'red', 'nir']
GROUPS = [(str(a), a, a) for a in range(15)] + [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('≥10', 10, 99), ('all', 0, 99)]
OUT = N.OUT / 'cross_scene'
SEED = 20261010
CAL = np.sqrt(2)


def wmean(x, w):
    ok = np.isfinite(x) & (w > 0)
    return float(np.sum(x[ok] * w[ok]) / w[ok].sum()) if ok.any() else np.nan


def main(B=600, clearance='1p5m'):
    out = OUT if clearance == '1p5m' else OUT / clearance
    out.mkdir(parents=True, exist_ok=True)
    rows, audit = [], []
    for tag in TAGS:
        t0 = time.time()
        d = N.load(tag, clearance)
        raw = d['tabs']['raw']; tidx, age, ct = d['tidx'], d['age'], d['ct']
        irn = d['native_sums'][2]; scope = cc.scope_filter(d['cells'], tag)[tidx]
        sel = d['sel']['all']
        nsame = d['sel']['same'].sum(1)
        pop = d['parent'] & scope & (irn >= 50) & (sel.sum(1) >= 2)
        ground = {k: cc.index_from_sums(*d['native_sums'])[k] for k in METRICS}

        def contrasts(wp=None):
            out = {}
            for k in METRICS:
                r, _ = E.A1.per_cell(d, 'all', k, wp)
                out[k] = (ground[k] - r['ic'] - r['dpre'], r['ic'], r['dpre'], r['raw'])
            return out
        pts = contrasts()
        valid = pop & np.isfinite(pts['ndvi'][0])
        st = np.where(ct.station.to_numpy() > 0, ct.station, 10 ** 9 + d['bt'])
        unique, inv = np.unique(st, return_inverse=True)
        masks = {g: (age >= lo) & (age <= hi) for g, lo, hi in GROUPS}
        boot = {(k, g): [] for k in METRICS for g in masks}
        rng = np.random.default_rng(SEED)
        for b in range(B):
            wp = es.multiplicities(d['blk_all'], rng)
            wt = np.bincount(rng.integers(0, len(unique), len(unique)), minlength=len(unique)).astype(float)[inv]
            pb = contrasts(wp)
            for k in METRICS:
                for g, m in masks.items():
                    boot[k, g].append(wmean(pb[k][0], wt * valid * m))
            if (b + 1) % 200 == 0:
                print('CROSS_SCENE', tag, b + 1, B, f'{time.time() - t0:.0f}s', flush=True)
        vm = d['valid']['same']
        orig = {k: E.A1.per_cell(d, 'same', k)[0] for k in METRICS}
        for k in METRICS:
            x = pts[k][0]; xo = ground[k] - orig[k]['ic'] - orig[k]['dpre']
            for g, m in masks.items():
                ok = valid & m & np.isfinite(x); nc = len(np.unique(st[ok]))
                est = wmean(x, ok.astype(float)); bb = np.asarray(boot[k, g]); bb = bb[np.isfinite(bb)]
                enough = nc >= J.MIN_CL and len(bb) > 50
                lo, hi = np.percentile(bb, [2.5, 97.5]) if enough else (np.nan, np.nan)
                se = float(bb.std()) if enough else np.nan
                tot = np.sqrt(se ** 2 + (CAL * J.SYS.get(k, 0.)) ** 2)
                common = ok & vm
                rows.append(dict(tag=tag, clearance=clearance, index=k, group=g, n_cells=int(ok.sum()), n_station_clusters=nc, estimate=est, joint_lo=lo, joint_hi=hi, joint_se=se,
                                 joint_plus_calibration_lo=est - 1.96 * tot if enough else np.nan, joint_plus_calibration_hi=est + 1.96 * tot if enough else np.nan,
                                 control_level=wmean(pts[k][1], ok.astype(float)), mixed=wmean(pts[k][3], ok.astype(float)),
                                 positive_share=float((x[ok] > 0).mean()) if ok.any() else np.nan, selected_km2=float(irn[ok].sum() * 4 / 1e6),
                                 n_stations=int(len(np.unique(ct.station.to_numpy()[ok][ct.station.to_numpy()[ok] > 0]))),
                                 common_with_original=int(common.sum()), cross_on_common=wmean(x, common.astype(float)), original_on_common=wmean(xo, common.astype(float))))
        u = pd.DataFrame(dict(cell=tidx, station=ct.station.to_numpy(), T=d['T'], age=age, x=ct.x.to_numpy(), y=ct.y.to_numpy(), valid_cross=valid,
                              valid_original=vm, n_controls=sel.sum(1), n_same_scene=nsame, native_footprints=irn))
        for k in METRICS:
            u[f'cross_{k}_dadj'] = pts[k][0]; u[f'cross_{k}_mixed'] = pts[k][3]
        u.to_csv(out / f'cells_{tag}.csv', index=False)
        audit.append(dict(tag=tag, units=int(valid.sum()), original=int(vm.sum()), overlap=int((valid & vm).sum()), parents=int(d['parent'].sum()),
                          median_controls=float(np.median(sel.sum(1)[valid])), selected_km2=float(irn[valid].sum() * 4 / 1e6), seconds=round(time.time() - t0)))
        print('CROSS_SCENE_DONE', audit[-1], flush=True)
    pd.DataFrame(rows).to_csv(out / 'effects.csv', index=False)
    (out / 'complete.json').write_text(json.dumps(dict(status='complete', B=B, seed=SEED, calibration_factor='sqrt(2) x frozen SYS',
        rule='original ten matched controls, any scene, >=2,000 valid >300 m ground pixels; >=2 controls; >=50 native footprints', runs=audit), indent=2), encoding='utf-8')
