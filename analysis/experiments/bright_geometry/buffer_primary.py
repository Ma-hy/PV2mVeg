"""Unmatched buffer baseline on exactly the primary units (v1.1).

Population: the primary cross-scene units (cross_scene/cells_<tag>.csv, valid_cross).
Comparison ground: every 150 m cell within R of the unit whose ground at a given distance
from the mapped arrays of the same acquisition supplies >=2,000 native 2 m pixels from the
unit's modal native image, subject to the land-cover quality rules of the control pool
(buffer_baseline.quality). Pre-construction Landsat difference against the same cells.
Main: ring 300-1,000 m, R = 5 km; 600 joint resamples with the primary seed, paired with
the primary matched design. Variants (2,000 whole-station resamples, ground fixed):
R = 3 or 8 km; ground >1 km within 8 km (with and without the pool's slope, elevation and
pre-construction support rules); whole cell >=500 m from arrays; pixel-count weights.
Outputs: analysis/buffer_primary/{effects.csv, variants.csv, cells_<tag>.csv, complete.json}
"""
# Mirrored from experiments/bright_geometry/buffer_primary.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import json
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import buffer_baseline as BB
import buffer_variants as BV

N, E, cc, es, J = BB.N, BB.E, BB.cc, BB.es, BB.J
TAGS = ['06', '08', '2026_08']
GROUPS = [(str(a), a, a) for a in range(15)] + [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('all', 0, 99)]
METRICS = BB.METRICS
OUT = N.OUT / 'buffer_primary'
SEED = 20261010


def matrix(d, cells, ok, radius, scenes):
    raw = d['tabs']['raw']; idx = np.flatnonzero(ok); xy = cells[['x', 'y']].values
    tree = cKDTree(xy[idx]); cs = raw.ctrl_scene.values
    rows, cols = [], []
    for i, (lst, s) in enumerate(zip(tree.query_ball_point(xy[d['tidx']], radius), scenes)):
        c = idx[lst]; c = c[cs[c] == s]
        rows.append(np.full(len(c), i)); cols.append(c)
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(d['tidx']), len(cells)))


def main(B=600):
    OUT.mkdir(parents=True, exist_ok=True)
    cells = pd.read_parquet(os.path.join(cc.TS, 'cells_pool.parquet'))
    okq = BB.quality(cells); supp = BV.support_rules(cells)
    eff, var, audit = [], [], []
    for tag in TAGS:
        t0 = time.time()
        d = N.load(tag, '1p5m'); raw = d['tabs']['raw']
        cs = pd.read_csv(N.OUT / 'cross_scene' / f'cells_{tag}.csv'); assert np.array_equal(cs.cell.values, d['tidx'])
        tidx, age, ct = d['tidx'], d['age'], d['ct']
        st = np.where(ct.station.to_numpy() > 0, ct.station, 10 ** 9 + d['bt'])
        unique, inv = np.unique(st, return_inverse=True)
        mode = d['native_mode']
        ground = {k: cc.index_from_sums(*d['native_sums'])[k] for k in METRICS}
        prim = {k: E.A1.per_cell(d, 'all', k)[0] for k in METRICS}
        pm = {k: ground[k] - prim[k]['ic'] - prim[k]['dpre'] for k in METRICS}
        vc = cs.valid_cross.values
        band_ok = raw.o300_1000_n.values >= BB.MIN_PIX
        specs = {'ring_R5km': ('o300_1000', okq & band_ok, 5000, None),
                 'ring_R3km': ('o300_1000', okq & band_ok, 3000, None),
                 'ring_R8km': ('o300_1000', okq & band_ok, 8000, None),
                 'beyond1km_R8km': ('ogt1000', okq & (raw.ogt1000_n.values >= BB.MIN_PIX), 8000, None),
                 'beyond1km_R8km_support': ('ogt1000', okq & supp & (raw.ogt1000_n.values >= BB.MIN_PIX), 8000, None),
                 'ring_cell500_R5km': ('o300_1000', okq & band_ok & (cells.dist_min.values >= 500), 5000, None),
                 'ring_R5km_pixel_weighted': ('o300_1000', okq & band_ok, 5000, 'pixels')}
        masks = {g: (age >= lo) & (age <= hi) for g, lo, hi in GROUPS}
        built = {}
        for name, (band, ok, radius, weight) in specs.items():
            W = matrix(d, cells, ok, radius, mode)
            if weight == 'pixels':
                W = W.multiply(raw[f'{band}_n'].values[None, :]).tocsr()
            bf = BB.Buffer(d, W, band); o = bf.compute()
            ncell = np.asarray((W > 0).sum(1)).ravel()
            x = {k: ground[k] - o[k]['ic'] - o[k]['dpre'] for k in METRICS}
            v = vc & (ncell >= BB.MIN_CELLS) & np.isfinite(x['ndvi'])
            built[name] = (bf, x, v, ncell, o)
            for k in ['ndvi', 'nirv', 'msavi']:
                bs = BB.station_boot(x[k], st, v.astype(float), {g: masks[g] for g in ['0–2', '3–5', '≥6', 'all']}, B=2000)
                for g in ['0–2', '3–5', '≥6', 'all']:
                    m = v & masks[g] & np.isfinite(x[k]); nc = len(np.unique(st[m])); bb = bs[g][np.isfinite(bs[g])]
                    lo, hi = np.percentile(bb, [2.5, 97.5]) if nc >= J.MIN_CL else (np.nan, np.nan)
                    var.append(dict(tag=tag, variant=name, index=k, group=g, estimate=BB.wmean(x[k], m.astype(float)), station_lo=lo, station_hi=hi,
                                    n_cells=int(m.sum()), primary_same_units=BB.wmean(pm[k], m.astype(float)), median_cells=float(np.median(ncell[m])) if m.any() else np.nan))
        bf, xb, vb, ncell, ob = built['ring_R5km']
        boot = {(s, k, g): [] for s in ['buffer', 'diff'] for k in METRICS for g in masks}
        rng = np.random.default_rng(SEED)
        for b in range(B):
            wp = es.multiplicities(d['blk_all'], rng)
            wt = np.bincount(rng.integers(0, len(unique), len(unique)), minlength=len(unique)).astype(float)[inv]
            obb = bf.compute(wp)
            for k in METRICS:
                rm, _ = E.A1.per_cell(d, 'all', k, wp)
                vm_ = ground[k] - rm['ic'] - rm['dpre']; vbk = ground[k] - obb[k]['ic'] - obb[k]['dpre']
                for g, m in masks.items():
                    boot['buffer', k, g].append(BB.wmean(vbk, wt * vb * m)); boot['diff', k, g].append(BB.wmean(vbk - vm_, wt * vb * m))
        for k in METRICS:
            for g, m in masks.items():
                ok = vb & m & np.isfinite(xb[k]); nc = len(np.unique(st[ok]))
                row = dict(tag=tag, index=k, group=g, n_cells=int(ok.sum()), n_station_clusters=nc,
                           n_stations=int(len(np.unique(ct.station.to_numpy()[ok][ct.station.to_numpy()[ok] > 0]))),
                           control_level=BB.wmean(ob[k]['ic'], ok.astype(float)), positive_share=float((xb[k][ok] > 0).mean()) if ok.any() else np.nan)
                for s, x in [('buffer', xb[k]), ('diff', xb[k] - pm[k])]:
                    est = BB.wmean(x, ok.astype(float)); bb = np.asarray(boot[s, k, g]); bb = bb[np.isfinite(bb)]
                    enough = nc >= J.MIN_CL and len(bb) > 50
                    lo, hi = np.percentile(bb, [2.5, 97.5]) if enough else (np.nan, np.nan); se = float(bb.std()) if enough else np.nan
                    row.update({f'{s}_estimate': est, f'{s}_joint_lo': lo, f'{s}_joint_hi': hi})
                    if s == 'buffer':
                        tot = np.sqrt(se ** 2 + J.SYS.get(k, 0.) ** 2)
                        row['buffer_joint_plus_calibration_lo'] = est - 1.96 * tot if enough else np.nan
                        row['buffer_joint_plus_calibration_hi'] = est + 1.96 * tot if enough else np.nan
                row['primary_same_units'] = BB.wmean(pm[k], ok.astype(float))
                eff.append(row)
        u = pd.DataFrame(dict(cell=tidx, station=ct.station.to_numpy(), age=age, x=ct.x.to_numpy(), y=ct.y.to_numpy(), valid_primary=vc, valid_buffer=vb,
                              n_buffer_cells=ncell, modal_scene=mode))
        for k in METRICS:
            u[f'primary_{k}_dadj'] = pm[k]; u[f'buffer_{k}_dadj'] = xb[k]
        u.to_csv(OUT / f'cells_{tag}.csv', index=False)
        audit.append(dict(tag=tag, primary=int(vc.sum()), buffer=int(vb.sum()), median_cells=float(np.median(ncell[vb])), seconds=round(time.time() - t0)))
        print('BUFFER_PRIMARY', audit[-1], flush=True)
    pd.DataFrame(eff).to_csv(OUT / 'effects.csv', index=False)
    pd.DataFrame(var).to_csv(OUT / 'variants.csv', index=False)
    (OUT / 'complete.json').write_text(json.dumps(dict(status='complete', B=B, seed=SEED, population='primary cross-scene units', runs=audit), indent=2), encoding='utf-8')
