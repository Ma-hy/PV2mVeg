"""Mechanism diagnostics under the primary (cross-scene) matched design (v1.1).

All use the original ten matched controls in any image (cross_scene.py) unless stated.
  spectral   four-band relative inter-row minus control reflectance (blue/green need the
             BJ3A1 visible calibration: units and controls in BJ3A1 images only)
  distance   adjusted NDVI in centre-distance bands 1-1.5, 1.5-2, 2-3, 3-4, 4-6 m
  aux        red-edge indices of the auxiliary 1.2 m imagery for the primary units
             (e2_spec_analysis with the primary eligibility; its matched controls are the
             original controls inside the auxiliary image, as before)
  optics     imaging-scenario artefacts recomputed with the primary control spectrum
             (mean of all usable original controls) for the units with forward optical
             models, and compared with their primary observed contrast
Outputs in analysis/cross_scene/: spectral.csv, distance.csv, aux/e2_indices.csv, optics.csv
"""
# Mirrored from experiments/bright_geometry/cross_scene_mechanisms.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import sys
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent))
import orientation_statistics as N

E = N.E; cc = E.cc; J = E.J
TAGS = ['06', '08', '2026_08']
GROUPS = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99), ('all', 0, 99)]
OUT = N.OUT / 'cross_scene'
BANDS4 = ['blue', 'green', 'red', 'nir']


def station_ci(x, clu, B=2000, seed=20261010):
    x = np.asarray(x, float); ok = np.isfinite(x); x, clu = x[ok], np.asarray(clu)[ok]
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0, 0
    g, inv = np.unique(clu, return_inverse=True)
    if len(g) < J.MIN_CL:
        return float(x.mean()), np.nan, np.nan, len(x), len(g)
    s = np.bincount(inv, x); c = np.bincount(inv)
    rng = np.random.default_rng(seed); idx = rng.integers(0, len(g), (B, len(g)))
    bs = s[idx].sum(1) / c[idx].sum(1)
    return float(x.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), len(x), len(g)


def base(tag, clearance='1p5m'):
    d = N.load(tag, clearance)
    cs = pd.read_csv((OUT if clearance == '1p5m' else OUT / clearance) / f'cells_{tag}.csv')
    assert np.array_equal(cs.cell.values, d['tidx'])
    ct = d['ct']
    clu = np.where(ct.station.to_numpy() > 0, ct.station, 10 ** 9 + d['bt'])
    return d, cs, clu


def spectral():
    import k_bg_greenness as K
    rows = []
    for tag in TAGS:
        d, cs, clu = base(tag)
        raw = pd.read_parquet(N.OUT / 'measurements' / tag / 'subsets.parquet')
        z = raw[raw.subset == 'ir_1p5m'].set_index('cell').reindex(d['tidx'])
        B = cc.read_bands(cc.prod_path(tag, '30m_bg'))
        ci, cn = K.cell_vals(B, cc.CTRL_BANDS)
        sel = d['sel']['all'] & (cn[np.clip(d['M'], 0, None)] >= 2000) & (d['M'] >= 0)
        ctrl = {b: cc.ctrl_mean(ci[b], d['M'], sel)[0] for b in BANDS4}
        ref = d['tabs']['raw'].ir_scene.values[d['tidx']]
        sensor = {s['index']: s['sensor'] for s in cc.bands_json(tag)['scenes']}
        a1 = np.array([sensor.get(int(s)) == 'BJ3A1' for s in ref])
        ok = cs.valid_cross.values & (z.n_bg.values >= 50) & (sel.sum(1) >= 2) & a1
        tr = {b: z[f'bg_{b}'].values / np.maximum(z.n_bg.values, 1) for b in BANDS4}
        for g, lo, hi in GROUPS:
            m = ok & (d['age'] >= lo) & (d['age'] <= hi)
            for b in BANDS4:
                diff = tr[b] - ctrl[b]
                est, l, h, n, nc = station_ci(diff[m], clu[m])
                c0 = float(np.nanmean(ctrl[b][m]))
                rows.append(dict(tag=tag, group=g, band=b, estimate=est, lo=l, hi=h, control=c0, relative=est / c0, n=n, clusters=nc))
        print('SPECTRAL', tag, int(ok.sum()), flush=True)
    pd.DataFrame(rows).to_csv(OUT / 'spectral.csv', index=False)


def distance():
    rows = []
    bands = ['d1_1.5', 'd1.5_2', 'd2_3', 'd3_4', 'd4_6']
    for tag in TAGS:
        d, cs, clu = base(tag)
        raw = pd.read_parquet(N.OUT / 'measurements' / tag / 'subsets.parquet')
        r, _ = E.A1.per_cell(d, 'all', 'ndvi')
        for name in bands:
            z = raw[raw.subset == name].set_index('cell').reindex(d['tidx'])
            idx = cc.index_from_sums(z.s_red.values, z.s_nir.values, z.n.values)['ndvi']
            x = idx - r['ic'] - r['dpre']
            ok = cs.valid_cross.values & (z.n.values >= 50) & np.isfinite(x)
            for g, lo, hi in GROUPS:
                m = ok & (d['age'] >= lo) & (d['age'] <= hi)
                est, l, h, n, nc = station_ci(x[m], clu[m])
                rows.append(dict(tag=tag, band=name, d_lo=float(name[1:].split('_')[0]), d_hi=float(name.split('_')[1]), group=g, estimate=est, lo=l, hi=h, n=n, clusters=nc))
        print('DISTANCE', tag, flush=True)
    pd.DataFrame(rows).to_csv(OUT / 'distance.csv', index=False)


def optics():
    """Re-aggregate the frozen forward-model weights with the primary control spectra."""
    sys.path.insert(0, str(N.ROOT if hasattr(N, 'ROOT') else Path(__file__).resolve().parents[2]) )
    import a5_model_set as A5
    comp = pd.read_csv(N.OUT / 'optical_scenario_composition.csv', dtype={'tag': str})
    key = ['clearance', 'tag', 'kernel', 'endmember', 'shadow', 'age']
    fam = {}
    for _, r in comp.iterrows():
        k = tuple(r[c] for c in key)
        fam.setdefault(k, set()).add(r.scenario_set)
    rows = []
    for clearance in ['1m', '1p5m', '2m']:
        root = N.folder(clearance) / 'adjacency/ew_only'
        sp = pd.read_csv(root / 'a4_shadow_params_camera.csv', dtype={'tag': str})
        for tag in TAGS:
            d, cs, clu = base(tag, clearance)
            G = pd.read_parquet(root / f'a4_{tag}_cells.parquet')
            same = {k: E.A1.per_cell(d, 'same', k)[0] for k in ['red', 'nir', 'ndvi']}
            allc = {k: E.A1.per_cell(d, 'all', k)[0] for k in ['red', 'nir']}
            T = pd.DataFrame(dict(cell=d['tidx'], age=d['age'], grp=clu, same_red_ic=same['red']['ic'], same_nir_ic=same['nir']['ic'],
                                  same_ndvi_dadj=d['valid']['same'] * 0 + np.where(d['valid']['same'], same['ndvi']['it'] - same['ndvi']['ic'] - same['ndvi']['dpre'], np.nan),
                                  all_red_ic=allc['red']['ic'], all_nir_ic=allc['nir']['ic'], cross_ndvi_dadj=np.where(cs.valid_cross, cs.cross_ndvi_dadj, np.nan)))
            D = T.merge(G, on='cell').reset_index(drop=True)
            sr = D.camera.map(sp[sp.tag == tag].set_index('camera').s_ratio)
            variants = sorted(set(c.split('_')[1] for c in D.columns if c.startswith('w_') and c.endswith('_red_core')))
            med = None
            for design, gr, gn, obs_col in [('same', 'same_red_ic', 'same_nir_ic', 'same_ndvi_dadj'), ('primary', 'all_red_ic', 'all_nir_ic', 'cross_ndvi_dadj')]:
                g_r, g_n = D[gr], D[gn]
                def pend(var, pk):
                    out = {}
                    for b, gk, pc in (('red', g_r, 'pr'), ('nir', g_n, 'pn')):
                        wi = D[f'w_{var}_{b}_{pk}']
                        p = (D[f'{pc}_{pk}'] - gk * (1 - wi)) / wi
                        p = p.where(D[f'n_{pk}'] >= 10)
                        out[b] = p.fillna(p.groupby(D.camera).transform('median'))
                    return out
                for var in variants:
                    ends = {pk: pend(var, pk) for pk in ('p10', 'p05', 'p15')}
                    b0 = ends['p10']
                    for nm, fr, fn in (('x0.8', .8, .8), ('x1.2', 1.2, 1.2), ('r+n−', 1.2, .8), ('r−n+', .8, 1.2)):
                        ends[nm] = {'red': b0['red'] * fr, 'nir': b0['nir'] * fn}
                    for en, p in ends.items():
                        rr = g_r + (p['red'] - g_r) * D[f'w_{var}_red_core']
                        nn = g_n + (p['nir'] - g_n) * D[f'w_{var}_nir_core']
                        art = (nn - rr) / (nn + rr) - (g_n - g_r) / (g_n + g_r)
                        for an, lo, hi in [('≥6', 6, 99), ('0–2', 0, 2), ('3–5', 3, 5)]:
                            s = (D.age >= lo) & (D.age <= hi) & np.isfinite(art) & np.isfinite(D[obs_col])
                            rows.append(dict(design=design, clearance=clearance, tag=tag, kernel=var, endmember=en, shadow='none', age=an,
                                             art=float(art[s].mean()), obs=float(D[obs_col][s].mean()), n=int(s.sum())))
                    if f'ws_{var}_red_core' in D.columns:
                        for mode, fnir in (('shadow', 1.0), ('shadow_blue', .85)):
                            s_r, s_n = sr * g_r, sr * fnir * g_n
                            pp = {}
                            for b, gk, pc, sv in (('red', g_r, 'pr', s_r), ('nir', g_n, 'pn', s_n)):
                                wb = D[f'wb_{var}_{b}_p10ns']; ws = D[f'ws_{var}_{b}_p10ns']
                                p = (D[f'{pc}_p10ns'] - gk * (1 - wb - ws) - sv * ws) / wb
                                p = p.where(D.n_p10ns >= 10)
                                pp[b] = p.fillna(p.groupby(D.camera).transform('median'))
                            rr = g_r + (pp['red'] - g_r) * D[f'wb_{var}_red_core'] + (s_r - g_r) * D[f'ws_{var}_red_core']
                            nn = g_n + (pp['nir'] - g_n) * D[f'wb_{var}_nir_core'] + (s_n - g_n) * D[f'ws_{var}_nir_core']
                            art = (nn - rr) / (nn + rr) - (g_n - g_r) / (g_n + g_r)
                            for an, lo, hi in [('≥6', 6, 99), ('0–2', 0, 2), ('3–5', 3, 5)]:
                                s = (D.age >= lo) & (D.age <= hi) & np.isfinite(art) & np.isfinite(D[obs_col])
                                rows.append(dict(design=design, clearance=clearance, tag=tag, kernel=var, endmember='p10ns', shadow=mode, age=an,
                                                 art=float(art[s].mean()), obs=float(D[obs_col][s].mean()), n=int(s.sum())))
            print('OPTICS', clearance, tag, len(D), flush=True)
    R_ = pd.DataFrame(rows)
    R_['sets'] = [';'.join(sorted(fam.get((r.clearance, r.tag, r.kernel, r.endmember, r.shadow, r.age), set()))) for r in R_.itertuples()]
    R_.to_csv(OUT / 'optics.csv', index=False)
    # validation: the 'same' design reproduces the frozen composition artefacts
    chk = R_[R_.design == 'same'].merge(comp[comp.scenario_set == 'all_common'][key + ['art', 'obs']], on=key, suffixes=('', '_frozen'))
    err = float((chk.art - chk.art_frozen).abs().max()); oerr = float((chk.obs - chk.obs_frozen).abs().max())
    print('OPTICS_VALIDATION matched rows', len(chk), 'max |art diff|', err, 'max |obs diff|', oerr, flush=True)
    (OUT / 'optics_validation.json').write_text(json.dumps(dict(rows=len(chk), max_art_diff=err, max_obs_diff=oerr), indent=2))
