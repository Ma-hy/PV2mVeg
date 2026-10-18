"""A2 forward geometry per treated unit: within the 150 m unit window of ts, the 0.5 m mask of the month and the
estimated per-band kernels give, for each 2 m pixel, the "kernel-weighted panel fraction" w_b (= the weight of the
panels in the observed value of the pixel); these are aggregated (for a3: expected artefact and correction), with the
profiles and rule variants needed by C/D/E.

Per unit:
  core       = class 3 (inter-row core), QA passed, pixel centre inside the unit (the rule of the 2 m product -> 30 m
               aggregation; n is checked against ir_n of a1)
  w̄_b        = mean w_b of the core pixels (per kernel scheme x band); also split by d 2-3 m / >=3 m and by the side
               (N/S/E/W) of the nearest panel
  panel end-member = observed mean p_obs and its w̄_int (per kernel scheme) of the 2 m pixels inside the unit that are
               fully panel with all 4x4 sub-pixels >=1 m from the ground (row centre)
  D profile  = inter-row pixels (class 3, and class-2 panel-edge pixels with gap <= row pitch) by d (0.5 m bins,
               0-6 m) x direction of the nearest panel: Σred, Σnir, Σw
  C profile  = same-scene 0.5 m panchromatic by signed sub-pixel distance to the edge (0.25 m bins, -3...+6 m,
               negative = inside the panel) x direction: ΣPAN
  E variants = mask eroded/dilated by 0.5 m (4-neighbourhood, 1 sub-pixel), 2 m centre distance recomputed; in the
               inter-row zone (original class 3, or class 2 by the gap rule) with d' >= 2.0/2.5/3.0 m: Σred, Σnir, n, Σw
The kernel schemes come from data/eco/adjacency/k2_kernels.csv (a2 only reads them, it does not select them).
Usage: python eco/adjacency/a2_cells.py <tag> [--limit n] [--workers 8]
Output: data/eco/adjacency/a2_<tag>_{cells.parquet, prof.parquet, pan.parquet, evar.parquet}
"""
# Mirrored from eco/adjacency/a2_cells.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import adjc as A  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from pyproj import Transformer  # noqa: E402
from scipy.ndimage import distance_transform_edt, binary_erosion, binary_dilation  # noqa: E402
from scipy.signal import fftconvolve  # noqa: E402

CELLS_CSV = {'06': 'a1_06_det_thr50_cells.csv', '08': 'a1_08_det_thr50_BJ3A1_cells.csv', '2026_08': 'a1_2026_08_det_thr50_cells.csv'}
MARG2 = 14            # margin in 2 m pixels (28 m)
DB = np.arange(0, 6.01, 0.5)
PB = np.arange(-3, 6.001, 0.25)
THRS = (2.0, 2.5, 3.0)
CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)


def camera_of(sid, sensor):
    if sensor == 'BJ3A1':
        return 'MS1' if '_MS1_' in sid else 'MS2'
    return sensor


def load_kernels(tag):
    K = pd.read_csv(os.path.join(A.OUT, 'k2_kernels.csv'), dtype={'tag': str})
    K = K[(K.tag == tag) & (K.variant != 'npk')]
    out = {}
    for (var, cam, band), r in K.set_index(['variant', 'camera', 'band']).iterrows():
        out[(var, cam, band)] = A.param_kernel([r.dx, r.dy, r.sx, r.sy, r.al, r.s2])
    pn = os.path.join(A.OUT, 'k2_np_kernels.npz')
    if os.path.exists(pn):
        Z = np.load(pn)
        for k in Z.files:
            t, cam, band = k.split('|')
            if t != tag:
                continue
            if band == 'nir#r':
                out[('npkr', cam, 'nir')] = Z[k]
            else:
                out[('npk', cam, band)] = Z[k]
                if band == 'red':
                    out[('npkr', cam, 'red')] = Z[k]
    return out


def sample_w(Mf, k, ro, co, h, w):
    """Mf on the fine grid (row/column 0 of the 2 m window (0,0) are fine-grid row ro, column co); returns w of shape (h,w) = Σ k[a,b] Mf[4I+ro+a−KM, 4J+co+b−KM]."""
    F = fftconvolve(Mf, k[::-1, ::-1], mode='full')
    o = A.KM + 3
    return F[ro + o:ro + o + 4 * h:4, co + o:co + o + 4 * w:4][:h, :w]


def centre_dist(pan, ro, co, h, w):
    """d exactly as in the product: minimum distance from the 2×2 0.5 m sub-pixels around the 2 m cell centre to the nearest panel sub-pixel centre, minus 0.25 m (grid2m algorithm).
    Note: the mask is aligned with whole 2 m cells; for row edges parallel to the grid, d = distance from the cell centre to the mask edge − 0.25 m (values 0.75, 1.25, 1.75, 2.25...)."""
    if not pan.any():
        return np.full((h, w), 200.0)
    D = distance_transform_edt(~pan) * 0.5
    best = np.full((h, w), np.inf)
    for dr in (1, 2):
        for dc in (1, 2):
            best = np.minimum(best, D[ro + dr::4, co + dc::4][:h, :w])
    return best - 0.25


def process(tag, rows, kernels, scen):
    c = A.pcfg(tag)
    with rasterio.open(os.path.join(A.BJ, tag, 'class_2m.tif')) as d:
        tg, crs = d.transform, d.crs
    with rasterio.open(c['mask']) as d:
        tk = d.transform
    rok, cok = A.lattice_offset(tg, tk)
    tr_u2g = Transformer.from_crs('EPSG:32647', crs, always_xy=True)
    tr_g2u = Transformer.from_crs(crs, 'EPSG:32647', always_xy=True)
    variants = sorted(set(k[0] for k in kernels))
    out_c, out_p, out_pan, out_e = [], [], [], []
    pan_cache = {}
    for _, cr in rows.iterrows():
        xs = cr.x + np.array([-75, 75, 75, -75]); ys = cr.y + np.array([75, 75, -75, -75])
        gx, gy = tr_u2g.transform(xs, ys)
        c0 = int(np.floor((min(gx) - tg.c) / 2)) - MARG2
        r0 = int(np.floor((tg.f - max(gy)) / 2)) - MARG2
        c1 = int(np.ceil((max(gx) - tg.c) / 2)) + MARG2
        r1 = int(np.ceil((tg.f - min(gy)) / 2)) + MARG2
        h, w = r1 - r0, c1 - c0
        rd = lambda p, f=0: A.read_win(p, r0, c0, h, w, fill=f)  # noqa: E731
        cls = rd(os.path.join(A.BJ, tag, 'class_2m.tif'))
        qa = rd(os.path.join(A.BJ, tag, 'qa_2m.tif'))
        sc = rd(os.path.join(A.BJ, tag, 'scene_2m.tif'))
        red = rd(os.path.join(A.BJ, tag, 'red_2m.tif'), -32768).astype(np.float64) / 1e4
        nir = rd(os.path.join(A.BJ, tag, 'nir_2m.tif'), -32768).astype(np.float64) / 1e4
        dctr = rd(os.path.join(A.G2, tag, 'dist_ctr_2m.tif'), np.nan)
        dirc = rd(os.path.join(A.G2, tag, 'dir_2m.tif'), 255)
        gns = rd(os.path.join(A.G2, tag, 'gap_ns_2m.tif'), np.nan)
        gew = rd(os.path.join(A.G2, tag, 'gap_ew_2m.tif'), np.nan)
        # fine-grid mask (2 m window (0,0) at fine-grid (4·MF/4 ...))
        MF = 4 * 12   # fine-grid margin (sub-pixels) = 24 m
        fr0, fc0 = 4 * r0 + rok - MF, 4 * c0 + cok - MF
        Mk = A.read_win(c['mask'], fr0, fc0, 4 * h + 2 * MF, 4 * w + 2 * MF, fill=255)
        pan = Mk == 1
        Mf = pan.astype(np.float64)
        ro = co = MF
        # pixel centres -> UTM, inside the unit
        jj, ii = np.meshgrid(np.arange(w), np.arange(h))
        px = tg.c + (c0 + jj + 0.5) * 2.0
        py = tg.f - (r0 + ii + 0.5) * 2.0
        ux, uy = tr_g2u.transform(px.ravel(), py.ravel())
        inc = ((np.abs(np.asarray(ux) - cr.x) < 75) & (np.abs(np.asarray(uy) - cr.y) < 75)).reshape(h, w)
        ok = ((qa & A.QA_EX) == 0) & (sc > 0) & (red > -1) & (nir > -1)
        d = dctr - getattr(A, 'CENTER_DISTANCE_OFFSET', 0.25)
        core = (cls == 3) & ok & inc
        gap = np.where(cr.orient == 2, gew, gns)
        irow = ((cls == 3) | ((cls == 2) & (gap <= cr.pitch_m))) & ok & inc
        # scene / camera
        sid = scen[int(cr.ir_scene) - 1]['sid']; sens = scen[int(cr.ir_scene) - 1]['sensor']
        cam = camera_of(sid, sens)
        rec = dict(cell=int(cr.cell), n_core=int(core.sum()), ir_n=float(cr.ir_n), camera=cam,
                   sr_core=float(red[core].sum()), sn_core=float(nir[core].sum()))
        for nm, sel in (('d23', core & (d < 3)), ('d3p', core & (d >= 3)),
                        ('sN', core & (dirc == 1)), ('sS', core & (dirc == 3)), ('sE', core & (dirc == 2)), ('sW', core & (dirc == 4))):
            rec[f'n_{nm}'] = int(sel.sum()); rec[f'sr_{nm}'] = float(red[sel].sum()); rec[f'sn_{nm}'] = float(nir[sel].sum())
        # pixels inside panels
        din = distance_transform_edt(pan) * 0.5   # from panel sub-pixel centres to the nearest non-panel sub-pixel centre
        blk_in = din[ro:ro + 4 * h, co:co + 4 * w].reshape(h, 4, w, 4).min((1, 3))
        cov = pan[ro:ro + 4 * h, co:co + 4 * w].reshape(h, 4, w, 4).mean((1, 3))
        pint = (cov == 1) & (blk_in >= 1.0) & ok & inc
        rec['n_pint'] = int(pint.sum())
        rec['pr_obs'] = float(red[pint].mean()) if pint.sum() >= 10 else np.nan
        rec['pn_obs'] = float(nir[pint].mean()) if pint.sum() >= 10 else np.nan
        # kernels
        W = {}
        for var in variants:
            for b in ('red', 'nir'):
                kk = kernels.get((var, cam, b), kernels.get((var, 'all', b)))
                if kk is None:
                    continue
                Wb = sample_w(Mf, kk, ro, co, h, w)
                W[(var, b)] = Wb
                for nm, sel in (('core', core), ('d23', core & (d < 3)), ('d3p', core & (d >= 3)), ('sN', core & (dirc == 1)),
                                ('sS', core & (dirc == 3)), ('sE', core & (dirc == 2)), ('sW', core & (dirc == 4)), ('pint', pint)):
                    rec[f'w_{var}_{b}_{nm}'] = float(Wb[sel].mean()) if sel.any() else np.nan
        out_c.append(rec)
        # D profile (main scheme)
        mainv = 'main' if ('main', 'red') in W else variants[0]
        dbin = np.digitize(d, DB) - 1
        for side in (1, 2, 3, 4):
            for bi in range(len(DB) - 1):
                sel = irow & (dirc == side) & (dbin == bi)
                n = int(sel.sum())
                if n == 0:
                    continue
                out_p.append(dict(cell=int(cr.cell), side=side, dlo=DB[bi], n=n, sr=float(red[sel].sum()), sn=float(nir[sel].sum()),
                                  swr=float(W[(mainv, 'red')][sel].sum()), swn=float(W[(mainv, 'nir')][sel].sum()),
                                  sndvi=float(A.ndvi(red[sel], nir[sel]).sum())))
        # E variants
        for mv, pm in (('m0', pan), ('ero', binary_erosion(pan, CROSS)), ('dil', binary_dilation(pan, CROSS))):
            dv = centre_dist(pm, ro, co, h, w)
            covv = pm[ro:ro + 4 * h, co:co + 4 * w].reshape(h, 4, w, 4).any((1, 3))
            if mv == 'm0':
                rec['d_check_max'] = float(np.nanmax(np.abs(dv - d)[core])) if core.any() else np.nan
            for thr in THRS:
                sel = irow & ~covv & (dv >= thr)
                e = dict(cell=int(cr.cell), mvar=mv, thr=thr, n=int(sel.sum()), sr=float(red[sel].sum()), sn=float(nir[sel].sum()))
                for b in ('red', 'nir'):
                    e[f'sw_{b}'] = float(W[(mainv, b)][sel].sum())
                out_e.append(e)
        # C panchromatic profile (modal scene of the unit)
        try:
            msp = scen[int(cr.ir_scene) - 1]['path']
            pp = A.pan_path(msp)
            if pp not in pan_cache:
                with rasterio.open(pp) as dpn:
                    pan_cache.clear()
                    pan_cache[pp] = dpn.transform
            tp = pan_cache[pp]
            pro, pco = A.lattice_offset(tk, tp, f=1)     # row/column of mask (0,0) on the PAN grid
            P = A.read_win(pp, fr0 + pro, fc0 + pco, Mk.shape[0], Mk.shape[1], fill=0).astype(np.float64)
            src = A.read_win(c['source'], fr0, fc0, Mk.shape[0], Mk.shape[1], fill=255)
            dout, (ir_, ic_) = distance_transform_edt(~pan, return_indices=True)
            sd = np.where(pan, -(din / 0.5 - 0.5) * 0.5, (dout - 0.5) * 0.5)
            rr_, cc_ = np.indices(pan.shape)
            # direction of the nearest panel (inside panels: of the nearest ground)
            _, (jr, jc) = distance_transform_edt(pan, return_indices=True)
            dr = np.where(pan, jr - rr_, ir_ - rr_)
            dcl = np.where(pan, jc - cc_, ic_ - cc_)
            side = np.where(np.abs(dr) >= np.abs(dcl), np.where(dr < 0, 1, 3), np.where(dcl > 0, 2, 4))
            # only the inter-row zone and panels (2 m parent pixel class 1 / inter-row), from this scene, inside the unit
            par = np.zeros(pan.shape, bool)
            par_ = (irow | ((cls == 1) & inc)) & ok
            par[ro:ro + 4 * h, co:co + 4 * w] = np.repeat(np.repeat(par_, 4, 0), 4, 1)
            sel = par & (src == int(cr.ir_scene) - 1) & (P > 0)
            pb = np.digitize(sd, PB) - 1
            for s_ in (1, 2, 3, 4):
                m_ = sel & (side == s_) & (pb >= 0) & (pb < len(PB) - 1)
                if m_.sum() == 0:
                    continue
                cnt = np.bincount(pb[m_], minlength=len(PB) - 1)
                sm = np.bincount(pb[m_], P[m_], minlength=len(PB) - 1)
                for bi in np.nonzero(cnt)[0]:
                    out_pan.append(dict(cell=int(cr.cell), side=s_, dlo=PB[bi], n=int(cnt[bi]), span=float(sm[bi])))
        except Exception as ex:  # noqa: BLE001
            A.log(f'PAN 失败 {tag} {cr.cell}: {ex}', 'log_a2.txt')
    return out_c, out_p, out_pan, out_e
