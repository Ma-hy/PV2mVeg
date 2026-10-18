"""E3 0.5 m panchromatic texture (after the scope change: Beijing-3 L4 only, east-west fixed arrays, BJ3A1, three acquisitions).

Hypothesis: sparse clumped vegetation / litter patches raise the sub-metre brightness variation; uniformly darker or wetter soil changes only the mean, not the relative variation.
Sampling windows (identical shapes for inter-row ground and controls):
  pix  = the 4×4 0.5 m sub-pixels of one 2 m pixel: cv16 = SD/mean; hp = SD of (panchromatic − 3×3 mean)/mean (high-pass <=1.5 m)
  strip= a 1 m × 6 m strip centred on the pixel along the east-west (row) direction (the two middle sub-pixel rows × the 12 columns of this pixel and its east and west neighbours),
         SD/mean of the means of its six 1 m blocks (local contrast at 1–3 m); the two east-west neighbouring 2 m pixels must be of the class of this pixel (inter-row core or pure control ground),
         so inter-row and control strips have the same shape and orientation, and inter-row strips stay within the core band (>=2 m from panels, without panel edges and shadow).
Treated units: units of the a1 cells table with valid_same and orient = 1 (east-west); pixels = inter-row core (class 3, QA passed, source scene = modal scene of the unit, centre inside the unit).
Controls: the same-scene matched controls of these units (ctrl_scene = ir_scene, ctrl_n >=2000, the same rule as a1), pixels = class 5 with dist_to_array >=300 m.
Discrete structures (mounting posts, cable trenches, combiner boxes): the main statistic is the median over the pixels of a unit (robust to a few outliers), the mean is also kept; pixels whose brightness departs >40% from the unit median are removed.
A sample of <=150 pixels per unit (mean and the three textures) is also kept for e7 ("brightness matching", to exclude the bias from the larger relative noise of dark surfaces).
Output: data/eco/adjacency/ew_only/e3_<tag>_{cells,pixels,pairs}.parquet
Usage: python eco/adjacency/e3_pan_texture.py <tag> [--workers 8]
"""
# Mirrored from eco/adjacency/e3_pan_texture.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
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
from scipy.ndimage import uniform_filter  # noqa: E402

EW = os.path.join(A.OUT, 'ew_only')
os.makedirs(EW, exist_ok=True)
CELLS_CSV = {'06': 'a1_06_det_thr50_cells.csv', '08': 'a1_08_det_thr50_BJ3A1_cells.csv', '2026_08': 'a1_2026_08_det_thr50_cells.csv'}
NSAMP = 150


def one(tag, cell, x, y, scene, kind, scen, tg, crs, tk, rok, cok, src_path, rng):
    tr_u2g = Transformer.from_crs('EPSG:32647', crs, always_xy=True)
    tr_g2u = Transformer.from_crs(crs, 'EPSG:32647', always_xy=True)
    xs = x + np.array([-75, 75, 75, -75]); ys = y + np.array([75, 75, -75, -75])
    gx, gy = tr_u2g.transform(xs, ys)
    c0 = int(np.floor((min(gx) - tg.c) / 2)) - 3; r0 = int(np.floor((tg.f - max(gy)) / 2)) - 3
    c1 = int(np.ceil((max(gx) - tg.c) / 2)) + 3; r1 = int(np.ceil((tg.f - min(gy)) / 2)) + 3
    h, w = r1 - r0, c1 - c0
    cls = A.read_win(os.path.join(A.BJ, tag, 'class_2m.tif'), r0, c0, h, w)
    qa = A.read_win(os.path.join(A.BJ, tag, 'qa_2m.tif'), r0, c0, h, w)
    sc = A.read_win(os.path.join(A.BJ, tag, 'scene_2m.tif'), r0, c0, h, w)
    jj, ii = np.meshgrid(np.arange(w), np.arange(h))
    ux, uy = tr_g2u.transform((tg.c + (c0 + jj + 0.5) * 2.0).ravel(), (tg.f - (r0 + ii + 0.5) * 2.0).ravel())
    inc = ((np.abs(np.asarray(ux) - x) < 75) & (np.abs(np.asarray(uy) - y) < 75)).reshape(h, w)
    ok = ((qa & A.QA_EX) == 0) & (sc == scene)
    if kind == 'treated':
        cand = ok & (cls == 3)
    else:
        dta = A.read_win(os.path.join(A.G2, tag, 'dist_to_array_2m.tif'), r0, c0, h, w, fill=0)
        cand = ok & (cls == 5) & (dta >= 300)
    pp = A.pan_path(scen[scene - 1]['path'])
    with rasterio.open(pp) as d:
        tp = d.transform
    pro, pco = A.lattice_offset(tk, tp, f=1)
    fr0, fc0 = 4 * r0 + rok, 4 * c0 + cok
    P = A.read_win(pp, fr0 + pro, fc0 + pco, 4 * h, 4 * w, fill=0).astype(np.float64)
    src = A.read_win(src_path, fr0, fc0, 4 * h, 4 * w, fill=255)
    hpf = P - uniform_filter(P, 3)
    blk = lambda Z: Z.reshape(h, 4, w, 4).transpose(0, 2, 1, 3).reshape(h, w, 16)  # noqa: E731
    Pb, Hb, Sb = blk(P), blk(hpf), blk(src)
    good = cand & (Sb == scene - 1).all(2) & (Pb > 0).all(2)
    gl = np.zeros_like(good); gl[:, 1:] = good[:, :-1]
    gr = np.zeros_like(good); gr[:, :-1] = good[:, 1:]
    good = good & gl & gr
    good[0, :] = good[-1, :] = False; good[:, 0] = good[:, -1] = False
    sel = good & inc
    if sel.sum() < 20:
        return None, None
    I, J = np.nonzero(sel)
    mu = Pb[I, J].mean(1)
    med = np.median(mu)
    keep = np.abs(mu / med - 1) <= 0.4
    I, J, mu = I[keep], J[keep], mu[keep]
    if len(I) < 20:
        return None, None
    cv = Pb[I, J].std(1) / mu
    hp = Hb[I, J].std(1) / mu
    rows = 4 * I[:, None] + np.array([1, 2])[None, :]
    cols = 4 * (J[:, None] - 1) + np.arange(12)[None, :]
    S = P[rows[:, :, None], cols[:, None, :]]                  # (n, 2, 12)
    bm = S.reshape(len(I), 2, 6, 2).mean((1, 3))               # (n, 6) blocks of 1 m
    st = bm.std(1) / bm.mean(1)
    rec = dict(cell=int(cell), kind=kind, scene=int(scene), n=int(len(I)), pan_mean=float(mu.mean()),
               cv16=float(np.median(cv)), hp=float(np.median(hp)), strip=float(np.median(st)),
               cv16_mean=float(cv.mean()), hp_mean=float(hp.mean()), strip_mean=float(st.mean()))
    k = rng.choice(len(I), min(NSAMP, len(I)), replace=False)
    pix = pd.DataFrame(dict(cell=int(cell), kind=kind, scene=int(scene), pan=mu[k], cv16=cv[k], hp=hp[k], strip=st[k]))
    return rec, pix


def work(tag, jobs, scen):
    c = A.pcfg(tag)
    with rasterio.open(os.path.join(A.BJ, tag, 'class_2m.tif')) as d:
        tg, crs = d.transform, d.crs
    with rasterio.open(c['mask']) as d:
        tk = d.transform
    rok, cok = A.lattice_offset(tg, tk)
    rng = np.random.default_rng(len(jobs) * 7 + (jobs[0][0] if jobs else 0))
    out, pix = [], []
    for j in jobs:
        try:
            r, p = one(tag, *j, scen, tg, crs, tk, rok, cok, c['source'], rng)
        except Exception as ex:  # noqa: BLE001
            A.log(f'{tag} {j[0]} 失败 {ex}', 'log_e3.txt'); r = None
        if r:
            out.append(r); pix.append(p)
    return out, pix
