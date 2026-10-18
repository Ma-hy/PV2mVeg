"""P1 0.5 m panchromatic inter-row profiles: extraction per treated unit / per control unit (once per acquisition).

Treated units (main a1 rule, see pp.py): within the 150 m unit window, with the 0.5 m mask of the acquisition and the panchromatic band of the modal scene of the unit (ir_scene, = the scene of a1):
  - row band = panel sub-pixels of the mask (the east-west mask includes the shadow strip along the north edge of the modules, treated as module by decision of the authors);
  - distance to the row band d = Euclidean distance from the sub-pixel centre to the nearest row-band sub-pixel edge (EDT, (dist_px − 0.5) × 0.5 m); 0.5 m bins, 12 bins over 0–6 m;
  - inter-row sub-pixels = not panel, source scene = ir_scene, panchromatic >0, parent 2 m pixel QA passed with centre inside the unit, parent 2 m class ∈ {1 panel, 2 panel edge, 3 inter-row core},
    with panels met both northward and southward (across the east-west rows) within one local row pitch (without a no-image area in between),
    clearance gap = number of non-panel sub-pixels between the two panels × 0.5 m <= row pitch pitch_m of the unit (the inter-row criterion of the 2 m product, applied per 0.5 m column);
  - side: southern panel nearer -> N (the pixel is on the north side of the row to its south, next to the north edge of the modules / shadow); northern panel nearer -> S (south side of a row); equal -> C;
  - row middle = |northward distance − southward distance| <= 1 sub-pixel (0.5 m);
  - 0.5 m pixels of the 2 m core pixels = parent 2 m class 3 (inter-row core, >=2 m from panels), QA passed, centre inside the unit: those of its 16 sub-pixels with source scene = ir_scene;
  - geometry: for each inter-row segment one centre sub-pixel, recording the clearance (sub-pixels); for each row-band segment its top sub-pixel, recording the north-south row-band width (sub-pixels).
Control units (same-scene matched controls, (unit, scene) deduplicated): parent 2 m class 5 (pure ground outside arrays), >=300 m from the array area, QA passed, source scene = that scene, centre inside the unit;
  sub-pixels not panel, source scene = that scene, panchromatic >0. The rule of the a1 control pixels (out_300_1000 + out_gt1000).
Statistics: sum and count of panchromatic DN (main statistic = mean); also the median of the 2 m pixel means (all 16 sub-pixels valid) and the "trimmed mean within 40% of the median" (the rule of the e6 texture analysis, for comparison).
Panchromatic DN are not radiometrically calibrated; ratios are taken within one scene only (treated / same-scene control).
Output: data/eco/panprof/p1_<tag>_{treated,control,prof,geom}.parquet
Usage: python eco/panprof/p1_extract.py <tag> [--workers 8] [--limit n]
"""
# Mirrored from eco/panprof/p1_extract.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pp  # noqa: E402
A = pp.A
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from pyproj import Transformer  # noqa: E402
from scipy.ndimage import distance_transform_edt  # noqa: E402

MG = 10          # margin of the 2 m window (20 m)
GMAX = 40        # largest recorded clearance / row-band width (sub-pixels)


def up(a):
    return np.repeat(np.repeat(a, 4, 0), 4, 1)


def scan_ns(pan, nod):
    """Per column: row distance from each sub-pixel northward / southward to the nearest panel sub-pixel (dn, ds), and whether a no-image area (nodata) lies in between."""
    H, W = pan.shape
    r = np.arange(H)[:, None] * np.ones((1, W), np.int64)
    big = 10 ** 6
    pu = np.maximum.accumulate(np.where(pan, r, -big), axis=0)
    nu = np.maximum.accumulate(np.where(nod, r, -big), axis=0)
    pd_ = np.minimum.accumulate(np.where(pan, r, big)[::-1], axis=0)[::-1]
    nd_ = np.minimum.accumulate(np.where(nod, r, big)[::-1], axis=0)[::-1]
    dn = r - pu
    ds = pd_ - r
    okn = (pu >= 0) & (pu > nu)
    oks = (pd_ < big) & (pd_ < nd_)
    return dn, ds, okn & oks


def window(tg, tr_u2g, x, y):
    xs = x + np.array([-75, 75, 75, -75]); ys = y + np.array([75, 75, -75, -75])
    gx, gy = tr_u2g.transform(xs, ys)
    c0 = int(np.floor((min(gx) - tg.c) / 2)) - MG; r0 = int(np.floor((tg.f - max(gy)) / 2)) - MG
    c1 = int(np.ceil((max(gx) - tg.c) / 2)) + MG; r1 = int(np.ceil((tg.f - min(gy)) / 2)) + MG
    return r0, c0, r1 - r0, c1 - c0


def pix_stats(Pb, valid16):
    """Median and trimmed mean (within 40% of the median) of the 2 m pixel means (all 16 sub-pixels valid), and the pixel count."""
    if valid16.sum() < 5:
        return np.nan, np.nan, int(valid16.sum())
    mu = Pb[valid16].mean(1)
    med = float(np.median(mu))
    k = np.abs(mu / med - 1) <= 0.4
    return med, float(mu[k].mean()), int(valid16.sum())


def geom_of(cell, pan, nod, irr, dn, ds, gap, inc):
    """Geometry: for each inter-row segment (per 0.5 m column) one centre sub-pixel with its clearance (sub-pixels); for each row-band segment its northern sub-pixel with the north-south row-band width (sub-pixels).
    Only segments whose parent 2 m pixel centre is inside the unit."""
    inc_f = up(inc)
    seg = irr & ((dn - ds) == 0) | (irr & ((dn - ds) == 1))
    gv = gap[seg & inc_f]
    gh = np.bincount(np.clip(gv, 0, GMAX), minlength=GMAX + 1)
    # row band: top of a panel segment (north end, non-panel and non-no-image above) -> length of the continuous panel southward; the lower end must be non-panel and non-no-image
    H_, W_ = pan.shape
    r = np.arange(H_)[:, None] * np.ones((1, W_), np.int64)
    big = 10 ** 6
    nxt = np.minimum.accumulate(np.where(~pan, r, big)[::-1], axis=0)[::-1]
    top = pan.copy(); top[0, :] = False
    top[1:, :] &= ~pan[:-1, :] & ~nod[:-1, :]
    L = nxt - r
    endok = np.zeros_like(top)
    rows_, cols_ = np.nonzero(top)
    ends = nxt[rows_, cols_]
    v_ = ends < H_
    ok_end = np.zeros(len(rows_), bool)
    ok_end[v_] = ~nod[ends[v_], cols_[v_]]
    endok[rows_, cols_] = ok_end
    bv = L[top & endok & inc_f]
    bh = np.bincount(np.clip(bv, 0, GMAX), minlength=GMAX + 1)
    geom = dict(cell=int(cell), **{f'g{k}': int(v) for k, v in enumerate(gh)}, **{f'b{k}': int(v) for k, v in enumerate(bh)})
    return geom


def one(tag, job, ctx):
    cell, x, y, scene, kind, pitch = job
    tg, crs, tk, rok, cok, cfgc, scen, tr_u2g, tr_g2u = ctx
    r0, c0, h, w = window(tg, tr_u2g, x, y)
    rd = lambda p, f=0: A.read_win(p, r0, c0, h, w, fill=f)  # noqa: E731
    cls = rd(os.path.join(A.BJ, tag, 'class_2m.tif'))
    qa = rd(os.path.join(A.BJ, tag, 'qa_2m.tif'))
    sc = rd(os.path.join(A.BJ, tag, 'scene_2m.tif'))
    jj, ii = np.meshgrid(np.arange(w), np.arange(h))
    ux, uy = tr_g2u.transform((tg.c + (c0 + jj + 0.5) * 2.0).ravel(), (tg.f - (r0 + ii + 0.5) * 2.0).ravel())
    inc = ((np.abs(np.asarray(ux) - x) < 75) & (np.abs(np.asarray(uy) - y) < 75)).reshape(h, w)
    ok2 = ((qa & A.QA_EX) == 0) & (sc > 0) & inc
    fr0, fc0 = 4 * r0 + rok, 4 * c0 + cok
    Mk = A.read_win(cfgc['mask'], fr0, fc0, 4 * h, 4 * w, fill=255)
    pan = Mk == 1
    nod = Mk == 255
    if kind == 'geom':      # geometry only (all 4759 east-west treated units; no panchromatic read, any scene)
        dn, ds, both = scan_ns(pan, nod)
        gap = dn + ds - 1
        irr = up(ok2 & np.isin(cls, (1, 2, 3))) & ~pan & ~nod & both & (gap * 0.5 <= pitch)
        return None, [], geom_of(cell, pan, nod, irr, dn, ds, gap, inc)
    src = A.read_win(cfgc['source'], fr0, fc0, 4 * h, 4 * w, fill=255)
    pp_ = A.pan_path(scen[scene - 1]['path'])
    with rasterio.open(pp_) as d:
        tp = d.transform
    pro, pco = A.lattice_offset(tk, tp, f=1)
    P = A.read_win(pp_, fr0 + pro, fc0 + pco, 4 * h, 4 * w, fill=0).astype(np.float64)
    good = (src == scene - 1) & (P > 0) & ~pan & ~nod
    blk = lambda Z: Z.reshape(h, 4, w, 4).transpose(0, 2, 1, 3).reshape(h, w, 16)  # noqa: E731
    if kind == 'control':
        dta = rd(os.path.join(A.G2, tag, 'dist_to_array_2m.tif'), 0)
        par = ok2 & (cls == 5) & (dta >= 300) & (sc == scene)
        sel = up(par) & good
        n = int(sel.sum())
        if n < 1000:
            return None, [], None
        Pb, Gb = blk(P), blk(good)
        med, tm, n2 = pix_stats(Pb, par & Gb.all(2))
        return dict(ctrl=int(cell), scene=int(scene), n=n, s=float(P[sel].sum()), ss=float((P[sel] ** 2).sum()),
                    med2m=med, tmean2m=tm, n2m=n2), [], None
    # ---------------- treated units
    dn, ds, both = scan_ns(pan, nod)
    gap = dn + ds - 1                                  # non-panel sub-pixels between the panels on both sides
    gap_ok = both & (gap * 0.5 <= pitch)
    par_ir = ok2 & np.isin(cls, (1, 2, 3))
    irr = up(par_ir) & good & gap_ok
    D = (distance_transform_edt(~pan) - 0.5) * 0.5
    b = np.floor(D / pp.DBIN).astype(int)
    side = np.where(ds < dn, 0, np.where(dn < ds, 1, 2))   # 0 N, 1 S, 2 C
    prof = []
    for si in range(3):
        m_ = irr & (side == si) & (b >= 0) & (b < pp.NBIN)
        if not m_.any():
            continue
        cnt = np.bincount(b[m_], minlength=pp.NBIN)
        sm = np.bincount(b[m_], P[m_], minlength=pp.NBIN)
        for bi in np.nonzero(cnt)[0]:
            prof.append(dict(cell=int(cell), side=pp.SIDES[si], bin=int(bi), n=int(cnt[bi]), s=float(sm[bi])))
    cen = irr & (np.abs(dn - ds) <= 1)
    par_core = ok2 & (cls == 3)
    core = up(par_core) & good
    Pb, Gb = blk(P), blk(good)
    med, tm, n2 = pix_stats(Pb, par_core & Gb.all(2))
    rec = dict(cell=int(cell), scene=int(scene), n_ir=int(irr.sum()), s_ir=float(P[irr].sum()),
               n_cen=int(cen.sum()), s_cen=float(P[cen].sum()),
               n_core=int(core.sum()), s_core=float(P[core].sum()), n_core2m=int(par_core.sum()),
               core_med2m=med, core_tmean2m=tm, core_n2m=n2,
               n_ir_core=int((irr & core).sum()))
    # inter-row sub-pixels >=2 m from panels (d >= 2.0): the "core" at sub-pixel level, for comparison with the 2 m core
    m2 = irr & (D >= 2.0)
    rec['n_d2'] = int(m2.sum()); rec['s_d2'] = float(P[m2].sum())
    return rec, prof, geom_of(cell, pan, nod, irr, dn, ds, gap, inc)


def work(tag, jobs, scen):
    c = A.pcfg(tag)
    with rasterio.open(os.path.join(A.BJ, tag, 'class_2m.tif')) as d:
        tg, crs = d.transform, d.crs
    with rasterio.open(c['mask']) as d:
        tk = d.transform
    rok, cok = A.lattice_offset(tg, tk)
    tr_u2g = Transformer.from_crs('EPSG:32647', crs, always_xy=True)
    tr_g2u = Transformer.from_crs(crs, 'EPSG:32647', always_xy=True)
    ctx = (tg, crs, tk, rok, cok, c, scen, tr_u2g, tr_g2u)
    recs, ctrls, profs, geoms = [], [], [], []
    for j in jobs:
        try:
            r, p, g = one(tag, j, ctx)
        except Exception as ex:  # noqa: BLE001
            pp.log(f'{tag} {j[4]} {j[0]} 失败 {ex!r}')
            continue
        if j[4] == 'geom':
            if g is not None:
                geoms.append(g)
            continue
        if r is None:
            continue
        if j[4] == 'control':
            ctrls.append(r)
        else:
            recs.append(r); profs += p; geoms.append(g)
    return recs, ctrls, profs, geoms
