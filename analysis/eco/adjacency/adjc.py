"""Shared functions of eco/adjacency: paths, scene tables, window reading, kernel (PSF) regression, mask convolution.

Every script starts with `import adjc` (this module adds src/, eco/, eco/bj3ndvi and eco/chrono to the path and imports talatan).
Reads only outputs/, data/raw/ and the folder of the 1.2 m imagery beside the project; writes its results to data/eco/adjacency/.
"""
# Mirrored from eco/adjacency/adjc.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import time

ROOT = str(PROJECT)  # project data root (pv2mveg configuration)
for p in (os.path.join(CODE, 'eco'), os.path.join(CODE, 'eco', 'bj3ndvi'),
          os.path.join(CODE, 'eco', 'chrono')):
    if p not in sys.path:
        sys.path.insert(0, p)
import pv2mveg  # noqa: F401,E402  (PROJ/GDAL data of the rasterio wheel, no WMI)

import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from rasterio.windows import Window  # noqa: E402

OUT = os.path.join(ROOT, 'data', 'eco', 'adjacency')
FIG = os.path.join(OUT, 'figs')
for _d in (OUT, FIG):
    os.makedirs(_d, exist_ok=True)
BJ = os.path.join(ROOT, 'data', 'eco', 'bj3ndvi')
G2 = os.path.join(ROOT, 'data', 'eco', 'grid2m')
CHR = os.path.join(ROOT, 'data', 'eco', 'chrono')
FONT_PATH = FONT_FILE
TAGS = ['06', '08', '2026_08']
LABEL = {'06': '2025-06', '08': '2025-08', '2026_08': '2026-08'}
QA_EX = 1 | 2 | 8 | 32 | 64 | 128

# kernel support: 0.5 m sub-pixel offsets relative to the 2 m pixel centre, k ∈ [−KM, KM+3] -> centre offset 0.5k − 0.75 m (−4.75...+4.75 m)
KM = 8
KN = 2 * KM + 4
OFFS = 0.5 * (np.arange(-KM, KM + 4)) - 0.75          # sub-pixel centres relative to the 2 m pixel centre (east / south positive)


# exclusion lists of anomalous arrays (shared with eco/ts/ew0_cells.py, eco/chrono/cc.py): data/eco/ts/ew/exclude_cells.txt (all acquisitions) ∪ exclude_cells_<tag>.txt (that acquisition only).
# the ew_only scripts (e3, e4, e6, k3, a4, a5) also remove these units from the east-west units (orient = 1); without the files, or with empty files, the behaviour is unchanged.
EXCL_DIR = os.path.join(ROOT, 'data', 'eco', 'ts', 'ew')


def log(msg, fn=None):
    s = time.strftime('%H:%M:%S') + ' ' + str(msg)
    print(s, flush=True)
    if fn:
        with open(os.path.join(OUT, fn), 'a', encoding='utf-8') as f:
            f.write(s + '\n')


def pcfg(tag):
    import cfg
    return cfg.pcfg(tag)


def scenes(tag):
    import cfg
    return cfg.scenes(cfg.pcfg(tag))


def pan_path(ms_path):
    b = os.path.basename(ms_path)
    b2 = b.replace('_MS1_', '_PAN1_').replace('_MS2_', '_PAN2_').replace('_MS_', '_PAN_').replace('_ORTHO_MS', '_ORTHO_PAN')
    p = os.path.join(os.path.dirname(ms_path), b2)
    assert os.path.exists(p), p
    return p


# ------------------------------------------------------------------ indices
def ndvi(r, n):
    with np.errstate(invalid='ignore', divide='ignore'):
        return (n - r) / (n + r)


def nirv(r, n):
    return ndvi(r, n) * n


def msavi(r, n):
    return (2 * n + 1 - np.sqrt(np.maximum((2 * n + 1) ** 2 - 8 * (n - r), 0))) / 2


IDXF = {'ndvi': ndvi, 'nirv': nirv, 'msavi': msavi}


# ------------------------------------------------------------------ geometry
def lattice_offset(t_coarse, t_fine, f=4):
    """Integer row/column offset of the upper-left corner of the coarse grid (2 m) on the fine grid (0.5 m): coarse pixel (I,J) covers the fine rows f*I+ro ... f*I+ro+f−1."""
    co = (t_coarse.c - t_fine.c) / t_fine.a
    ro = (t_coarse.f - t_fine.f) / t_fine.e
    assert abs(co - round(co)) < 1e-6 and abs(ro - round(ro)) < 1e-6, (co, ro)
    return int(round(ro)), int(round(co))


def read_win(path, row0, col0, h, w, band=1, fill=0):
    with rasterio.open(path) as d:
        return d.read(band, window=Window(col0, row0, w, h), boundless=True, fill_value=fill)


def _lap_pen():
    """2-D Laplacian (second-difference) smoothness penalty DᵀD on the KN×KN grid (penalises curvature; does not shrink the kernel towards 0)."""
    n = KN
    D1 = np.zeros((n - 2, n))
    for i in range(n - 2):
        D1[i, i:i + 3] = [1, -2, 1]
    I = np.eye(n)
    Dy = np.kron(D1, I)
    Dx = np.kron(I, D1)
    P = Dy.T @ Dy + Dx.T @ Dx
    return P / np.trace(P) * (n * n)


PEN = _lap_pen()


def _B(o, s):
    from scipy.special import ndtr
    s = max(s, 0.03)
    return ndtr((o + 1.0) / s) - ndtr((o - 1.0) / s)


def kern1d(o, d, s, al, s2):
    """1-D kernel: 2 m pixel box (half-width 1 m) ⊛ [(1+α)G(σ) − α G(σ2)], centre shifted by d (m); evaluated at the sub-pixel centres o (not normalised)."""
    return (1 + al) * _B(o - d, s) - al * _B(o - d, s2)


def param_kernel(th):
    """th = (dx east, dy south, σx, σy, α, σ2); returns the normalised KN×KN kernel (rows = southward offset, columns = eastward offset)."""
    dx, dy, sx, sy, al, s2 = th
    kx = kern1d(OFFS, dx, sx, al, s2)
    ky = kern1d(OFFS, dy, sy, al, s2)
    k = np.outer(ky, kx)
    return k / k.sum()


def setup_fonts():
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    import matplotlib.pyplot as plt
    font_manager.fontManager.addfont(FONT_PATH)
    name = font_manager.FontProperties(fname=FONT_PATH).get_name()
    plt.rcParams['font.family'] = name
    plt.rcParams['axes.unicode_minus'] = False
    return plt


def savefig(fig, name, dpi=150):
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    name_ = font_manager.FontProperties(fname=FONT_PATH).get_name()
    for ax in fig.get_axes():
        for t in [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels():
            t.set_fontfamily(name_)
        leg = ax.get_legend()
        if leg:
            for t in leg.get_texts():
                t.set_fontfamily(name_)
    fig.savefig(os.path.join(FIG, name), dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def boot_ci(x, grp, B=1000, seed=0, w=None):
    """Mean and 95% interval by bootstrap over groups (station clusters); x: unit values, units weighted equally."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    x, grp = x[ok], np.asarray(grp)[ok]
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0, 0
    g, inv = np.unique(grp, return_inverse=True)
    s = np.bincount(inv, x)
    c = np.bincount(inv)
    rng = np.random.default_rng(seed)
    G = len(g)
    idx = rng.integers(0, G, (B, G))
    bs = s[idx].sum(1) / c[idx].sum(1)
    return float(x.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), len(x), G
