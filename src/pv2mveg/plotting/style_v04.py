"""Style of the earlier figure code that draws the redrawn supplementary Figs. S1 and S2 (radiometric calibration,
unmixing tests), at the v2.0 sizes.

The figures are composed at 190 mm; 8.65 pt equals 7.5 pt in the 165 mm figures of the v1.2 family (printed 7.2 pt
at 158.75 mm). Colours are colour-blind safe (Okabe-Ito / Tol based):
  indices  : NDVI blue, NIRv vermillion, MSAVI bluish green
  periods  : 2025-06 orange, 2025-08 (BJ3A1) reddish purple, 2026-08 black
  bands    : red band red, NIR band purple
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pv2mveg.config import P  # noqa: E402
from pv2mveg.plotting import profile as _v2  # noqa: E402
from pv2mveg.plotting.fonts import FONT_NAME  # noqa: E402

FIGDATA = str(P.figdata)
os.makedirs(FIGDATA, exist_ok=True)

V2 = _v2.V2
FS = 8.65          # base font size (pt) at the composed width of 190 mm
FS_LABEL = 10      # panel labels
MM = 1 / 25.4
SINGLE = 90     # mm
DOUBLE = 190    # mm

plt.rcParams.update({
    'font.family': FONT_NAME,
    'font.sans-serif': [FONT_NAME],
    'font.serif': [FONT_NAME],
    'font.size': FS,
    'axes.titlesize': FS,
    'axes.labelsize': FS,
    'xtick.labelsize': FS,
    'ytick.labelsize': FS,
    'legend.fontsize': FS,
    'figure.titlesize': FS,
    'axes.unicode_minus': True,
    'mathtext.fontset': 'custom',
    'mathtext.rm': FONT_NAME,
    'mathtext.it': FONT_NAME,
    'mathtext.bf': FONT_NAME,
    'mathtext.default': 'regular',
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'lines.linewidth': 1.1,
    'lines.markersize': 3.5,
    'patch.linewidth': 0.5,
    'legend.frameon': False,
    'legend.handlelength': 1.8,
    'legend.borderaxespad': 0.3,
    'legend.labelspacing': 0.3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.dpi': 600,
    'figure.dpi': 150,
    'hatch.linewidth': 0.5,
})
# axis labels and titles: the equivalent of 8.0 and 8.2 pt in a 165 mm figure
plt.rcParams.update({'axes.labelsize': 9.2, 'axes.titlesize': 9.45, 'text.color': _v2.BLACK, 'axes.labelcolor': _v2.BLACK,
                     'xtick.labelcolor': _v2.BLACK, 'ytick.labelcolor': _v2.BLACK, 'savefig.dpi': _v2.DPI})

# ---------------------------------------------------------------- colours
C_INDEX = {'ndvi': '#0072B2', 'nirv': '#D55E00', 'msavi': '#009E73'}
L_INDEX = {'ndvi': 'NDVI', 'nirv': 'NIRv', 'msavi': 'MSAVI', 'red': 'Red', 'nir': 'NIR'}
C_PERIOD = {'2025-06': '#E69F00', '2025-08': '#CC79A7', '2026-08': '#000000'}
M_PERIOD = {'2025-06': 'o', '2025-08': 's', '2026-08': 'D'}
L_PERIOD = {'2025-06': 'Jun 2025', '2025-08': 'Aug 2025 (BJ3A1)', '2026-08': 'Aug 2026'}
C_SENSOR = {'Landsat': '#332288', 'S2': '#56B4E9'}
C_RAW = '#8C8C8C'
C_BAND = {'red': '#CC3311', 'nir': '#6A3D9A'}
C_CONSTR = '#EFEFEF'


def mm(x):
    return x * MM


def new_fig(w_mm, h_mm, **kw):
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM), **kw)
    return fig


def panel_label(ax, s, x=-0.02, y=1.0, ha='right', va='bottom', **kw):
    """(a), (b)... at upper-left outside the axes."""
    ax.text(x, y, f'({s})', transform=ax.transAxes, ha=ha, va=va, fontsize=FS_LABEL,
            **kw)


def zero_line(ax, **kw):
    ax.axhline(0, color='k', lw=0.5, zorder=1, **kw)


def apply_font(fig):
    """Force the figure font on every text artist right before saving."""
    for t in fig.findobj(matplotlib.text.Text):
        t.set_fontfamily(FONT_NAME)


def save(fig, name):
    """Save under the published name through the v2.0 pipeline (profile.save; no role normalisation)."""
    return _v2.save(fig, name, apply_font, 'style')


def fmt(v, nd=3, sign=True):
    if v is None or not np.isfinite(v):
        return '—'
    s = f'{v:+.{nd}f}' if sign else f'{v:.{nd}f}'
    return s.replace('-', '−')
