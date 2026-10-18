"""House style of the v1.2 figure family (main Figs. 1, 3-8 and supplementary Figs. S4-S12).

Figure scripts start with `from pv2mveg.plotting.style import *`, which also provides numpy (np), pandas (pd),
matplotlib.pyplot (plt), json, Path and the matplotlib classes used by the scripts, the paths of the results
(CS12, BP12, SEC12, TS, FIGDATA) and the colours, sizes and helpers below.

Figures are composed at 165 mm width; `save()` hands them to the v2.0 save pipeline (profile.py: 500 dpi, black
text, role-normalised type sizes, printed floor of 6.5 pt, contrast check). The figure font is set in the
configuration (PV2MVEG_FONT); the published figures used times+simsun.ttf (a Times New Roman / SimSun font that is
not redistributed here). Without it, STIX (bundled with matplotlib) is used and the layout differs slightly.
The font has no bold face; emphasis uses size, colour and a thin stroke ("fake bold").
"""
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm, TwoSlopeNorm, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon, Circle, FancyBboxPatch, Patch
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.ticker import FuncFormatter, MultipleLocator, MaxNLocator
from pyproj import Transformer

from pv2mveg.config import P, TAGS
from pv2mveg.plotting import profile as _v2
from pv2mveg.plotting.fonts import FONT_NAME, FONT_FILE
from pv2mveg.plotting.house import house  # noqa: F401  (re-exported for the scripts)

# v2.0 profile: black text, larger type, 500 dpi
V2 = _v2.V2
V2SIZES = _v2.SIZES
v2pick = _v2.v2pick
v2fix = _v2.v2fix

# result folders (config.P; the same names in the repository and the project layout)
V12 = P.v12
CS12 = P.cs12                 # primary cross-scene design (1 m and 2 m in subfolders), lanes, spectral, distance, optics, aux
BP12 = P.bp12                 # unmatched buffer baseline and its variants
SEC12 = P.sec12               # supplementary diagnostics of the primary design
TS = P.ts                     # Landsat unit series and the mixed-pixel event study
FIGDATA = P.figdata           # plotted values written by the figure scripts
FIGDATA.mkdir(parents=True, exist_ok=True)
VERSION = 'v2.0'

MM = 1 / 25.4
WIDTH = 165  # mm

INK = '#1F2933'
SUB = '#52606D'
MUTED = '#7B8794'
FAINT = '#CBD2D9'
GRID = '#E4E7EB'
PAPERBG = '#FFFFFF'
if V2:  # every text colour is pure black (lines that used these colours become black as well)
    INK = SUB = MUTED = _v2.BLACK

# semantic colours ---------------------------------------------------------
GROUND = '#0E8A6A'      # resolved inter-row ground
GROUND_L = '#7CC9B2'
MIXED = '#4F6D8F'       # mixed Landsat pixels
MIXED_L = '#AFC3D8'
LANE = '#C98A1B'        # service lanes
LANE_L = '#EBC98A'
TREAT = '#1F6FD1'
CTRL = '#E8792B'
TREAT_GLOW = '#3FD9FF'
CTRL_GLOW = '#FFB13B'
BUFFER_GLOW = '#C9B8FF'
DATE = {'06': '#0090D0', '08': '#F0592A', '2026_08': '#8E3BD8'}
DATE_MARK = {'06': 'o', '08': 's', '2026_08': 'D'}
DATE_LABEL = {'06': '2025-06', '08': '2025-08', '2026_08': '2026-08'}
DATE_SHORT = DATE_LABEL
AGE = {'0–2': '#F6B800', '3–5': '#14A38B', '≥6': '#3B2C96'}
INDEXC = {'ndvi': '#1F6FD1', 'nirv': '#E4572E', 'msavi': '#12A086'}
BANDC = {'red': '#D62839', 'nir': '#7B2CBF'}
NAME = {'ndvi': 'NDVI', 'nirv': 'NIRv', 'msavi': 'MSAVI', 'red': 'Red', 'nir': 'NIR'}
GROUPS = ['0–2', '3–5', '≥6']

# diverging map of the adjusted NDVI contrast: earth to teal
DIVERGE = LinearSegmentedColormap.from_list('earth_teal', ['#7A3B07', '#C4762B', '#E9C27D', '#F6F1E4', '#9ED3C2', '#2E9C83', '#00594A'])
# installation year: deep blue (older) to warm red (recent)
YEAR_COLORS = ['#1B2A57', '#22427A', '#1F6698', '#18879E', '#21A08A', '#62B56A', '#A9C24A', '#E3C23A', '#F2A33A', '#EE7A3A', '#DE5143', '#C2334F', '#9E1F57', '#7A1458', '#57104F']
YEAR_CMAP = ListedColormap(YEAR_COLORS)
YEAR_NORM = BoundaryNorm(np.arange(2011.5, 2027.5), len(YEAR_COLORS))
GEO = Transformer.from_crs(32647, 4326, always_xy=True)

plt.rcParams.update({
    'font.family': FONT_NAME, 'font.size': 7.2, 'axes.titlesize': 7.8, 'axes.labelsize': 7.4,
    'xtick.labelsize': 6.8, 'ytick.labelsize': 6.8, 'legend.fontsize': 6.8,
    'axes.edgecolor': '#8A96A3', 'axes.linewidth': 0.55, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': '#5F6B7A', 'ytick.color': '#5F6B7A', 'xtick.major.width': .5, 'ytick.major.width': .5,
    'xtick.major.size': 2.4, 'ytick.major.size': 2.4, 'xtick.minor.size': 1.3, 'ytick.minor.size': 1.3,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': True,
    'lines.linewidth': 1.0, 'lines.markersize': 3.4, 'legend.frameon': False, 'legend.handlelength': 1.6,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 600, 'figure.dpi': 150,
    'mathtext.fontset': 'custom', 'mathtext.rm': FONT_NAME, 'mathtext.it': FONT_NAME, 'mathtext.default': 'regular',
    'axes.titlepad': 4, 'axes.labelpad': 2.5, 'hatch.linewidth': .4,
})
if V2:  # base sizes of the v2.0 figures (font sizes at the composed width of 165 mm); tick lines and spines stay grey
    plt.rcParams.update({
        'font.size': V2SIZES['tick'], 'axes.titlesize': V2SIZES['title'], 'axes.labelsize': V2SIZES['label'],
        'xtick.labelsize': V2SIZES['tick'], 'ytick.labelsize': V2SIZES['tick'], 'legend.fontsize': V2SIZES['legend'],
        'legend.title_fontsize': V2SIZES['legend'], 'xtick.labelcolor': _v2.BLACK, 'ytick.labelcolor': _v2.BLACK,
        'savefig.dpi': _v2.DPI,
    })

BOLD = [pe.withStroke(linewidth=.45, foreground=INK)]


def bold(color=INK, w=.45):
    return [pe.withStroke(linewidth=w, foreground=color)]


def halo(w=1.6, color='white', alpha=.9):
    return [pe.withStroke(linewidth=w, foreground=color, alpha=alpha)]


def figure(height_mm, width_mm=WIDTH):
    return plt.figure(figsize=(width_mm * MM, height_mm * MM))


def letter(ax, s, x=-0.02, y=1.0, fig=None, size=9.5, color=INK, ha='right'):
    """Panel letter placed just outside the upper-left corner of the axes box."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=size, color=color, ha=ha, va='bottom',
            path_effects=bold(color, .55), zorder=50)


def title(ax, s, sub=None, loc='left', size=7.8, pad=3.5):
    t = ax.set_title(s, loc=loc, fontsize=size, color=INK, pad=pad + (6.5 if sub else 0), path_effects=bold(INK, .25))
    if sub:
        ax.text(0 if loc == 'left' else .5, 1.0, sub, transform=ax.transAxes, fontsize=6.5, color=SUB,
                ha='left' if loc == 'left' else 'center', va='bottom')
    return t


def tidy(ax, grid='y', zero=True, zero_label=None):
    if zero:
        ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=.45, zorder=0)
        ax.set_axisbelow(True)
    ax.tick_params(length=2.4, width=.5, pad=1.6)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#8A96A3')


def signed(v, nd=3):
    return f'{v:+.{nd}f}'.replace('-', '−')


def _font(fig):
    for item in fig.findobj(matplotlib.text.Text):
        item.set_fontfamily(FONT_NAME)


def save(fig, name, normalize=None):
    """Save under the published name (e.g. 'Fig3') through the v2.0 pipeline (profile.save)."""
    return _v2.save(fig, name, _font, 'v12', normalize_roles=normalize)


# ------------------------------------------------------------------ maps (restricted-data figures)
# Drawing defaults of the published map figures. In the research code, the true-colour background and the native
# example came from earlier map modules that changed matplotlib defaults when first imported; the published Figs. 1,
# 5, 6 and S8 were drawn under these defaults from that point on. They are applied at the same point here:
# MAP_RC at the first true-colour background of a figure (Figs. 5, 6, S8) or before Fig. 1 is created, AGE_RC
# (Fig. 1 only) after it.
MAP_RC = {'font.size': 9, 'axes.spines.top': True, 'axes.spines.right': True, 'axes.linewidth': .65, 'pdf.fonttype': 42, 'savefig.dpi': 400}
AGE_RC = {'font.size': 16, 'axes.spines.top': True, 'axes.spines.right': True, 'axes.unicode_minus': False, 'pdf.fonttype': 42}
_MAP_RC_APPLIED = []


def map_defaults():
    """Apply MAP_RC once per process."""
    if not _MAP_RC_APPLIED:
        plt.rcParams.update(MAP_RC)
        _MAP_RC_APPLIED.append(True)


def lonlat_frame(ax, ext, step=.1, labels=True, nd=1, color='#5F6B7A', size=6.4, sides='lb'):
    ax.set_xlim(ext[:2]); ax.set_ylim(ext[2:])
    ax.set_aspect(1 / np.cos(np.deg2rad(np.mean(ext[2:]))))
    ax.xaxis.set_major_locator(MultipleLocator(step)); ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.{nd}f}°E'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.{nd}f}°N'))
    ax.tick_params(labelsize=size, length=2, width=.45, colors=color, pad=1.2,
                   labelbottom='b' in sides, labelleft='l' in sides, labeltop='t' in sides, labelright='r' in sides,
                   bottom='b' in sides, left='l' in sides, top='t' in sides, right='r' in sides)
    if V2:  # black labels, grey tick marks
        ax.tick_params(labelcolor=_v2.BLACK)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color('#6B7785'); s.set_linewidth(.5)
    if not labels:
        ax.set_xticklabels([]); ax.set_yticklabels([])


def scalebar(ax, km, x=.05, y=.06, color='white', text_color=None, lw=1.4, size=6.3, ext=None, segments=2):
    """Alternating black/white scale bar in a lon/lat axes."""
    ext = ext or (*ax.get_xlim(), *ax.get_ylim())
    lat = np.mean(ext[2:])
    deg = km / (111.32 * np.cos(np.deg2rad(lat)))
    x0 = ext[0] + x * (ext[1] - ext[0]); y0 = ext[2] + y * (ext[3] - ext[2])
    seg = deg / segments
    for i in range(segments):
        c = color if i % 2 == 0 else ('#1F2933' if color == 'white' else '#9AA5B1')
        ax.plot([x0 + i * seg, x0 + (i + 1) * seg], [y0, y0], color=c, lw=lw, solid_capstyle='butt', zorder=30)
    ax.plot([x0, x0 + deg], [y0, y0], color='#1F2933' if color == 'white' else 'white', lw=lw + .9, solid_capstyle='butt', zorder=29, alpha=.7)
    ax.text(x0 + deg / 2, y0 + .012 * (ext[3] - ext[2]), f'{km:g} km', ha='center', va='bottom', fontsize=size,
            color=text_color or color, zorder=31, path_effects=halo(1.2, '#1F2933' if color == 'white' else 'white', .6))


def north(ax, x=.06, y=.86, size=6.5, color='white'):
    ax.annotate('N', xy=(x, y + .07), xytext=(x, y), xycoords='axes fraction', textcoords='axes fraction', ha='center', va='center',
                fontsize=size, color=color, arrowprops=dict(arrowstyle='-|>,head_width=.25,head_length=.5', color=color, lw=.8),
                path_effects=halo(1.1, '#1F2933' if color == 'white' else 'white', .55), zorder=31)


def north_arrow(ax, x=.05, y=.82, h=.09, color=INK, fill='white', size=6.4):
    """Half-filled north needle with N above it, at axes fraction (x, y); h is its height in axes fraction."""
    ax.apply_aspect()
    bb = ax.get_window_extent()
    w = h * bb.height / bb.width / 1.6
    nax = ax.inset_axes([x - w / 2, y, w, h], zorder=40)
    nax.set_xlim(-1, 1); nax.set_ylim(0, 3.2); nax.axis('off')
    nax.add_patch(Polygon([[0, 2.15], [-.62, 0], [0, .5]], closed=True, fc=color, ec=color, lw=.5))
    nax.add_patch(Polygon([[0, 2.15], [.62, 0], [0, .5]], closed=True, fc=fill, ec=color, lw=.5))
    nax.text(0, 2.3, 'N', ha='center', va='bottom', fontsize=size, color=color)
    return nax


def cell_polygons(x, y, half=75.):
    xx = np.asarray(x)[:, None] + np.array([-half, half, half, -half])[None, :]
    yy = np.asarray(y)[:, None] + np.array([-half, -half, half, half])[None, :]
    lon, lat = GEO.transform(xx, yy)
    return np.stack([lon, lat], axis=-1)


def true_color_tag(ext, width, tag):
    """Original true colour of one acquisition (bands 3/2/1, scene-wise 1-99% stretch); display only.

    Read from the display cache of the project (P.rgb_cache, one npz per extent, width and acquisition); the cache
    is built from the Beijing-3 multispectral scenes, which are restricted (see README.md)."""
    import hashlib
    P.need_project('A Beijing-3 true-colour background')
    map_defaults()
    height = int(width * (ext[3] - ext[2]) / ((ext[1] - ext[0]) * np.cos(np.deg2rad(36))))
    spec = dict(ext=[float(v) for v in ext], width=width, height=height, tag=tag, version=3)
    key = hashlib.sha256(json.dumps(spec).encode()).hexdigest()[:12]
    cache = P.rgb_cache / f'rgb_{tag}_{key}.npz'
    if not cache.exists():
        from pv2mveg.geometry.truecolor import build_true_color
        rgb = build_true_color(ext, width, height, tag)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, rgb=rgb)
        return rgb
    return np.load(cache)['rgb']


def rgb_cache(ext, width):
    """2025-08 true colour with gaps filled from 2025-06, then 2026-08 (display only)."""
    base = true_color_tag(ext, width, '08')
    fills = {}
    for tag in ['06', '2026_08']:
        gap = (base == 232).all(axis=2)
        if not gap.any():
            break
        other = true_color_tag(ext, width, tag)
        ok = gap & ~(other == 232).all(axis=2)
        base = base.copy(); base[ok] = other[ok]
        fills[tag] = int(ok.sum())
    return base, dict(fills=fills)


def stations():
    """Hand-drawn station polygons (restricted), in lon/lat."""
    import geopandas as gpd
    P.need_project('The hand-drawn station polygons')
    return gpd.read_file(P.stations).to_crs(4326)


def polygon_paths(gdf):
    out = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        parts = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
        for p in parts:
            out.append(np.asarray(p.exterior.coords))
    return out


def dark_basemap(rgb, strength=.30, tint=(0.03, 0.07, 0.12), desat=.75):
    """Night-mode version of a true-colour image for luminous overlays."""
    a = rgb.astype(float) / 255.
    grey = a.mean(axis=2, keepdims=True)
    a = grey * desat + a * (1 - desat)
    a = a * strength + np.array(tint)[None, None, :]
    return np.clip(a, 0, 1)


def light_basemap(rgb, lift=.62, desat=.85):
    a = rgb.astype(float) / 255.
    grey = a.mean(axis=2, keepdims=True)
    a = grey * desat + a * (1 - desat)
    return np.clip(lift + (1 - lift) * a, 0, 1)


def bezier(p0, p1, bend=.22, n=24):
    """Quadratic Bezier arcs between point arrays p0 and p1, bulging to the left of p0->p1."""
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    mid = (p0 + p1) / 2
    d = p1 - p0
    perp = np.stack([-d[:, 1], d[:, 0]], 1)
    ctrl = mid + bend * perp
    t = np.linspace(0, 1, n)[None, :, None]
    return (1 - t) ** 2 * p0[:, None, :] + 2 * (1 - t) * t * ctrl[:, None, :] + t ** 2 * p1[:, None, :]


def gradient_arcs(ax, arcs, c0, c1, lw=.3, alpha=.08, z=5):
    """Draw arcs (m, n, 2) as segments coloured from c0 to c1 along each arc."""
    m, n, _ = arcs.shape
    seg = np.stack([arcs[:, :-1, :], arcs[:, 1:, :]], axis=2).reshape(-1, 2, 2)
    t = np.tile(np.linspace(0, 1, n - 1), m)
    c0 = np.array(to_rgb(c0)); c1 = np.array(to_rgb(c1))
    cols = np.concatenate([(1 - t)[:, None] * c0 + t[:, None] * c1, np.full((len(t), 1), alpha)], 1)
    lc = LineCollection(seg, colors=cols, linewidths=lw, zorder=z, rasterized=True, capstyle='round')
    ax.add_collection(lc)
    return lc


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float), encoding='utf-8')
