"""Shared helpers of the supplementary diagnostic figures (Figs. S5, S6 and S9-S11).

Inputs are the supplementary diagnostics of the primary design (v12/cross_scene/secondary) and the primary
mechanism files (v12/cross_scene). Plotted values go to figdata/supp_diag/.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *

CSD = CS12
SECD = CSD / 'secondary'
OUTD = FIGDATA / 'supp_diag'; OUTD.mkdir(parents=True, exist_ok=True)
POL = ['1m', '1p5m', '2m']
CLR = {'1m': '#8FD0BB', '1p5m': '#0E8A6A', '2m': '#16435A'}
CLAB = {'1m': '≥1 m', '1p5m': '≥1.5 m (primary)', '2m': '≥2 m'}
CSHORT = {'1m': '≥1 m', '1p5m': '≥1.5 m', '2m': '≥2 m'}
PAIRS = [('06', '08'), ('08', '2026_08')]
PAIR_LAB = {('06', '08'): '2025-08 − 2025-06', ('08', '2026_08'): '2026-08 − 2025-08'}
QCOL = {'common_change': GROUND, 'own_change_same_units': MIXED, 'support_component': LANE}
QLAB = {'common_change': 'Common footprints', 'own_change_same_units': 'Own footprints, same units', 'support_component': 'Own − common'}


def rd(p):
    return pd.read_csv(p, dtype={'tag': str, 'group': str, 'a': str, 'b': str, 'date': str, 'age_group': str})


def keep(name, d):
    d.to_csv(OUTD / f'{name}.csv', index=False)


def head(ax, s, sub=None, x=0.0, y=1.0, size=7.3):
    # the heading 12 pt above the axes when a (larger) subtitle follows; sub='' keeps the heading at the height of its row
    ax.annotate(s, (x, y), xycoords='axes fraction', xytext=(0, 12 if sub is not None else 3), textcoords='offset points',
                fontsize=size, va='bottom', ha='left', path_effects=bold(INK, .2))
    if sub:
        ax.annotate(sub, (x, y), xycoords='axes fraction', xytext=(0, 2.2), textcoords='offset points',
                    fontsize=6.0, va='bottom', ha='left', color=SUB)


def pr(ax, x, est, lo, hi, color, marker='o', filled=True, ms=3.6, lw=.9, horizontal=False, z=3, label=None, mew=.9):
    est, lo, hi = (np.atleast_1d(np.asarray(v, float)) for v in (est, lo, hi))
    x = np.atleast_1d(np.asarray(x, float))
    ok = np.isfinite(lo) & np.isfinite(hi)
    err = np.vstack([np.where(ok, np.maximum(est - lo, 0), 0), np.where(ok, np.maximum(hi - est, 0), 0)])
    kw = dict(fmt=marker, ms=ms, color=color, ecolor=color, elinewidth=lw, capsize=1.3, capthick=lw * .8, zorder=z,
              mfc=color if filled else 'white', mec=color, mew=mew, label=label)
    if horizontal:
        ax.errorbar(est, x, xerr=err, **kw)
    else:
        ax.errorbar(x, est, yerr=err, **kw)


def legend_row(fig, handles, y, x=.5, ncol=None, size=6.6):
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(x, y), ncol=ncol or len(handles), frameon=False,
               fontsize=size, handletextpad=.4, columnspacing=1.3, borderaxespad=0)


def H(color, label, marker='o', filled=True, ls='none', lw=1.0, ms=4):
    return Line2D([], [], color=color, marker=marker, ls=ls, lw=lw, ms=ms, mfc=color if filled else 'white', mec=color, label=label)


def effects(p):
    return rd(CSD / 'effects.csv' if p == '1p5m' else CSD / p / 'effects.csv')


def cells(p, tag):
    return pd.read_csv((CSD if p == '1p5m' else CSD / p) / f'cells_{tag}.csv')
