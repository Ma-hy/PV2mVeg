"""Fig. S7: single-year cohorts of 2025-06 (a) and 2026-08 (b), drawn as in Fig. 4a.

Inputs and plotted values as Fig. 4 (_cohorts.data). Output: FigS7.png/.pdf
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from _cohorts import data, surface_panel


def main():
    rows, lanes, curves = data()
    fig = figure(70)
    for i, tag in enumerate(['06', '2026_08']):
        ax = fig.add_axes([.07 + i * .49, .2, .41, .68])
        surface_panel(ax, 'ndvi', rows, lanes, curves, tag=tag, big=False)
        ax.set_ylim(-.152, .16); ax.set_xlim(-4.6, 14.6); ax.set_xticks(range(-4, 15, 2)); ax.set_ylabel('ΔNDVI'); ax.set_xlabel('Years since installation, e  ·  array age')
        letter(ax, 'ab'[i], x=-.15, y=1.03, ha='left')
        ax.text(-.08, 1.03, f'{DATE_LABEL[tag]}: single-year cohorts', transform=ax.transAxes, fontsize=7.3, va='bottom', path_effects=bold(INK, .18))
    handles = [Line2D([], [], marker='o', ls='', ms=4, mfc=GROUND, mec='white', label='inter-row ground, 2 m'),
               Line2D([], [], marker='s', ls='', ms=3.2, mfc='white', mec=MIXED, label='mixed 30 m, same units'),
               Line2D([], [], color=MIXED, lw=1.2, label='mixed 30 m, units installed by 2018'),
               Line2D([], [], color=LANE, lw=1.4, label='service lanes')]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6.0, bbox_to_anchor=(.5, .0))
    save(fig, 'FigS7')


if __name__ == '__main__':
    main()
