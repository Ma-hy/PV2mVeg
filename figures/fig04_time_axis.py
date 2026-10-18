"""Fig. 4: one time axis, two satellite records, three surfaces.

(a) NDVI, 2025-08: inter-row ground of single-year cohorts (2 m, cross-section), the mixed Landsat composite of the
    same units in 2025 (cross-section), the mixed Landsat record of a fixed set of units followed since installation
    (longitudinal), and service lanes by age group.
(b) The age-group step in the three acquisitions (inter-row ground).
(c, d) NIRv and MSAVI with the same elements as (a).

Inputs (repository data): v12/cross_scene/{effects,lanes,cells_*}.csv, the station/block labels of
v12/cross_scene/secondary/1p5m/chrono/a1_*_det_thr50_cells.csv and the event-study curves (ts/ew).
Output: Fig4.png/.pdf; plotted values figdata/fig3_{cohorts,lanes,longitudinal}.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from _cohorts import data, surface_panel


def main():
    rows, lanes, curves = data()
    fig = figure(140)
    axa = fig.add_axes([.075, .585, .555, .36])
    axb = fig.add_axes([.725, .585, .255, .36])
    axc = fig.add_axes([.075, .145, .405, .31])
    axd = fig.add_axes([.575, .145, .405, .31])
    # ------------------------------------------------ (a)
    R = surface_panel(axa, 'ndvi', rows, lanes, curves, big=True)
    axa.set_ylim(-.152, .158)
    axa.set_ylabel('ΔNDVI')
    axa.set_xlabel('Years since detected installation, e (line)  ·  years since construction (2 m cohorts)')
    # direct labels
    axa.text(13.6, .153, 'inter-row ground, 2 m\n2025-08 cohorts', color=GROUND, fontsize=6.4, ha='right', va='top', linespacing=1.0, path_effects=halo(1.4))
    sq0 = R.set_index('age').loc[0]
    axa.annotate('mixed 30 m pixels,\nsame units in 2025', xy=(.2, sq0.mixed), xytext=(-4.3, -.075), fontsize=6.2, color=MIXED, ha='left', va='center',
                 linespacing=1.0, arrowprops=dict(arrowstyle='-', color=MIXED, lw=.5, shrinkB=2.5), path_effects=halo(1.4))
    axa.annotate('mixed 30 m pixels of 2,304 units (93 stations)\nfollowed since installation (Landsat)', xy=(3.0, curves['ndvi'].set_index('e').loc[3, 'est']), xytext=(-4.3, .095),
                 fontsize=6.2, color=MIXED, ha='left', va='center', linespacing=1.0, arrowprops=dict(arrowstyle='-', color=MIXED, lw=.5, shrinkB=1.5), path_effects=halo(1.4))
    axa.text(13.6, -.006, 'service lanes', color=LANE, fontsize=6.2, ha='right', va='top', path_effects=halo(1.4))
    # above the zero line it labels: below it the larger text would touch the curve at e = -1
    axa.text(-4.4, .004, 'matched-control\nexpectation', color=MUTED, fontsize=5.6, ha='left', va='bottom', linespacing=1.0)
    cn = curves['ndvi'].set_index('e')
    axa.annotate('e ≥ 9: only early cohorts remain\n(1,969 → 914 units)', xy=(10.5, (cn.loc[10, 'est'] + cn.loc[11, 'est']) / 2), xytext=(9.6, -.058), fontsize=5.5, color=MUTED,
                 ha='center', va='top', linespacing=1.0, arrowprops=dict(arrowstyle='-', color='#B0B8C1', lw=.45))
    for x, lab in [(1.6, '0–2 y'), (4, '3–5 y'), (7.0, '≥6 y')]:
        axa.text(x, .153, lab, ha='center', va='top', fontsize=6.0, color=MUTED)
    letter(axa, 'a', x=-.115, y=1.02, ha='left')
    axa.text(-.07, 1.02, 'One time axis, two satellite records, three surfaces (NDVI)', transform=axa.transAxes, fontsize=7.7, va='bottom', path_effects=bold(INK, .2))
    # ------------------------------------------------ (b)
    for j, t in enumerate(TAGS):
        q = rows[(rows.tag == t) & (rows['index'] == 'ndvi') & rows.age_group.isin(GROUPS)].set_index('age_group').loc[GROUPS]
        x = np.arange(3) + (j - 1) * .24
        axb.errorbar(x, q.ground, yerr=[q.ground - q.ground_lo, q.ground_hi - q.ground], fmt=DATE_MARK[t], ms=4.3, mfc=DATE[t], mec='white', mew=.5,
                     ecolor=DATE[t], elinewidth=1.1, capsize=1.6, lw=0, zorder=4, label=DATE_LABEL[t])
    axb.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    axb.set_xticks(range(3)); axb.set_xticklabels(['0–2 y', '3–5 y', '≥6 y'])
    axb.set_xlim(-.55, 2.55); axb.set_ylim(-.022, .105)
    axb.set_ylabel('ΔNDVI, inter-row ground')
    axb.grid(axis='y', color=GRID, lw=.45); axb.set_axisbelow(True); axb.tick_params(length=2.4, width=.5, pad=1.6)
    axb.legend(loc='upper left', fontsize=6.1, handletextpad=.25, borderaxespad=.1, labelspacing=.25, handlelength=1.0)
    letter(axb, 'b', x=-.27, y=1.02, ha='left')
    axb.text(-.17, 1.02, 'The same step in three acquisitions', transform=axb.transAxes, fontsize=7.7, va='bottom', path_effects=bold(INK, .2))
    # ------------------------------------------------ (c, d)
    for ax, k, lab, ylim in [(axc, 'nirv', 'c', (-.045, .036)), (axd, 'msavi', 'd', (-.088, .07))]:
        surface_panel(ax, k, rows, lanes, curves, big=False)
        ax.set_ylim(*ylim)
        ax.set_ylabel(f'Δ{NAME[k]}')
        ax.set_xlabel('Years since detected installation, e  ·  years since construction')
        # at the title size the heading takes two lines; the letter stands at the first
        letter(ax, lab, x=-.15, y=1.111, ha='left')
        ax.text(-.09, 1.03, f'{NAME[k]}: mixed pixels turn negative,' + '\n' + 'inter-row ground stays positive', transform=ax.transAxes, fontsize=7.2,
                va='bottom', path_effects=bold(INK, .18), linespacing=1.05)
    handles = [Line2D([], [], marker='o', ls='', ms=4, mfc=GROUND, mec='white', label='Inter-row ground, 2 m (2025-08 cohorts; joint + calibration 95%)'),
               Line2D([], [], marker='o', ls='', ms=4, mfc='white', mec=GROUND, label='Cohort with <6 stations (no interval)'),
               Line2D([], [], marker='s', ls='', ms=3.2, mfc='white', mec=MIXED, label='Mixed 30 m, same units in 2025'),
               Line2D([], [], color=MIXED, lw=1.2, label='Mixed 30 m, units installed by 2018, complete to e = 8 (2 km block 95%)'),
               Line2D([], [], color=LANE, lw=1.4, label='Service lanes by age group')]
    # two centred rows at 7 pt: the three surfaces, then the open cohorts and the fixed-unit record
    kw = dict(loc='lower center', fontsize=7.0, handletextpad=.35, columnspacing=1.4, handlelength=1.5, borderaxespad=0)
    fig.legend(handles=[handles[i] for i in (0, 2, 4)], bbox_to_anchor=(.5, .031), ncol=3, **kw)
    fig.legend(handles=[handles[i] for i in (1, 3)], bbox_to_anchor=(.5, .006), ncol=2, **kw)
    save(fig, 'Fig4')
    out = rows[rows.age_group.isin(GROUPS)].merge(lanes, on=['tag', 'index', 'age_group'])
    print(out[['tag', 'index', 'age_group', 'ground', 'mixed', 'lane']].round(4).to_string())


if __name__ == '__main__':
    main()
