"""Fig. 3: the mixed-pixel (Landsat, 30 m) event study.

(a) NDVI of the original cohort around the detected installation year (matched difference, regression
    counterfactual, units installed by 2018, spatial placebo), with the pre-period inset, the same-target
    estimator comparison and the unit counts; (b, c) NIRv and MSAVI; (d) red and NIR reflectance effects.

Inputs (repository data):
  ts/ew/es_raw_sum_{ndvi,nirv,msavi}_det_treated/curves.csv   event-study curves
  ts/ew/placebo_sum_ndvi_det/space_placebo_curve_null.csv      NDVI spatial placebo
  legacy/fig2_estimators_x1.csv, fig2_counts.csv, fig2d_band_relative.csv
                                                               estimator comparison, counts, band effects
Output: Fig3.png/.pdf
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *

EW = P.es
FD = P.legacy
PLACEBO = '#E6EAEE'
CONSTR = '#F2F4F6'


def curves(k, estimator='diff', subset='all'):
    d = pd.read_csv(EW / f'es_raw_sum_{k}_det_treated' / 'curves.csv')
    d = d[(d.estimator == estimator) & (d.weight == 'cell') & (d.subset == subset)]
    return d.sort_values('e')


def event(ax, k, placebo=None, big=False):
    col = INDEXC[k]
    d = curves(k); r = curves(k, 'reg'); b = curves(k, subset='balanced T≤2018')
    ax.axvspan(-1.5, .5, color=CONSTR, lw=0, zorder=0)
    if placebo is not None:
        ax.fill_between(placebo.e, placebo['null_p2.5'], placebo['null_p97.5'], color=PLACEBO, lw=0, zorder=1)
    ax.axhline(0, color='#8A96A3', lw=.6, zorder=1)
    ax.fill_between(d.e, d.lo, d.hi, color=col, alpha=.2, lw=0, zorder=2)
    ax.plot(r.e, r.est, color=col, lw=.8, ls=(0, (3.5, 1.6)), zorder=3)
    ax.plot(b.e, b.est, color=INK, lw=.8, ls=(0, (1, 1.3)), zorder=3)
    ax.plot(d.e, d.est, color=col, lw=1.3 if big else 1.1, zorder=4)
    ax.scatter(d.e, d.est, s=9 if big else 6, color=col, edgecolors='white', linewidths=.4, zorder=5)
    ax.set_xlim(-12.6, 12.6); ax.set_xticks(range(-12, 13, 4)); ax.set_xticks(range(-12, 13, 2), minor=True)
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.3, pad=1.4)
    return d


def inset_pre(ax, k, rect, placebo=None):
    d = curves(k); p = d[d.e <= -2]
    ins = ax.inset_axes(rect)
    if placebo is not None:
        nn = placebo[placebo.e <= -2]
        ins.fill_between(nn.e, nn['null_p2.5'], nn['null_p97.5'], color=PLACEBO, lw=0)
    ins.axhline(0, color='#8A96A3', lw=.45)
    ins.fill_between(p.e, p.lo, p.hi, color=INDEXC[k], alpha=.2, lw=0)
    ins.plot(p.e, p.est, color=INDEXC[k], lw=.8); ins.scatter(p.e, p.est, s=4, color=INDEXC[k], edgecolors='white', linewidths=.3, zorder=3)
    lim = max(abs(p.lo).max(), abs(p.hi).max()) * 1.15
    ins.set_ylim(-lim, lim); ins.set_xlim(-12.4, -1.6)
    ins.set_xticks([-12, -8, -4]); ins.set_yticks([0])
    ins.tick_params(labelsize=5.4, length=1.4, pad=.8)
    ins.set_facecolor('white')
    for s in ins.spines.values():
        s.set_color('#B0B8C1'); s.set_linewidth(.4)
    if placebo is None:  # the small insets of (b, c) are narrower than their label; it goes above them
        ins.text(0, 1.04, f'e ≤ −2,\n±{lim:.3f}'.replace('-', '−'), transform=ins.transAxes, fontsize=5.4, color=SUB, va='bottom', linespacing=1.05)
    else:
        ins.text(.03, .96, f'e ≤ −2, ±{lim:.3f}'.replace('-', '−'), transform=ins.transAxes, fontsize=5.4, color=SUB, va='top')
    return ins


def strip(ax, k, x1, show_labels=False):
    s = x1[x1['index'] == k]
    spec = [('matched', 'matched controls', 'o', INDEXC[k], INDEXC[k]), ('never', 'never-treated DiD', 's', '#6B7785', '#6B7785'), ('nyt', 'not-yet-treated DiD', 'D', 'white', INK)]
    ax.axhline(0, color='#8A96A3', lw=.6, zorder=1)
    for gx, tgt in [(0, 'e=-1'), (1, 'plateau e=6-10')]:
        for j, (est, lab, mk, fc, ec) in enumerate(spec):
            r = s[(s.target == tgt) & (s.estimator == est)].iloc[0]
            x = gx + (j - 1) * .26
            ax.errorbar(x, r.est, yerr=[[r.est - r.lo], [r.hi - r.est]], fmt=mk, ms=3.2, mfc=fc, mec=ec, mew=.8, ecolor=ec if est != 'nyt' else INK,
                        elinewidth=.8, capsize=0, zorder=3, label=lab if gx == 0 else None)
    ax.axvline(.5, color='#D5DBE1', lw=.5, zorder=0)
    ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels(['−1', '6–10'], fontsize=5.9)
    ax.tick_params(axis='y', labelleft=False, length=0); ax.tick_params(axis='x', length=2, pad=1.2)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True)
    return ax


def main():
    x1 = pd.read_csv(FD / 'fig2_estimators_x1.csv')
    placebo = pd.read_csv(EW / 'placebo_sum_ndvi_det' / 'space_placebo_curve_null.csv')
    cnt = pd.read_csv(FD / 'fig2_counts.csv')
    band = pd.read_csv(FD / 'fig2d_band_relative.csv')
    fig = figure(124)
    # panel rectangles were laid out for a 118 mm figure: 6 mm more height for a two-row key at the bottom,
    # the upper block 2.4 mm further up; panels keep their size in mm
    pos = lambda r, up: [r[0], (r[1] * 118 + up) / 124, r[2], r[3] * 118 / 124]
    axa = fig.add_axes(pos([.07, .565, .735, .37], 8.4)); sa = fig.add_axes(pos([.815, .565, .17, .37], 8.4), sharey=axa)
    axn = fig.add_axes(pos([.07, .505, .735, .045], 8.4), sharex=axa)
    axb = fig.add_axes(pos([.07, .095, .205, .32], 6)); sb = fig.add_axes(pos([.28, .095, .055, .32], 6), sharey=axb)
    axc = fig.add_axes(pos([.41, .095, .205, .32], 6)); sc = fig.add_axes(pos([.62, .095, .055, .32], 6), sharey=axc)
    axd = fig.add_axes(pos([.775, .095, .21, .32], 6))
    # (a)
    event(axa, 'ndvi', placebo, big=True)
    axa.set_ylim(-.04, .056); axa.set_ylabel('ΔNDVI, treated − control')
    axa.tick_params(labelbottom=False)
    inset_pre(axa, 'ndvi', [.05, .6, .3, .34], placebo)
    axa.text(-.5, .053, 'construction', ha='center', va='top', fontsize=5.8, color=MUTED)
    d = curves('ndvi')
    axa.annotate(f'plateau e ≥ 6: {signed(0.031, 3)}', xy=(9, d.set_index("e").loc[9, "est"]), xytext=(4.6, .049), fontsize=6.0, color=INDEXC['ndvi'],
                 arrowprops=dict(arrowstyle='-', color=INDEXC['ndvi'], lw=.5), ha='left', va='center')
    axa.annotate('first panel-covered\nmidsummer (e = 0)', xy=(0, d.set_index('e').loc[0, 'est']), xytext=(2.2, -.031), fontsize=5.9, color=SUB,
                 arrowprops=dict(arrowstyle='-', color='#9AA5B1', lw=.5), ha='left', va='center', linespacing=1.0)
    strip(sa, 'ndvi', x1)
    sa.legend(loc='upper center', bbox_to_anchor=(.5, 1.0), fontsize=5.8, handletextpad=.2, borderaxespad=.1, labelspacing=.25, ncol=1, frameon=True, framealpha=.9, edgecolor='none')
    sa.set_xticklabels(['e = −1', 'e = 6–10'], fontsize=6.0)
    sa.text(.5, -.165, 'same target, common support', transform=sa.transAxes, ha='center', va='top', fontsize=5.6, color=MUTED)
    # counts sparkline
    c = cnt.sort_values('e')
    axn.bar(c.e, c.n_units, width=.72, color='#C5CED8', edgecolor='none')
    axn.set_ylim(0, 5200); axn.set_yticks([]); axn.spines['left'].set_visible(False)
    for e in [0, 6, 12]:
        r = c[c.e == e].iloc[0]
        if e == 0:  # inside the tall bar, clear of the axis of (a) above
            axn.text(e, r.n_units - 250, f'{int(r.n_units):,} · {int(r.n_stations)} st.', ha='center', va='top', fontsize=5.3, color=SUB)
        else:
            axn.text(e, r.n_units + 150, f'{int(r.n_units):,} · {int(r.n_stations)} st.', ha='center', va='bottom', fontsize=5.3, color=SUB)
    axn.set_xlabel('Years since installation, e'); axn.tick_params(length=2.3, pad=1.2)
    axn.set_ylabel('units', rotation=0, ha='right', va='center', fontsize=5.8, color=MUTED, labelpad=2)
    # (b, c)
    for ax, s_, k, ylim in [(axb, sb, 'nirv', (-.031, .006)), (axc, sc, 'msavi', (-.055, .012))]:
        event(ax, k)
        ax.set_ylim(*ylim); ax.set_ylabel(f'Δ{NAME[k]}, treated − control'); ax.set_xlabel('Years since installation, e')
        strip(s_, k, x1)
        v2fix(s_, 6.8)  # the narrow strips keep '−1' and '6–10' apart from the main axis' '12'
        ins = inset_pre(ax, k, [.06, .07, .36, .28])
    # (d) red and NIR
    for b_, lab in [('red', 'red'), ('nir', 'NIR')]:
        z = band[band.band == b_].sort_values('e')
        axd.fill_between(z.e, z.lo_pct, z.hi_pct, color=BANDC[b_], alpha=.18, lw=0)
        axd.plot(z.e, z.est_pct, color=BANDC[b_], lw=1.1)
        axd.scatter(z.e, z.est_pct, s=5, color=BANDC[b_], edgecolors='white', linewidths=.3, zorder=3)
        last = z[z.e >= 6].est_pct.mean()
        axd.text(12.4, last + (-2.6 if b_ == 'red' else 2.2), f'{lab} {last:.0f}%'.replace('-', '−'), color=BANDC[b_], fontsize=6.0, ha='right', va='center')
    axd.axvspan(-1.5, .5, color=CONSTR, lw=0, zorder=0); axd.axhline(0, color='#8A96A3', lw=.6)
    axd.set_xlim(-12.6, 12.6); axd.set_xticks([-12, -6, 0, 6, 12]); axd.set_ylim(-42, 8)
    axd.set_ylabel('Reflectance effect (% of control)'); axd.set_xlabel('Years since installation, e')
    axd.grid(axis='y', color=GRID, lw=.4); axd.set_axisbelow(True); axd.tick_params(length=2.3, pad=1.4)
    axd.annotate('red rises at\ne = −1\n(cleared\nground)', xy=(-1, 4.6), xytext=(-11.8, -15),
                 fontsize=5.6, color=SUB, va='center', linespacing=1.0, arrowprops=dict(arrowstyle='-', color='#9AA5B1', lw=.5))
    # labels
    for ax, l, t, dx in [(axa, 'a', 'NDVI: mixed pixels green after installation', -.065), (axb, 'b', 'NIRv stays below the counterfactual', -.23),
                         (axc, 'c', 'MSAVI stays below the counterfactual', -.23), (axd, 'd', 'Red and NIR both fall by a third', -.24)]:
        letter(ax, l, x=dx, y=1.02, ha='left')
        ax.text(dx + (.035 if ax is axa else .11), 1.02, t, transform=ax.transAxes, fontsize=7.3, va='bottom', path_effects=bold(INK, .18))
    handles = [Line2D([], [], color=INDEXC['ndvi'], lw=1.2, marker='o', ms=2.8, mec='white', label='matched difference (2 km block 95%)'),
               Line2D([], [], color=INDEXC['ndvi'], lw=.8, ls=(0, (3.5, 1.6)), label='regression counterfactual'),
               Line2D([], [], color=INK, lw=.8, ls=(0, (1, 1.3)), label='units installed by 2018'),
               Patch(facecolor=PLACEBO, label='spatial placebo 95% (NDVI)'),
               Patch(facecolor=CONSTR, edgecolor='#D5DBE1', lw=.4, label='e = −1 and 0')]
    # two rows of three (legend columns fill top-down: reorder so the rows read left to right)
    fig.legend(handles=[handles[i] for i in (0, 3, 1, 4, 2)], loc='lower center', bbox_to_anchor=(.5, .0), ncol=3, fontsize=7.5,
               handletextpad=.35, columnspacing=1.6, handlelength=1.8, labelspacing=.3)
    save(fig, 'Fig3')


if __name__ == '__main__':
    main()
