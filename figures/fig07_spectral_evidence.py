"""Fig. 7: spectral, red-edge, distance and clearance evidence on the inter-row signal.

(a) relative inter-row minus control reflectance in the four calibrated 2 m bands (BJ3A1 scenes with visible
    calibration), >=6 and 0-2 years, three acquisitions;
(b) red-edge index differences in the auxiliary 1.2 m six-band imagery (>=6 years);
(c) the contrast in centre-distance bands across the gap, 2025-08;
(d) the >=6-year NDVI contrast under the three clearance rules.

Inputs (repository data): v12/cross_scene/{spectral,distance,effects}.csv, 1m/ and 2m/effects.csv,
red_edge/e2_indices.csv (folder aux/ in the project data)
Output: Fig7.png/.pdf; plotted values figdata/fig6{a,b,c,d}_*.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *

BANDS = [('blue', 485, 450, 520, '#3B6FD8'), ('green', 555, 520, 590, '#2FA35A'), ('red', 660, 630, 690, '#D93A3A'), ('nir', 830, 770, 890, '#8A1C1C')]


CS = CS12


def panel_a(ax):
    sp = pd.read_csv(CS / 'spectral.csv', dtype={'tag': str})
    rows = []
    for tag in TAGS:
        for g, ls, mfc, lw in [('≥6', '-', None, 1.3), ('0–2', (0, (2.4, 1.6)), 'white', .9)]:
            z = sp[(sp.tag == tag) & (sp.group == g)].set_index('band')
            y = [z.loc[b, 'relative'] * 100 for b, *_ in BANDS]
            x = [c for _, c, *_ in BANDS]
            ax.plot(x, y, color=DATE[tag], lw=lw, ls=ls, zorder=3, alpha=1 if g == '≥6' else .75)
            ax.scatter(x, y, s=14 if g == '≥6' else 10, marker=DATE_MARK[tag], facecolor=DATE[tag] if mfc is None else 'white', edgecolor=DATE[tag], linewidths=.8, zorder=4)
            rows.append(dict(tag=tag, group=g, **{b: v for (b, *_), v in zip(BANDS, y)}))
    for b, c, lo, hi, col in BANDS:
        ax.axvspan(lo, hi, ymin=0, ymax=.075, color=col, alpha=.92, lw=0, zorder=2)
        ax.text(c, -22.6, b if b != 'nir' else 'NIR', ha='center', va='center', fontsize=5.8, color='white', zorder=5)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
    ax.set_xlim(440, 900); ax.set_ylim(-24, 14)
    ax.set_xticks([500, 600, 700, 800, 900]); ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Inter-row − control reflectance (%)')
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    r08 = [r for r in rows if r['tag'] == '08' and r['group'] == '≥6'][0]
    ax.annotate('visible darkening,\nstrongest in red', xy=(660, r08['red']), xytext=(705, -19.2), fontsize=5.8, color=SUB, ha='left', va='center', linespacing=1.0,
                arrowprops=dict(arrowstyle='-', color='#9AA5B1', lw=.5))
    ax.text(822, -12.0, 'NIR almost\nunchanged', fontsize=5.8, color=SUB, ha='left', va='center', linespacing=1.0)  # inside the axes
    hh = [Line2D([], [], color=DATE[t], marker=DATE_MARK[t], ms=3.4, lw=1.2, label=DATE_SHORT[t]) for t in TAGS]
    hh += [Line2D([], [], color=INK, lw=1.2, label='≥6 y'), Line2D([], [], color=INK, lw=.9, ls=(0, (2.4, 1.6)), marker='o', mfc='white', ms=3, label='0–2 y')]
    ax.legend(handles=hh, loc='upper center', fontsize=5.8, ncol=5, handletextpad=.25, columnspacing=.6, borderaxespad=.15, labelspacing=.25, handlelength=1.5, bbox_to_anchor=(.5, 1.0))
    pd.DataFrame(rows).to_csv(FIGDATA / 'fig6a_relative_reflectance.csv', index=False)


def panel_b(ax):
    e2 = pd.read_csv(P.aux12 / 'e2_indices.csv', dtype={'date': str})
    q = e2[(e2.age == '≥6') & (e2.ctrl == 'matched')].copy()
    q['doy'] = pd.to_datetime(q.date).dt.dayofyear
    out = []
    for k, col, mk, lab, off in [('ndre', '#9B2C5E', 'o', 'NDRE', -2.2), ('cire', '#E07B39', 'D', 'CI$_{re}$', 2.2)]:
        z = q[q['index'] == k].sort_values('doy')
        ax.plot(z.doy + off, z.est, color=col, lw=.7, alpha=.5, zorder=2)
        ax.errorbar(z.doy + off, z.est, yerr=[z.est - z.lo, z.hi - z.est], fmt=mk, ms=3.8, mfc=col, mec='white', mew=.4, ecolor=col, elinewidth=1.0, capsize=1.3, zorder=3, label=lab)
        out.append(z)
    pd.concat(out).to_csv(FIGDATA / 'fig6b_red_edge.csv', index=False)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
    ax.set_xticks([148, 179, 202, 236]); ax.set_xticklabels(['05-28', '06-28', '07-21', '08-24'])
    ax.set_xlim(138, 246); ax.set_ylim(-.004, .06)
    ax.set_ylabel('Index difference\nvs matched controls')  # two lines fit the axes height
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    ax.legend(loc='upper left', fontsize=6.0, handletextpad=.2, borderaxespad=.2, labelspacing=.3)
    ax.text(.98, .97, f'≥6-year units, {int(q[q["index"] == "ndre"].n_cells.min())}–{int(q[q["index"] == "ndre"].n_cells.max())} per date', transform=ax.transAxes, ha='right',
            va='top', fontsize=5.7, color=MUTED)  # top right, clear of the NDRE series


def panel_c(ax):
    dd = pd.read_csv(CS / 'distance.csv', dtype={'tag': str})
    q = dd[(dd.tag == '08') & dd.group.isin(GROUPS) & (dd.n >= 30) & (dd.clusters >= 6)]
    q = q[q.n >= .5 * q.groupby('group').n.transform('max')]  # bands sampled by at least half of each age group
    ax.axvspan(-.75, 0, color='#2E3B48', lw=0, zorder=1)
    ax.text(-.375, .118, 'row\nband', color='white', fontsize=5.7, ha='center', va='top', linespacing=1.0, zorder=6)
    for j, (g, col) in enumerate([('0–2', AGE['0–2']), ('3–5', AGE['3–5']), ('≥6', AGE['≥6'])]):
        z = q[q.group == g].sort_values('d_lo')
        x = (z.d_lo + z.d_hi) / 2 + (j - 1) * .05
        ax.errorbar(x, z.estimate, yerr=[z.estimate - z.lo, z.hi - z.estimate], fmt='o-', color=col, ms=3.6, mfc=col, mec='white', mew=.4, lw=1.2,
                    elinewidth=.9, capsize=1.2, zorder=3, label=f'{g} y')
    q.to_csv(FIGDATA / 'fig6c_distance_08.csv', index=False)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
    ax.set_xlim(-.75, 3.1); ax.set_ylim(-.004, .12)
    for xb in [1, 1.5, 2, 3]:
        ax.axvline(xb, color='#E4E7EB', lw=.5, zorder=0)
    ax.text(.98, .84, 'bands 1–1.5, 1.5–2 and 2–3 m', transform=ax.transAxes, fontsize=5.6, color=MUTED, va='top', ha='right')
    ax.set_xticks([0, 1, 1.5, 2, 3]); ax.set_xticklabels(['0', '1', '1.5', '2', '3'])
    ax.set_xlabel('Distance of footprint centre from the row band (m)'); ax.set_ylabel(f'ΔNDVI, {DATE_LABEL["08"]}')
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    ax.legend(loc='upper right', fontsize=5.9, handletextpad=.3, borderaxespad=.2, labelspacing=.25, ncol=3, columnspacing=.8)


def panel_d(ax):
    rows = []
    for j, tag in enumerate(TAGS):
        e = []
        for p in ['1m', '1p5m', '2m']:
            f = CS / 'effects.csv' if p == '1p5m' else CS / p / 'effects.csv'
            z = pd.read_csv(f, dtype={'tag': str, 'group': str})
            e.append(z[(z.tag == tag) & (z['index'] == 'ndvi') & (z.group == '≥6')].iloc[0])
        e = pd.DataFrame(e)
        x = np.arange(3) + (j - 1) * .17
        ax.plot(x, e.estimate, color=DATE[tag], lw=.9, alpha=.6, zorder=2)
        for i, (p, r) in enumerate(zip(['1m', '1p5m', '2m'], e.itertuples())):
            ax.errorbar(x[i], r.estimate, yerr=[[r.estimate - r.joint_plus_calibration_lo], [r.joint_plus_calibration_hi - r.estimate]], fmt=DATE_MARK[tag], ms=4,
                        mfc=DATE[tag] if p == '1p5m' else 'white', mec=DATE[tag], mew=.9, ecolor=DATE[tag], elinewidth=.9, capsize=1.3, zorder=3)
            rows.append(dict(tag=tag, clearance=p, estimate=r.estimate, lo=r.joint_plus_calibration_lo, hi=r.joint_plus_calibration_hi, n=int(r.n_cells)))
    pd.DataFrame(rows).to_csv(FIGDATA / 'fig6d_clearance.csv', index=False)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
    ax.set_xticks(range(3)); ax.set_xticklabels(['≥1 m', '≥1.5 m\nprimary', '≥2 m'])
    ax.set_xlim(-.5, 2.5); ax.set_ylim(-.004, .1)
    ax.set_xlabel('Centre clearance rule'); ax.set_ylabel('ΔNDVI, ≥6-year arrays')
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    hh = [Line2D([], [], color=DATE[t], marker=DATE_MARK[t], ms=3.6, lw=.9, label=DATE_SHORT[t]) for t in TAGS]
    ax.legend(handles=hh, loc='lower left', fontsize=5.9, handletextpad=.3, borderaxespad=.2, labelspacing=.25, ncol=3, columnspacing=.8)


def main():
    fig = figure(110)
    axa = fig.add_axes([.075, .585, .4, .35]); axb = fig.add_axes([.585, .585, .4, .35])
    axc = fig.add_axes([.075, .115, .4, .35]); axd = fig.add_axes([.585, .115, .4, .35])
    panel_a(axa); panel_b(axb); panel_c(axc); panel_d(axd)
    for ax, l, t in [(axa, 'a', 'Spectral shape of the contrast (2 m)'), (axb, 'b', 'Red-edge contrasts (1.2 m, 2025)'),
                     (axc, 'c', 'The whole gap is greener, not only its edges'), (axd, 'd', 'The contrast under three clearance rules')]:
        letter(ax, l, x=-.16, y=1.03, ha='left')
        ax.text(-.09, 1.03, t, transform=ax.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .18))
    save(fig, 'Fig7')


if __name__ == '__main__':
    main()
