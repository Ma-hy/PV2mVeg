"""Fig. S11: imaging-scenario and boundary sensitivity of the >=6-year contrast.

(a-c) retained imaging scenarios: simulated range and observed mean at the three clearances;
(d-f) observed minus simulated NDVI for the retained, with-sharpening and all candidate scenario sets (>=1.5 m);
(g-i) paired NDVI change under 0.5 m perturbations of the mapped row-band boundary.
Inputs (repository data): v12/cross_scene/optics.csv, v12/cross_scene/secondary/measurements/<tag>/boundary_paired.csv
Output: FigS11.png/.pdf; plotted values figdata/supp_diag/figS22_optics_boundary.csv
"""
from _diag import *

SETS = [('restricted_common', 'Retained'), ('with_sharpening_common', 'With sharpening'), ('all_common', 'All candidates')]
SETC = ['#0E8A6A', '#C98A1B', '#8E3BD8']
BOUND = ['wid-0.5', 'wid+0.5', 'shiftN0.5', 'shiftS0.5', 'Nedge-0.5', 'Nedge+0.5', 'Sedge-0.5', 'Sedge+0.5']
BLAB = ['Band 0.5 m narrower', 'Band 0.5 m wider', 'Band shifted 0.5 m north', 'Band shifted 0.5 m south',
        'North edge 0.5 m in', 'North edge 0.5 m out', 'South edge 0.5 m in', 'South edge 0.5 m out']


def envelope(o, p, tag, key):
    q = o[(o.design == 'primary') & (o.clearance == p) & (o.tag == tag) & (o.age == '≥6') & o.sets.str.contains(key, na=False)]
    obs = q.obs.iloc[0]
    return dict(models=len(q), art_min=q.art.min(), art_max=q.art.max(), obs=obs, net_min=(obs - q.art).min(), net_max=(obs - q.art).max(), n=int(q.n.min()))


def main():
    o = rd(CSD / 'optics.csv')
    fig = figure(150)
    W, X0, DX = .235, .205, .275
    A = [fig.add_axes([X0 + j * DX, .725, W, .19]) for j in range(3)]
    B = [fig.add_axes([X0 + j * DX, .43, W, .19]) for j in range(3)]
    C = [fig.add_axes([X0 + j * DX, .06, W, .265]) for j in range(3)]  # top clear of the labels of (d-f)
    # the subtitle of the first column describes the row; the repeats in the other columns are left out, so that
    # the first subtitle can run under the next heading
    rep = lambda j, s: s if not j else ''
    rows = []
    for j, tag in enumerate(TAGS):
        ax = A[j]
        for i, p in enumerate(POL):
            e = envelope(o, p, tag, 'restricted_common'); y = 2 - i
            ax.plot([e['art_min'], e['art_max']], [y, y], color=CLR[p], lw=5.2, solid_capstyle='butt', alpha=.85, zorder=2)
            ax.plot(e['obs'], y, 'o', ms=4.2, color=INK, mec='white', mew=.6, zorder=4)
            rows.append(dict(panel='a-c', tag=tag, clearance=p, set='retained', **e))
        tidy(ax, grid='x', zero=False); ax.axvline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
        ax.set_xlim(-.03, .09); ax.set_ylim(-.6, 2.6); ax.set_yticks([2, 1, 0])
        ax.set_yticklabels([CSHORT[p] for p in POL] if j == 0 else []); ax.tick_params(axis='y', length=0)
        ax.set_xticks([0, .03, .06]); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, 2) if v else '0'))
        ax.set_xlabel('ΔNDVI', labelpad=1)
        head(ax, f'{"abc"[j]}   {DATE_LABEL[tag]}', rep(j, 'retained scenarios: simulated range, observed mean' if not j else 'retained scenarios'))
        ax = B[j]
        for i, (key, lab) in enumerate(SETS):
            e = envelope(o, '1p5m', tag, key); y = 2 - i
            ax.plot([e['net_min'], e['net_max']], [y, y], color=SETC[i], lw=5.2, solid_capstyle='butt', alpha=.85, zorder=2)
            # above the bar, from its start: right of the bar the label would leave the panel
            ax.text(max(e['net_min'], -.058), y + .3, f'{signed(e["net_min"])} to {signed(e["net_max"])}', fontsize=5.4, color=SETC[i], va='bottom',
                    path_effects=halo(1.4), zorder=5)
            rows.append(dict(panel='d-f', tag=tag, clearance='1p5m', set=lab, **e))
        tidy(ax, grid='x', zero=False); ax.axvline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
        ax.axvspan(-.06, 0, color='#A85B62', alpha=.06, lw=0)
        ax.set_xlim(-.06, .14); ax.set_ylim(-.6, 2.6); ax.set_yticks([2, 1, 0])
        ax.set_yticklabels([f'{lab} ({envelope(o, "1p5m", tag, key)["models"]})' for key, lab in SETS] if j == 0 else []); ax.tick_params(axis='y', length=0)
        ax.set_xticks([-.05, 0, .05, .1]); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, 2) if v else '0'))
        ax.set_xlabel('Observed − simulated NDVI', labelpad=1)
        head(ax, f'{"def"[j]}   {DATE_LABEL[tag]}', rep(j, '≥1.5 m; scenario counts in brackets' if not j else '≥1.5 m'))
        ax = C[j]
        b = rd(SECD / f'measurements/{tag}/boundary_paired.csv')
        b = b[(b['index'] == 'ndvi') & (b.group == '≥6')].set_index('subset').reindex(['boundary_' + v for v in BOUND])
        for i, (name, r) in enumerate(b.iterrows()):
            y = 7 - i
            pr(ax, y, r.estimate, r.lo, r.hi, '#6B5B2E', horizontal=True, ms=3.4)
            ax.text(.0145, y, f'{int(r.n):,}', ha='right', va='center', fontsize=5.3, color=SUB)
            rows.append(dict(panel='g-i', tag=tag, subset=name, estimate=r.estimate, lo=r.lo, hi=r.hi, n=int(r.n)))
        for y in [1.5, 3.5, 5.5]:
            ax.axhline(y, color=GRID, lw=.5)
        tidy(ax, grid='x', zero=False); ax.axvline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
        ax.set_xlim(-.0125, .0148); ax.set_ylim(-.6, 7.9); ax.set_yticks(np.arange(7, -1, -1))
        ax.set_yticklabels(BLAB if j == 0 else []); ax.tick_params(axis='y', length=0)
        ax.set_xticks([-.01, 0, .01]); ax.xaxis.set_minor_locator(MultipleLocator(.005)); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, 2) if v else '0'))
        ax.text(.0145, 7.42, 'units', ha='right', va='bottom', fontsize=5.3, color=SUB)  # inside the axes, under the subtitle
        ax.set_xlabel('Paired NDVI change', labelpad=1)
        head(ax, f'{"ghi"[j]}   {DATE_LABEL[tag]}', rep(j, '0.5 m boundary perturbations, arrays ≥6 years' if not j else '0.5 m boundary perturbations'))
    keep('figS22_optics_boundary', pd.DataFrame(rows))
    legend_row(fig, [Line2D([], [], color=CLR[p], lw=4.5, label=f'{CSHORT[p]} simulated range') for p in POL] + [H(INK, 'observed mean')], .993, size=6.3)
    save(fig, 'FigS11')


if __name__ == '__main__':
    main()
