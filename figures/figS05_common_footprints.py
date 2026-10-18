"""Fig. S5: change on common ground footprints between acquisitions, by clearance rule and age group.

Input (repository data): v12/cross_scene/secondary/crossperiod/{06_08,08_2026_08}/summary.csv
Output: FigS5.png/.pdf; plotted values figdata/supp_diag/figS8_common_by_clearance.csv
"""
from _diag import *


def common_panels(fig, axs, kind):
    rows = []
    for j, (ax, pair) in enumerate(zip(axs, PAIRS)):
        d = rd(SECD / f'crossperiod/{pair[0]}_{pair[1]}/summary.csv')
        if kind == 'clearance':
            series = [(p, d[(d.clearance == p) & (d['index'] == 'ndvi') & (d.quantity == 'common_change')], CLR[p]) for p in POL]
        else:
            series = [(q, d[(d.clearance == '1p5m') & (d['index'] == 'ndvi') & (d.quantity == q)], QCOL[q]) for q in QCOL]
        for i, (name, q, col) in enumerate(series):
            q = q.set_index('group').reindex(GROUPS)
            x = np.arange(3) + (i - 1) * .2
            pr(ax, x, q.estimate, q.lo, q.hi, col, ms=3.8)
            for xx, n in zip(x, q.n):
                if np.isfinite(n):
                    ax.text(xx, -.02, f'{int(n):,}', transform=ax.get_xaxis_transform(), rotation=90, fontsize=4.9, color=col, ha='center', va='top')
            rows.append(q.reset_index().assign(pair=PAIR_LAB[pair], series=name))
        tidy(ax); ax.set_xlim(-.5, 2.5); ax.set_xticks(range(3)); ax.set_xticklabels([f'{g} years' for g in GROUPS])
        ax.tick_params(axis='x', pad=21)  # below the unit counts
        ax.set_xlabel('Age at the first acquisition', labelpad=1)
        if j == 0:
            ax.set_ylabel('NDVI change between acquisitions')
        head(ax, f'{"ab"[j]}   {PAIR_LAB[pair]}')
    return pd.concat(rows)


def main():
    fig = figure(70)
    axs = [fig.add_axes([.085 + i * .47, .25, .4, .6]) for i in range(2)]
    rows = common_panels(fig, axs, 'clearance')
    keep('figS8_common_by_clearance', rows)
    legend_row(fig, [H(CLR[p], CLAB[p]) for p in POL], .985)
    save(fig, 'FigS5')


if __name__ == '__main__':
    main()
