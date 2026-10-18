"""Fig. S6: mixed minus inter-row contrast by the mapped row-band fraction of the 30 m pixel.

Input (repository data): v12/cross_scene/secondary/1p5m/statistics/mixed_bias.csv
Output: FigS6.png/.pdf; plotted values figdata/supp_diag/figS23_mixed_minus_interrow.csv
"""
from _diag import *


def main():
    d = rd(SECD / '1p5m/statistics/mixed_bias.csv')
    d = d[d['index'].isin(['ndvi', 'nirv', 'msavi']) & d.pf_lo.round(1).isin([.3, .4, .5])]
    fig = figure(66)
    axs = [fig.add_axes([.075 + i * .318, .27, .262, .58]) for i in range(3)]
    for i, (ax, k) in enumerate(zip(axs, ['ndvi', 'nirv', 'msavi'])):
        for j, tag in enumerate(TAGS):
            q = d[(d.tag == tag) & (d['index'] == k)].sort_values('pf_lo')
            x = np.arange(len(q)) + (j - 1) * .16
            ax.plot(x, q.estimate, color=DATE[tag], lw=.7, alpha=.7, zorder=2)
            pr(ax, x, q.estimate, q.lo, q.hi, DATE[tag], marker=DATE_MARK[tag], ms=3.4)
        tidy(ax); ax.set_xticks(range(3)); ax.set_xticklabels(['30–40%', '40–50%', '50–60%']); ax.set_xlim(-.5, 2.5)
        ax.set_xlabel('Mapped row-band fraction' + '\n' + 'of the 30 m pixel')  # two lines, apart from the neighbours
        if i == 0:
            ax.set_ylabel('Mixed − inter-row Δ')
        head(ax, f'{"abc"[i]}   {NAME[k]}')
    keep('figS23_mixed_minus_interrow', d)
    cnt = d[d['index'] == 'ndvi']
    txt = '   '.join(f'{DATE_SHORT[t]}: ' + ' / '.join(f'{int(n):,}' for n in cnt[cnt.tag == t].sort_values('pf_lo').n) for t in TAGS)
    fig.text(.075, .06, 'Units per bin (30–40% / 40–50% / 50–60%)   ' + txt, fontsize=6.0, color=SUB)
    legend_row(fig, [H(DATE[t], DATE_LABEL[t], marker=DATE_MARK[t]) for t in TAGS], .985)
    save(fig, 'FigS6')


if __name__ == '__main__':
    main()
