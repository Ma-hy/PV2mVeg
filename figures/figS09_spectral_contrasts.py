"""Fig. S9: visible and near-infrared spectral contrasts of the inter-row ground (>=1.5 m).

(a-c) band contrasts as a share of the control reflectance, (d-f) visible-index differences, 0-2 and >=6 years.
Input (repository data): v12/cross_scene/secondary/measurements/<tag>/subset_summary.csv
Output: FigS9.png/.pdf; plotted values figdata/supp_diag/figS21_spectral.csv
"""
from _diag import *


def main():
    fig = figure(98)
    rows = []
    bands = ['blue', 'green', 'red', 'nir']; blab = ['Blue', 'Green', 'Red', 'NIR']
    vis = ['grvi', 'vari', 'exg']; vlab = ['GRVI', 'VARI', 'ExG']
    top = [fig.add_axes([.075 + i * .315, .575, .262, .304]) for i in range(3)]  # top 1.6 mm below the key
    bot = [fig.add_axes([.075 + i * .315, .1, .262, .33]) for i in range(3)]
    for i, tag in enumerate(TAGS):
        s = rd(SECD / f'measurements/{tag}/subset_summary.csv')
        s = s[(s.subset == 'ir_1p5m') & (s.quantity == 'direct_spectral_difference')]
        ax = top[i]
        for j, g in enumerate(['0–2', '≥6']):
            q = s[(s.group == g) & s['index'].isin(bands)].set_index('index').reindex(bands)
            x = np.arange(4) + (j - .5) * .12
            rel = q.estimate / q.control * 100; lo = q.lo / q.control * 100; hi = q.hi / q.control * 100
            ax.plot(x, rel, color=AGE[g], lw=.8, alpha=.75, zorder=2)
            pr(ax, x, rel, lo, hi, AGE[g], ms=3.6)
            rows.append(q.reset_index().assign(tag=tag, panel='bands'))
        tidy(ax); ax.set_xticks(range(4)); ax.set_xticklabels(blab); ax.set_xlim(-.5, 3.5)
        if i == 0:
            ax.set_ylabel('Inter-row − control (% of control)')
        n6 = int(s[(s.group == '≥6') & (s['index'] == 'red')].n.iloc[0])
        head(ax, f'{"abc"[i]}   {DATE_LABEL[tag]}', f'≥6 years: {n6:,} units with visible calibration')
        ax = bot[i]
        for j, g in enumerate(['0–2', '≥6']):
            q = s[(s.group == g) & s['index'].isin(vis)].set_index('index').reindex(vis)
            pr(ax, np.arange(3) + (j - .5) * .16, q.estimate, q.lo, q.hi, AGE[g], ms=3.6)
            rows.append(q.reset_index().assign(tag=tag, panel='visible'))
        tidy(ax); ax.set_xticks(range(3)); ax.set_xticklabels(vlab); ax.set_xlim(-.5, 2.5)
        if i == 0:
            ax.set_ylabel('Visible-index difference')
        head(ax, f'{"def"[i]}   {DATE_LABEL[tag]}')
    for grp in (top, bot):
        lim = [min(a.get_ylim()[0] for a in grp), max(a.get_ylim()[1] for a in grp)]
        for a in grp:
            a.set_ylim(lim)
    keep('figS21_spectral', pd.concat(rows))
    legend_row(fig, [H(AGE['0–2'], '0–2 years'), H(AGE['≥6'], '≥6 years')], .995)
    save(fig, 'FigS9')


if __name__ == '__main__':
    main()
