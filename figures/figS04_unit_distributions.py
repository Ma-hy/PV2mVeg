"""Fig. S4: unit-level distributions behind the age-group means of the primary design.

Columns are acquisitions, rows ΔNDVI, ΔNIRv and ΔMSAVI; every point is one valid 150 m unit (deterministic jitter),
with a kernel density outline, the interquartile box, the median and the unit-weighted mean (diamond).
Inputs (repository data): v12/cross_scene/cells_{06,08,2026_08}.csv, v12/cross_scene/effects.csv
Output: FigS4.png/.pdf
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *


def main():
    from scipy.stats import gaussian_kde
    E = {t: pd.read_csv(CS12 / f'cells_{t}.csv').query('valid_cross') for t in TAGS}
    eff = pd.read_csv(CS12 / 'effects.csv', dtype={'tag': str, 'group': str})
    fig = figure(150)
    gs = fig.add_gridspec(3, 3, left=.08, right=.985, bottom=.06, top=.95, hspace=.42, wspace=.28)
    rng = np.random.default_rng(7)
    for i, k in enumerate(['ndvi', 'nirv', 'msavi']):
        for j, t in enumerate(TAGS):
            ax = fig.add_subplot(gs[i, j]); d = E[t]
            g = np.where(d.age <= 2, '0–2', np.where(d.age <= 5, '3–5', '≥6'))
            for x0, grp in enumerate(GROUPS):
                y = d[f'cross_{k}_dadj'].values[g == grp]; y = y[np.isfinite(y)]
                col = AGE[grp]
                if len(y) > 20:
                    yy = np.linspace(np.percentile(y, .5), np.percentile(y, 99.5), 160); kd = gaussian_kde(y)(yy)
                    ax.fill_betweenx(yy, x0 - .38 * kd / kd.max(), x0 + .38 * kd / kd.max(), color=col, alpha=.22, lw=.6, edgecolor=col)
                ax.scatter(x0 + rng.uniform(-.16, .16, len(y)), y, s=.6, color=col, alpha=.35, lw=0, rasterized=True)
                q1, md, q3 = np.percentile(y, [25, 50, 75])
                ax.add_patch(Rectangle((x0 - .07, q1), .14, q3 - q1, fc='white', ec=col, lw=.8, zorder=4))
                ax.plot([x0 - .07, x0 + .07], [md, md], color=col, lw=1.2, zorder=5)
                r = eff[(eff.tag == t) & (eff['index'] == k) & (eff.group == grp)].iloc[0]
                ax.plot(x0, r.estimate, 'D', ms=3.2, mfc=INK, mec='white', mew=.4, zorder=6)
                ax.text(x0, 1.0, f'n {int(r.n_cells):,}', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.3, color=SUB)
            ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)))
            ax.set_xticks(range(3)); ax.set_xticklabels(['0–2 y', '3–5 y', '≥6 y'], fontsize=6.2)
            lim = {'ndvi': (-.15, .22), 'nirv': (-.04, .055), 'msavi': (-.08, .11)}[k]; ax.set_ylim(*lim)
            ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2, pad=1.2, labelsize=6.0)
            if j == 0: ax.set_ylabel(f'Δ{NAME[k]}')
            if i == 0: ax.set_title(DATE_LABEL[t], fontsize=7.2, pad=9)
    save(fig, 'FigS4')


if __name__ == '__main__':
    main()
