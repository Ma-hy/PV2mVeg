"""Fig. S12: mixed-pixel rainfall diagnostic.

Change of the post-installation NDVI and NIRv contrasts of mixed Landsat pixels per standard deviation of
growing-season precipitation: the observed coefficient beta_obs with its year-jackknife 95% interval
(beta_obs +- 1.96 x se_year_jackknife) and the reflectance-level panel-dilution null null_beta_mid
(rows outcome == 'raw' and e_min == 1 of the frozen rainfall table).

Input (repository data): frozen/rainfall_frozen.csv
Output: FigS12.png/.pdf; plotted values figdata/figS12_rainfall.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *

SRC = P.frozen / 'rainfall_frozen.csv'
NULL = '#D1495B'   # reference expectation without an ecological response (as the 100% line of Fig. 8)


def main():
    raw = pd.read_csv(SRC)
    q = raw[(raw.outcome == 'raw') & (raw.e_min == 1)].set_index('index').reindex(['ndvi', 'nirv'])
    assert len(q) == 2 and q.notna().beta_obs.all(), q
    half = 1.96 * q.se_year_jackknife
    out = pd.DataFrame(dict(index=q.index, beta_obs=q.beta_obs.values, lo=(q.beta_obs - half).values, hi=(q.beta_obs + half).values,
                            se_year_jackknife=q.se_year_jackknife.values, null_beta_mid=q.null_beta_mid.values, n_years=q.n_years.values))
    out.to_csv(FIGDATA / 'figS12_rainfall.csv', index=False)

    fig = figure(66)
    ax = fig.add_axes([.105, .12, .865, .73])
    x = np.arange(2)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    ax.plot(x, q.null_beta_mid, ls='', marker='_', ms=22, mew=2.4, color=NULL, zorder=3)
    eb = ax.errorbar(x, q.beta_obs, yerr=half, fmt='o', ms=5, mfc=MIXED, mec='white', mew=.6, ecolor=MIXED,
                     elinewidth=1.3, capsize=3, capthick=1.1, zorder=4, label='Observed; year-jackknife 95% interval')
    ax.set_xticks(x); ax.set_xticklabels(['NDVI', 'NIRv'])
    ax.set_xlim(-.6, 1.6)
    ax.set_ylabel('Change per precipitation SD')
    ax.grid(axis='y', color=GRID, lw=.45, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=2.4, width=.5, pad=1.6)
    handles = [eb, Line2D([], [], color=NULL, ls='', marker='_', ms=14, mew=2.4, label='Panel-dilution null')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .985), ncol=2, frameon=False, handletextpad=.5, columnspacing=2.2)
    save(fig, 'FigS12')
    print(out.round(5).to_string())


if __name__ == '__main__':
    main()
