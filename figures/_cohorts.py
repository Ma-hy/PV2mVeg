"""Shared data and drawing of Fig. 4 and Fig. S7: single-year cohorts of the inter-row ground, the mixed Landsat
pixels of the same units, the mixed record of units followed since installation, and the service lanes.

Mixed-pixel cross-sections are resampled by whole station (2,000 draws, controls fixed); the inter-row values and
their joint sampling-plus-calibration intervals are read from cross_scene/effects.csv.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *

INDEX = ['ndvi', 'nirv', 'msavi']
SEED = 20261011


def station_boot(values, stations, B=2000, seed=SEED):
    values = np.asarray(values, float); stations = np.asarray(stations)
    ok = np.isfinite(values); values, stations = values[ok], stations[ok]
    ids = np.unique(stations)
    if len(ids) < 6:
        return float(values.mean()), np.nan, np.nan
    rng = np.random.default_rng(seed)
    sums = np.array([values[stations == s].sum() for s in ids]); cnt = np.array([(stations == s).sum() for s in ids])
    draws = rng.integers(0, len(ids), size=(B, len(ids)))
    means = sums[draws].sum(1) / cnt[draws].sum(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(values.mean()), float(lo), float(hi)


def group_of(age):
    return np.where(age <= 2, '0–2', np.where(age <= 5, '3–5', '≥6'))


def data():
    """Primary population: all v1.2 units with >=50 footprints and >=2 original controls in any scene."""
    CS = CS12
    eff = pd.read_csv(CS / 'effects.csv', dtype={'tag': str, 'group': str})
    rows = []
    for t in TAGS:
        d = pd.read_csv(CS / f'cells_{t}.csv').query('valid_cross')
        blk = pd.read_csv(SEC12 / '1p5m/chrono' / f'a1_{t}_det_thr50_cells.csv')[['cell', 'blk']]
        d = d.merge(blk, on='cell', how='left')
        st = np.where(d.station > 0, d.station, 10 ** 9 + d.blk).astype(np.int64)
        for k in INDEX:
            for g in [str(a_) for a_ in range(15)] + GROUPS:
                m = (group_of(d.age.values) == g) if g in GROUPS else (d.age.values == int(g))
                if not m.any():
                    continue
                e = eff[(eff.tag == t) & (eff['index'] == k) & (eff.group == g)].iloc[0]
                assert int(e.n_cells) == int(m.sum()) and abs(e.estimate - d.loc[m, f'cross_{k}_dadj'].mean()) < 1e-10, (t, k, g)
                mm, lo, hi = station_boot(d.loc[m, f'cross_{k}_mixed'].values, st[m])
                rows.append(dict(tag=t, index=k, age_group=g, n=int(m.sum()), stations=int(e.n_stations), ground=e.estimate,
                                 ground_lo=e.joint_plus_calibration_lo if e.n_stations >= 6 else np.nan,
                                 ground_hi=e.joint_plus_calibration_hi if e.n_stations >= 6 else np.nan, mixed=mm, mixed_lo=lo, mixed_hi=hi))
    rows = pd.DataFrame(rows)
    ln = pd.read_csv(CS / 'lanes.csv', dtype={'tag': str, 'group': str})
    lanes = ln[ln.group.isin(GROUPS)].rename(columns={'group': 'age_group'})[['tag', 'index', 'age_group', 'lane', 'lane_lo', 'lane_hi', 'lane_n', 'lane_clusters']]
    curves = {}
    for k in INDEX:
        c = pd.read_csv(P.es / f'es_raw_sum_{k}_det_treated' / 'curves.csv')
        c = c[(c.estimator == 'diff') & (c.weight == 'cell') & (c.subset == 'balanced T≤2018')].sort_values('e')
        curves[k] = c
    fixed = curves['ndvi']
    assert (fixed[fixed.e.between(0, 8)].n_cells == 2304).all() and fixed[fixed.e == 9].n_cells.iloc[0] < 2304
    rows.to_csv(FIGDATA / 'fig3_cohorts.csv', index=False); lanes.to_csv(FIGDATA / 'fig3_lanes.csv', index=False)
    pd.concat([c.assign(index=k) for k, c in curves.items()]).to_csv(FIGDATA / 'fig3_longitudinal.csv', index=False)
    return rows, lanes, curves


def surface_panel(ax, k, rows, lanes, curves, tag='08', big=True):
    ms = 4.2 if big else 3.3
    ax.axvspan(-1.5, .5, color='#EEF1F4', lw=0, zorder=0)
    ax.text(-.5, .99, 'construction', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=5.8 if big else 5.4, color=MUTED, rotation=90 if not big else 0)
    for lo, hi in [(2.5, 2.5), (5.5, 5.5)]:
        ax.axvline(lo, color='#DDE2E7', lw=.6, zorder=0)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    # longitudinal mixed record of fixed units
    c = curves[k]
    a = c[c.e.between(-4, 8)]; b = c[c.e.between(8, 12)]
    ax.fill_between(a.e, a.lo, a.hi, color=MIXED, alpha=.16, lw=0, zorder=2)
    ax.plot(a.e, a.est, color=MIXED, lw=1.25 if big else 1.05, zorder=3, solid_capstyle='round')
    ax.fill_between(b.e, b.lo, b.hi, color=MIXED, alpha=.07, lw=0, zorder=2)
    ax.plot(b.e, b.est, color=MIXED, lw=.9, ls=(0, (2.2, 1.6)), zorder=3, alpha=.8)
    # lanes, age groups
    L = lanes[(lanes.tag == tag) & (lanes['index'] == k)].set_index('age_group')
    for g, (x0, x1) in zip(GROUPS, [(-.35, 2.35), (2.65, 5.35), (5.65, 13.4)]):
        r = L.loc[g]
        ax.fill_between([x0, x1], [r.lane_lo] * 2, [r.lane_hi] * 2, color=LANE, alpha=.17, lw=0, zorder=2)
        ax.plot([x0, x1], [r.lane] * 2, color=LANE, lw=1.4, solid_capstyle='butt', zorder=3)
    # cross-sections of the acquisition
    R = rows[(rows.tag == tag) & (rows['index'] == k) & ~rows.age_group.isin(GROUPS)].copy()
    R['age'] = R.age_group.astype(int); R = R.sort_values('age')
    ok = R.mixed_lo.notna()
    ax.errorbar(R.age[ok] + .2, R.mixed[ok], yerr=[R.mixed[ok] - R.mixed_lo[ok], R.mixed_hi[ok] - R.mixed[ok]], fmt='s', ms=ms * .78,
                mfc='white', mec=MIXED, mew=.8, ecolor=MIXED, elinewidth=.6, capsize=0, zorder=4, lw=0)
    ax.plot(R.age[~ok] + .2, R.mixed[~ok], 's', ms=ms * .78, mfc='white', mec=MIXED, mew=.6, alpha=.6, zorder=4)
    ok = R.ground_lo.notna()
    ax.errorbar(R.age[ok] - .12, R.ground[ok], yerr=[R.ground[ok] - R.ground_lo[ok], R.ground_hi[ok] - R.ground[ok]], fmt='o', ms=ms,
                mfc=GROUND, mec='white', mew=.45, ecolor=GROUND, elinewidth=.9, capsize=1.3 if big else 1.0, zorder=6, lw=0)
    ax.plot(R.age[~ok] - .12, R.ground[~ok], 'o', ms=ms, mfc='white', mec=GROUND, mew=.9, zorder=6)
    ax.set_xlim(-4.6, 13.8)
    ax.set_xticks(range(-4, 14, 2))
    ax.grid(axis='y', color=GRID, lw=.45, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=2.4, width=.5, pad=1.6)
    return R
