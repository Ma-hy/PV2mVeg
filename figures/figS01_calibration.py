"""Fig. S1: Beijing-3 radiometric traceability and per-scene cross-calibration to Sentinel-2.   [restricted inputs]

Historical figure: drawn by the figure code of an earlier release on its original inputs; only the text on the
figure follows the v2.0 house style (pv2mveg.plotting.historical).

(a)-(d) blue, green, red, NIR: Level-4 DN vs top-of-atmosphere reflectance computed from Level-1 (gain/bias/ESUN, RPC
    orthorectification), 10 m cells of one MS1 SC_002 scene per acquisition; random subsample of the cells used in
    the fit, OLS line over the 1-99th percentile of TOA; coefficients from the result table.
(e) Per-scene Huber calibration coefficients S2 = a + b*DN (red, NIR), 40 scenes.
(f) Per-scene NDVI bias against Sentinel-2 (pure-ground 10 m cells): raw L4 (DN/1e4) vs calibrated (2 km
    checkerboard cross-validation), grouped by acquisition x satellite.  (g) Per-scene NDVI RMSE, same grouping.
(h) Per-scene cross-validated reflectance RMSE of the four calibrated bands (blue and green: BJ3A1 scenes only);
    bars = cell-weighted pooled RMSE of each group.
(i) 2025-08 -> 2026-08 NDVI change over stable ground, calibrated BJ-3 vs Sentinel-2: 10 m cells (density) and
    1 km block means (points).

Inputs (project layout only; derived from the Beijing-3 scenes): bj3cal/trace_l4_vs_toa_s2.csv,
cache/trace_*_L1.npz, coef_per_scene.csv, ndvi_cv_metrics.csv, s2_bg/bg_cal_metrics.csv, yoy_2025_2026_stable.csv,
cache/yoy_pairs.npz. Output: FigS1.png/.pdf; figdata/figS1fg_scene_ndvi_metrics.csv, figS1h_band_cv_rmse.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style_v04 import *
from pv2mveg.plotting import historical
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, LogNorm

P.need_project('Fig. S1')
BJ = str(P.bj3cal)
# Beijing-3 satellites (Tol "bright", colour-blind safe; distinct from the period colours C_PERIOD)
C_SAT = {'BJ3A1': '#4477AA', 'BJ3A4': '#EE6677', 'BJ3N3': '#228833'}
M_SAT = {'BJ3A1': 'o', 'BJ3A4': '^', 'BJ3N3': 'v'}
C_B4 = {'blue': '#4477AA', 'green': '#228833', 'red': C_BAND['red'], 'nir': C_BAND['nir']}
M_B4 = {'blue': 'o', 'green': '^', 'red': 's', 'nir': 'D'}
L_B4 = {'blue': 'Blue', 'green': 'Green', 'red': 'Red', 'nir': 'NIR'}
CASES = [('2025-06 L1', '2025-06'), ('2025-08 L1', '2025-08'), ('2026-08 L1', '2026-08')]
GROUPS = [('2025-06', 'BJ3A1'), ('2025-08', 'BJ3A1'), ('2025-08', 'BJ3A4'), ('2025-08', 'BJ3N3'),
          ('2026-08', 'BJ3A1')]
GLAB = ['A1', 'A1', 'A4', 'N3', 'A1']
TAG = {'2025-06': '06', '2025-08': '08', '2026-08': '2026_08'}
LIMS = {'blue': ((0.08, 0.21), (0, 2400)), 'green': ((0.08, 0.24), (300, 2800)),
        'red': ((0.05, 0.30), (200, 3400)), 'nir': ((0.12, 0.36), (1000, 3900))}


def trace_panel(ax, fit, band, rng):
    for case, per in CASES:
        z = np.load(os.path.join(BJ, 'cache', f'trace_{case.replace(" ", "_")}.npz'))
        ok = z['ok10']
        x = z[f'toa_{band}'][ok]
        y = z[f'l4_{band}'][ok]
        i = rng.choice(x.size, 2500, replace=False)
        ax.scatter(x[i], y[i], s=0.6, color=C_PERIOD[per], alpha=0.25, lw=0, rasterized=True, zorder=2)
    for k, (case, per) in enumerate(CASES):
        r = fit[(fit.case == case) & (fit.scale == 10) & (fit.pair == 'L4~TOA') & (fit.band == band)].iloc[0]
        g = np.linspace(r.x_p01, r.x_p99, 20)
        ax.plot(g, r.a + r.b * g, color=C_PERIOD[per], lw=1.0, zorder=3)
        ax.text(0.03, 0.97 - k * 0.09, f'$R^2$ = {r.r2:.3f}', transform=ax.transAxes, ha='left', va='top',
                color=C_PERIOD[per], fontsize=FS)
    ax.set_xlabel(f'{L_B4[band]} TOA reflectance')
    ax.set_ylabel(f'{L_B4[band]} Level-4 value (DN)')
    ax.set_xlim(*LIMS[band][0])
    ax.set_ylim(*LIMS[band][1])


def group_axis(ax, y_text=-0.13):
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GLAB, fontsize=FS)
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    tr = ax.get_xaxis_transform()
    for xc, x0, x1, t in ((0, None, None, 'Jun 25'), (2, 0.75, 3.25, 'Aug 2025'), (4, None, None, 'Aug 26')):
        ax.text(xc, y_text, t, transform=tr, ha='center', va='top', fontsize=FS)
        if x0 is not None:
            ax.plot([x0, x1], [y_text + 0.01, y_text + 0.01], transform=tr, color='0.4', lw=0.5, clip_on=False)


def main():
    rng = np.random.default_rng(7)
    fit = pd.read_csv(os.path.join(BJ, 'trace_l4_vs_toa_s2.csv'))
    coef = pd.read_csv(os.path.join(BJ, 'coef_per_scene.csv'))
    cv = pd.read_csv(os.path.join(BJ, 'ndvi_cv_metrics.csv'))
    bg = pd.read_csv(os.path.join(BJ, 's2_bg', 'bg_cal_metrics.csv'))
    yoy = pd.read_csv(os.path.join(BJ, 'yoy_2025_2026_stable.csv'))

    fig = new_fig(DOUBLE, 178)
    g1 = fig.add_gridspec(1, 4, left=0.075, right=0.985, bottom=0.735, top=0.975, wspace=0.50)
    g2 = fig.add_gridspec(1, 3, left=0.075, right=0.985, bottom=0.415, top=0.635, wspace=0.42)
    g3 = fig.add_gridspec(1, 2, left=0.075, right=0.985, bottom=0.075, top=0.315, wspace=0.30,
                          width_ratios=[1.9, 1.0])

    # (a)-(d) traceability
    for k, band in enumerate(('blue', 'green', 'red', 'nir')):
        ax = fig.add_subplot(g1[0, k])
        trace_panel(ax, fit, band, rng)
        panel_label(ax, 'abcd'[k], x=-0.27)
    handles = [Line2D([], [], color=C_PERIOD[p], marker='o', ms=2.5, lw=1.0, label=L_PERIOD[p].replace(' (BJ3A1)', ''))
               for _, p in CASES]
    fig.axes[3].legend(handles=handles, loc='lower right', handlelength=1.2, fontsize=FS, borderaxespad=0.1,
                       labelspacing=0.2)

    # (e) coefficients
    ax = fig.add_subplot(g2[0, 0])
    h = coef[coef.method == 'huber']
    for band, mfc in (('red', None), ('nir', 'none')):
        for sat in C_SAT:
            s = h[(h.band == band) & (h.sensor == sat)]
            ax.scatter(s.b * 1e4, s.a, s=9, marker=M_SAT[sat], facecolor=C_SAT[sat] if mfc is None else 'none',
                       edgecolor=C_SAT[sat], lw=0.6, zorder=3)
    ax.set_xlabel('Gain b × 10$^4$ (reflectance per 10$^4$ DN)')
    ax.set_ylabel('Offset a (reflectance)')
    zero_line(ax)
    hs = [Line2D([], [], ls='', marker=M_SAT[s], color=C_SAT[s], ms=3, label=s) for s in C_SAT]
    hs += [Line2D([], [], ls='', marker='s', color='0.3', ms=3, label='Red'),
           Line2D([], [], ls='', marker='s', mfc='none', color='0.3', ms=3, label='NIR')]
    ax.legend(handles=hs, loc='upper right', ncol=1, fontsize=FS, handletextpad=0.2, labelspacing=0.2)
    panel_label(ax, 'e', x=-0.2)

    # (f)(g) per-scene NDVI bias and RMSE vs S2
    m = cv[cv.main == True]  # noqa: E712
    out = []
    for k, (stat, lab) in enumerate((('bias', 'NDVI bias vs S2'), ('rmse', 'NDVI RMSE vs S2'))):
        ax = fig.add_subplot(g2[0, k + 1])
        for j, (per, sat) in enumerate(GROUPS):
            for sub, dx, filled in (('raw_dn_ndvi_vs_s2', -0.17, False), ('all_pure_cv', 0.17, True)):
                v = m[(m.period == per) & (m.sensor == sat) & (m.subset == sub)][stat].values
                xx = j + dx + rng.uniform(-0.07, 0.07, v.size)
                ax.scatter(xx, v, s=8, marker=M_SAT[sat], facecolor=C_SAT[sat] if filled else 'none',
                           edgecolor=C_SAT[sat], lw=0.6, zorder=3)
                ax.hlines(v.mean(), j + dx - 0.13, j + dx + 0.13, color='k', lw=0.9, zorder=4)
                if k == 0:
                    out.append(dict(period=per, sensor=sat, subset=sub, n_scenes=v.size,
                                    bias_mean=v.mean(), bias_min=v.min(), bias_max=v.max()))
                else:
                    out[[i for i, o in enumerate(out) if o['period'] == per and o['sensor'] == sat
                         and o['subset'] == sub][0]].update(rmse_mean=v.mean(), rmse_min=v.min(), rmse_max=v.max())
        group_axis(ax)
        ax.set_ylabel(lab)
        if k == 0:
            zero_line(ax)
            ax.set_ylim(-0.005, 0.062)
            hs = [Line2D([], [], ls='', marker='o', mfc='none', mec='0.3', ms=3, label='Raw L4 (DN/10$^4$)'),
                  Line2D([], [], ls='', marker='o', color='0.3', ms=3, label='Calibrated (spatial CV)'),
                  Line2D([], [], color='k', lw=0.9, label='Mean of scenes')]
            ax.legend(handles=hs, loc='upper right', fontsize=FS, labelspacing=0.2, borderaxespad=0.1)
        else:
            ax.set_ylim(0, 0.06)
        panel_label(ax, 'fg'[k], x=-0.25)
    pd.DataFrame(out).to_csv(os.path.join(FIGDATA, 'figS1fg_scene_ndvi_metrics.csv'), index=False)

    # (h) per-scene reflectance CV RMSE, four bands
    ax = fig.add_subplot(g3[0, 0])
    boff = dict(zip(('blue', 'green', 'red', 'nir'), np.linspace(-0.27, 0.27, 4)))
    rows = []
    for j, (per, sat) in enumerate(GROUPS):
        for band in ('blue', 'green', 'red', 'nir'):
            if band in ('blue', 'green'):
                if sat != 'BJ3A1':
                    continue
                s = bg[(bg.tag.astype(str) == TAG[per]) & (bg.band == band) & bg.sid.str.startswith(sat)]
                v, n = s.cv_rmse.values, s.n.values
            else:
                s = m[(m.period == per) & (m.sensor == sat) & (m.subset == f'{band}_refl_cv')]
                v, n = s.rmse.values, s.n.values
            pooled = float(np.sqrt(np.sum(n * v ** 2) / np.sum(n)))
            x = j + boff[band]
            ax.scatter(x + rng.uniform(-0.04, 0.04, v.size), v, s=7, marker=M_B4[band], facecolor='none',
                       edgecolor=C_B4[band], lw=0.6, zorder=3)
            ax.hlines(pooled, x - 0.08, x + 0.08, color=C_B4[band], lw=1.4, zorder=4)
            rows.append(dict(period=per, sensor=sat, band=band, n_scenes=v.size, rmse_pooled=pooled,
                             rmse_min=v.min(), rmse_max=v.max()))
    pd.DataFrame(rows).to_csv(os.path.join(FIGDATA, 'figS1h_band_cv_rmse.csv'), index=False)
    group_axis(ax)
    ax.set_ylim(0, 0.0225)
    ax.set_ylabel('Reflectance RMSE vs S2')
    hs = [Line2D([], [], ls='', marker=M_B4[b], mfc='none', mec=C_B4[b], ms=3.2, label=L_B4[b])
          for b in ('blue', 'green', 'red', 'nir')]
    hs += [Line2D([], [], color='0.3', lw=1.4, label='Pooled over scenes')]
    ax.legend(handles=hs, loc='upper left', ncol=5, fontsize=FS, handletextpad=0.2, columnspacing=0.9,
              borderaxespad=0.1)
    panel_label(ax, 'h', x=-0.085)

    # (i) year-on-year change over stable ground
    ax = fig.add_subplot(g3[0, 1])
    z = np.load(os.path.join(BJ, 'cache', 'yoy_pairs.npz'))
    cmap = LinearSegmentedColormap.from_list('g', ['#F0F0F0', '#202020'])
    ax.hexbin(z['d_s2'], z['d_bj'], gridsize=70, extent=(-0.15, 0.08, -0.15, 0.08), mincnt=20,
              cmap=cmap, norm=LogNorm(), lw=0, rasterized=True, zorder=1)
    ax.plot([-0.15, 0.08], [-0.15, 0.08], color='k', lw=0.5, ls='--', zorder=2)
    ax.scatter(z['ms'], z['mb'], s=3, facecolor='#D55E00', edgecolor='none', lw=0, zorder=3, alpha=0.9)
    ax.set_xlim(-0.15, 0.08); ax.set_ylim(-0.15, 0.08)
    ax.set_aspect('equal')
    ax.set_xticks([-0.15, -0.1, -0.05, 0, 0.05])
    ax.set_xlabel('S2 ΔNDVI, Aug 26 − Aug 25')
    ax.set_ylabel('BJ-3 ΔNDVI (calibrated)')
    r10 = yoy[(yoy.method == 'calibrated_cv') & (yoy.scale == '10m')].iloc[0]
    r1k = yoy[(yoy.method == 'calibrated_cv') & (yoy.scale == '1km_block')].iloc[0]
    ax.text(0.98, 0.03, f'10 m: RMSE {r10.err_rmse:.3f}\n'
            f'1 km: RMSE {r1k.err_rmse:.4f}\n'
            f'$R^2$ = {r1k.err_r2:.2f}, bias {fmt(r1k.err_bias, 4)}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=FS, zorder=5,
            bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85))
    hs = [Line2D([], [], ls='', marker='h', color='0.45', ms=4, label='10 m cells'),
          Line2D([], [], ls='', marker='o', color='#D55E00', ms=2.5, label='1 km blocks')]
    ax.legend(handles=hs, loc='upper left', fontsize=FS, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.2)
    panel_label(ax, 'i', x=-0.32)

    log = historical.finish(fig, [], FS_LABEL)
    save(fig, 'FigS1')
    for a, b in sorted(set(log)):
        print('   text:', repr(a), '->', repr(b))


if __name__ == '__main__':
    main()
