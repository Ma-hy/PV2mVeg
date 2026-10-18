"""Fig. S2: failure test of per-scene panel end-member estimation and unmixing (Sentinel-2, east-west arrays).

Historical figure: drawn by the figure code of an earlier release on its original inputs (not recomputed for the
primary design); only the text on the figure follows the v2.0 house style (pv2mveg.plotting.historical).

(a) red, (b) NIR panel end-member per test combination (ordered by Sentinel-2 date): scene-internal fixed-effect (FE)
    estimate, FE minus the offset calibrated on the other combinations (leave-one-out), and the reference-implied
    ("oracle") end-member; all on the same reference footprint (east-west group).
(c) FE - oracle offset per combination (red, NIR) with the 7-combination mean +/- SD.
(d)-(i) Bias and RMSE of the 150 m block inter-row ground estimate against the reference for raw Sentinel-2, FE,
    FE calibrated (leave-one-out) and oracle unmixing, NDVI / NIRv / MSAVI, east-west array blocks (f >= 0.05).

Inputs (repository data): purify/transfer_endmembers.csv, fe_calibration.json, metrics_main_*.csv (raw, FE =
fe300_orient, oracle = oracle_cv; 150 m, ground g0 for reference A and gO for reference B), transfer_metrics.csv
(fe_cal). Output: FigS2.png/.pdf; plotted values figdata/figS6c_fe_offset.csv, figS6_bias_rmse_150m.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style_v04 import *
from pv2mveg.plotting import historical
from pv2mveg.data import purify_tables
from matplotlib.lines import Line2D

CASES = ['B0528', 'A06', 'B0628', 'B0721', 'A08b', 'A08a', 'B0824', 'A26']   # by Sentinel-2 date
CDATE = {'B0528': '26 May', 'A06': '15 Jun', 'B0628': '28 Jun', 'B0721': '22 Jul', 'A08b': '14 Aug',
         'A08a': '17 Aug', 'B0824': '24 Aug', 'A26': '4 Aug 2026'}
MAIN = [c for c in CASES if c != 'B0528']          # 7 main combinations (B0528 hazy: shown in a-b only)
C_M = {'raw': C_RAW, 'fe': '#882255', 'fe_cal': '#44AA99', 'oracle': '#000000'}
M_M = {'raw': 'o', 'fe': 'v', 'fe_cal': 's', 'oracle': 'D'}
L_M = {'raw': 'Raw (mixed)', 'fe': 'FE end-member', 'fe_cal': 'FE − calibrated offset (LOO)',
       'oracle': 'Reference-implied (oracle)'}


def xlab(c):
    return f'{c} ({CDATE[c]})'


def main():
    te, cal, mm, tm = purify_tables()
    te = te[te.group == 'east-west'].set_index(['case', 'band'])

    fig = new_fig(DOUBLE, 140)
    gs_top = fig.add_gridspec(1, 3, left=0.07, right=0.985, bottom=0.66, top=0.975, wspace=0.30,
                              width_ratios=[1, 1, 0.8])
    gs_bot = fig.add_gridspec(1, 6, left=0.07, right=0.985, bottom=0.07, top=0.43, wspace=0.22)

    # (a)(b) end-members
    for k, bnd in enumerate(('red', 'nir')):
        ax = fig.add_subplot(gs_top[0, k])
        x = np.arange(len(CASES))
        for m, col, dx in (('fe', 'fe', -0.18), ('fe_cal', 'fe_cal', 0.0), ('oracle', 'orc_own', 0.18)):
            y = np.array([te.loc[(c, bnd), col] for c in CASES])
            mfc = C_M[m] if m != 'oracle' else 'white'
            ax.plot(x + dx, y, ls='', marker=M_M[m], ms=3.8, mfc=mfc, mec=C_M[m], mew=0.8, label=L_M[m], zorder=3)
        for i, c in enumerate(CASES):
            ax.plot([i - 0.18, i + 0.18], [te.loc[(c, bnd), 'fe'], te.loc[(c, bnd), 'orc_own']], color='0.6', lw=0.5,
                    zorder=1)
        ax.axvspan(-0.5, 0.5, color=C_CONSTR, lw=0, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels([xlab(c) for c in CASES], fontsize=FS, rotation=40, ha='right',
                           rotation_mode='anchor')
        ax.set_xlim(-0.5, len(CASES) - 0.5)
        ax.set_ylabel(f'Panel end-member, {L_INDEX[bnd]} reflectance')
        ax.set_ylim((0.03, 0.13) if bnd == 'red' else (0.06, 0.19))
        if k == 0:
            ax.legend(loc='upper right', fontsize=FS, handletextpad=0.2)
            ax.text(0, 0.03, 'haze', transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                    fontsize=FS, color='0.35')
        panel_label(ax, 'ab'[k], x=-0.13)

    # (c) FE - oracle offset
    ax = fig.add_subplot(gs_top[0, 2])
    zero_line(ax)
    rows = []
    for j, (bnd, key) in enumerate((('red', 'EW_red'), ('nir', 'EW_nir'))):
        pc = cal[key]['per_case']
        v = np.array([pc[c] for c in MAIN])
        xx = j + np.linspace(-0.22, 0.22, len(MAIN))
        ax.scatter(xx, v, s=10, color=C_BAND[bnd], zorder=3, lw=0)
        ax.errorbar(j + 0.36, cal[key]['offset'], yerr=cal[key]['sd_between_cases'], fmt='D', ms=4, color='k',
                    mfc=C_BAND[bnd], mec='k', mew=0.5, elinewidth=0.8, capsize=1.5, zorder=4)
        ax.text(j + 0.36, cal[key]['offset'] + cal[key]['sd_between_cases'] + 0.003,
                f"{fmt(cal[key]['offset'], 4)}\n± {cal[key]['sd_between_cases']:.4f}", ha='center', va='bottom',
                fontsize=FS)
        for c in MAIN:
            rows.append(dict(band=bnd, case=c, fe_minus_oracle=pc[c],
                             check_transfer=te.loc[(c, bnd), 'fe'] - te.loc[(c, bnd), 'orc_own']))
    pd.DataFrame(rows).to_csv(os.path.join(FIGDATA, 'figS6c_fe_offset.csv'), index=False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Red', 'NIR'])
    ax.set_xlim(-0.5, 1.7)
    ax.set_ylim(-0.02, 0.075)
    ax.set_ylabel('FE − oracle end-member')
    panel_label(ax, 'c', x=-0.2)

    # (d)-(i) bias / RMSE at 150 m
    mm = mm[(mm.unit_m == 150) & (mm.stratum == 'orient') & (mm.group == 'east-west')
            & (((mm.ref == 'A') & (mm.ground == 'g0')) | ((mm.ref == 'B') & (mm.ground == 'gO')))]
    tm = tm[(tm.orient == 'east-west') & (tm.variant == 'fe_cal')]
    tab = []
    for c in MAIN:
        for idx in ('ndvi', 'nirv', 'msavi'):
            fe = mm[(mm.case == c) & (mm['index'] == idx) & (mm.config == 'fe300_orient')].iloc[0]
            orc = mm[(mm.case == c) & (mm['index'] == idx) & (mm.config == 'oracle_cv')].iloc[0]
            fc = tm[(tm.case == c) & (tm['index'] == idx)].iloc[0]
            for m, b, r in (('raw', fe.bias_raw, fe.rmse_raw), ('fe', fe.bias_pur, fe.rmse_pur),
                            ('fe_cal', fc.bias_pur, fc.rmse_pur), ('oracle', orc.bias_pur, orc.rmse_pur)):
                tab.append(dict(case=c, index=idx, method=m, bias=b, rmse=r, n_blocks=int(fe.n)))
    tab = pd.DataFrame(tab)
    tab.to_csv(os.path.join(FIGDATA, 'figS6_bias_rmse_150m.csv'), index=False)
    y = -np.arange(len(MAIN))
    dy = {'raw': 0.27, 'fe': 0.09, 'fe_cal': -0.09, 'oracle': -0.27}
    labs = iter('defghi')
    TICKS = {('ndvi', 'bias'): [-0.08, -0.04, 0, 0.04], ('ndvi', 'rmse'): [0, 0.05, 0.1],
             ('nirv', 'bias'): [-0.03, 0, 0.03], ('nirv', 'rmse'): [0, 0.02, 0.04],
             ('msavi', 'bias'): [-0.06, -0.03, 0, 0.03], ('msavi', 'rmse'): [0, 0.04, 0.08]}
    lims = {('ndvi', 'bias'): (-0.105, 0.05), ('ndvi', 'rmse'): (0, 0.125), ('nirv', 'bias'): (-0.04, 0.024),
            ('nirv', 'rmse'): (0, 0.042), ('msavi', 'bias'): (-0.078, 0.035), ('msavi', 'rmse'): (0, 0.085)}
    for k, idx in enumerate(('ndvi', 'nirv', 'msavi')):
        for s, stat in enumerate(('bias', 'rmse')):
            ax = fig.add_subplot(gs_bot[0, 2 * k + s])
            if stat == 'bias':
                ax.axvline(0, color='k', lw=0.5)
            for i in range(len(MAIN)):
                if i % 2 == 0:
                    ax.axhspan(y[i] - 0.5, y[i] + 0.5, color='#F5F5F5', lw=0, zorder=0)
            lo, hi = lims[(idx, stat)]
            for m in ('raw', 'fe', 'fe_cal', 'oracle'):
                d = tab[(tab['index'] == idx) & (tab.method == m)].set_index('case').reindex(MAIN)
                v = d[stat].values
                clip = v > hi
                vv = np.where(clip, hi * 0.985, v)
                mfc = C_M[m] if m != 'oracle' else 'white'
                ax.plot(vv, y + dy[m], ls='', marker=M_M[m], ms=3.3, mfc=mfc, mec=C_M[m], mew=0.7, zorder=3,
                        label=L_M[m])
                for i in np.nonzero(clip)[0]:
                    ax.annotate(f'{v[i]:.3f}', xy=(hi, y[i] + dy[m]), xytext=(-2, 0), textcoords='offset points',
                                ha='right', va='center', fontsize=FS, color=C_M[m])
            ax.set_xlim(lo, hi)
            ax.set_ylim(y[-1] - 0.5, 0.5)
            ax.set_yticks(y)
            ax.set_yticklabels(MAIN if k == 0 and s == 0 else [], fontsize=FS)
            ax.tick_params(axis='y', length=0)
            ax.set_xlabel(f'{L_INDEX[idx]} {"bias" if stat == "bias" else "RMSE"}', fontsize=FS)
            ax.set_xticks(TICKS[(idx, stat)])
            # on the two bias axes with four ticks, label every other tick so they do not touch
            skip = {('ndvi', 'bias'): -0.04, ('msavi', 'bias'): -0.03}.get((idx, stat))
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
                lambda v, _, sk=skip: '' if sk is not None and abs(v - sk) < 1e-9 else f'{v:g}'.replace('-', '−')))
            ax.tick_params(axis='x', labelsize=FS)
            panel_label(ax, next(labs), x=-0.03 if (k or s) else -0.30)
    fig.axes[-6].legend(loc='upper center', bbox_to_anchor=(3.5, 1.17), ncol=4, fontsize=FS,
                        handletextpad=0.2, columnspacing=1.2)

    log = historical.finish(fig, [], FS_LABEL, rotate_first_two=True)
    save(fig, 'FigS2')
    for a, b in sorted(set(log)):
        print('   text:', repr(a), '->', repr(b))


if __name__ == '__main__':
    main()
