"""E2 (spectral test), 1.2 m six bands: spectral difference of the inter-row ground (>=2.1 m from panels) relative to controls, chlorophyll / visible / soil-line indices, seasonality, two-component decomposition.

Units and controls as in b2_analysis.py, restricted to east-west fixed arrays (cells_pool.orient = 1): ts units (array age = 2025 − detection year) with >=20 inter-row pixels >=2.1 m from panels (clouds and cloud shadows removed by e1);
  local   pure-ground units in the image >=500 m from the nearest panel (>=5000 pixels of 1.2 m), >=3 within 3 km; control spectrum = equal-weight mean of the mean spectra of the control units;
  matched the ts matched controls inside the image (>=2). The Landsat baseline is not subtracted (it has red/NIR only, and b2 shows differences <=0.01 with and without it).
Indices (from mean reflectances): NDVI, NDRE=(N−RE)/(N+RE), CIre=N/RE−1, MTCI=(N−RE)/(RE−R), GRVI=(G−R)/(G+R), VARI=(G−R)/(G+R−B),
  ExG=2g−r−b (chromatic coordinates), CI=1−(R−B)/(R+B) (crust index of Karnieli 1997: lower when blue is low relative to red), PVI (perpendicular distance to the soil line, see below).
Soil line: per date, a Theil–Sen line N = a + bR through the mean red/NIR of all local control units (main); sensitivity with the 25% of control units of lowest NDVI ("bare-soil" line).
Two-component decomposition: group-mean spectral difference Δρ_b = α(V_b − S_b) + κ S_b, V = mean spectrum of the 1% highest-NDVI pure-ground pixels outside arrays in the same scene (the "green vegetation direction" of the scene),
  S = mean control spectrum of the group; κ<0 is uniform (proportional) darkening. A darkening-only model (κ only) is also fitted; residuals are compared.
Intervals: station-clustered bootstrap, 95%.
Output: data/eco/adjacency/spectral/e2_{spectrum,indices,season,decomp,soilline}.csv; figures figs/e2_spectrum.png, figs/e2_season.png
"""
# Mirrored from eco/adjacency/e2_spec_analysis.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import adjc as A  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import cc  # noqa: E402

SP = os.path.join(A.OUT, 'spectral')
DATES = ['20250528', '20250628', '20250721', '20250824']
BANDS = ['deepblue', 'blue', 'green', 'red', 'rededge', 'nir']
WL = [436, 483, 562, 654, 723, 828]
AGEG = [('0–2', 0, 2), ('3–5', 3, 5), ('≥6', 6, 99)]
IDXN = ['ndvi', 'ndre', 'cire', 'mtci', 'grvi', 'vari', 'exg', 'ci', 'pvi', 'pvi_bare']


def indices(S, soil):
    """S: dict band -> array (mean reflectance). soil: (a, b, a_bare, b_bare)."""
    B_, G, R, RE, N = S['blue'], S['green'], S['red'], S['rededge'], S['nir']
    with np.errstate(invalid='ignore', divide='ignore'):
        tot = R + G + B_
        out = dict(ndvi=(N - R) / (N + R), ndre=(N - RE) / (N + RE), cire=N / RE - 1, mtci=(N - RE) / (RE - R),
                   grvi=(G - R) / (G + R), vari=(G - R) / (G + R - B_), exg=(2 * G - R - B_) / tot, ci=1 - (R - B_) / (R + B_),
                   pvi=(N - soil[0] - soil[1] * R) / np.sqrt(1 + soil[1] ** 2),
                   pvi_bare=(N - soil[2] - soil[3] * R) / np.sqrt(1 + soil[3] ** 2))
    return out


def theil_sen(x, y):
    from scipy.stats import theilslopes
    b, a, _, _ = theilslopes(y, x)
    return a, b


def main():
    cells = pd.read_parquet(os.path.join(A.ROOT, 'data', 'eco', 'ts', 'cells_pool.parquet'), columns=['cell', 'x', 'y', 'dist_min', 'station', 'orient'])
    tidx, Td, M = cc.matched('det')
    blk_all = cc.es.block_id(cells.x.values, cells.y.values)
    veg = pd.read_csv(os.path.join(SP, 'e1_veg.csv'), dtype={'date': str})
    rows_sp, rows_ix, rows_dc, rows_sl, percell = [], [], [], [], []
    from scipy.spatial import cKDTree
    for date in DATES:
        P = pd.read_parquet(os.path.join(SP, f'e1_{date}_prof.parquet'))
        Cc = pd.read_parquet(os.path.join(SP, f'e1_{date}_ctrl.parquet'))
        Cc = Cc[Cc.n >= 5000].merge(cells, on='cell')
        for b in BANDS:
            Cc[b] = Cc[b] / Cc.n
        # soil line
        a0, b0 = theil_sen(Cc.red.values, Cc.nir.values)
        nd = (Cc.nir - Cc.red) / (Cc.nir + Cc.red)
        bare = Cc[nd <= nd.quantile(0.25)]
        a1, b1 = theil_sen(bare.red.values, bare.nir.values)
        soil = (a0, b0, a1, b1)
        rows_sl.append(dict(date=date, a_all=a0, b_all=b0, a_bare=a1, b_bare=b1, n_all=len(Cc), n_bare=len(bare)))
        Cl = Cc[Cc.dist_min >= 500].reset_index(drop=True)
        CI = indices({b: Cl[b].values for b in BANDS}, soil)
        CI_m = indices({b: Cc[b].values for b in BANDS}, soil)
        # treated units
        Q = P[P.sdlo >= globals().get('CLEARANCE_M', 2.1) - 1e-6].groupby('cell')[['n'] + BANDS].sum()
        Q = Q[Q.n >= globals().get('MIN_PIXELS', 20)]
        for b in BANDS:
            Q[b] = Q[b] / Q.n
        t = pd.DataFrame(dict(cell=tidx, T=Td, row=np.arange(len(tidx))))
        t['age'] = 2025 - t['T']
        t = t[t.age >= 0].merge(cells, on='cell').merge(Q, left_on='cell', right_index=True)
        t = t[t.orient == 1].reset_index(drop=True)      # scope: east-west fixed arrays only
        if 'ELIGIBLE_CELLS' in globals():
            t = t[t.cell.isin(ELIGIBLE_CELLS[date])].reset_index(drop=True)
        t['clu'] = np.where(t.station > 0, t.station, 10 ** 9 + blk_all[t.cell.values])
        TI = indices({b: t[b].values for b in BANDS}, soil)
        # local
        kt = cKDTree(Cl[['x', 'y']].values)
        nb = kt.query_ball_point(t[['x', 'y']].values, 3000)
        okl = np.array([len(j) >= 3 for j in nb])
        # matched
        Mm = M[t.row.values]
        pos = pd.Series(np.arange(len(Cc)), index=Cc.cell.values)
        inB = np.isin(Mm, Cc.cell.values) & (Mm >= 0)
        okm = inB.sum(1) >= 2
        D = {}
        for cv, okc in (('local', okl), ('matched', okm)):
            for b in BANDS:
                if cv == 'local':
                    cvals = np.array([Cl[b].values[j].mean() if len(j) >= 3 else np.nan for j in nb])
                else:
                    V = Cc[b].values
                    cvals = np.array([V[pos[m[i]].values].mean() if i.sum() >= 2 else np.nan for m, i in zip(Mm, inB)])
                D[(cv, 'band', b)] = t[b].values - cvals
                D[(cv, 'rel', b)] = t[b].values / cvals - 1
            for k in IDXN:
                src = CI if cv == 'local' else CI_m
                if cv == 'local':
                    cvals = np.array([np.nanmean(src[k][j]) if len(j) >= 3 else np.nan for j in nb])
                else:
                    cvals = np.array([np.nanmean(src[k][pos[m[i]].values]) if i.sum() >= 2 else np.nan for m, i in zip(Mm, inB)])
                D[(cv, 'idx', k)] = TI[k] - cvals
                D[(cv, 'ctrl', k)] = cvals
                D[(cv, 'trt', k)] = np.where(okc, TI[k], np.nan)
            D[(cv, 'ok', '')] = okc
        # summary
        for cv in ('local', 'matched'):
            okc = D[(cv, 'ok', '')]
            for an, lo, hi in AGEG:
                s = okc & (t.age >= lo).values & (t.age <= hi).values
                for b in BANDS:
                    for q in ('band', 'rel'):
                        m, l, h, n, g = A.boot_ci(D[(cv, q, b)][s], t.clu.values[s], B=1000, seed=11)
                        rows_sp.append(dict(date=date, ctrl=cv, age=an, band=b, wl=WL[BANDS.index(b)], qty=q, est=m, lo=l, hi=h, n_cells=n, n_clu=g))
                    # absolute spectra of controls and treated units
                    rows_sp.append(dict(date=date, ctrl=cv, age=an, band=b, wl=WL[BANDS.index(b)], qty='trt_mean', est=float(np.nanmean(t[b].values[s]))))
                for k in IDXN:
                    m, l, h, n, g = A.boot_ci(D[(cv, 'idx', k)][s], t.clu.values[s], B=1000, seed=12)
                    rows_ix.append(dict(date=date, ctrl=cv, age=an, index=k, est=m, lo=l, hi=h, n_cells=n, n_clu=g,
                                        ctrl_mean=float(np.nanmean(D[(cv, 'ctrl', k)][s])), trt_mean=float(np.nanmean(D[(cv, 'trt', k)][s]))))
                # two-component decomposition (group means)
                if s.sum() >= 10:
                    dr = np.array([np.nanmean(D[(cv, 'band', b)][s]) for b in BANDS])
                    Sg = np.array([np.nanmean(t[b].values[s] - D[(cv, 'band', b)][s]) for b in BANDS])
                    Vv = veg[(veg.date == date) & (veg.kind == 'veg_top1pct')][BANDS].values[0]
                    X2 = np.c_[Vv - Sg, Sg]
                    co, *_ = np.linalg.lstsq(X2, dr, rcond=None)
                    res2 = dr - X2 @ co
                    k1 = np.dot(Sg, dr) / np.dot(Sg, Sg)
                    res1 = dr - k1 * Sg
                    rows_dc.append(dict(date=date, ctrl=cv, age=an, n_cells=int(s.sum()), alpha_veg=co[0], kappa_dark=co[1],
                                        rms_2comp=float(np.sqrt(np.mean(res2 ** 2))), rms_darkonly=float(np.sqrt(np.mean(res1 ** 2))), kappa_darkonly=k1,
                                        rms_obs=float(np.sqrt(np.mean(dr ** 2))), veg_ndvi=float((Vv[5] - Vv[3]) / (Vv[5] + Vv[3])),
                                        **{f'd_{b}': dr[i] for i, b in enumerate(BANDS)}, **{f'res2_{b}': res2[i] for i, b in enumerate(BANDS)},
                                        **{f'res1_{b}': res1[i] for i, b in enumerate(BANDS)}, **{f'S_{b}': Sg[i] for i, b in enumerate(BANDS)},
                                        **{f'V_{b}': Vv[i] for i, b in enumerate(BANDS)}))
        pc = pd.DataFrame({'cell': t.cell.values, 'age': t.age.values, 'clu': t.clu.values})
        for k in IDXN + ['ndvi']:
            pc[f'd_local_{k}'] = D[('local', 'idx', k)]
            pc[f'trt_{k}'] = TI[k]
            pc[f'ctrl_local_{k}'] = D[('local', 'ctrl', k)]
        for b in BANDS:
            pc[f'rel_local_{b}'] = D[('local', 'rel', b)]
        pc['date'] = date
        percell.append(pc)
        A.log(f'{date}: 处理单元 {len(t)}（local {okl.sum()}，matched {okm.sum()}）；土壤线 N = {a0:.3f} + {b0:.3f} R（裸土 {a1:.3f} + {b1:.3f} R）', 'log_e2.txt')
    pd.DataFrame(rows_sp).to_csv(os.path.join(SP, 'e2_spectrum.csv'), index=False)
    pd.DataFrame(rows_ix).to_csv(os.path.join(SP, 'e2_indices.csv'), index=False)
    pd.DataFrame(rows_dc).to_csv(os.path.join(SP, 'e2_decomp.csv'), index=False)
    pd.DataFrame(rows_sl).to_csv(os.path.join(SP, 'e2_soilline.csv'), index=False)
    PC = pd.concat(percell)
    PC.to_parquet(os.path.join(SP, 'e2_percell.parquet'))
    # seasonality: the same units, with local controls on all four dates
    W = PC.pivot_table(index='cell', columns='date', values='d_local_ndvi')
    common = W.dropna().index
    rows_se = []
    for an, lo, hi in AGEG:
        for k in ('ndvi', 'ndre', 'cire', 'grvi', 'vari', 'pvi', 'pvi_bare', 'ci', 'mtci', 'exg'):
            for date in DATES:
                q = PC[(PC.date == date) & PC.cell.isin(common) & (PC.age >= lo) & (PC.age <= hi)]
                m, l, h, n, g = A.boot_ci(q[f'd_local_{k}'].values, q.clu.values, B=1000, seed=13)
                rows_se.append(dict(age=an, index=k, date=date, d=m, lo=l, hi=h, n_cells=n, n_clu=g,
                                    trt=float(q[f'trt_{k}'].mean()), ctrl=float(q[f'ctrl_local_{k}'].mean())))
    pd.DataFrame(rows_se).to_csv(os.path.join(SP, 'e2_season.csv'), index=False)
    figs()


def figs():
    plt = A.setup_fonts()
    S = pd.read_csv(os.path.join(SP, 'e2_spectrum.csv'), dtype={'date': str})
    DC = pd.read_csv(os.path.join(SP, 'e2_decomp.csv'), dtype={'date': str})
    fig, axs = plt.subplots(2, 4, figsize=(15, 7), sharex=True)
    for i, date in enumerate(DATES):
        for j, q in enumerate(('rel', 'band')):
            ax = axs[j, i]
            for an, c in (('0–2', 'C1'), ('3–5', 'C2'), ('≥6', 'C0')):
                g = S[(S.date == date) & (S.ctrl == 'local') & (S.age == an) & (S.qty == q)].sort_values('wl')
                ax.errorbar(g.wl, g.est, yerr=[g.est - g.lo, g.hi - g.est], marker='o', color=c, capsize=2, label=f'{an} 岁（{int(g.n_cells.iloc[0])} 单元）')
                if q == 'band' and an == '≥6':
                    d = DC[(DC.date == date) & (DC.ctrl == 'local') & (DC.age == an)].iloc[0]
                    Sg = np.array([d[f'S_{b}'] for b in BANDS]); Vv = np.array([d[f'V_{b}'] for b in BANDS])
                    ax.plot(WL, d.kappa_darkonly * Sg, 'k:', label='≥6 岁：纯同比例变暗拟合')
                    ax.plot(WL, d.alpha_veg * (Vv - Sg) + d.kappa_dark * Sg, 'k--', label='≥6 岁：植被 + 变暗拟合')
            ax.axhline(0, color='0.5', lw=0.6)
            ax.set_title(f'{date[:4]}-{date[4:6]}-{date[6:]}  ' + ('相对差 ρ行间/ρ对照 − 1' if q == 'rel' else '反射率差 ρ行间 − ρ对照'))
            if j == 1:
                ax.set_xlabel('波长 (nm)：深蓝、蓝、绿、红、红边、近红外')
    axs[0, 0].legend(fontsize=7)
    axs[1, 0].legend(fontsize=7)
    A.savefig(fig, 'e2_spectrum.png')
    SE = pd.read_csv(os.path.join(SP, 'e2_season.csv'), dtype={'date': str})
    ks = ['ndvi', 'ndre', 'cire', 'grvi', 'vari', 'pvi_bare']
    fig, axs = plt.subplots(2, len(ks), figsize=(3.2 * len(ks), 6.5))
    x = np.arange(4)
    for j, k in enumerate(ks):
        for an, c in (('0–2', 'C1'), ('≥6', 'C0')):
            g = SE[(SE['index'] == k) & (SE.age == an)].set_index('date').reindex(DATES)
            axs[0, j].plot(x, g.trt, '-o', color=c, label=f'行间 {an} 岁')
            axs[0, j].plot(x, g.ctrl, '--', color=c, alpha=0.6, label=f'对照（{an} 岁的）')
            axs[1, j].errorbar(x, g.d, yerr=[g.d - g.lo, g.hi - g.d], fmt='-o', color=c, capsize=2, label=f'{an} 岁')
        axs[1, j].axhline(0, color='0.5', lw=0.6)
        axs[0, j].set_title(k.upper()); axs[1, j].set_title(f'{k.upper()} 行间 − 对照')
        for a in axs[:, j]:
            a.set_xticks(x); a.set_xticklabels(['05-28*', '06-28', '07-21', '08-24'], fontsize=8)
    axs[0, 0].legend(fontsize=7); axs[1, 0].legend(fontsize=7)
    A.savefig(fig, 'e2_season.png')
