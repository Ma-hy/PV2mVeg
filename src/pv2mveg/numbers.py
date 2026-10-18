"""Every number quoted in the text, read from the result files (numbers/numbers.json).

Keys (one block per acquisition '06', '08', '2026_08', then study-wide blocks):
  <tag>.units, parents, units_original, units_phased, km2, footprints_m, stations     sample (Table 1, Section 3.2)
  <tag>.<index>_<group>                 ΔNDVI/ΔNIRv/ΔMSAVI/Δred/ΔNIR: estimate, joint + calibration 95% interval,
                                        units, station clusters, control level, mixed-pixel value, positive share
  <tag>.ndvi_≥6_1m / _2m                other clearance rules
  <tag>.lane_<group>                    service lanes
  <tag>.spec_<group>, dist_<group>      relative band contrasts, centre-distance bands
  <tag>.heldout                         held-out pre-construction check
  <tag>.optics_<set>                    imaging-scenario envelopes (Table 2)
  <tag>.buf_*                           unmatched buffer baseline, its variants and its agreement with the primary design
  aux                                   red-edge indices (1.2 m imagery)
  fig5                                  station-level summary of Fig. 6 (side output of the Fig. 6 script)
  fixed                                 mixed-pixel record of the units installed by 2018
  coverage, coverage_counts             coverage decomposition (side output of the Fig. S8 script)
  supp                                  supplementary diagnostics of the primary design
  units_v12, timing                     analysis units and the detected-minus-construction-year offsets
  origin, terrain                       contrasts by unit origin, terrain and greening of the comparison grounds
  optics_retained_max_any_clearance, coverage_v11
"""
import json
import numpy as np
import pandas as pd

from pv2mveg.config import P, TAGS
from pv2mveg import data as D


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float), encoding='utf-8')


def eff(clearance='1p5m'):
    p = P.cs12 / 'effects.csv' if clearance == '1p5m' else P.cs12 / clearance / 'effects.csv'
    return pd.read_csv(p, dtype={'tag': str, 'group': str})


def pick(e, tag, k, g):
    r = e[(e.tag == tag) & (e['index'] == k) & (e.group == g)].iloc[0]
    return r


def compute():
    """numbers.json as a dictionary. Needs the outputs of the tables step (P.tables/table2.csv, P.figdata/
    heldout_primary.csv) and the side outputs of the Fig. 6 and Fig. S8 scripts (P.summaries)."""
    CS, BP = P.cs12, P.bp12
    N = {}
    e = eff(); e1 = eff('1m'); e2 = eff('2m')
    for tag in TAGS:
        o = {}
        c = pd.read_csv(CS / f'cells_{tag}.csv'); v = c[c.valid_cross]
        o['units'] = int(len(v)); o['parents'] = int(len(c))
        o['units_original'] = int((v.origin == 'original').sum()); o['units_phased'] = int((v.origin == 'phased').sum())
        o['km2'] = float(v.native_footprints.sum() * 4 / 1e6); o['footprints_m'] = float(v.native_footprints.sum() / 1e6)
        o['stations'] = int(v[v.station > 0].station.nunique())
        for k in ['ndvi', 'nirv', 'msavi', 'red', 'nir']:
            for g in ['0–2', '3–5', '≥6', 'all'] + [str(a) for a in range(15)]:
                q = e[(e.tag == tag) & (e['index'] == k) & (e.group == g)]
                if len(q) == 0 or q.n_cells.iloc[0] == 0:
                    continue
                r = q.iloc[0]
                o[f'{k}_{g}'] = dict(est=r.estimate, lo=r.joint_plus_calibration_lo, hi=r.joint_plus_calibration_hi, n=int(r.n_cells), st=int(r.n_stations),
                                     ctrl=r.control_level, mixed=r.mixed, pos=r.positive_share)
        for lab, ee in [('1m', e1), ('2m', e2)]:
            r = pick(ee, tag, 'ndvi', '≥6'); o[f'ndvi_≥6_{lab}'] = dict(est=r.estimate, lo=r.joint_plus_calibration_lo, hi=r.joint_plus_calibration_hi, n=int(r.n_cells))
        ln = pd.read_csv(CS / 'lanes.csv', dtype={'tag': str, 'group': str})
        for g in ['0–2', '3–5', '≥6']:
            r = ln[(ln.tag == tag) & (ln['index'] == 'ndvi') & (ln.group == g)].iloc[0]
            o[f'lane_{g}'] = dict(est=r.lane, lo=r.lane_lo, hi=r.lane_hi, n=int(r.lane_n))
        sp = pd.read_csv(CS / 'spectral.csv', dtype={'tag': str})
        for g in ['≥6', '0–2']:
            z = sp[(sp.tag == tag) & (sp.group == g)].set_index('band')
            o[f'spec_{g}'] = {b: float(z.loc[b, 'relative'] * 100) for b in ['blue', 'green', 'red', 'nir']}
        dd = pd.read_csv(CS / 'distance.csv', dtype={'tag': str})
        for g in ['≥6', '0–2']:
            z = dd[(dd.tag == tag) & (dd.group == g) & dd.band.isin(['d1_1.5', 'd1.5_2', 'd2_3'])].set_index('band')
            o[f'dist_{g}'] = {b: float(z.loc[b, 'estimate']) for b in ['d1_1.5', 'd1.5_2', 'd2_3']}
        h = pd.read_csv(P.figdata / 'heldout_primary.csv', dtype={'tag': str})
        r = h[(h.tag == tag) & (h.group == 'all')].iloc[0]; o['heldout'] = dict(est=r['mean'], lo=r.lo, hi=r.hi)
        t4 = pd.read_csv(P.tables / 'table2.csv', dtype={'tag': str})
        for s_ in ['retained', 'with sharpening', 'all']:
            r = t4[(t4.tag == tag) & (t4.set == s_)].iloc[0]
            o[f'optics_{s_}'] = dict(art_min=r.art_min, art_max=r.art_max, share=r.share_max, net_min=r.net_min, net_max=r.net_max, obs=r.obs, n=int(r.n))
        if (BP / 'effects.csv').exists():
            bp = pd.read_csv(BP / 'effects.csv', dtype={'tag': str, 'group': str})
            for k in ['ndvi', 'nirv', 'msavi']:
                for g in ['0–2', '3–5', '≥6', 'all']:
                    r = bp[(bp.tag == tag) & (bp['index'] == k) & (bp.group == g)].iloc[0]
                    o[f'buf_{k}_{g}'] = dict(est=r.buffer_estimate, lo=r.buffer_joint_plus_calibration_lo, hi=r.buffer_joint_plus_calibration_hi, n=int(r.n_cells),
                                             d=r.diff_estimate, dlo=r.diff_joint_lo, dhi=r.diff_joint_hi, prim=r.primary_same_units)
            vv = pd.read_csv(BP / 'variants.csv', dtype={'tag': str, 'group': str})
            z = vv[(vv.tag == tag) & (vv['index'] == 'ndvi') & (vv.group == '≥6')]
            o['buf_variants_ndvi_≥6'] = {r.variant: dict(est=r.estimate, lo=r.station_lo, hi=r.station_hi, n=int(r.n_cells), prim=r.primary_same_units) for r in z.itertuples()}
            u = pd.read_csv(BP / f'cells_{tag}.csv'); u = u[u.valid_buffer]
            o['buf_r'] = float(np.corrcoef(u.primary_ndvi_dadj, u.buffer_ndvi_dadj)[0, 1]); o['buf_meandiff'] = float((u.buffer_ndvi_dadj - u.primary_ndvi_dadj).mean())
        N[tag] = o
    aux = pd.read_csv(P.aux12 / 'e2_indices.csv', dtype={'date': str}); q = aux[(aux.age == '≥6') & (aux.ctrl == 'matched')]
    N['aux'] = {k: dict(min=float(q[q['index'] == k].est.min()), max=float(q[q['index'] == k].est.max()), lo_min=float(q[q['index'] == k].lo.min()),
                        n_min=int(q[q['index'] == k].n_cells.min()), n_max=int(q[q['index'] == k].n_cells.max())) for k in ['ndre', 'cire', 'ndvi']}
    f5 = json.loads((P.summaries / 'fig5_summary.json').read_text(encoding='utf-8')); N['fig5'] = f5
    lc = pd.read_csv(P.es / 'es_raw_sum_ndvi_det_treated' / 'curves.csv')
    lc = lc[(lc.estimator == 'diff') & (lc.weight == 'cell') & (lc.subset == 'balanced T≤2018')].set_index('e')
    N['fixed'] = dict(e0=lc.loc[0, 'est'], e1=lc.loc[1, 'est'], min18=lc.loc[1:8, 'est'].min(), max18=lc.loc[1:8, 'est'].max(),
                      min1012=lc.loc[10:12, 'est'].min(), max1012=lc.loc[10:12, 'est'].max(), n10=int(lc.loc[10, 'n_cells']), n12=int(lc.loc[12, 'n_cells']))
    cov = pd.read_csv(P.summaries / 'figS27_coverage.csv')
    N['coverage'] = {rg: (z.cat.value_counts() / len(z)).to_dict() | {'n': int(len(z))} for rg, z in cov.groupby('region')}
    N['coverage_counts'] = cov.cat.value_counts().to_dict()
    # ---- supplementary diagnostics of the primary design
    SECD = P.sec12; sup = {}
    rd = lambda q: pd.read_csv(q, dtype={'tag': str, 'group': str})
    for a, b in [('06', '08'), ('08', '2026_08')]:
        d = rd(SECD / f'crossperiod/{a}_{b}/summary.csv')
        for g in ['0–2', '3–5', '≥6']:
            for qq in ['common_change', 'own_change_same_units', 'support_component']:
                r = d[(d.clearance == '1p5m') & (d['index'] == 'ndvi') & (d.quantity == qq) & (d.group == g)].iloc[0]
                sup[f'{qq}_{a}_{b}_{g}'] = dict(est=r.estimate, lo=r.lo, hi=r.hi, n=int(r.n))
    com = rd(SECD / 'threshold_common.csv')
    for t in TAGS:
        for q in ['1m', '1p5m', '2m']:
            r = com[(com.tag == t) & (com.clearance == q) & (com.group == '≥6')].iloc[0]
            sup[f'thr_common_{t}_{q}'] = dict(est=r.estimate, lo=r.joint_plus_calibration_lo, hi=r.joint_plus_calibration_hi, n=int(r.n_cells))
    rob = rd(SECD / '1p5m/statistics/robustness.csv'); rob = rob[(rob.group == '≥6') & (rob['index'] == 'ndvi')]
    prim = {t: N[t]['ndvi_≥6']['est'] for t in TAGS}
    sup['minpix_maxdiff'] = float(max(abs(r.estimate - prim[r.tag]) for r in rob[rob.scenario.str.startswith('minimum_')].itertuples()))
    for sc in ['screen_all', 'screen_targeted', 'screen_dark']:
        z = rob[rob.scenario == sc + '_paired_change']
        sup[sc] = dict(min=float(z.estimate.min()), max=float(z.estimate.max()))
    lo = rd(SECD / '1p5m/statistics/leave_station_out.csv'); lo = lo[lo['index'] == 'ndvi']
    sup['loso_maxdiff'] = float(max(abs(lo[lo.tag == t].estimate - prim[t]).max() for t in TAGS))
    sup['loso'] = {t: dict(min=float(lo[lo.tag == t].estimate.min()), max=float(lo[lo.tag == t].estimate.max())) for t in TAGS}
    bd = pd.concat([rd(SECD / f'measurements/{t}/boundary_paired.csv') for t in TAGS]); bd = bd[(bd['index'] == 'ndvi') & (bd.group == '≥6')]
    sup['boundary'] = dict(min=float(bd.estimate.min()), max=float(bd.estimate.max()))
    w = pd.concat([rd(SECD / f'measurements/{t}/within_unit.csv') for t in TAGS]); w = w[(w.group == '≥6') & (w.a == 'ir_1p5m') & (w.b == 'lane_1p5m')]
    sup['lane_minus_ir'] = {r.tag: float(r.estimate) for r in w.itertuples()}
    pl = pd.read_parquet(SECD / '1p5m/statistics/placebo_pairs.parquet'); pl['tag'] = pl.tag.astype(str)
    sup['loco_sd'] = {t: float(pl[pl.tag == t].groupby('ctrl').value.mean().std()) for t in TAGS}
    reg = rd(SECD / '1p5m/statistics/regression.csv'); reg = reg[(reg['index'] == 'ndvi') & (reg.term == 'age')]
    sup['age_slope'] = {f'{r.tag}_{r.model}': dict(est=r.estimate, lo=r.lo, hi=r.hi) for r in reg.itertuples()}
    fb = rd(SECD / 'temporal/forward_backward.csv'); fb = fb[(fb.clearance == '1p5m') & (fb.group == 'all') & (fb['index'] == 'ndvi')]
    sup['forward'] = {f'{r.sensor}_{r.direction}': dict(bias=r.bias, rmse=r.rmse, null=r.null_rmse, skill=r.skill, n=int(r.n)) for r in fb.itertuples()}
    ax = []
    for q in ['1m', '1p5m', '2m']:
        a = pd.read_csv(P.sec_aux12 / q / 'e2_indices.csv', dtype={'date': str}); a = a[(a.age == '≥6') & (a.ctrl == 'matched') & a['index'].isin(['ndre', 'cire'])]
        ax.append(a.assign(clearance=q))
    ax = pd.concat(ax)
    sup['aux_all_positive'] = bool((ax.lo > 0).all()); sup['aux_lo_min'] = float(ax.lo.min())
    h = rd(SECD / '1p5m/statistics/holdout.csv'); h = h[(h.controls == 'primary') & (h.group == 'all')]
    sup['heldout_slope'] = {t: float(h[(h.tag == t) & (h.quantity == 'heldout_slope')].estimate.iloc[0]) for t in TAGS}
    N['supp'] = sup
    # ---- v1.2: units, timing, origin, terrain
    U = D.units()
    dd = U[U.origin == 0].timing_offset.values
    N['units_v12'] = dict(total=int(len(U)), original=int((U.origin == 0).sum()), phased=int((U.origin == 1).sum()), grid_km2=float(len(U) * .0225),
                          original_cohort=4759)
    N['timing'] = dict(same=float((dd == 0).mean()), later1=float((dd == 1).mean()), earlier1=float((dd == -1).mean()), far=float((abs(dd) >= 2).mean()), n=int(len(dd)))
    og = pd.read_csv(CS / 'origin.csv', dtype={'tag': str, 'group': str})
    N['origin'] = {f'{r.tag}_{r.origin}_{r["index"]}_{r.group}': dict(est=r.estimate, lo=r.station_lo, hi=r.station_hi, n=int(r.n_cells), clusters=int(r.clusters),
                   v11=(float(r.v11_detected_year) if 'v11_detected_year' in og and pd.notna(r.v11_detected_year) else None)) for _, r in og.iterrows()}
    N['terrain'] = json.loads((BP / 'terrain.json').read_text(encoding='utf-8'))
    oo = pd.read_csv(CS / 'optics.csv', dtype={'tag': str}); oo = oo[(oo.design == 'primary') & (oo.age == '≥6') & oo.sets.str.contains('restricted_common', na=False)]
    N['optics_retained_max_any_clearance'] = float(max(q.art.max() / q.obs.iloc[0] for _, q in oo.groupby(['clearance', 'tag'])))
    c11 = pd.read_csv(P.legacy_v110 / 'figS27_coverage.csv')
    N['coverage_v11'] = {rg.split(' ')[0]: float((z.cat == 'primary sample, August 2025').mean()) for rg, z in c11.groupby('region')}
    return N


def write():
    P.mkdirs()
    N = compute()
    dump(P.numbers / 'numbers.json', N)
    return N
