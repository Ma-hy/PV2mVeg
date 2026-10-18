"""Supplementary tables recomputed from the result files and compared cell by cell with the published supplement.

Covered (every printed estimate, interval bound and unit count; tolerance 5.1e-5, i.e. the printed rounding):
  S1   units of 2025-08 by age group and satellite (pv2mveg.tables.table_s1; exact)
  S5   matching balance and held-out pre-construction checks (secondary/<clearance>/statistics/balance, holdout)
  S8   inter-row contrasts by age group (cross_scene/effects.csv)
  S9   >=6-year ΔNDVI for unit subsets and alternative unmatched baselines (cross_scene/origin.csv,
       buffer_primary/effects.csv and variants.csv)
  S10  changes on common ground footprints between acquisitions (secondary/crossperiod)
  S11  single-year cohorts (cross_scene/effects.csv)
  S12  red-edge indices at three clearances (secondary/red_edge/<clearance>/e2_indices.csv; aux/ in the project data)
  S13  clearance, footprint-minimum and screen sensitivity, within-unit contrasts and row-band perturbations
       (effects of 1, 1.5 and 2 m, secondary threshold_common, robustness, within_unit, boundary_paired)
Not recomputed here: S2-S3 (blind interpretation of earlier mask designs), S4 (rules), S6-S7 (identification tests
and robust inference of the mixed-pixel event study, frozen outputs of the event-study release; see docs/PIPELINE.md).
"""
import re

import numpy as np
import pandas as pd

from pv2mveg.config import P, TAGS

LAB = {'06': '2025-06', '08': '2025-08', '2026_08': '2026-08'}
TOL = 5.1e-5


def num(x):
    return float(x.replace('−', '-').replace(',', '').replace('%', ''))


def nums(cell):
    """All signed or unsigned decimal numbers and integers of a cell, in order (thousand separators allowed)."""
    return [num(v) for v in re.findall(r'[+−-]?\d[\d,]*\.?\d*', cell)]


def rd(p):
    return pd.read_csv(p, dtype={'tag': str, 'group': str, 'a': str, 'b': str, 'date': str, 'age_group': str})


def eff(p='1p5m'):
    return rd(P.cs12 / 'effects.csv' if p == '1p5m' else P.cs12 / p / 'effects.csv')


class Checker:
    def __init__(self):
        self.ok, self.bad = 0, []

    def close(self, where, got, want, tol=TOL):
        got, want = list(got), list(want)
        if len(got) != len(want) or not np.allclose(got, want, atol=tol, rtol=0, equal_nan=True):
            self.bad.append(f'{where}: published {got} vs recomputed {[round(float(w), 6) for w in want]}')
        else:
            self.ok += 1

    def equal(self, where, got, want):
        if got != want:
            self.bad.append(f'{where}: published {got!r} vs recomputed {want!r}')
        else:
            self.ok += 1


def _rows(blocks, name):
    return blocks[name][1][1:]       # without the header row


def s5(C, B):
    rows = {r[0]: r for r in _rows(B, 'Table S5')}
    b = rd(P.sec12 / '1p5m/statistics/balance.csv'); b = b[b.controls == 'primary']
    h = rd(P.sec12 / '1p5m/statistics/holdout.csv'); h = h[h.controls == 'primary']
    for t in TAGS:
        r = rows[LAB[t]]; z = b[b.tag == t]
        C.close(f'S5 {t} SMD', nums(r[1]) + nums(r[2]),
                [np.average(abs(z.smd_before), weights=z.n), np.average(abs(z.smd_after), weights=z.n), abs(z.smd_after).max()])
        for col, q, g in [(3, 'heldout_mean', 'all'), (4, 'heldout_slope', 'all'), (5, 'heldout_mean', '≥6')]:
            x = h[(h.tag == t) & (h.quantity == q) & (h.group == g)].iloc[0]
            C.close(f'S5 {t} {q} {g}', nums(r[col]), [x.estimate, x.lo, x.hi, x.n])
    # caption: SMD after matching and held-out means at 1 m and 2 m
    cap = B['Table S5'][0]
    smd, held = [], []
    for q in ['1m', '2m']:
        b = rd(P.sec12 / q / 'statistics/balance.csv'); b = b[b.controls == 'primary']
        h = rd(P.sec12 / q / 'statistics/holdout.csv'); h = h[(h.controls == 'primary') & (h.group == 'all') & (h.quantity == 'heldout_mean')]
        for t in TAGS:
            z = b[b.tag == t]; smd.append(np.average(abs(z.smd_after), weights=z.n)); held.append(float(h[h.tag == t].estimate.iloc[0]))
    m = re.search(r'SMD after matching was (\d+\.\d+)–(\d+\.\d+) and the held-out means ([+−]?\d+\.\d+) to ([+−]?\d+\.\d+)', cap)
    C.close('S5 caption 1 m/2 m', [num(v) for v in m.groups()], [min(smd), max(smd), min(held), max(held)])


def s8(C, B):
    e = eff()
    rows = _rows(B, 'Table S8')
    group = None
    vals = {}
    for r in rows:
        if r[0].startswith('**Arrays aged'):
            group = re.search(r'aged (.+) years', r[0]).group(1)
            continue
        vals[(group, r[0])] = r[1:]
    for g in ['0–2', '3–5', '≥6']:
        for j, t in enumerate(TAGS):
            nd = e[(e.tag == t) & (e.group == g) & (e['index'] == 'ndvi')].iloc[0]
            C.close(f'S8 {g} {t} units', nums(vals[(g, 'Units (stations)')][j]), [nd.n_cells, nd.n_stations])
            C.close(f'S8 {g} {t} positive', nums(vals[(g, 'Positive units, ΔNDVI (%)')][j]), [nd.positive_share * 100], tol=.051)
            for k, lab in [('ndvi', 'ΔNDVI'), ('nirv', 'ΔNIRv'), ('msavi', 'ΔMSAVI'), ('red', 'ΔRed reflectance'), ('nir', 'ΔNIR reflectance')]:
                if (g, lab) not in vals:
                    continue
                x = e[(e.tag == t) & (e.group == g) & (e['index'] == k)].iloc[0]
                C.close(f'S8 {g} {t} {k}', nums(vals[(g, lab)][j]), [x.estimate, x.joint_plus_calibration_lo, x.joint_plus_calibration_hi])


def s9(C, B):
    rows = {r[0]: r[1:] for r in _rows(B, 'Table S9') if not r[0].startswith('**')}
    e = eff(); og = rd(P.cs12 / 'origin.csv')
    bp = rd(P.bp12 / 'effects.csv'); vv = rd(P.bp12 / 'variants.csv')
    VAR = {'Within 3 km': 'ring_R3km', 'Within 8 km': 'ring_R8km', 'Whole cell ≥500 m from the arrays': 'ring_cell500_R5km',
           'Pixel-count weights': 'ring_R5km_pixel_weighted', 'Ground >1 km from the arrays, within 8 km': 'beyond1km_R8km',
           '>1 km, with the terrain and pre-construction support of the matched controls': 'beyond1km_R8km_support'}
    for j, t in enumerate(TAGS):
        x = e[(e.tag == t) & (e.group == '≥6') & (e['index'] == 'ndvi')].iloc[0]
        C.close(f'S9 {t} primary', nums(rows['Matched controls (primary)'][j]), [x.estimate, x.joint_plus_calibration_lo, x.joint_plus_calibration_hi, x.n_cells])
        for lab, org in [('Original cohort only', 'original'), ('Phased arrays only', 'phased')]:
            z = og[(og.tag == t) & (og.origin == org) & (og['index'] == 'ndvi') & (og.group == '≥6')].iloc[0]
            C.close(f'S9 {t} {org}', nums(rows[lab][j]), [z.estimate, z.station_lo, z.station_hi, z.n_cells])
        z = og[(og.tag == t) & (og.origin == 'original') & (og['index'] == 'ndvi') & (og.group == '≥6')].iloc[0]
        C.close(f'S9 {t} detected year', nums(rows['Original cohort, timed by the detected year'][j]), [z.v11_detected_year, z.v11_units])
        b = bp[(bp.tag == t) & (bp['index'] == 'ndvi') & (bp.group == '≥6')].iloc[0]
        C.close(f'S9 {t} buffer main', nums(rows['Main definition'][j]),
                [b.buffer_estimate, b.buffer_joint_plus_calibration_lo, b.buffer_joint_plus_calibration_hi, b.n_cells])
        z = vv[(vv.tag == t) & (vv['index'] == 'ndvi') & (vv.group == '≥6')].set_index('variant')
        for lab, var in VAR.items():
            if var in z.index:
                q = z.loc[var]
                C.close(f'S9 {t} {var}', nums(rows[lab][j]), [q.estimate, q.station_lo, q.station_hi, q.n_cells])
            else:
                C.bad.append(f'S9 {t}: variant {var} not in variants.csv ({sorted(z.index)})')


def s10(C, B):
    rows = {r[0]: r[1:] for r in _rows(B, 'Table S10')}
    QL = {'Common footprints': 'common_change', 'Own footprints, same units': 'own_change_same_units', 'Own minus common': 'support_component'}
    for j, (a, b) in enumerate([('06', '08'), ('08', '2026_08')]):
        d = rd(P.sec12 / f'crossperiod/{a}_{b}/summary.csv')
        for lab, q in QL.items():
            z = d[(d.clearance == '1p5m') & (d['index'] == 'ndvi') & (d.quantity == q) & (d.group == '≥6')].iloc[0]
            C.close(f'S10 {a}-{b} {q}', nums(rows[lab][j]), [z.estimate, z.lo, z.hi])
            if q == 'common_change':
                C.close(f'S10 {a}-{b} units', nums(rows['Units'][j]), [z.n])


def s11(C, B):
    e = eff(); e = e[e['index'] == 'ndvi']
    for r in _rows(B, 'Table S11'):
        age = r[0]
        for j, t in enumerate(TAGS):
            n_cell, v_cell = r[1 + 2 * j], r[2 + 2 * j]
            z = e[(e.tag == t) & (e.group == age)]
            if nums(n_cell) == [0]:
                C.equal(f'S11 {t} {age}', len(z) == 0 or int(z.n_cells.iloc[0]) == 0, True)
                continue
            z = z.iloc[0]
            C.close(f'S11 {t} {age} n', nums(n_cell), [z.n_cells, z.n_stations])
            want = [z.estimate] + ([z.joint_plus_calibration_lo, z.joint_plus_calibration_hi] if z.n_stations >= 6 else [])
            C.close(f'S11 {t} {age} est', nums(v_cell), want)


def s12(C, B):
    AD = {'2025-05-28': '20250528', '2025-06-28': '20250628', '2025-07-21': '20250721', '2025-08-24': '20250824'}
    for r in _rows(B, 'Table S12'):
        for j, q in enumerate(['1m', '1p5m', '2m']):
            a = rd(P.sec_aux12 / q / 'e2_indices.csv')
            z = a[(a.date == AD[r[0]]) & (a['index'] == r[1].lower()) & (a.age == '≥6') & (a.ctrl == 'matched')].iloc[0]
            C.close(f'S12 {r[0]} {r[1]} {q}', nums(r[2 + j]), [z.est, z.lo, z.hi, z.n_cells])


def s13(C, B):
    rows = _rows(B, 'Table S13')
    sect, data = None, {}
    for r in rows:
        if r[0].startswith('**'):
            sect = r[0].strip('*'); continue
        key = next(k for k in ('Clearance, each', 'Clearance, units', 'Footprint minimum', 'Leave-one-station-out', 'Within-unit contrast', 'Row band')
                   if sect.startswith(k))
        data[(key, r[0])] = r[1:]
    eff_ = {q: eff(q) for q in ['1m', '1p5m', '2m']}
    com = rd(P.sec12 / 'threshold_common.csv')
    rob = rd(P.sec12 / '1p5m/statistics/robustness.csv'); rob = rob[(rob.group == '≥6') & (rob['index'] == 'ndvi')]
    lo = rd(P.sec12 / '1p5m/statistics/leave_station_out.csv'); lo = lo[lo['index'] == 'ndvi']
    CL = {'≥1 m': '1m', '≥1.5 m (primary)': '1p5m', '≥1.5 m': '1p5m', '≥2 m': '2m'}
    SC = {'Minimum 20 footprints': 'minimum_20_pixels', 'Minimum 100 footprints': 'minimum_100_pixels', 'Minimum 200 footprints': 'minimum_200_pixels',
          'All control-pixel screens: contrast': 'screen_all_effect', 'All control-pixel screens: paired change': 'screen_all_paired_change',
          'Targeted screen: paired change': 'screen_targeted_paired_change', 'Dark-pixel screen: paired change': 'screen_dark_paired_change'}
    BL = {'Band 0.5 m narrower': 'wid-0.5', 'Band 0.5 m wider': 'wid+0.5', 'Band shifted 0.5 m north': 'shiftN0.5', 'Band shifted 0.5 m south': 'shiftS0.5',
          'North edge 0.5 m in': 'Nedge-0.5', 'North edge 0.5 m out': 'Nedge+0.5', 'South edge 0.5 m in': 'Sedge-0.5', 'South edge 0.5 m out': 'Sedge+0.5'}
    for (sec, lab), cells in data.items():
        for j, t in enumerate(TAGS):
            c = cells[j]
            if sec == 'Clearance, each':
                z = eff_[CL[lab]]; z = z[(z.tag == t) & (z['index'] == 'ndvi') & (z.group == '≥6')].iloc[0]
                C.close(f'S13 own {lab} {t}', nums(c), [z.estimate, z.joint_plus_calibration_lo, z.joint_plus_calibration_hi, z.n_cells])
            elif sec == 'Clearance, units':
                z = com[(com.tag == t) & (com.clearance == CL[lab]) & (com.group == '≥6')].iloc[0]
                C.close(f'S13 common {lab} {t}', nums(c), [z.estimate, z.joint_plus_calibration_lo, z.joint_plus_calibration_hi])
            elif sec == 'Footprint minimum':
                z = rob[(rob.tag == t) & (rob.scenario == SC[lab])].iloc[0]
                C.close(f'S13 {SC[lab]} {t}', nums(c), [z.estimate, z.lo, z.hi, z.n])
            elif sec == 'Leave-one-station-out':
                z = lo[lo.tag == t].estimate
                C.close(f'S13 leave-one-station-out {lab} {t}', nums(c), [z.min(), z.max()])
            elif sec == 'Within-unit contrast':
                w = rd(P.sec12 / f'measurements/{t}/within_unit.csv')
                a, b = ('ir_1p5m', 'lane_1p5m') if lab.startswith('Lanes') else ('d1.5_2', 'd2_3')
                z = w[(w.group == '≥6') & (((w.a == a) & (w.b == b)) | ((w.a == b) & (w.b == a)))]
                if len(z) != 1:
                    C.bad.append(f'S13 within-unit {lab} {t}: {len(z)} rows ({sorted(set(zip(w.a, w.b)))})'); continue
                z = z.iloc[0]
                C.close(f'S13 within {lab} {t}', nums(c), [z.estimate, z.lo, z.hi, z.n])
            elif sec == 'Row band':
                b = rd(P.sec12 / f'measurements/{t}/boundary_paired.csv')
                z = b[(b.subset == 'boundary_' + BL[lab]) & (b['index'] == 'ndvi') & (b.group == '≥6')].iloc[0]
                C.close(f'S13 {BL[lab]} {t}', nums(c), [z.estimate, z.lo, z.hi, z.n])
            else:
                C.bad.append(f'S13: unknown section {sec!r}')


def check(blocks):
    C = Checker()
    from pv2mveg.tables import table_s1
    from pv2mveg.verify import _rows as md_rows
    _, md = table_s1()
    C.equal('S1 cells', md_rows(md), blocks['Table S1'][1])
    for f in (s5, s8, s9, s10, s11, s12, s13):
        try:
            f(C, blocks)
        except Exception as ex:  # a missing row or file is a failed check, reported with the others
            C.bad.append(f'{f.__name__}: {type(ex).__name__}: {ex}')
    return dict(status='all cells agree' if not C.bad else 'differences', checked=C.ok, failed=C.bad[:40], n_failed=len(C.bad))
