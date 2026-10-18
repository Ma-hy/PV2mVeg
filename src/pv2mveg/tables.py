"""Main Tables 1-2 and supplementary Table S1 of the manuscript, computed from the result files.

Table 1   analysis units and selected native ground by acquisition, clearance rule and age group
          (v12/cross_scene/cells_*.csv, 1m/ and 2m/ subfolders)
Table 2   imaging-scenario artefacts for arrays aged >=6 years at >=1.5 m (v12/cross_scene/optics.csv)
Table S1  units of the primary design in 2025-08 by age group and satellite of the modal native image
          (v12/buffer_primary/cells_08.csv)
held-out  held-out pre-construction check of the primary units (v12/cross_scene/secondary/1p5m/statistics/
          holdout.csv), a source of numbers.json

Outputs (P.tables): table1.csv/.md, table2.csv/.md, tableS1.csv/.md; P.figdata/heldout_primary.csv.
The Markdown files reproduce the published tables cell by cell (titles included).
"""
import numpy as np
import pandas as pd

from pv2mveg.config import P, TAGS
from pv2mveg.plotting.house import house

LAB = {'06': '2025-06', '08': '2025-08', '2026_08': '2026-08'}
UNITS = 7544


def f4(v):
    return '—' if not np.isfinite(v) else f'{v:+.4f}'.replace('-', '−')


def iv(lo, hi):
    return f'[{f4(lo)}, {f4(hi)}]'


def effects(clearance):
    p = P.cs12 / 'effects.csv' if clearance == '1p5m' else P.cs12 / clearance / 'effects.csv'
    return pd.read_csv(p, dtype={'tag': str, 'group': str})


def cells(tag, clearance):
    p = P.cs12 / f'cells_{tag}.csv' if clearance == '1p5m' else P.cs12 / clearance / f'cells_{tag}.csv'
    return pd.read_csv(p)


def md_table(caption, header, rows, widths):
    """Pipe table; the dash counts of the separator line set the relative column widths of the Word table."""
    assert all(len(r) == len(header) for r in rows) and len(widths) == len(header)
    sep = '|' + '|'.join('-' * w for w in widths) + '|'
    return '\n'.join([caption, '', '| ' + ' | '.join(header) + ' |', sep] + ['| ' + ' | '.join(r) + ' |' for r in rows]) + '\n'


def table1():
    out, rows = [], []
    for t in TAGS:
        for p, lab in [('1m', '1'), ('1p5m', '1.5'), ('2m', '2')]:
            c = cells(t, p); v = c[c.valid_cross]
            a = v.age.values
            area = v.native_footprints.sum() * 4 / 1e4
            rows.append([LAB[t], lab, f'{len(v):,}', f'{int((a <= 2).sum()):,}', f'{int(((a >= 3) & (a <= 5)).sum()):,}', f'{int((a >= 6).sum()):,}', f'{area:,.1f}'])
            out.append(dict(tag=t, clearance=p, units=len(v), g02=int((a <= 2).sum()), g35=int(((a >= 3) & (a <= 5)).sum()), g6=int((a >= 6).sum()), area_ha=area))
    head = ['Acquisition', 'Clearance (m)', 'All ages', '0–2 years', '3–5 years', '≥6 years', 'Selected area (ha)']
    cap = (f'**Table 1.** Analysis units and selected native ground of the primary design by acquisition, clearance and array age (years since construction). '
           f'Units are 150 m × 150 m cells of east–west fixed arrays ({UNITS:,} in the sampling frame); a unit enters an acquisition with at least 50 '
           'complete 2 m ground footprints and two usable matched controls; area is the sum of selected footprints.')
    return pd.DataFrame(out), house(md_table(cap, head, rows, [14, 12, 10, 10, 10, 10, 14]))


def table2():
    o = pd.read_csv(P.cs12 / 'optics.csv', dtype={'tag': str})
    o = o[(o.design == 'primary') & (o.clearance == '1p5m') & (o.age == '≥6')]
    out, rows, units = [], [], {}
    for t in TAGS:
        z = o[o.tag == t]
        for lab, key in [('retained', 'restricted_common'), ('with sharpening', 'with_sharpening_common'), ('all', 'all_common')]:
            q = z[z.sets.str.contains(key, na=False)]
            obs = q.obs.iloc[0]
            n = f'{int(q.n.min()):,}'
            units.setdefault(LAB[t], n); assert units[LAB[t]] == n
            rows.append([f'{LAB[t]} ({f4(obs)})', lab[0].upper() + lab[1:], f'{len(q)}', f4(q.art.min()), f4(q.art.max()), f'{q.art.max() / obs * 100:.0f}%',
                         f4((obs - q.art).min()), f4((obs - q.art).max())])
            out.append(dict(tag=t, set=lab, models=len(q), art_min=q.art.min(), art_max=q.art.max(), share_max=q.art.max() / obs, net_min=(obs - q.art).min(), net_max=(obs - q.art).max(), obs=obs, n=int(q.n.min())))
    u = list(units.values())
    cap = ('**Table 2.** Imaging-scenario artefacts in ΔNDVI for arrays aged ≥6 years at ≥1.5 m clearance (original-cohort units with forward optical '
           f'models: {u[0]}, {u[1]} and {u[2]} units; Section 3.6); the observed ΔNDVI of the same units is in parentheses. “Retained” excludes sharpening and '
           'kernels rejected in the held-out test; “with sharpening” adds the partly rejected sharpening kernel; “all” adds the 48 rejected kernels. '
           'Ranges are scenario envelopes, not confidence intervals; residual, observed minus simulated.')
    head = ['Acquisition (observed)', 'Scenario set', 'Scenarios', 'Artefact min.', 'Artefact max.', 'Largest share', 'Residual min.', 'Residual max.']
    return pd.DataFrame(out), house(md_table(cap, head, rows, [17, 15, 9, 10, 10, 9, 10, 10]))


def table_s1():
    """Table S1: 2025-08 units of the primary design by age group and satellite of the modal native image.

    Native scenes 1-7 of 2025-08 are BJ3A1, 8-9 BJ3A4 and 10-13 BJ3N3; 2025-06 and 2026-08 are BJ3A1 alone."""
    c08 = pd.read_csv(P.bp12 / 'cells_08.csv')
    v08 = c08[c08.valid_primary].copy()
    v08['sat'] = pd.cut(v08.modal_scene, [0, 7, 9, 13], labels=['BJ3A1', 'BJ3A4', 'BJ3N3'])
    v08['grp'] = pd.cut(v08.age, [-1, 2, 5, 99], labels=['0–2 years', '3–5 years', '≥6 years'])
    assert v08.sat.notna().all() and v08.grp.notna().all()
    ct = pd.crosstab(v08.grp, v08.sat)
    rows = [[g] + [f'{int(ct.loc[g, s]):,}' for s in ['BJ3A1', 'BJ3A4', 'BJ3N3']] + [f'{int(ct.loc[g].sum()):,}'] for g in ct.index]
    rows += [['All ages'] + [f'{int(ct[s].sum()):,}' for s in ['BJ3A1', 'BJ3A4', 'BJ3N3']] + [f'{len(v08):,}']]
    cap = ('**Table S1.** Units of the primary design in 2025-08 by age group and by the satellite of their modal native image, with the number of '
           'scenes and the imaging date. The 2025-06 and 2026-08 acquisitions were imaged by BJ3A1 alone.')
    head = ['Age group', 'BJ3A1 (7 scenes, 2025-08-17)', 'BJ3A4 (2 scenes, 2025-08-13)', 'BJ3N3 (4 scenes, 2025-08-14)', 'All']
    sep = '|' + '|'.join('-' * w for w in [14, 16, 16, 16, 8]) + '|'
    md = '\n'.join([cap, '', '| ' + ' | '.join(head) + ' |', sep] + ['| ' + ' | '.join(r) + ' |' for r in rows])
    tab = ct.copy(); tab.columns = [str(c) for c in tab.columns]; tab['All'] = tab.sum(axis=1)
    return tab.reset_index().rename(columns={'grp': 'age_group'}), house(md)


def heldout():
    """Held-out pre-construction check of the primary units.

    The Landsat NDVI difference between each unit and its usable matched controls in e = -5..-2, relative to the mean
    difference in the matching window (e <= -6), averaged over the four years."""
    h = pd.read_csv(P.sec12 / '1p5m/statistics/holdout.csv', dtype={'tag': str, 'group': str})
    h = h[(h.controls == 'primary') & (h.quantity == 'heldout_mean') & h.group.isin(['all', '≥6'])]
    rows = [dict(tag=t, group=g, mean=float(r.estimate), lo=float(r.lo), hi=float(r.hi), n=int(r.n))
            for t in TAGS for g in ['all', '≥6'] for r in h[(h.tag == t) & (h.group == g)].itertuples()]
    return pd.DataFrame(rows)


def write_all():
    P.mkdirs()
    t1, m1 = table1(); t1.to_csv(P.tables / 'table1.csv', index=False); (P.tables / 'table1.md').write_text(m1, encoding='utf-8')
    t2, m2 = table2(); t2.to_csv(P.tables / 'table2.csv', index=False); (P.tables / 'table2.md').write_text(m2, encoding='utf-8')
    s1, ms1 = table_s1(); s1.to_csv(P.tables / 'tableS1.csv', index=False); (P.tables / 'tableS1.md').write_text(ms1 + '\n', encoding='utf-8')
    h = heldout(); h.to_csv(P.figdata / 'heldout_primary.csv', index=False)
    return dict(table1=t1, table2=t2, tableS1=s1, heldout=h)
