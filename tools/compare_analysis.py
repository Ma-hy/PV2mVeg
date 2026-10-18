"""Compare a re-run of the analysis with the archived results of the release, file by file.

usage: python tools/compare_analysis.py [--new DIR] [--ref DIR] [subfolder ...]

Defaults: --new <output>/analysis/v12 (written by scripts/10-13), --ref the archived v1.2 analysis folder of the
project data (<project_root>/data/eco/releases/fixed_ew_v2_orientation/analysis/v12; restricted). Every file of
the reference must exist in the re-run and be byte-identical; tables (CSV, Parquet) that differ in bytes are
compared value by value, so a report says which columns differ and by how much. The stage records (complete.json)
hold run times and dates; they are compared without those fields. Files of the re-run that the reference does not
hold are listed.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pv2mveg.config import P  # noqa: E402

RECORDS = {'complete.json'}          # stage records: counts and checks, with run times and dates
RUN_TIME = {'seconds', 'time'}       # fields of the stage records that depend on the run


def without_run_time(x):
    if isinstance(x, dict):
        return {k: without_run_time(v) for k, v in x.items() if k not in RUN_TIME}
    if isinstance(x, list):
        return [without_run_time(v) for v in x]
    return x


def table(p):
    return pd.read_parquet(p) if p.suffix == '.parquet' else pd.read_csv(p, low_memory=False)


def compare_tables(a, b):
    if list(a.columns) != list(b.columns):
        return f'columns differ: {sorted(set(a.columns) ^ set(b.columns))[:8]}'
    if len(a) != len(b):
        return f'rows {len(a)} vs {len(b)}'
    worst, cols = 0.0, []
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            xv, yv = x.to_numpy(float), y.to_numpy(float)
            d = np.where(np.isnan(xv) & np.isnan(yv), 0, np.abs(xv - yv))
            if np.isnan(d).any():
                cols.append(c + ' (missing values differ)')
                continue
            m = float(d.max()) if len(d) else 0.0
            if m > 0:
                cols.append(c)
                worst = max(worst, m)
        elif not (x.astype(str).values == y.astype(str).values).all():
            cols.append(c + ' (text)')
    return 'equal values, bytes differ' if not cols else f'max abs diff {worst:.3g} in {cols[:6]}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', type=Path, default=P.out / 'analysis' / 'v12')
    ap.add_argument('--ref', type=Path, default=None)
    ap.add_argument('subfolders', nargs='*')
    a = ap.parse_args()
    if a.ref is None:
        P.need_project('The comparison with the archived analysis')
        a.ref = P.v12
    groups, records = {}, 0
    for sub in a.subfolders or ['']:
        for p in sorted((a.ref / sub).rglob('*')):
            if not p.is_file():
                continue
            rel = p.relative_to(a.ref)
            q = a.new / rel
            g = groups.setdefault(str(rel.parent).replace('\\', '/'), dict(identical=0, missing=[], differing=[]))
            if not q.exists():
                g['missing'].append(p.name)
            elif p.read_bytes() == q.read_bytes():
                g['identical'] += 1
            elif p.name in RECORDS:
                same = without_run_time(json.loads(p.read_text(encoding='utf-8'))) ==                     without_run_time(json.loads(q.read_text(encoding='utf-8')))
                if same:
                    g['identical'] += 1
                    records += 1
                else:
                    g['differing'].append(f'{p.name}: differs beyond run times')
            elif p.suffix in ('.csv', '.parquet'):
                try:
                    g['differing'].append(f'{p.name}: {compare_tables(table(q), table(p))}')
                except Exception as ex:          # unreadable table: report, do not stop
                    g['differing'].append(f'{p.name}: {type(ex).__name__} {ex}')
            else:
                g['differing'].append(f'{p.name}: bytes differ')
    extra = sorted(str(q.relative_to(a.new)) for sub in (a.subfolders or ['']) if (a.new / sub).exists()
                   for q in (a.new / sub).rglob('*') if q.is_file() and not (a.ref / q.relative_to(a.new)).exists())
    tot = dict(identical=0, missing=0, differing=0)
    for k, g in sorted(groups.items()):
        tot['identical'] += g['identical']
        tot['missing'] += len(g['missing'])
        tot['differing'] += len(g['differing'])
        if g['missing'] or g['differing']:
            print(f'{k}: identical {g["identical"]}, missing {len(g["missing"])}, differing {len(g["differing"])}')
            for o in (g['missing'][:5] + g['differing'][:12]):
                print('    ', o)
    print(f'stage records equal apart from run times and dates: {records} (counted as identical)')
    print(f'files of the re-run that the reference does not hold: {len(extra)}' + ''.join(f'\n     {e}' for e in extra[:8]))
    print('TOTAL', tot, '->', 'IDENTICAL' if not tot['missing'] and not tot['differing'] else 'DIFFERENCES')
    return 0 if not tot['missing'] and not tot['differing'] else 1


if __name__ == '__main__':
    sys.exit(main())
