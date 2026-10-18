"""Write the de-identified repository data (data/derived) from the project data tree.

Run by the authors with the project layout (PV2MVEG_LAYOUT=project, PV2MVEG_PROJECT_ROOT=<project>):

    python tools/export_derived.py [--dest data/derived]

What is removed or changed (see data/README.md):
  * unit coordinates and grid positions: x, y, lon, lat (and the full Landsat/grid unit tables that hold them);
  * construction and detection years of units: T, T_hand, hand_year, det_year (array age, the acquisition year
    minus the construction year, is kept: every result is stratified by it); the timing statistic of the text keeps
    only the detected-minus-construction-year offset of each original-cohort unit;
  * station labels: station_label is removed;
  * identifiers are replaced by order-preserving pseudonyms (dense ranks over all exported tables, one map per
    kind): cell and ctrl (150 m grid cells; ranks over every cell that appears), station (hand-drawn station
    polygons, positive values; non-positive values are kept) and blk (2 km control blocks). Order is preserved so
    that every sort, grouping and resampling in the figure and number code is unchanged; the maps are not kept;
  * columns that the figure and number code does not read are dropped from the large unit tables of the
    supplementary diagnostics (a1_*_det_thr50_cells.csv keeps cell and blk).

CSV files are processed as text: kept fields are copied verbatim, so every value parses exactly as in the source.
The manifest (data/derived/MANIFEST.json) lists each file with its source, SHA-256 of source and copy, rows and the
removed and pseudonymised columns.
"""
import argparse
import csv
import hashlib
import io
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pv2mveg.config import P, TAGS, REPO  # noqa: E402
from pv2mveg import data as D  # noqa: E402

DROP_UNITS = ['T', 'x', 'y', 'hand_year', 'station_label', 'det_year']    # cross_scene cells tables
DROP_BUFFER = ['x', 'y', 'hand_year', 'station_label', 'det_year']        # buffer_primary cells tables
PURIFY_FILTERS = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_rows(path):
    with open(path, newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def write_rows(path, header, rows):
    """CSV with LF line ends on every platform (the checksums of MANIFEST.json hold after a git checkout)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator='\n')
    w.writerow(header); w.writerows(rows)
    path.write_bytes(buf.getvalue().encode('utf-8'))


def copy_text(src, dst):
    """Copy a text file with LF line ends."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(Path(src).read_bytes().replace(b'\r\n', b'\n'))


class Pseudonyms:
    """Order-preserving dense ranks (1, 2, ...) over every value seen in the first pass."""

    def __init__(self):
        self.values = {'cell': set(), 'station': set(), 'blk': set()}
        self.maps = {}

    def see(self, kind, tokens):
        for t in tokens:
            if t == '':
                continue
            v = int(float(t))
            if kind == 'station' and v <= 0:
                continue
            self.values[kind].add(v)

    def freeze(self):
        for k, vals in self.values.items():
            self.maps[k] = {v: i + 1 for i, v in enumerate(sorted(vals))}

    def token(self, kind, t):
        """Pseudonym of one CSV token, written in the token's own number format ('12' or '12.0')."""
        if t == '':
            return t
        v = int(float(t))
        assert float(t) == v, (kind, t)
        if kind == 'station' and v <= 0:
            return t
        s = str(self.maps[kind][v])
        return s + '.0' if '.' in t else s

    def array(self, kind, a):
        m = self.maps[kind]
        return np.array([m[int(v)] for v in a], dtype=np.int64)


KIND = {'cell': 'cell', 'ctrl': 'cell', 'station': 'station', 'blk': 'blk'}


def text_job(src, dst, drop=(), keep=None, ids=(), row_filter=None, translate=None):
    return dict(kind='csv', src=Path(src), dst=Path(dst), drop=list(drop), keep=keep, ids=list(ids), row_filter=row_filter, translate=translate or {})


def plan(dest):
    """All exported files: (job, ...). Paths: source in the project layout, destination in the repository layout."""
    v12, rv12 = P.v12, dest / 'v12'
    J = []
    for t in TAGS:
        for sub in ['', '1m/', '2m/']:
            J.append(text_job(P.cs12 / f'{sub}cells_{t}.csv', rv12 / f'cross_scene/{sub}cells_{t}.csv', drop=DROP_UNITS, ids=['cell', 'station']))
        J.append(text_job(P.bp12 / f'cells_{t}.csv', rv12 / f'buffer_primary/cells_{t}.csv', drop=DROP_BUFFER, ids=['cell', 'station']))
        J.append(text_job(P.sec12 / f'1p5m/chrono/a1_{t}_det_thr50_cells.csv', rv12 / f'cross_scene/secondary/1p5m/chrono/a1_{t}_det_thr50_cells.csv',
                          keep=['cell', 'blk'], ids=['cell', 'blk']))
        for n in ['boundary_paired', 'within_unit', 'subset_summary']:
            J.append(text_job(P.sec12 / f'measurements/{t}/{n}.csv', rv12 / f'cross_scene/secondary/measurements/{t}/{n}.csv'))
    for n in ['effects', 'lanes', 'optics', 'spectral', 'distance', 'origin', '1m/effects', '2m/effects']:
        J.append(text_job(P.cs12 / f'{n}.csv', rv12 / f'cross_scene/{n}.csv'))
    # 'aux' is a reserved device name on Windows: the red-edge folders are called red_edge in the repository
    J.append(text_job(P.aux12 / 'e2_indices.csv', rv12 / 'cross_scene/red_edge/e2_indices.csv'))
    for q in ['1m', '1p5m', '2m']:
        J.append(text_job(P.sec_aux12 / q / 'e2_indices.csv', rv12 / f'cross_scene/secondary/red_edge/{q}/e2_indices.csv'))
    for n in ['effects', 'variants']:
        J.append(text_job(P.bp12 / f'{n}.csv', rv12 / f'buffer_primary/{n}.csv'))
    J.append(dict(kind='copy', src=P.bp12 / 'terrain.json', dst=rv12 / 'buffer_primary/terrain.json'))
    for n in ['threshold_common', '1p5m/statistics/holdout', '1m/statistics/holdout', '2m/statistics/holdout', '1p5m/statistics/balance',
              '1m/statistics/balance', '2m/statistics/balance', '1p5m/statistics/leave_station_out',
              '1p5m/statistics/regression', '1p5m/statistics/robustness', '1p5m/statistics/mixed_bias',
              'crossperiod/06_08/summary', 'crossperiod/08_2026_08/summary', 'temporal/forward_backward']:
        J.append(text_job(P.sec12 / f'{n}.csv', rv12 / f'cross_scene/secondary/{n}.csv'))
    J.append(dict(kind='parquet', src=P.sec12 / '1p5m/statistics/placebo_pairs.parquet', dst=rv12 / 'cross_scene/secondary/1p5m/statistics/placebo_pairs.parquet',
                  ids=['cell', 'ctrl']))
    for k in ['ndvi', 'nirv', 'msavi']:
        J.append(text_job(P.es / f'es_raw_sum_{k}_det_treated/curves.csv', dest / f'ts/ew/es_raw_sum_{k}_det_treated/curves.csv'))
    J.append(text_job(P.es / 'placebo_sum_ndvi_det/space_placebo_curve_null.csv', dest / 'ts/ew/placebo_sum_ndvi_det/space_placebo_curve_null.csv'))
    for n in ['fig2_counts', 'fig2_estimators_x1', 'fig2d_band_relative']:
        J.append(text_job(P.legacy / f'{n}.csv', dest / f'legacy/{n}.csv'))
    J.append(text_job(P.legacy_v110 / 'figS27_coverage.csv', dest / 'legacy/v110/figS27_coverage.csv', keep=['region', 'cat']))
    J.append(text_job(P.frozen / 'rainfall_frozen.csv', dest / 'frozen/rainfall_frozen.csv'))
    # side outputs of the restricted map figures (Fig. 6, Fig. S8), read by the numbers
    J.append(dict(kind='copy', src=P.figdata / 'fig5_summary.json', dst=dest / 'summaries/fig5_summary.json'))
    J.append(text_job(P.figdata / 'figS27_coverage.csv', dest / 'summaries/figS27_coverage.csv', keep=['region', 'cat']))
    # historical unmixing tests (Fig. S2): the rows the figure reads, labels in English
    J.append(text_job(P.purify / 'transfer_endmembers.csv', dest / 'purify/transfer_endmembers.csv', row_filter=lambda r: r['group'] == '东西向',
                      translate={'group': D.PURIFY_LABELS}))
    J.append(text_job(P.purify / 'transfer_metrics.csv', dest / 'purify/transfer_metrics.csv',
                      row_filter=lambda r: r['orient'] == '东西向' and r['variant'] == 'fe_cal', translate={'orient': D.PURIFY_LABELS}))
    for f in sorted(P.purify.glob('metrics_main_*.csv')):
        J.append(text_job(f, dest / f'purify/{f.name}', translate={'group': D.PURIFY_LABELS},
                          row_filter=lambda r: (float(r['unit_m']) == 150 and r['stratum'] == 'orient' and r['group'] == '东西向'
                                                and ((r['ref'] == 'A' and r['ground'] == 'g0') or (r['ref'] == 'B' and r['ground'] == 'gO')))))
    J.append(dict(kind='copy', src=P.purify / 'fe_calibration.json', dst=dest / 'purify/fe_calibration.json'))
    J.append(dict(kind='units', src=v12 / 'inputs', dst=dest / 'units/units.csv'))
    return J


def run(dest):
    P.need_project('tools/export_derived.py')
    dest = Path(dest)
    jobs = plan(dest)
    ps = Pseudonyms()
    # first pass: every identifier that will be written
    for j in jobs:
        if j['kind'] == 'csv' and j['ids']:
            h, rows = read_rows(j['src'])
            for c in j['ids']:
                ps.see(KIND[c], (r[h.index(c)] for r in rows))
        elif j['kind'] == 'parquet':
            d = pd.read_parquet(j['src'])
            for c in j['ids']:
                ps.see(KIND[c], d[c].astype(str))
    U = D.units()
    ps.see('cell', U.cell.astype(str))
    ps.freeze()
    manifest = {}
    for j in jobs:
        src, dst = j['src'], j['dst']
        try:
            rec = dict(source=src.relative_to(P.project).as_posix())
        except ValueError:  # side outputs of the restricted figures, written into the output folder
            rec = dict(source='<output>/' + src.relative_to(P.out).as_posix())
        if j['kind'] == 'copy':
            copy_text(src, dst)
        elif j['kind'] == 'csv':
            h, rows = read_rows(src)
            if j['row_filter'] is not None:
                rows = [r for r in rows if j['row_filter'](dict(zip(h, r)))]
            for c, m in j['translate'].items():
                k = h.index(c)
                for r in rows:
                    r[k] = m.get(r[k], r[k])
            for c in j['ids']:
                k = h.index(c)
                for r in rows:
                    r[k] = ps.token(KIND[c], r[k])
            cols = [c for c in h if c not in j['drop']] if j['keep'] is None else [c for c in h if c in j['keep']]
            idx = [h.index(c) for c in cols]
            write_rows(dst, cols, [[r[i] for i in idx] for r in rows])
            rec.update(rows=len(rows), removed=[c for c in h if c not in cols], pseudonymised=j['ids'],
                       filtered=j['row_filter'] is not None, translated=sorted(j['translate']))
        elif j['kind'] == 'parquet':
            d = pd.read_parquet(src)
            for c in j['ids']:
                d[c] = ps.array(KIND[c], d[c].values).astype(d[c].dtype)
            dst.parent.mkdir(parents=True, exist_ok=True); d.to_parquet(dst, index=False)
            rec.update(rows=len(d), pseudonymised=j['ids'])
        elif j['kind'] == 'units':
            u = U.copy(); u['cell'] = ps.array('cell', u.cell.values)
            dst.parent.mkdir(parents=True, exist_ok=True); u.to_csv(dst, index=False, lineterminator='\n')
            rec.update(source='data/eco/releases/fixed_ew_v2_orientation/analysis/v12/inputs/{matched_parent_hand.npz, matched_ew_hy.npz, cells_expansion.parquet}',
                       rows=len(u), derived='cell, origin (0 original cohort, 1 phased), timing_offset (det_year - hand_year, original cohort)',
                       pseudonymised=['cell'])
            src = None
        if src is not None:
            rec['source_sha256'] = sha(src)
        rec['sha256'] = sha(dst)
        rec['bytes'] = dst.stat().st_size
        manifest[dst.relative_to(dest).as_posix()] = rec
    (dest / 'MANIFEST.json').write_bytes(json.dumps(dict(
        description='De-identified derived data of the PV2mVeg repository, written by tools/export_derived.py from the project data.',
        pseudonyms={k: len(v) for k, v in ps.maps.items()}, files=manifest), indent=1, ensure_ascii=False).encode('utf-8'))
    total = sum(r['bytes'] for r in manifest.values())
    print(f'EXPORTED {len(manifest)} files, {total / 1e6:.1f} MB into {dest}; pseudonyms {dict((k, len(v)) for k, v in ps.maps.items())}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dest', default=str(REPO / 'data' / 'derived'))
    a = ap.parse_args()
    run(a.dest)
