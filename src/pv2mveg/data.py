"""Readers shared by figures, numbers and the export of the repository data.

Each reader returns the same table in both layouts (config.P): in the project layout it reads the research data
tree, in the repository layout the de-identified copy written by tools/export_derived.py from the same reader.
"""
import glob
import os

import numpy as np
import pandas as pd

from pv2mveg.config import P

# labels of the historical unmixing tables (written in Chinese by the original pipeline)
PURIFY_LABELS = {'东西向': 'east-west', '南北向': 'north-south', '全部(f≥0.05)': 'all (f≥0.05)', '全部': 'all'}


def _label(d, cols):
    for c in cols:
        if c in d:
            d[c] = d[c].replace(PURIFY_LABELS)
    return d


def merged():
    """(tidx, T, M, origin) of the 7,544 v1.2 units in cell order (restricted: needs the project layout).

    tidx unit cells; T construction year; M the ten matched control cells (-1 = none); origin 0 original cohort,
    1 phased units. Frozen inputs: v12/inputs/matched_parent_hand.npz and matched_ew_hy.npz."""
    P.need_project('The matched sets and construction years of the units')
    a = np.load(P.v12 / 'inputs/matched_parent_hand.npz'); b = np.load(P.v12 / 'inputs/matched_ew_hy.npz')
    tidx = np.concatenate([a['tidx'], b['tidx']]).astype(np.int64)
    T = np.concatenate([a['T'], b['T']]).astype(int)
    M = np.concatenate([a['M'], b['M']]).astype(np.int64)
    origin = np.r_[np.zeros(len(a['tidx']), int), np.ones(len(b['tidx']), int)]
    assert len(np.unique(tidx)) == len(tidx)
    o = np.argsort(tidx, kind='stable')
    return tidx[o], T[o], M[o], origin[o]


def cells_pool(ew=False):
    """The 150 m grid table with coordinates, panel fractions, detection and land-cover screens (restricted).
    ew=True: the east-west variant with the orientation shares and exclusion flags (ts/ew/cells_pool.parquet)."""
    P.need_project('The 150 m unit table with coordinates')
    return pd.read_parquet((P.ts / 'ew' if ew else P.ts) / 'cells_pool.parquet')


def units():
    """The 7,544 east-west analysis units of v1.2 in cell order.

    cell            unit identifier (row-major index of the 150 m grid in the project data; an order-preserving
                    pseudonym in the repository copy)
    origin          0 = original cohort (4,730 units with one construction year), 1 = phased east-west units (2,814)
    timing_offset   Landsat-detected installation year minus construction year (original cohort; empty for the
                    phased units, whose pixels do not share one detected year)
    """
    if P.layout == 'project':
        a = np.load(P.v12 / 'inputs/matched_parent_hand.npz'); b = np.load(P.v12 / 'inputs/matched_ew_hy.npz')
        tidx = np.concatenate([a['tidx'], b['tidx']]).astype(np.int64)
        origin = np.r_[np.zeros(len(a['tidx']), int), np.ones(len(b['tidx']), int)]
        assert len(np.unique(tidx)) == len(tidx)
        o = np.argsort(tidx, kind='stable')
        tidx, origin = tidx[o], origin[o]
        u = pd.read_parquet(P.v12 / 'inputs/cells_expansion.parquet').set_index('cell')
        off = (u.loc[tidx, 'det_year'] - u.loc[tidx, 'hand_year']).values.astype(float)
        return pd.DataFrame(dict(cell=tidx, origin=origin, timing_offset=np.where(origin == 0, off, np.nan)))
    return pd.read_csv(P.units / 'units.csv')


def purify_tables():
    """Tables of the panel end-member and unmixing tests (Fig. S2), east-west rows used by the figure.

    Returns (transfer_endmembers, fe_calibration_10m, metrics_main, transfer_metrics); metrics_main concatenates the
    metrics_main_*.csv files in name order."""
    import json
    te = _label(pd.read_csv(os.path.join(P.purify, 'transfer_endmembers.csv')), ['group'])
    cal = json.load(open(os.path.join(P.purify, 'fe_calibration.json'), encoding='utf-8'))['10m']
    files = sorted(glob.glob(os.path.join(P.purify, 'metrics_main_*.csv')))
    mm = _label(pd.concat([pd.read_csv(f) for f in files]), ['group'])
    tm = _label(pd.read_csv(os.path.join(P.purify, 'transfer_metrics.csv')), ['orient'])
    return te, cal, mm, tm


def purify_subset(te, mm, tm):
    """The rows of the unmixing tables that Fig. S2 uses (the repository copy keeps only these)."""
    te = te[te.group == 'east-west']
    mm = mm[(mm.unit_m == 150) & (mm.stratum == 'orient') & (mm.group == 'east-west')
            & (((mm.ref == 'A') & (mm.ground == 'g0')) | ((mm.ref == 'B') & (mm.ground == 'gO')))]
    tm = tm[(tm.orient == 'east-west') & (tm.variant == 'fe_cal')]
    return te, mm, tm
