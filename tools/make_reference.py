"""Collect the reference values of the published manuscript (v2.0) into reference/.

Run once by the authors on the research repository at the release commit (read only):

    python tools/make_reference.py --paper <project_root>/paper_eco

reference/numbers.json            every number of the text (paper_eco/figdata/v120/numbers.json)
reference/tables/                 Tables 1-2 as published (paper_eco/tables/v20: Markdown and CSV)
reference/supplement_tables.md    the pipe tables of the supplement (paper_eco/supplement/v20_parts/S1-S7.md),
                                  with their titles
reference/figures.json            per figure: pixel size and SHA-256 of the RGBA pixels of the published PNG
                                  (paper_eco/figures/v20), so that a re-render can be compared pixel by pixel
reference/figdata/                plotted-value tables written by the figure scripts of the release
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
FIGURES = ['Fig1', 'Fig2', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'Fig7', 'Fig8'] + [f'FigS{i}' for i in range(1, 13)]
FIGDATA_V120 = ['fig1_counts.json', 'fig3_cohorts.csv', 'fig3_lanes.csv', 'fig3_longitudinal.csv', 'fig4_summary.json', 'fig4e_designs.csv',
                'fig5_summary.json', 'fig6a_relative_reflectance.csv', 'fig6b_red_edge.csv', 'fig6c_distance_08.csv', 'fig6d_clearance.csv',
                'fig7_clearance.csv', 'fig7_scenario_shares.csv', 'fig7_summary.csv', 'heldout_primary.csv',
                'supp_diag/figS8_common_by_clearance.csv', 'supp_diag/figS14_distance.csv', 'supp_diag/figS21_spectral.csv',
                'supp_diag/figS22_optics_boundary.csv', 'supp_diag/figS23_mixed_minus_interrow.csv']
FIGDATA_OLD = ['figS1fg_scene_ndvi_metrics.csv', 'figS1h_band_cv_rmse.csv', 'figS6c_fe_offset.csv', 'figS6_bias_rmse_150m.csv']


def pixel_sha(path):
    Image.MAX_IMAGE_PIXELS = None
    a = np.ascontiguousarray(np.asarray(Image.open(path).convert('RGBA')))
    return a.shape[1], a.shape[0], hashlib.sha256(a.tobytes()).hexdigest()


def copy_lf(src, dst):
    """Copy a text file with LF line ends (the reference is compared byte for byte on every platform)."""
    Path(dst).write_bytes(Path(src).read_bytes().replace(b'\r\n', b'\n'))


def write_lf(dst, text):
    Path(dst).write_bytes(text.encode('utf-8'))


def tables_md(text):
    """Title lines and pipe tables of a Markdown text, in order."""
    out = []
    for block in re.split(r'\n\s*\n', text):
        lines = block.strip().split('\n')
        if lines[0].startswith('**Table S') or (lines[0].startswith('|') and out):
            out.append(block.strip())
    return out


def main(paper):
    ref = REPO / 'reference'
    (ref / 'tables').mkdir(parents=True, exist_ok=True)
    (ref / 'figdata' / 'supp_diag').mkdir(parents=True, exist_ok=True)
    copy_lf(paper / 'figdata/v120/numbers.json', ref / 'numbers.json')
    for n in ['table1.md', 'table1.csv', 'table2.md', 'table2.csv']:
        copy_lf(paper / 'tables/v20' / n, ref / 'tables' / n)
    blocks = []
    for i in range(1, 8):
        blocks += tables_md((paper / f'supplement/v20_parts/S{i}.md').read_text(encoding='utf-8'))
    write_lf(ref / 'supplement_tables.md', '\n\n'.join(blocks) + '\n')
    figs = {}
    for name in FIGURES:
        p = paper / 'figures/v20' / f'{name}.png'
        w, h, s = pixel_sha(p)
        figs[name] = dict(width=w, height=h, rgba_sha256=s, png_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    write_lf(ref / 'figures.json', json.dumps(dict(
        description='Published v2.0 figures (PNG, 500 dpi): pixel size and SHA-256 of the RGBA pixel array (row-major, 8 bit).',
        figures=figs), indent=1))
    for n in FIGDATA_V120:
        copy_lf(paper / 'figdata/v120' / n, ref / 'figdata' / n)
    for n in FIGDATA_OLD:
        copy_lf(paper / 'figdata' / n, ref / 'figdata' / n)
    print('REFERENCE', len(figs), 'figures,', len(blocks), 'supplement blocks')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--paper', required=True, help='paper_eco folder of the research repository at the release commit')
    a = ap.parse_args()
    main(Path(a.paper))
