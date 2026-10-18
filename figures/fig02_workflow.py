"""Fig. 2: workflow of the study (data, processing, design, evidence).

Four layers, each labelled with the sections that describe it. The boxes of the data, processing and design layers
sit on one four-column grid (a box spans one or two columns), so that every connection is a straight vertical or
horizontal arrow; the contrasts reach the five evidence boxes through one distribution line (trunk, horizontal bus and
vertical drops). Solid arrows carry data and results, dashed arrows lead to the diagnostics. The two baselines are
drawn as parallel alternatives in a dashed group, the matched controls next to the event study whose matching
selected them. Small italic labels at the lower right of the boxes name the figures and tables that report each step.
The three analysed surfaces keep the colours of Figs 4 and S6 (GROUND, MIXED, LANE); everything else is monochrome.

Input: numbers/numbers.json (scripts/30_numbers.py): the unit count and the >=6-year ΔNDVI printed on the figure.
Output: Fig2.png/.pdf
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

NUMBERS = P.numbers / 'numbers.json'
if not NUMBERS.exists():
    sys.exit(f'{NUMBERS} is missing: run scripts/30_numbers.py first')
N = json.loads(NUMBERS.read_text(encoding='utf-8'))
UV = N['units_v12']
g6 = [N[t]['ndvi_≥6']['est'] for t in TAGS]


def s(v, nd=4):
    return f'{v:+.{nd}f}'.replace('-', '−')


W, H = 165.0, 142.0                          # mm
INKC = '#000000'
EDGE = '#3E4C59'
FILL = '#F4F6F8'
LAYER = '#E9EDF1'
FS, FS_HEAD, FS_TAG, FS_LAYER = 7.0, 7.6, 6.9, 7.8
DASH = (0, (3, 2))

# grid: four columns for the upper three layers, five for the evidence layer (mm)
L, R = 24.0, 163.5
G4 = 6.0
W4 = (R - L - 3 * G4) / 4
COL = [(L + i * (W4 + G4), L + i * (W4 + G4) + W4) for i in range(4)]       # A, B, C, D
G5 = 2.5
W5 = (R - L - 4 * G5) / 5
EV = [(L + i * (W5 + G5), L + i * (W5 + G5) + W5) for i in range(5)]
cx = lambda c: (c[0] + c[1]) / 2
span = lambda a, b: (COL[a][0], COL[b][1])

# layers (label, sections, band y0, band y1) and the box rows inside them
LAYERS = [('Data', 'Section 2.2', 117.5, 141.0), ('Processing', 'Sections 3.1–3.2', 88.0, 112.5),
          ('Design', 'Sections 3.3–3.5', 38.5, 83.0), ('Evidence', 'Sections 3.6, 4', 1.0, 33.5)]
DATA_Y, PROC_Y, DES_Y, EV_Y, RES_Y = (119.5, 139.0), (90.0, 110.5), (40.5, 81.0), (12.0, 31.5), (2.0, 9.0)
BUS_Y = 36.0


def box(ax, x0, x1, y0, y1, head, body, tag=None, fill=FILL, edge=EDGE, lw=.6, chips=None):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle='round,pad=0,rounding_size=1.2',
                                fc=fill, ec=edge, lw=lw, zorder=2))
    c = (x0 + x1) / 2
    top = y1 - 2.2
    if head:
        ax.text(c, top, head, ha='center', va='top', fontsize=FS_HEAD, color=INKC, zorder=4, linespacing=1.05,
                path_effects=bold(INKC, .35))
    nl = head.count('\n') + 1 if head else 0
    if body:
        ax.text(c, top - 3.3 * nl - .4, body, ha='center', va='top', fontsize=FS, color=INKC, zorder=4, linespacing=1.12)
    if chips:
        n = len(chips); cw = 4.0; x = c - (n * cw + (n - 1) * 1.4) / 2
        for col in chips:
            ax.add_patch(Rectangle((x, y0 + 6.2), cw, 2.2, fc=col, ec='none', zorder=4))
            x += cw + 1.4
    if tag:
        ax.text(x1 - 1.0, y0 + 1.0, tag, ha='right', va='bottom', fontsize=FS_TAG, color=INKC, style='italic', zorder=4)
    return dict(x0=x0, x1=x1, y0=y0, y1=y1, cx=c, cy=(y0 + y1) / 2)


def arrow(ax, p, q, dashed=False, lw=.8):
    """A straight arrow; every call below is vertical or horizontal."""
    assert abs(p[0] - q[0]) < 1e-9 or abs(p[1] - q[1]) < 1e-9, (p, q)
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=7, lw=lw, color=EDGE, zorder=3,
                                 linestyle=DASH if dashed else '-', shrinkA=0, shrinkB=0))


def line(ax, xs, ys, dashed=False, lw=.8):
    ax.plot(xs, ys, color=EDGE, lw=lw, ls=DASH if dashed else '-', zorder=3, solid_capstyle='butt')


def main():
    fig = figure(H, W)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
    for lab, sec, y0, y1 in LAYERS:
        ax.add_patch(Rectangle((0, y0), W, y1 - y0, fc=LAYER, ec='none', zorder=0))
        ax.text(1.6, (y0 + y1) / 2 + 1.6, lab, ha='left', va='center', fontsize=FS_LAYER, color=INKC, path_effects=bold(INKC, .35))
        ax.text(1.6, (y0 + y1) / 2 - 2.0, sec, ha='left', va='center', fontsize=FS_TAG, color=INKC)

    # data: Beijing-3 over the two mapping columns, the time series and the auxiliary data over the timing box
    d1 = box(ax, *span(0, 1), *DATA_Y, 'Beijing-3 imagery',
             '0.5 m panchromatic, 2 m multispectral\n2025-06, 2025-08, 2026-08 (40 scenes)\ncalibrated to Sentinel-2', 'Fig. S1; Table S1')
    d2 = box(ax, *COL[2], *DATA_Y, 'Satellite time series', 'Landsat 2000–2026,\n674 scenes;\nSentinel-2 2019–2026,\n354 scenes')
    d3 = box(ax, *COL[3], *DATA_Y, 'Auxiliary data', '1.2 m red-edge imagery\nterrain, land cover\nAgERA5 rainfall\nstation polygons')
    # processing
    p1 = box(ax, *COL[0], *PROC_Y, 'Row-band mapping', 'U-Net and geometric rules\non the panchromatic band;\neligibility layer', 'Tables S2–S4')
    p2 = box(ax, *COL[1], *PROC_Y, 'Inter-row sampling', 'complete 2 m footprints,\ncentre ≥1.5 m (≥1, ≥2 m)\nfrom the row band; lanes', 'Fig. 1; Table 1')
    p3 = box(ax, *span(2, 3), *PROC_Y, 'Installation timing and units',
             'construction year from the station polygons,\ndetected year from the Landsat step;\n'
             f'{UV["total"]:,} east–west units of 150 m × 150 m', 'Fig. S3')
    # design: the contrasts (two columns), the two baselines stacked as parallel alternatives in a dashed group, and
    # the event study, whose matching selected the matched controls, beside them
    g3 = box(ax, *span(0, 1), *DES_Y, 'Baseline-adjusted contrasts',
             'Δ = inter-row ground − baseline ground\n− pre-construction Landsat difference\n\n'
             'ΔNDVI, ΔNIRv and ΔMSAVI of inter-row\nground, mixed pixels and lanes\n\n'
             'equal weight for every 150 m unit;\n600 joint resamples + calibration error', 'Table S8')
    gx0, gx1 = COL[2][0] - 1.0, COL[2][1] + 1.0
    ax.add_patch(FancyBboxPatch((gx0, DES_Y[0]), gx1 - gx0, DES_Y[1] - DES_Y[0], boxstyle='round,pad=0,rounding_size=1.2',
                                fc='#FFFFFF', ec=EDGE, lw=.6, ls=DASH, zorder=1))
    ax.text(COL[2][0] + .3, DES_Y[1] - .6, 'Alternative baselines', ha='left', va='top', fontsize=FS_TAG, color=INKC, style='italic', zorder=4)
    gm = box(ax, *COL[2], 59.2, 77.6, 'Matched controls', 'ten per unit outside\nthe park, matched on\npre-construction history', 'Fig. 5; Table S5')
    gu = box(ax, *COL[2], 41.5, 58.0, 'Unmatched buffer', 'cells 300–1,000 m from\nthe arrays, within 5 km', 'Table S9')
    g4 = box(ax, *COL[3], *DES_Y, 'Event study on\nmixed pixels',
             '4,759 original-cohort\nunits, detected year;\nLandsat 2000–2026,\nSentinel-2 check;\nmatched and DiD\nestimators', 'Fig. 3; Tables S6–S7')
    # evidence and result
    e = [box(ax, *EV[0], *EV_Y, 'Three surfaces', 'inter-row ground\nmixed pixel, lane', 'Fig. 4', chips=[GROUND, MIXED, LANE]),
         box(ax, *EV[1], *EV_Y, 'Age and space', 'single-year cohorts\nstation means', 'Figs 4, 6'),
         box(ax, *EV[2], *EV_Y, 'Spectra and\nred edge', 'band contrasts\nNDRE, CIre', 'Fig. 7a–b'),
         box(ax, *EV[3], *EV_Y, 'Distance and\nclearance', 'distance bands,\nclearance, boundaries', 'Fig. 7c–d; Table S13'),
         box(ax, *EV[4], *EV_Y, 'Imaging\nscenarios', 'simulated artefact\nshare of the contrast', 'Fig. 8; Table 2')]
    res = box(ax, L, R, *RES_Y, '', '', fill='#FFFFFF', edge=GROUND, lw=.9)
    ax.text(res['cx'], res['cy'] + 1.6, 'Inter-row ground of arrays aged ≥6 years is greener than comparable steppe',
            ha='center', va='center', fontsize=FS_HEAD, color=INKC, zorder=4, path_effects=bold(INKC, .35))
    ax.text(res['cx'], res['cy'] - 1.7, 'ΔNDVI ' + ', '.join(s(v, 3) for v in g6[:2]) + f' and {s(g6[2], 3)} in 2025-06, 2025-08 and 2026-08',
            ha='center', va='center', fontsize=FS, color=INKC, zorder=4)

    # data -> processing (one arrow per column) and Sentinel-2 -> Beijing-3 calibration
    arrow(ax, (cx(COL[0]), d1['y0']), (cx(COL[0]), p1['y1']))
    arrow(ax, (cx(COL[1]), d1['y0']), (cx(COL[1]), p2['y1']))
    arrow(ax, (cx(COL[2]), d2['y0']), (cx(COL[2]), p3['y1']))
    arrow(ax, (cx(COL[3]), d3['y0']), (cx(COL[3]), p3['y1']))
    arrow(ax, (d2['x0'], d2['cy']), (d1['x1'], d2['cy']))
    # processing
    arrow(ax, (p1['x1'], p1['cy']), (p2['x0'], p1['cy']))
    arrow(ax, (cx(COL[1]), p2['y0']), (cx(COL[1]), g3['y1']))
    ax_in = COL[2][1] - 4.0                                  # into the group, clear of its label
    arrow(ax, (ax_in, p3['y0']), (ax_in, DES_Y[1]))
    arrow(ax, (cx(COL[3]), p3['y0']), (cx(COL[3]), g4['y1']))
    # design: event-study matching -> matched controls; each baseline -> contrasts
    arrow(ax, (g4['x0'], gm['cy']), (gm['x1'], gm['cy']))
    arrow(ax, (gm['x0'], gm['cy']), (g3['x1'], gm['cy']))
    arrow(ax, (gu['x0'], gu['cy']), (g3['x1'], gu['cy']))
    # contrasts -> evidence through one distribution line: trunk, bus, drops (dashed to the diagnostics)
    tx = g3['cx']
    line(ax, [tx, tx], [g3['y0'], BUS_Y])
    xs = [b['cx'] for b in e]
    line(ax, [xs[0], xs[1]], [BUS_Y, BUS_Y])
    line(ax, [xs[1], xs[-1]], [BUS_Y, BUS_Y], dashed=True)
    assert xs[0] < tx < xs[1]
    for i, b in enumerate(e):
        arrow(ax, (b['cx'], BUS_Y), (b['cx'], b['y1']), dashed=i >= 2)
    # evidence -> result
    for b in e:
        arrow(ax, (b['cx'], b['y0']), (b['cx'], res['y1']))
    save(fig, 'Fig2')


if __name__ == '__main__':
    main()
