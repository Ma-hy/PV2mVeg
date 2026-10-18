"""Fig. S3: installation year detected from Landsat compared with the station label year.   [restricted inputs]

Historical figure: drawn by the figure code of an earlier release on its original inputs; only the text on the
figure follows the v2.0 house style (pv2mveg.plotting.historical).

(a) Pixel-level distribution of detected minus label year by label year (east-west array pixels); (b) station-level
modal detected minus label year (frozen detection summaries, legacy/figS2a_*.csv and figS2b_*.csv);
(c) Landsat midsummer true-colour composites (ls_annual/ls_<year>.tif, bands sum_blue/green/red) for four stations,
T-2..T+2 around the station's modal detected year T, with the station outline and the cumulative share of the
station's array pixels detected by each year (ts/detect_30m.tif).

Inputs (project layout only): legacy/figS2a_detect_minus_label.csv, figS2b_station_mode_diff.csv,
ts/detect_station.csv, ts/detect_30m.tif, grid/panel_frac_{06,08}_30m.tif, grid/station_id_30m.tif,
ls_annual/ls_<year>.tif (public Landsat, derived), the labelled station outlines (SelectedStations.shp).
Output: FigS3.png/.pdf; figdata/figS3_selected_examples.csv, figS3_station_example_audit.json
"""
import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style_v04 import *
from pv2mveg.plotting import historical
from pv2mveg.plotting import fonts
from pv2mveg.plotting import profile as _v2

import rasterio
from rasterio.features import rasterize
from rasterio.windows import Window, from_bounds
import geopandas as gpd
from matplotlib.patches import Rectangle

P.need_project('Fig. S3')
# drawing defaults of the earlier example code that drew this figure (font size 9 pt, boxed axes)
plt.rcParams.update({'font.family': fonts.FONT_NAME, 'font.size': 9, 'axes.spines.top': True, 'axes.spines.right': True,
                     'pdf.fonttype': 42, 'axes.unicode_minus': True})
C_DIFF = {'share_m2minus': ('≤ −2', '#2166AC'), 'share_m1': ('−1', '#92C5DE'), 'share_0': ('0', '#D9D9D9'),
          'share_p1': ('+1', '#F4A582'), 'share_p2plus': ('≥ +2', '#B2182B')}
C_DET = '#D55E00'
ECO = P.eco


def _font(fig):
    for txt in fig.findobj(matplotlib.text.Text):
        txt.set_fontproperties(fonts.properties(txt.get_fontsize()))


def scale_bar(ax, length, label, loc=(0.05, 0.05), height_frac=0.012, color='k', fontsize=FS, bg=True,
              segments=2):
    """Alternating black/white scale bar in data units (metres). loc = lower-left in axes fraction."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    w = x1 - x0
    h = y1 - y0
    bx = x0 + loc[0] * w
    by = y0 + loc[1] * h
    bh = height_frac * h
    seg = length / segments
    if bg:
        ax.add_patch(Rectangle((bx - 0.015 * w, by - 0.025 * h), length + 0.03 * w, bh + 0.085 * h,
                               facecolor='white', edgecolor='none', alpha=0.8, zorder=20))
    for i in range(segments):
        ax.add_patch(Rectangle((bx + i * seg, by), seg, bh, facecolor=color if i % 2 == 0 else 'white',
                               edgecolor=color, lw=0.5, zorder=21))
    ax.text(bx + length / 2, by + bh + 0.008 * h, label, ha='center', va='bottom', fontsize=fontsize, zorder=22)


def stretch(imgs):
    allv = np.concatenate([i[np.isfinite(i).all(-1)] for i in imgs])
    lo = np.percentile(allv, 1, axis=0)
    hi = np.percentile(allv, 99, axis=0)
    out = []
    for i in imgs:
        v = np.clip((i - lo) / (hi - lo), 0, 1) ** 0.85
        v[~np.isfinite(i).all(-1)] = 1
        out.append(v)
    return out


def station_cube(sid, geom, years):
    # an integer window and its actual bounds, rather than a fractional requested bounding box on the rounded read
    x0, y0, x1, y1 = geom.bounds; cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    half = np.ceil((max(x1 - x0, y1 - y0) / 2 * 1.25 + 60) / 30) * 30
    requested = (cx - half, cy - half, cx + half, cy + half)
    path = ECO / 'ls_annual' / f'ls_{years[0]}.tif'
    with rasterio.open(path) as ds:
        raw = from_bounds(*requested, ds.transform)
        c0, r0 = int(np.floor(raw.col_off)), int(np.floor(raw.row_off))
        size = int(np.ceil(max(raw.width, raw.height)))
        window = Window(c0, r0, size, size); bounds = ds.window_bounds(window); transform = ds.window_transform(window)
        reference = (ds.crs, ds.transform, ds.width, ds.height)
    images = []
    for year in years:
        with rasterio.open(ECO / 'ls_annual' / f'ls_{year}.tif') as ds:
            assert (ds.crs, ds.transform, ds.width, ds.height) == reference
            assert list(ds.descriptions[:3]) == ['sum_blue', 'sum_green', 'sum_red']
            a = ds.read([3, 2, 1], window=window, boundless=True, fill_value=-32768).astype(float)
        a[a == -32768] = np.nan; images.append(np.moveaxis(a / 1e4, 0, -1))

    def read(p):
        with rasterio.open(p) as ds:
            assert (ds.crs, ds.transform, ds.width, ds.height) == reference
            return np.nan_to_num(ds.read(1, window=window, boundless=True, fill_value=0))
    detected = read(ECO / 'ts/detect_30m.tif')
    fmax = np.maximum(read(ECO / 'grid/panel_frac_06_30m.tif'), read(ECO / 'grid/panel_frac_08_30m.tif'))
    station = read(ECO / 'grid/station_id_30m.tif')
    polygon = rasterize([(geom, 1)], out_shape=station.shape, transform=transform).astype(bool)
    array = (fmax > .2) & (station == sid)
    meta = dict(station=sid, crs=str(reference[0]), transform=list(reference[1])[:6],
                requested_bounds=requested, plotted_bounds=bounds, integer_window=[c0, r0, size, size],
                station_pixels=int((station == sid).sum()), polygon_pixels=int(polygon.sum()),
                polygon_station_agreement=int(((station == sid) & polygon).sum()), array_pixels=int(array.sum()),
                array_fraction_of_polygon=float(array.sum() / polygon.sum()),
                temporal_grids_identical=True)
    return images, bounds, detected, array, meta


def main():
    # panels a/b are frozen detection summaries, not re-estimated statistics
    W, H = 190, 208
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    def axmm(l, b, w, h): return fig.add_axes([l / W, b / H, w / W, h / H])
    ta = pd.read_csv(P.legacy / 'figS2a_detect_minus_label.csv')
    ta = ta[ta.build_year != 'All arrays (both orientations)'].copy()
    y = np.arange(len(ta))[::-1].astype(float); y[ta.build_year == 'All'] = -1
    ax = axmm(16, 158, 100, 40); left = np.zeros(len(ta))
    for key, (label, color) in C_DIFF.items():
        ax.barh(y, ta[key] * 100, left=left, color=color, height=.75, edgecolor='white', lw=.3, label=label)
        left += ta[key].values * 100
    for yy, n in zip(y, ta.n_px): ax.text(101, yy, f'{n:,}', ha='left', va='center', fontsize=8)
    ax.text(101, y.max() + .9, 'Pixels', fontsize=8)
    ax.set(yticks=y, yticklabels=ta.build_year, xlim=(0, 100), ylim=(-1.6, y.max() + .6), xlabel='Share of detected array pixels (%)', ylabel='Label year')
    ax.axhline(-.25, color='k', lw=.4)
    ax.legend(title='Detected − label year', ncol=5, loc='lower left', bbox_to_anchor=(0, 1.0), fontsize=8, title_fontsize=8, handlelength=1, columnspacing=.8, handletextpad=.3)
    panel_label(ax, 'a', x=-.085, y=1.13)
    for spine in ['top', 'right']: ax.spines[spine].set_visible(False)
    tb = pd.read_csv(P.legacy / 'figS2b_station_mode_diff.csv', index_col=0).iloc[:, 0]
    ax = axmm(140, 158, 46, 40)
    ax.bar(np.arange(5), tb.values, color=[v[1] for v in C_DIFF.values()], edgecolor='k', lw=.3, width=.75)
    for j, v in enumerate(tb.values): ax.text(j, v + 2, str(v), ha='center', fontsize=8)
    ax.set(xticks=range(5), xticklabels=['≤ −2', '−1', '0', '+1', '≥ +2'], xlabel='Modal detected − label year', ylabel='Stations', ylim=(0, tb.max() * 1.15))
    ax.text(.98, .97, f'n = {int(tb.sum())} stations', transform=ax.transAxes, ha='right', va='top', fontsize=8)
    panel_label(ax, 'b', x=-.30, y=1.13)
    for spine in ['top', 'right']: ax.spines[spine].set_visible(False)
    station_table = pd.read_csv(ECO / 'ts/detect_station.csv').set_index('station')
    shapes = gpd.read_file(P.selected_stations); assert shapes.crs.to_epsg() == 32647
    records = []; diagnostics = []
    for r, sid in enumerate([4, 85, 94, 138]):
        T = int(station_table.loc[sid, 'detect_mode']); L = int(station_table.loc[sid, 'build_year']); years = list(range(T - 2, T + 3))
        geom = shapes.geometry.iloc[sid - 1]; assert int(str(shapes.O_Name.iloc[sid - 1])[:4]) == L
        imgs, b, dy, arr, meta = station_cube(sid, geom, years); diagnostics.append(meta); imgs = stretch(imgs)
        yb = 140 - (r + 1) * 30 - r * 6
        for c, year in enumerate(years):
            ax = axmm(30 + c * 31.5, yb, 30, 30)
            ax.imshow(imgs[c], extent=(b[0], b[2], b[1], b[3]), interpolation='nearest')
            for poly in getattr(geom, 'geoms', [geom]): ax.plot(*poly.exterior.xy, color='#F0E442', lw=.6)
            ax.set(xlim=(b[0], b[2]), ylim=(b[1], b[3]), xticks=[], yticks=[])
            for spine in ax.spines.values(): spine.set_visible(True); spine.set_linewidth(1.6 if year == T else .5); spine.set_color(C_DET if year == T else 'k')
            ax.set_title(f'{year}' + (' (label)' if year == L else ''), fontsize=8, pad=1.5, color=C_DET if year == T else 'k')
            share = float((arr & (dy > 0) & (dy <= year)).sum() / arr.sum())
            ax.text(.03, .03, f'{share * 100:.0f}%', transform=ax.transAxes, fontsize=8, bbox=dict(boxstyle='square,pad=.15', fc='white', ec='none', alpha=.85))
            if c == 0:
                length = 500 if b[2] - b[0] < 2500 else 1000
                scale_bar(ax, length, f'{length} m' if length < 1000 else '1 km', loc=(.5, .04), segments=1, bg=True)
            records.append(dict(station=sid, label_year=L, detect_mode=T, year=year, cum_share_detected=share, n_array_px=int(arr.sum())))
        fig.text(28 / W, (yb + 15) / H, f'Station {sid}\nLabel {L}\nDetected {T}\n({station_table.loc[sid, "mode_share"] * 100:.0f}% of\ndetected px)', ha='right', va='center', fontsize=8)
    fig.text(14 / W, 145 / H, '(c)', fontsize=9)
    log = historical.finish(fig, [], FS_LABEL, tight=True)
    _v2.save(fig, 'FigS3', _font, 'rev', bbox_inches='tight')
    for a, b_ in sorted(set(log)):
        print('   text:', repr(a), '->', repr(b_))
    # the questioned example of the earlier release stays in the diagnostic record
    sid = 30; T = int(station_table.loc[sid, 'detect_mode'])
    _, _, _, _, old_meta = station_cube(sid, shapes.geometry.iloc[sid - 1], list(range(T - 2, T + 3)))
    old_meta.update(display_disposition='Replaced by Station 85 for clearer station/array correspondence; label and detection remain frozen.')
    diagnostics.append(old_meta)
    pd.DataFrame(records).to_csv(Path(FIGDATA) / 'figS3_selected_examples.csv', index=False)
    (Path(FIGDATA) / 'figS3_station_example_audit.json').write_text(json.dumps(diagnostics, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
