"""Fig. 1: the landscape, one 150 m unit, and what is measured between the rows.   [restricted inputs]

(a) Beijing-3 true colour (2025-08) with the 7,544 east-west analysis units coloured by construction year and the
    hand-drawn station outlines; the white square locates (b) and is joined to it.
(b) The 150 m unit at H086 in native 0.5 m panchromatic imagery (2026-08): mapped row bands and complete native 2 m
    ground footprints nested by centre clearance.
(c) Cross-section across three rows, drawn to scale from the median geometry (10 m pitch, ~5.2 m row band, 1.5 m
    centre clearance): what the 0.5 m map, the native 2 m footprints and a 30 m Landsat pixel each see.

Inputs (project layout only): unit coordinates (ts/cells_pool.parquet), construction years (v12/inputs), the
true-colour display cache (work/maps_v110), station polygons, the 0.5 m panchromatic scene, source and row-band
rasters and the 2 m footprint layers of the row-band release.
Output: Fig1.png/.pdf; figdata/fig1_counts.json
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from pv2mveg import data as D
from pv2mveg import geometry as G
from matplotlib.patches import ConnectionPatch

EX_TAG, EX_CELL = '2026_08', 52538
SEL = {'1m': '#F6A8B8', '1p5m': '#E23D5B', '2m': '#7E0E2C'}
BAND_AMBER = '#F2A31A'
# 30 m and 10 m pixel squares of (b): tones that read on the image and in the key on white
PIX30, PIX10 = '#E0A800', '#0096B4'


def parent_units():
    P.need_project('Fig. 1')
    # drawing defaults of the published figure (see style.MAP_RC): the map and age modules of the research code
    # were imported here
    map_defaults(); plt.rcParams.update(AGE_RC)
    tidx, T, _, origin = D.merged()
    t = D.cells_pool().set_index('cell', drop=False).loc[tidx].copy()
    t['T'] = T; t['origin'] = np.where(origin == 1, 'phased', 'original')
    t['lon'], t['lat'] = GEO.transform(t.x.to_numpy(), t.y.to_numpy())
    assert len(t) == 7544
    return t


def panel_a(ax, t):
    ext = (float(t.lon.min() - .016), float(t.lon.max() + .016), float(t.lat.min() - .013), float(t.lat.max() + .013))
    rgb, meta = rgb_cache(ext, 2200)
    a = rgb.astype(float) / 255.
    a = np.clip((a - .04) * 1.08, 0, 1) ** .95
    a[(rgb == 232).all(axis=2)] = 1.0
    ax.imshow(a, extent=ext, origin='upper', zorder=0, interpolation='bilinear')
    ax.add_collection(PolyCollection(cell_polygons(t.x, t.y), array=t['T'].to_numpy().astype(float), cmap=YEAR_CMAP, norm=YEAR_NORM,
                                     edgecolor='none', alpha=.86, zorder=2, rasterized=True))
    st = stations()
    ax.add_collection(PolyCollection(polygon_paths(st), facecolor='none', edgecolor='white', linewidths=.35, alpha=.75, zorder=3))
    lonlat_frame(ax, ext, step=.1)
    # the panel corners lie outside the image footprint (white): dark scale bar and north needle
    scalebar(ax, 10, x=.05, y=.07, color=INK, text_color=INK, ext=ext)
    north_arrow(ax, x=.06, y=.80)
    sm = plt.cm.ScalarMappable(norm=YEAR_NORM, cmap=YEAR_CMAP)
    # no shaded boxes: the colour key is drawn in black in the white corner outside the image footprint (checked
    # below against the image); the unit count and grid area are given in the caption and the text
    cax = ax.inset_axes([.665, .085, .30, .026])
    cb = ax.figure.colorbar(sm, cax=cax, orientation='horizontal')
    cb.set_ticks([2012, 2016, 2020, 2024]); cb.set_ticklabels(['≤2012', '2016', '2020', '2024'])
    cb.ax.tick_params(labelsize=5.8, length=1.6, pad=1, colors=INK); cb.outline.set_edgecolor(INK); cb.outline.set_linewidth(.4)
    cax.set_title('construction year', fontsize=5.9, color=INK, pad=2.5)
    cax.set_zorder(10)
    white = (rgb == 232).all(axis=2)
    nr, nc = white.shape
    x0, x1, y0, y1 = .645, .985, .0, .17                       # key with its title and tick labels (axes fraction)
    under = white[int((1 - y1) * nr):int((1 - y0) * nr), int(x0 * nc):int(x1 * nc)]
    assert under.mean() > .97, f'colour key of (a) overlaps the image ({under.mean():.3f} white)'
    p = t.set_index('cell').loc[EX_CELL]
    return ext, (p.lon, p.lat)


def panel_b(ax):
    pan, mask, inside, sel = G.native_example(EX_TAG, EX_CELL)
    ext = (0, 192, 0, 192)
    ax.imshow(pan, extent=ext, cmap='gray', vmin=0, vmax=255, origin='upper', zorder=0, interpolation='nearest')
    rb = np.zeros((*mask.shape, 4)); rb[mask == 1] = [*to_rgb(BAND_AMBER), .5]
    ax.imshow(rb, extent=ext, origin='upper', interpolation='nearest', zorder=1)
    for p in ['1m', '1p5m', '2m']:
        m = sel[p]
        rgba = np.zeros((*m.shape, 4)); rgba[m] = [*to_rgb(SEL[p]), .97]
        ax.imshow(rgba, extent=ext, origin='upper', interpolation='nearest', zorder=2)
    ax.contour(inside.astype(float), levels=[.5], colors='white', linewidths=.9, extent=ext, origin='upper', zorder=4)
    counts = {p: int((sel[p] & inside).sum()) for p in sel}
    # 30 m and 10 m supports for scale (aligned to the window origin; illustrative); named in the key below
    ax.add_patch(Rectangle((8, 150), 30, 30, fill=False, ec=PIX30, lw=1.3, ls=(0, (3, 1.5)), zorder=5))
    ax.add_patch(Rectangle((8, 150), 10, 10, fill=False, ec=PIX10, lw=1.2, zorder=5))
    ax.set_xlim(0, 192); ax.set_ylim(0, 192); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True); s.set_color('#6B7785'); s.set_linewidth(.5)
    return counts


def panel_c(ax):
    """To-scale cross-section (metres), south (left) to north (right), across three rows.

    Typical geometry: pitch 10.1 m, row band 5.2 m (module 3.2 m in plan + 2.0 m dark strip), low/high module edges
    0.8/3.0 m; the native 2 m lattice phase is fixed, so the three gaps show how the same rule selects footprints at
    different clearances.
    """
    pitch, band, low, high, proj, origin = 10.1, 5.2, .8, 3.0, 3.2, -1.8
    rows = [0., pitch, 2 * pitch]
    x_lo, x_hi = -4.0, 28.4
    y_hi = 5.8  # sky above the drawing for the footprint key
    ax.set_xlim(x_lo, x_hi); ax.set_ylim(-4.1, y_hi); ax.set_aspect('equal')
    ax.axis('off')
    sky = np.linspace(0, 1, 50)[:, None] * np.ones((1, 10))
    ax.imshow(sky, extent=(x_lo, x_hi, 0, y_hi), cmap=LinearSegmentedColormap.from_list('sky', ['#EAF1F7', '#FFFFFF']), origin='lower', aspect='auto', zorder=0)
    ax.add_patch(Rectangle((x_lo, -.32), x_hi - x_lo, .32, fc='#D9CBAA', ec='none', zorder=1))
    ax.plot([x_lo, x_hi], [0, 0], color='#8C7B5A', lw=.8, zorder=2)
    for k in range(4):
        x0 = 4.4 + k * 1.25
        ax.annotate('', xy=(x0 + 1.15, 2.85), xytext=(x0, 3.85), arrowprops=dict(arrowstyle='-|>,head_width=.12,head_length=.3', color='#F2B705', lw=.6, alpha=.8), zorder=3)
    ax.text(4.2, 4.05, 'sun from the south', fontsize=5.7, color='#A87A06', va='center', ha='left')
    for xx in np.linspace(-3.8, 28.2, 67):
        for dx_, h_ in [(-.05, .15), (0, .21), (.05, .14)]:
            ax.plot([xx, xx + dx_ * 2.4], [0, h_], color='#A3B57A', lw=.42, alpha=.8, zorder=2)
    for x0 in rows:
        ax.add_patch(Rectangle((x0 + proj, -.02), band - proj, .07, fc='#30353C', ec='none', alpha=.92, zorder=3))
        ax.add_patch(Polygon([[x0 + proj + .05, 0], [x0 + band, 0], [x0 + proj, high]], closed=True, fc='#5B6470', alpha=.11, ec='none', zorder=2))
        for xp in [x0 + .45, x0 + 2.7]:
            hp = low + (high - low) * (xp - x0) / proj
            ax.plot([xp, xp], [0, hp - .05], color='#7D8794', lw=1.0, zorder=3)
        ax.add_patch(Polygon([[x0, low], [x0 + proj, high], [x0 + proj + .07, high - .1], [x0 + .07, low - .1]], closed=True, fc='#1E3A5F', ec='#0F2238', lw=.4, zorder=4))
        ax.plot([x0 + .05, x0 + proj], [low + .03, high + .02], color='#86B8EA', lw=.6, zorder=5)
    # 0.5 m row-band map
    y0m, hm = -1.0, .45
    for j in range(int((x_hi - x_lo) / .5)):
        a = x_lo + j * .5
        inb = any((a + .25 >= r0) and (a + .25 < r0 + band) for r0 in [-pitch] + rows)
        ax.add_patch(Rectangle((a, y0m), .5, hm, fc=BAND_AMBER if inb else '#F3EBDC', ec='white', lw=.25, zorder=3))
    # black on the amber band (white had a contrast of about 2:1)
    ax.text(band / 2, y0m + hm / 2, 'row band ≈ 5.2 m', color=INK, fontsize=5.6, ha='center', va='center', zorder=5)
    ax.text(band + (pitch - band) / 2, y0m + hm / 2, 'inter-row ground', color=INK, fontsize=5.6, ha='center', va='center', zorder=5)
    # native 2 m footprints
    y0f, hf = -1.95, .65
    for j in range(-3, 20):
        a = origin + j * 2; b = a + 2; c = (a + b) / 2
        if b <= x_lo or a >= x_hi:
            continue
        on_ground = any(a >= r0 + band and b <= r0 + pitch for r0 in [-pitch] + rows)
        dist = min(min(abs(c - (r0 + band)), abs(c - (r0 + pitch))) for r0 in [-pitch] + rows)
        fc = SEL['2m'] if on_ground and dist >= 2 else SEL['1p5m'] if on_ground and dist >= 1.5 else SEL['1m'] if on_ground and dist >= 1 else '#E9ECEF'
        ax.add_patch(Rectangle((max(a, x_lo), y0f), min(b, x_hi) - max(a, x_lo), hf, fc=fc, ec='white', lw=.6, zorder=3))
        if on_ground and dist >= 1:
            ax.scatter([c], [y0f + hf * .5], s=2.2, color='white', zorder=4, linewidths=0)
    # the primary rule in the second gap: a footprint is kept when its centre (white dot) lies at least 1.5 m from
    # both row-band edges, i.e. between the dashed lines
    g0, g1 = rows[1] + band, rows[2]
    top, yd = y0m + hm, y0f - .32
    for xx in (g0, g1):
        ax.plot([xx, xx], [yd - .12, top], color='#9C7A2E', lw=.55, zorder=4)
    for xx in (g0 + 1.5, g1 - 1.5):
        ax.plot([xx, xx], [yd - .12, top], color='#8E1230', lw=.8, ls=(0, (3, 1.5)), zorder=4)
    for xa, xb in ((g0, g0 + 1.5), (g1 - 1.5, g1)):
        ax.annotate('', xy=(xb, yd), xytext=(xa, yd), arrowprops=dict(arrowstyle='<->,head_width=.1,head_length=.2', color='#8E1230', lw=.6, shrinkA=0, shrinkB=0), zorder=5)
        ax.text((xa + xb) / 2, yd - .1, '1.5 m', ha='center', va='top', fontsize=5.4, color='#8E1230')
    ax.text(g0 - .3, yd, 'primary rule: centre ≥1.5 m\nfrom the row band (dashed)', ha='right', va='center', fontsize=5.4, color='#8E1230', linespacing=1.05)
    # 30 m Landsat pixel (exactly 30 m)
    lx0 = -1.6
    ax.add_patch(Rectangle((lx0, -3.55), 30, .42, fc=MIXED_L, ec='none', alpha=.6, zorder=2))
    ax.annotate('', xy=(lx0, -3.0), xytext=(lx0 + 30, -3.0), arrowprops=dict(arrowstyle='|-|,widthA=.3,widthB=.3', color=MIXED, lw=.8))
    ax.text(lx0 + 15, -3.34, 'one 30 m Landsat pixel mixes three rows of modules, dark strips, posts and ground', ha='center', va='center', fontsize=5.9, color=INK, zorder=3)
    # row labels on the left
    for y, lab in [(y0m + hm / 2, '0.5 m row-band map'), (y0f + hf / 2, 'native 2 m footprints'), (-3.34, 'mixed 30 m pixel')]:
        v2fix(ax.text(x_lo - .25, y, lab, ha='right', va='center', fontsize=6.0, color=INK), 6.8)
    ax.annotate('module', xy=(1.4, low + (high - low) * 1.4 / proj + .1), xytext=(-.9, 3.3), fontsize=5.9, color='#1E3A5F',
                arrowprops=dict(arrowstyle='-', color='#1E3A5F', lw=.5), ha='center')
    ax.annotate('dark strip', xy=(rows[1] + 4.3, .05), xytext=(rows[1] + 5.5, 1.25), fontsize=5.9, color=INK,
                arrowprops=dict(arrowstyle='-', color='#30353C', lw=.5), ha='left')
    ax.annotate('', xy=(rows[2], 3.9), xytext=(rows[1], 3.9), arrowprops=dict(arrowstyle='<->,head_width=.1,head_length=.25', color=SUB, lw=.5))
    ax.text((rows[1] + rows[2]) / 2, 4.0, 'pitch ≈ 10 m', ha='center', va='bottom', fontsize=5.8, color=SUB)
    hh = [Patch(fc=SEL['2m'], label='≥2 m'), Patch(fc=SEL['1p5m'], label='1.5–2 m'), Patch(fc=SEL['1m'], label='1–1.5 m'), Patch(fc='#E9ECEF', label='excluded')]
    # one row above the drawing, clear of the pitch arrow and the third module
    ax.legend(handles=hh, loc='upper right', bbox_to_anchor=(1.0, 1.0), ncol=4, fontsize=7.0, handlelength=1.1, handletextpad=.3, columnspacing=.9,
              title='2 m footprint centre from the row band', title_fontsize=7.0, frameon=False, borderaxespad=0, alignment='right')


def main():
    t = parent_units()
    fig = figure(150)
    axa = fig.add_axes([.06, .425, .585, .53])
    # (b) leaves room below for its keys; (c) starts far enough right for its row labels
    axb = fig.add_axes([.69, .548, .295, .39])
    axc = fig.add_axes([.15, .0, .84, .36])
    ext, (lon0, lat0) = panel_a(axa, t)
    counts = panel_b(axb)
    panel_c(axc)
    # locator square and zoom connectors
    dlon = .0035 / np.cos(np.deg2rad(lat0)) * 1.25; dlat = .0035 * 1.25
    axa.add_patch(Rectangle((lon0 - dlon, lat0 - dlat), 2 * dlon, 2 * dlat, fill=False, ec='white', lw=1.0, zorder=12))
    for (xa, ya), (xb, yb) in [((lon0 + dlon, lat0 + dlat), (0, 1)), ((lon0 + dlon, lat0 - dlat), (0, 0))]:
        cp = ConnectionPatch(xyA=(xa, ya), coordsA=axa.transData, xyB=(xb, yb), coordsB=axb.transAxes, color='#4A5866', lw=.55, alpha=.85, zorder=20)
        fig.add_artist(cp)
    letter(axa, 'a', x=-.075, y=1.005, ha='left')
    axa.text(-.03, 1.005, 'Talatan: east–west fixed arrays by construction year', transform=axa.transAxes, fontsize=7.6, va='bottom', path_effects=bold(INK, .2))
    letter(axb, 'b', x=-.03, y=1.025, ha='left')
    axb.text(.07, 1.025, 'One 150 m unit (H086, 2026-08)', transform=axb.transAxes, fontsize=7.3, va='bottom', path_effects=bold(INK, .2))
    hb = [Patch(fc=SEL['1m'], label='1–1.5 m'), Patch(fc=SEL['1p5m'], label='1.5–2 m'), Patch(fc=SEL['2m'], label='≥2 m')]
    # the keys below (b) are set at 7 pt, the long title on two lines and the second key in one column
    lg = axb.legend(handles=hb, loc='upper left', bbox_to_anchor=(-.012, -.01), ncol=3, fontsize=7.0, handlelength=1.1, handletextpad=.35, columnspacing=1.1,
                    title='2 m ground footprints,\ncentre from the row band', title_fontsize=7.0, frameon=False, alignment='left', borderaxespad=0)
    axb.add_artist(lg)
    ub = Line2D([], [], color='white', lw=1.0, path_effects=[pe.Stroke(linewidth=2.2, foreground='#6B7785'), pe.Normal()], label='unit boundary')
    pix = [Patch(fill=False, ec=PIX30, lw=1.3, ls=(0, (3, 1.5)), label='30 m Landsat pixel'),
           Patch(fill=False, ec=PIX10, lw=1.2, label='10 m Sentinel-2 pixel')]
    lg2 = axb.legend(handles=[Patch(fc=(*to_rgb(BAND_AMBER), .55), ec='none', label='row band (module + dark strip)'), ub] + pix, loc='upper left', bbox_to_anchor=(-.012, -.15),
                     ncol=1, fontsize=7.0, handlelength=1.1, handletextpad=.35, columnspacing=1.1, frameon=False, borderaxespad=0, labelspacing=.3)
    tb = axb.text(.0, -.235, f'footprints kept in this unit:\n≥1 m {counts["1m"]:,} · ≥1.5 m {counts["1p5m"]:,} (primary) · ≥2 m {counts["2m"]:,}',
                  transform=axb.transAxes, fontsize=5.8, color=SUB, va='top', linespacing=1.15)
    # stack the two keys and the counts under (b) with 0.8 mm gaps
    r = fig.canvas.get_renderer(); fig.canvas.draw()
    gap = .8 / (axb.get_position().height * 150)
    below = lambda a: axb.transAxes.inverted().transform((0, a.get_window_extent(r).y0))[1] - gap
    lg2.set_bbox_to_anchor((-.012, below(lg)), transform=axb.transAxes); fig.canvas.draw()
    tb.set_y(below(lg2))
    # the letter and title of (c) keep their place on the page although the axes start further right
    letter(axc, 'c', x=-.1482, y=.985, ha='left')
    axc.text(-.1173, .985, 'What is measured between the rows (to scale, typical geometry)', transform=axc.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    save(fig, 'Fig1')
    dump(FIGDATA / 'fig1_counts.json', dict(example=dict(tag=EX_TAG, cell=EX_CELL, **counts), parents=len(t)))
    print(counts)


if __name__ == '__main__':
    main()
