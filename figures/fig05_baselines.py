"""Fig. 5: matched-control and unmatched-buffer baselines.   [restricted inputs]

(a) Every usable matched pair of 2025-08 drawn as a curved link from the treated unit (cyan) to its control (amber)
    on a night-mode image.
(b, c) Two units seen through a 5 km lens: their matched controls and the unselected buffer ground (300-1,000 m from
    arrays, same scene, within 5 km).
(d) Unit-level ΔNDVI under the two baselines, 2025-08.
(e) Age-group means under both baselines in the three acquisitions.
Examples are chosen by geometry only: the 10th-percentile matched-set spread overall, and within >=6-year units of
the north-eastern cluster.

Inputs (project layout only): unit coordinates and screens (ts/cells_pool.parquet), the Landsat unit series
(ts/cell_y.npz) for the buffer membership, the matched sets (v12/inputs), the control-ground table of 2025-08
(analysis/1p5m/chrono/celltab_08_raw.parquet), the true-colour display cache; v12/buffer_primary and
v12/cross_scene results.
Output: Fig5.png/.pdf; figdata/fig4e_designs.csv, fig4_summary.json
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from pv2mveg import data as D

NIGHT = '#0A1520'
CAPTIONS = []  # captions under the two lenses (the buffer-ground note is placed under the second)


def load():
    P.need_project('Fig. 5')
    cells = D.cells_pool().set_index('cell', drop=False)
    cells['lon'], cells['lat'] = GEO.transform(cells.x.to_numpy(), cells.y.to_numpy())
    units = {t: pd.read_csv(BP12 / f'cells_{t}.csv') for t in TAGS}
    tidx_, _, M_, _ = D.merged(); z = dict(tidx=tidx_, M=M_)
    cs = pd.read_csv(CS12 / 'cells_08.csv')
    keep = set(cs.cell[cs.valid_cross])
    raw = pd.read_parquet(P.analysis / '1p5m/chrono/celltab_08_raw.parquet')
    okc = raw.ctrl_n.values >= 2000
    rows = [(t, m) for t, mm in zip(z['tidx'], z['M']) if t in keep for m in mm if m >= 0 and okc[m]]
    pairs = pd.DataFrame(rows, columns=['cell', 'ctrl']); pairs['tag'] = '08'
    units['08s'] = cs
    return cells, pairs, units


def choose_examples(cells, pairs, units):
    q = pairs[pairs.tag == '08']
    u = units['08'].set_index('cell')
    def pick(mask):
        scores = []
        for cell, g in q.groupby('cell'):
            if len(g) < 7 or not mask(cell):
                continue
            pts = cells.loc[np.r_[cell, g.ctrl.to_numpy()]]
            scores.append((max(pts.x.max() - pts.x.min(), pts.y.max() - pts.y.min()), int(cell)))
        scores.sort()
        return scores[len(scores) // 10][1]
    # 60587: 10th-percentile matched-set spread; 52539: 10th-percentile older (>=6 y) north-eastern unit whose
    # 5 km lens lies inside the imagery
    first = 60587 if 60587 in set(q.cell) else pick(lambda c: True)
    second = 52539 if 52539 in set(q.cell) else pick(lambda c: (u.loc[c, 'age'] >= 6) and (cells.loc[c, 'lon'] > 100.5) and (c != first))
    return [first, second]


def buffer_members(cells, unit_row, tag='08', radius=5000.):
    """Recompute the main-variant buffer cells of one unit (as in the buffer baseline)."""
    Y = np.load(TS / 'cell_y.npz')
    nvalid = np.isfinite(Y['sum_ndvi']).sum(1)
    c = cells.sort_index()
    q = (c.in_aoi.values & (nvalid >= 20) & (c.station_any.values == 0) & (c.strict_n.values == 0) & (c.wc40.values < .05) & (c.wc50.values < .05)
         & (c.wc80.values < .01) & (c.wc90.values < .05) & (c.wc10.values < .1) & (c.water_frac_max.values < .05))
    raw = pd.read_parquet(P.analysis / '1p5m/chrono' / f'celltab_{tag}_raw.parquet')
    ok = q & (raw.o300_1000_n.values >= 2000) & (raw.ctrl_scene.values == unit_row.reference_scene)
    d = np.hypot(c.x.values - unit_row.x, c.y.values - unit_row.y)
    idx = np.flatnonzero(ok & (d <= radius))
    assert len(idx) == unit_row.n_buffer_cells, (len(idx), unit_row.n_buffer_cells)
    return c.iloc[idx]


def nightify(rgb, strength, tint, desat):
    a = dark_basemap(rgb, strength=strength, tint=tint, desat=desat)
    nodata = (rgb == 232).all(axis=2)
    a[nodata] = np.array(to_rgb(NIGHT))
    return a


def night(ax, ext, width):
    rgb, _ = rgb_cache(ext, width)
    ax.imshow(nightify(rgb, .60, (.006, .012, .02), .12), extent=(ext[0], ext[1], ext[2], ext[3]), origin='upper', zorder=0, interpolation='bilinear')


def panel_a(ax, cells, pairs, units):
    q = pairs[pairs.tag == '08']
    t = cells.loc[np.unique(q.cell)]
    ctrl = cells.loc[np.unique(q.ctrl)]
    allc = cells.loc[np.unique(np.r_[t.cell.values, q.ctrl.values])]
    pad = .012
    ext = (allc.lon.min() - pad, allc.lon.max() + pad, allc.lat.min() - pad * .8, allc.lat.max() + pad * .8)
    night(ax, ext, 2200)
    p0 = cells.loc[q.cell, ['lon', 'lat']].to_numpy(); p1 = cells.loc[q.ctrl, ['lon', 'lat']].to_numpy()
    arcs = bezier(p0, p1, bend=.18, n=16)
    ax.add_collection(PolyCollection(cell_polygons(ctrl.x, ctrl.y), facecolor=CTRL_GLOW, edgecolor='none', alpha=.55, zorder=3, rasterized=True))
    ax.add_collection(PolyCollection(cell_polygons(t.x, t.y), facecolor=TREAT_GLOW, edgecolor='none', alpha=.55, zorder=3, rasterized=True))
    gradient_arcs(ax, arcs, TREAT_GLOW, CTRL_GLOW, lw=.9, alpha=.010, z=4)
    gradient_arcs(ax, arcs, TREAT_GLOW, CTRL_GLOW, lw=.2, alpha=.075, z=5)
    lonlat_frame(ax, ext, step=.1)
    ax.set_facecolor(NIGHT)
    scalebar(ax, 10, x=.05, y=.06, color='white', ext=ext)
    north_arrow(ax, x=.05, y=.80, color='white', fill='#0E1A26')
    # the counts on two lines: one line would not fit the box, which stays clear of the data
    ax.add_patch(FancyBboxPatch((.555, .01), .43, .125, boxstyle='round,pad=0.004,rounding_size=0.01', transform=ax.transAxes, fc=NIGHT, ec='none', alpha=.72, zorder=8))
    ax.text(.975, .098, f'{DATE_LABEL["08"]} · matched design', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.6, color='white', zorder=9)
    ax.text(.975, .058, f'{len(t):,} treated units · {len(ctrl):,} controls', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.1, color='#B8C7D6', zorder=9)
    ax.text(.975, .02, f'{len(q):,} links', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.1, color='#B8C7D6', zorder=9)
    ax.text(.025, .015, 'background: original true colour, darkened for contrast', transform=ax.transAxes, ha='left', va='bottom', fontsize=5.4, color='#AEBBC8', path_effects=halo(1.0, NIGHT, .6))
    handles = [Line2D([], [], marker='s', ls='', ms=4, mfc=TREAT_GLOW, mec='none', label='treated unit'),
               Line2D([], [], marker='s', ls='', ms=4, mfc=CTRL_GLOW, mec='none', label='matched control'),
               Line2D([], [], color=CTRL_GLOW, lw=1.0, label='matching link')]
    ax.legend(handles=handles, loc='upper right', fontsize=6.1, labelcolor='white', frameon=True, framealpha=.55, facecolor=NIGHT, edgecolor='none', borderpad=.5, handletextpad=.3, labelspacing=.3)
    return dict(treated=len(t), controls=len(ctrl), links=len(q), extent=ext)


def lens(ax, cells, pairs, units, cell, label):
    u = pd.read_csv(BP12 / 'cells_08.csv').set_index('cell').loc[cell].copy()
    u['reference_scene'] = u.modal_scene
    q = pairs[(pairs.tag == '08') & (pairs.cell == cell)]
    buf = buffer_members(cells, u)
    lon0, lat0 = GEO.transform(u.x, u.y)
    r_km = 5.0
    dlat = r_km / 110.6; dlon = r_km / (111.32 * np.cos(np.deg2rad(lat0)))
    ext = (lon0 - dlon * 1.04, lon0 + dlon * 1.04, lat0 - dlat * 1.04, lat0 + dlat * 1.04)
    rgb, _ = rgb_cache(ext, 900)
    im = ax.imshow(nightify(rgb, .74, (.004, .008, .012), .08), extent=ext, origin='upper', zorder=0, interpolation='bilinear')
    ax.set_xlim(ext[:2]); ax.set_ylim(ext[2:]); ax.set_aspect(1 / np.cos(np.deg2rad(lat0)))
    # circle in display-consistent geometry: an ellipse in lon/lat
    from matplotlib.patches import Ellipse
    clip = Ellipse((lon0, lat0), 2 * dlon, 2 * dlat, transform=ax.transData)
    im.set_clip_path(clip)
    bcol = PolyCollection(cell_polygons(buf.x, buf.y), facecolor=BUFFER_GLOW, edgecolor='none', alpha=.46, zorder=2, rasterized=True)
    ax.add_collection(bcol); bcol.set_clip_path(clip)
    ctrl = cells.loc[q.ctrl]
    p0 = np.repeat([[lon0, lat0]], len(ctrl), 0); p1 = ctrl[['lon', 'lat']].to_numpy()
    lc = gradient_arcs(ax, bezier(p0, p1, bend=.2, n=30), TREAT_GLOW, CTRL_GLOW, lw=.9, alpha=.9, z=4); lc.set_clip_path(clip)
    cc_ = PolyCollection(cell_polygons(ctrl.x, ctrl.y), facecolor=CTRL_GLOW, edgecolor='white', linewidths=.25, zorder=5); ax.add_collection(cc_)
    tc = PolyCollection(cell_polygons([u.x], [u.y]), facecolor=TREAT_GLOW, edgecolor='white', linewidths=.4, zorder=6); ax.add_collection(tc)
    ax.add_patch(Ellipse((lon0, lat0), 2 * dlon, 2 * dlat, transform=ax.transData, fill=False, ec='#DDE6EE', lw=.8, zorder=7))
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    cap = ax.text(.5, -.02, f'unit {cell} · {int(u.age)} y · {len(ctrl)} controls · {int(u.n_buffer_cells)} buffer cells', transform=ax.transAxes,
                  ha='center', va='top', fontsize=6.8, color=SUB)
    CAPTIONS.append(v2fix(cap, 6.8))
    # 1 km tick on the rim
    km = 1 / (111.32 * np.cos(np.deg2rad(lat0)))
    ax.plot([lon0 - km, lon0 + km], [lat0 - dlat * .86] * 2, color='white', lw=1.1, zorder=8, solid_capstyle='butt')
    ax.text(lon0, lat0 - dlat * .82, '2 km', color='white', fontsize=5.6, ha='center', va='bottom', zorder=8, path_effects=halo(1.0, NIGHT, .6))
    ax.text(lon0, lat0 + dlat * .9, 'r = 5 km', color='#C9D4DE', fontsize=5.4, ha='center', va='top', zorder=8, path_effects=halo(1.0, NIGHT, .6))
    return dict(cell=int(cell), age=int(u.age), controls=len(ctrl), buffer_cells=int(u.n_buffer_cells))


def panel_d(ax, units):
    bp = pd.read_csv(BP12 / 'cells_08.csv')
    u = bp[bp.valid_buffer].copy()
    u['rematched_ndvi_dadj'] = u.primary_ndvi_dadj; u['expanded_ndvi_dadj'] = u.buffer_ndvi_dadj
    g = np.where(u.age <= 2, '0–2', np.where(u.age <= 5, '3–5', '≥6'))
    for grp in ['0–2', '3–5', '≥6']:
        m = g == grp
        ax.scatter(u.rematched_ndvi_dadj[m], u.expanded_ndvi_dadj[m], s=2.0, c=AGE[grp], alpha=.5, linewidths=0, zorder=3, rasterized=True, label=f'{grp} y')
    lim = (-.16, .26)
    ax.plot(lim, lim, color=INK, lw=.7, ls=(0, (3, 2)), zorder=2)
    ax.axhline(0, color='#B0B8C1', lw=.5, zorder=1); ax.axvline(0, color='#B0B8C1', lw=.5, zorder=1)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
    r = np.corrcoef(u.rematched_ndvi_dadj, u.expanded_ndvi_dadj)[0, 1]
    dmean = (u.expanded_ndvi_dadj - u.rematched_ndvi_dadj).mean()
    ax.text(.04, .97, f'{DATE_LABEL["08"]} · {len(u):,} units\nr = {r:.2f}' + '\n' + f'mean difference {signed(dmean, 4)}', transform=ax.transAxes, va='top', fontsize=6.1,
            color=INK, linespacing=1.15, path_effects=halo(1.4), zorder=6)
    ax.set_xlabel('ΔNDVI, matched controls'); ax.set_ylabel('ΔNDVI, unmatched buffer')
    ax.legend(loc='lower right', fontsize=6.0, markerscale=2.2, handletextpad=.1, borderaxespad=.1, labelspacing=.2, title='age', title_fontsize=6.0)
    ax.grid(color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    ax.set_xticks([-.1, 0, .1, .2]); ax.set_yticks([-.1, 0, .1, .2])
    return dict(r=float(r), mean_diff=float(dmean), n=int(len(u)))


def panel_e(ax):
    cs = pd.read_csv(CS12 / 'effects.csv', dtype={'tag': str, 'group': str}); cs = cs[cs['index'] == 'ndvi']
    bx = pd.read_csv(BP12 / 'effects.csv', dtype={'tag': str, 'group': str}); bx = bx[bx['index'] == 'ndvi']
    designs = [('matched controls', cs, 'estimate', 'joint_plus_calibration_lo', 'joint_plus_calibration_hi', 'o', 'fill'),
               ('unmatched buffer ground', bx, 'buffer_estimate', 'buffer_joint_plus_calibration_lo', 'buffer_joint_plus_calibration_hi', 'D', 'open')]
    xt, xl, rows = [], [], []
    for j, t in enumerate(TAGS):
        for i, g in enumerate(GROUPS):
            x0 = j * 3.6 + i
            for k, (lab, df, ce, cl, ch, mk, fill) in enumerate(designs):
                r = df[(df.tag == t) & (df.group == g)].iloc[0]
                x = x0 + (k - .5) * .26
                alpha = .45 if fill == 'faint' else 1
                ax.errorbar(x, r[ce], yerr=[[r[ce] - r[cl]], [r[ch] - r[ce]]], fmt=mk, ms=3.2, mfc=DATE[t] if fill in ('fill', 'faint') else 'white', mec=DATE[t], mew=.8,
                            ecolor=DATE[t], elinewidth=.85, capsize=.9, zorder=4, alpha=alpha)
                rows.append(dict(tag=t, group=g, design=lab, estimate=r[ce], lo=r[cl], hi=r[ch], n=int(r.n_cells)))
            xt.append(x0); xl.append(g)
        ax.text(j * 3.6 + 1, .122, DATE_LABEL[t], ha='center', va='top', fontsize=6.5, color=DATE[t], path_effects=bold(DATE[t], .15))
        if j < 2:
            ax.axvline(j * 3.6 + 2.8, color='#E1E5EA', lw=.6, zorder=0)
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    ax.set_ylim(-.024, .14); ax.set_xlim(-.6, 2 * 3.6 + 2.6)
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=6.3); ax.set_xlabel('Age group')
    ax.set_ylabel('ΔNDVI')
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    handles = [Line2D([], [], marker=mk, ls='', ms=3.6, mfc=INK if fill in ('fill', 'faint') else 'white', mec=INK, mew=.8, alpha=.45 if fill == 'faint' else 1, label=lab)
               for k, (lab, _, _, _, _, mk, fill) in enumerate(designs)]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, 1.0), ncol=2, fontsize=6.0, handletextpad=.25, borderaxespad=.15, labelspacing=.2, columnspacing=1.6)
    pd.DataFrame(rows).to_csv(FIGDATA / 'fig4e_designs.csv', index=False)


def main():
    cells, pairs, units = load()
    ex = choose_examples(cells, pairs, units)
    fig = figure(158)
    axa = fig.add_axes([.055, .40, .63, .57])
    axb = fig.add_axes([.725, .705, .255, .245])
    axc = fig.add_axes([.725, .405, .255, .245])
    axd = fig.add_axes([.075, .07, .25, .24])
    axe = fig.add_axes([.43, .08, .555, .216])  # top low enough for the buffer note
    sa = panel_a(axa, cells, pairs, units)
    sb = lens(axb, cells, pairs, units, ex[0], '')
    sc = lens(axc, cells, pairs, units, ex[1], '')
    sd = panel_d(axd, units)
    panel_e(axe)
    letter(axa, 'a', x=-.075, y=1.005, ha='left')
    axa.text(-.03, 1.005, 'Every treated unit linked to its usable matched controls', transform=axa.transAxes, fontsize=7.6, va='bottom', path_effects=bold(INK, .2))
    for ax, l, txt in [(axb, 'b', 'Matched controls and buffer ground'), (axc, 'c', 'The same for an older array')]:
        letter(ax, l, x=-.10, y=1.02, ha='left')
        ax.text(.0, 1.02, txt, transform=ax.transAxes, fontsize=7.1, va='bottom', path_effects=bold(INK, .2))
    letter(axd, 'd', x=-.3, y=1.03, ha='left')
    axd.text(-.17, 1.03, 'Matched vs unmatched, unit by unit', transform=axd.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    letter(axe, 'e', x=-.115, y=1.06, ha='left')
    axe.text(-.06, 1.06, 'Matched and unmatched comparisons agree across dates and ages', transform=axe.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    sq = fig.text(.735, .372, '■', color=BUFFER_GLOW, fontsize=8.5, ha='left', va='center')
    nt = fig.text(.752, .372, 'buffer ground: 300–1,000 m from arrays,\nsame scene, within 5 km, no matching', fontsize=6.8, color=SUB, ha='left', va='center', linespacing=1.05)
    # directly under the caption of (c), 0.5 mm apart
    v2fix(nt, 6.8); fig.canvas.draw()
    y = fig.transFigure.inverted().transform((0, CAPTIONS[1].get_window_extent(fig.canvas.get_renderer()).y0))[1] - .5 / 158
    nt.set_va('top'); nt.set_y(y)
    sq.set_y(y - nt.get_window_extent(fig.canvas.get_renderer()).height / 2 / fig.bbox.height)
    save(fig, 'Fig5')
    dump(FIGDATA / 'fig4_summary.json', dict(panel_a=dict(treated=sa['treated'], controls=sa['controls'], links=sa['links']), examples=[sb, sc], agreement=sd))
    print(sa, sb, sc, sd)


if __name__ == '__main__':
    main()
