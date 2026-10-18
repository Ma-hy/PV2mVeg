"""Fig. 6: spatial pattern of the inter-row contrast and its station-level structure.   [restricted inputs]

(a) 2025-08 unit-level ΔNDVI (>=1.5 m) on a faded true-colour image, with the hand-drawn station polygons.
(b) Station means against the station's construction year.
(c, d) Station means in pairs of acquisitions (stations valid in both).
Station means use the hand-drawn station polygons that define the resampling clusters.

Inputs (project layout only): v12/cross_scene/cells_*.csv with unit coordinates and construction years, the station
polygons, the unit table (ts/cells_pool.parquet), the matched sets (v12/inputs) and the true-colour display cache.
Output: Fig6.png/.pdf; figdata/fig5_station_means.csv (restricted: station construction years) and
fig5_summary.json (aggregates; shipped in data/derived/summaries for the numbers)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from pv2mveg import data as D


def unit_tables():
    P.need_project('Fig. 6')
    out = {}
    for t in TAGS:
        d = pd.read_csv(CS12 / f'cells_{t}.csv').query('valid_cross').copy()
        d['same_ndvi_dadj'] = d.cross_ndvi_dadj
        d['lon'], d['lat'] = GEO.transform(d.x.to_numpy(), d.y.to_numpy())
        out[t] = d
    return out


def station_means(units):
    rows = []
    for t, d in units.items():
        z = d[d.station > 0]
        g = z.groupby('station')
        s = pd.DataFrame(dict(n=g.size(), T=g['T'].agg(lambda v: int(v.mode().iloc[0])), age=g.age.mean(), ndvi=g.same_ndvi_dadj.mean()))
        s['tag'] = t
        rows.append(s.reset_index())
    s = pd.concat(rows)
    s.to_csv(FIGDATA / 'fig5_station_means.csv', index=False)
    return s


def panel_a(ax, units):
    t = units['08']
    allu = pd.concat(units.values()).drop_duplicates('cell')
    pad = .01
    ext = (allu.lon.min() - pad, allu.lon.max() + pad, allu.lat.min() - pad * .8, allu.lat.max() + pad * .8)
    rgb, _ = rgb_cache(ext, 2200)
    base = light_basemap(rgb, lift=.34, desat=.62)
    nodata = (rgb == 232).all(axis=2)
    base[nodata] = 1.0
    ax.imshow(base, extent=ext, origin='upper', zorder=0, interpolation='bilinear')
    st = stations()
    import geopandas as gpd
    cp = D.cells_pool().set_index('cell')
    par = cp.loc[D.merged()[0], ['x', 'y']].copy(); par['lon'], par['lat'] = GEO.transform(par.x.to_numpy(), par.y.to_numpy())
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(par.lon, par.lat), crs=4326)
    inscope = np.zeros(len(st), bool); inscope[gpd.sjoin(st.reset_index(), pts, predicate='contains')['index'].unique()] = True
    paths_in = polygon_paths(st[inscope]); paths_out = polygon_paths(st[~inscope])
    ax.add_collection(PolyCollection(paths_out, facecolor='none', edgecolor='#7B8794', linewidths=.35, linestyles=(0, (2, 1.5)), zorder=3, alpha=.8))
    ax.add_collection(PolyCollection(paths_in, facecolor='none', edgecolor='#2F3B45', linewidths=.4, zorder=3, alpha=.9))
    paths = paths_in + paths_out
    norm = TwoSlopeNorm(vmin=-.12, vcenter=0, vmax=.15)
    cmap = LinearSegmentedColormap.from_list('vivid_earth_teal', ['#7F2704', '#D9480F', '#FDAE61', '#FFF7E6', '#6DD3B6', '#00977C', '#004D40'])
    coll = PolyCollection(cell_polygons(t.x, t.y), array=t.same_ndvi_dadj.to_numpy(), cmap=cmap, norm=norm, edgecolor='none', zorder=2, rasterized=True)
    ax.add_collection(coll)
    lonlat_frame(ax, ext, step=.1)
    scalebar(ax, 10, x=.05, y=.07, color='#1F2933', text_color='#1F2933', ext=ext)
    north_arrow(ax, x=.05, y=.80)
    cax = ax.inset_axes([.60, .085, .34, .03])
    cb = ax.figure.colorbar(coll, cax=cax, orientation='horizontal', extend='both')
    cb.set_ticks([-.1, 0, .1]); cb.ax.tick_params(labelsize=5.9, length=1.8, pad=1)
    cb.outline.set_linewidth(.4)
    cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, 1) if v else '0'))
    cax.set_title('ΔNDVI', fontsize=6.2, color=INK, pad=2.5)
    ax.add_patch(FancyBboxPatch((.575, .02), .395, .135, boxstyle='round,pad=0.006,rounding_size=0.01', transform=ax.transAxes, fc='white', ec='none', alpha=.8, zorder=10))
    cax.set_zorder(11)
    ax.text(.10, .975, f'{len(t):,} of the {len(par):,} east–west units ({len(t) / len(par):.0%})\nsolid: {int(inscope.sum())} hand-drawn stations with analysis units\ndashed: {int((~inscope).sum())} stations outside the east–west fixed scope', transform=ax.transAxes, ha='left', va='top', fontsize=6.0, color=INK, linespacing=1.15)
    return dict(units=len(t), polygons=len(paths), in_scope=int(inscope.sum()), out_of_scope=int((~inscope).sum()))


def panel_b(ax, s):
    rng = np.random.default_rng(5)
    for j, t in enumerate(TAGS):
        z = s[s.tag == t]
        x = z['T'] + (j - 1) * .24 + rng.uniform(-.07, .07, len(z))
        ax.scatter(x, z.ndvi, s=2.5 + .32 * z.n, c=DATE[t], marker=DATE_MARK[t], edgecolors='white', linewidths=.35, alpha=.85, zorder=3 + j, label=DATE_LABEL[t])
    ax.axhline(0, color='#8A96A3', lw=.6, ls=(0, (2.5, 2)), zorder=1)
    ax.set_xlim(2012.3, 2026.0); ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.set_ylim(-.11, .19)
    ax.set_xlabel('Station construction year')
    ax.set_ylabel('Station-mean ΔNDVI')
    ax.grid(axis='y', color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4)
    # two columns, so that the key stays in the empty space below zero
    lg = ax.legend(loc='lower left', fontsize=7.0, markerscale=.9, handletextpad=.2, borderaxespad=.15, labelspacing=.25, ncol=2,
                   bbox_to_anchor=(.0, .0), title='marker area ∝ units', title_fontsize=7.0, frameon=False, alignment='left', columnspacing=.8)
    lg.get_title().set_color(MUTED)
    for h in lg.legend_handles:
        h.set_sizes([16])


def panel_pair(ax, s, ta, tb):
    a = s[s.tag == ta].set_index('station'); b = s[s.tag == tb].set_index('station')
    c = a.index.intersection(b.index); a, b = a.loc[c], b.loc[c]
    grp = np.where(b.age >= 6, '≥6', np.where(b.age >= 3, '3–5', '0–2'))
    lim = (-.07, .165)
    ax.plot(lim, lim, color=INK, lw=.6, ls=(0, (3, 2)), zorder=1)
    for g in ['0–2', '3–5', '≥6']:
        m = grp == g
        ax.scatter(a.ndvi[m], b.ndvi[m], s=3 + .4 * np.minimum(a.n[m], b.n[m]), c=AGE[g], edgecolors='white', linewidths=.35, alpha=.88, zorder=3, label=f'{g} y')
    r = np.corrcoef(a.ndvi, b.ndvi)[0, 1]
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
    ax.set_xticks([0, .1]); ax.set_yticks([0, .1])
    ax.set_xticks([-.05, .05, .15], minor=True); ax.set_yticks([-.05, .05, .15], minor=True)
    ax.set_xlabel(DATE_LABEL[ta]); ax.set_ylabel(DATE_LABEL[tb])
    ax.grid(color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(length=2.2, pad=1.4, labelsize=6.2)
    return dict(stations=int(len(c)), r=float(r))


def main():
    units = unit_tables()
    s = station_means(units)
    fig = figure(88)
    axa = fig.add_axes([.06, .03, .56, .90])
    # (b) ends 3 mm short of the right edge (its last year label); (c, d) low enough that their two-line headings
    # and the colour key between them and (b) do not touch; (c, d) narrow enough that the y label of (d) stays out of (c)
    axb = fig.add_axes([.695, .628, .285, .317])
    axc = fig.add_axes([.695, .1, .118, .30])
    axd = fig.add_axes([.867, .1, .118, .30])
    sa = panel_a(axa, units)
    panel_b(axb, s)
    sc = panel_pair(axc, s, '06', '08')
    sd = panel_pair(axd, s, '08', '2026_08')
    fig.canvas.draw()
    pos = axa.get_position()
    fig.text(pos.x0 - .045, pos.y1 + .012, 'a', fontsize=9.5, path_effects=bold(INK, .55), va='bottom')
    fig.text(pos.x0 - .02, pos.y1 + .012, f'Inter-row ΔNDVI across the park, {DATE_LABEL["08"]}', fontsize=7.6, va='bottom', path_effects=bold(INK, .2))
    letter(axb, 'b', x=-.25, y=1.03, ha='left')
    axb.text(-.16, 1.03, 'Station means by construction year', transform=axb.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    # c and d: one title per panel; the shared colour key (station mean age) sits above both
    for ax, ch, x0, x1, head, st in [(axc, 'c', -.58, -.25, 'Between seasons', sc), (axd, 'd', -.36, -.03, 'Between years', sd)]:
        letter(ax, ch, x=x0, y=1.155, ha='left')
        ax.text(x1, 1.155, head, transform=ax.transAxes, fontsize=7.2, va='bottom', path_effects=bold(INK, .2))
        ax.text(x1, 1.025, f'{st["stations"]} stations · r = {st["r"]:.2f}', transform=ax.transAxes, fontsize=6.0, color=SUB, va='bottom')
    ky = .508
    fig.text(.695, ky, 'station mean age', color=MUTED, fontsize=6.0, va='center')
    for i, g in enumerate(['0–2', '3–5', '≥6']):  # the keys start after 'station mean age'
        fig.text(.812 + i * .058, ky, '●', color=AGE[g], fontsize=7, va='center')
        fig.text(.826 + i * .058, ky, f'{g} y', color=SUB, fontsize=6.0, va='center')
    save(fig, 'Fig6')
    dump(FIGDATA / 'fig5_summary.json', dict(panel_a=sa, c=sc, d=sd))
    print(sa, sc, sd)


if __name__ == '__main__':
    main()
