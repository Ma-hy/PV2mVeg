"""Fig. S8: where the analysis units are, and why the rest are not.   [restricted inputs]

(a) Every 150 m cell with mapped panels (maximum panel fraction > 0.2), one category each, on the true-colour image;
(b) shares of the categories north and south of 36.03°N.

Inputs (project layout only): the east-west unit table (ts/ew/cells_pool.parquet: coordinates, orientation shares,
exclusion flags), the matched sets (v12/inputs), v12/cross_scene/cells_08.csv, the station polygons and the
true-colour display cache.
Output: FigS8.png/.pdf; figdata/figS27_coverage.csv (cell categories with coordinates, restricted; the region and
category columns are shipped in data/derived/summaries for the numbers)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from pv2mveg import data as D

CATS = [('analysis unit with an estimate, August 2025', '#1F6FD1'),
        ('east–west analysis unit without an estimate', '#9EC5F0'),
        ('north–south trackers (outside scope)', '#E07B39'),
        ('removed by the anomaly register', '#4A4A4A'),
        ('other east–west cells (array edges, mixed rows, no single construction year)', '#B9B2A6')]


def map_frame(ax, ext, labels='lb'):
    lonlat_frame(ax, ext, step=.1, sides=labels)


def base(ax, ext):
    rgb, _ = rgb_cache(ext, 1600)
    b = light_basemap(rgb, lift=.42, desat=.7); b[(rgb == 232).all(axis=2)] = 1.0
    ax.imshow(b, extent=ext, origin='upper', zorder=0, interpolation='bilinear')
    st = stations()
    ax.add_collection(PolyCollection(polygon_paths(st), facecolor='none', edgecolor='#56636F', linewidths=.25, zorder=3, alpha=.7))


def coverage_cells():
    """Every 150 m cell with mapped panels (maximum panel fraction > 0.2), one category each."""
    P.need_project('Fig. S8')
    c = D.cells_pool(ew=True)
    par = set(D.merged()[0])
    prim = set(pd.read_csv(CS12 / 'cells_08.csv').query('valid_cross').cell)
    a = c[c.pfmax_mean > .2].copy()
    cat = np.full(len(a), CATS[4][0], object)
    ns = a.ns_share.values >= .5
    reg = (a.excluded_station | a.excluded_cell).values
    ispar = a.cell.isin(par).values; ispr = a.cell.isin(prim).values
    cat[reg & ~ispar] = CATS[3][0]
    cat[ns & ~ispar] = CATS[2][0]
    cat[ispar] = CATS[1][0]
    cat[ispr] = CATS[0][0]
    a['cat'] = cat
    a['lon'], a['lat'] = GEO.transform(a.x.to_numpy(), a.y.to_numpy())
    a['region'] = np.where(a.lat < 36.03, 'south of 36.03°N', 'north of 36.03°N')
    a[['cell', 'x', 'y', 'lon', 'lat', 'region', 'cat']].to_csv(FIGDATA / 'figS27_coverage.csv', index=False)
    return a


def main():
    a = coverage_cells()
    ext = (a.lon.min() - .012, a.lon.max() + .012, a.lat.min() - .01, a.lat.max() + .01)
    fig = figure(92)
    # the map starts far enough right for its latitude labels
    axa = fig.add_axes([.07, .1, .54, .82]); axb = fig.add_axes([.71, .2, .27, .62])
    base(axa, ext)
    for name, col in CATS:
        z = a[a.cat == name]
        axa.add_collection(PolyCollection(cell_polygons(z.x, z.y), facecolor=col, edgecolor='none', zorder=2, rasterized=True))
    axa.axhline(36.03, color=INK, lw=.6, ls=(0, (3, 2)), zorder=4)
    axa.text(ext[0] + .01, 36.032, '36.03°N', fontsize=5.8, color=INK, va='bottom')
    map_frame(axa, ext)
    # the key at 6.8 pt with its longest entry on two lines, so that it stays clear of the northern arrays
    wrap = lambda s: s.replace('mixed rows, no single', 'mixed rows,\nno single')
    hh = [Patch(fc=col, label=wrap(f'{house(name)} ({int((a.cat == name).sum()):,})')) for name, col in CATS]
    v2fix(axa.legend(handles=hh, loc='upper left', fontsize=6.8, frameon=True, framealpha=.88, edgecolor='none', handlelength=1.0, labelspacing=.3,
                     title=f'{len(a):,} cells of 150 m with mapped panels', title_fontsize=6.8), 6.8)
    regs = ['north of 36.03°N', 'south of 36.03°N']
    for j, rg in enumerate(regs):
        z = a[a.region == rg]; left = 0
        for name, col in CATS:
            n = int((z.cat == name).sum()); share = n / len(z)
            axb.barh(j, share * 100, left=left, color=col, height=.55, edgecolor='white', lw=.4)
            if share > .07:
                axb.text(left + share * 50, j, f'{share:.0%}', ha='center', va='center', fontsize=5.8, color='white' if col in ('#1F6FD1', '#4A4A4A', '#E07B39') else INK)
            left += share * 100
        axb.text(101, j, f'{len(z):,}', ha='left', va='center', fontsize=5.8, color=SUB)
    axb.set_yticks(range(2)); axb.set_yticklabels(['north', 'south']); axb.invert_yaxis()
    axb.set_xlim(0, 112); axb.set_xticks([0, 25, 50, 75, 100]); axb.set_xlabel('Share of cells with mapped panels (%)')
    axb.spines['left'].set_visible(False); axb.tick_params(axis='y', length=0)
    # letters and headings: that of (b) far enough left that the heading ends inside the page
    letter(axa, 'a', x=-.12, y=1.01, ha='left'); axa.text(-.068, 1.01, 'Where the analysis units are, and why the rest are not', transform=axa.transAxes, fontsize=7.3, va='bottom', path_effects=bold(INK, .18))
    letter(axb, 'b', x=-.24, y=1.04, ha='left'); axb.text(-.14, 1.04, 'What remains in the south is mostly trackers', transform=axb.transAxes, fontsize=7.3, va='bottom', path_effects=bold(INK, .18))
    save(fig, 'FigS8')
    print({rg: a[a.region == rg].cat.value_counts().to_dict() for rg in regs})


if __name__ == '__main__':
    main()
