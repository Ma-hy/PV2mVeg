"""Fig. 8: how much of the observed inter-row contrast can imaging alone reproduce?

Every optical scenario assigns the matched-control spectrum to array ground and simulates the NDVI contrast that
row-band geometry, panel end-members, the dark strip and a candidate sensor response would create. Dividing each
simulated contrast by the observed contrast on exactly the same units gives the share of the observation that
imaging alone could reproduce.
(a) every scenario as one dot, three acquisitions, >=6-year arrays, >=1.5 m;
(b) schematic of the forward model; (c) the retained-family share at 1, 1.5 and 2 m.

Input (repository data): v12/cross_scene/optics.csv
Output: Fig8.png/.pdf; plotted values figdata/fig7_{scenario_shares,summary,clearance}.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.plotting.style import *
from scipy.ndimage import gaussian_filter1d

FAM = {'retained': '#1F6FD1', 'sharpening': '#F2A007', 'rejected': '#9AA3AD'}
FAM_LABEL = {'retained': 'Retained (53): not rejected or untested',
             'sharpening': 'Sharpening (7): partly rejected',
             'rejected': 'Rejected in held-out test (48)'}
LIMIT = '#D1495B'


def scenarios():
    """Primary-design artefacts (cross_scene/optics.csv) with the frozen scenario families."""
    o = pd.read_csv(CS12 / 'optics.csv', dtype={'tag': str})
    a = o[(o.design == 'primary') & o.sets.str.contains('all_common', na=False)].copy()
    a['family'] = np.where(a.sets.str.contains('restricted_common'), 'retained', np.where(a.sets.str.contains('with_sharpening_common'), 'sharpening', 'rejected'))
    a = a[a.age == '≥6'].copy()
    a['share'] = a.art / a.obs * 100
    cnt = a[a.clearance == '1p5m'].groupby(['tag', 'family']).size().unstack()
    assert (cnt['retained'] == 53).all() and (cnt['sharpening'] == 7).all() and (cnt['rejected'] == 48).all(), cnt
    a.to_csv(FIGDATA / 'fig7_scenario_shares.csv', index=False)
    return a


def swarm(x, dx, dy, max_layers=12):
    """Deterministic beeswarm offsets for values x; dx, dy are the dot diameter in data units."""
    order = np.argsort(x)
    placed = []
    off = np.zeros(len(x))
    cands = [0] + [s * k for k in range(1, max_layers) for s in (1, -1)]
    for i in order:
        for c in cands:
            yo = c * dy
            ok = True
            for (xp, yp) in placed:
                if abs(xp - x[i]) < dx and abs(yp - yo) < dy * .98:
                    ok = False; break
            if ok:
                off[i] = yo; placed.append((x[i], yo)); break
    return off


def panel_a(ax, a):
    q = a[a.clearance == '1p5m']
    rows = {'06': 2.5, '08': 1.25, '2026_08': 0.0}
    xmin, xmax = -50, 185
    ax.set_xlim(xmin, xmax); ax.set_ylim(-.78, 3.12)
    ax.axvspan(0, 100, color='#FBF3EA', zorder=0, lw=0)
    ax.axvline(0, color='#8A96A3', lw=.6, zorder=1)
    ax.axvline(100, color=LIMIT, lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.text(101.5, 3.1, 'imaging alone reproduces\nthe observed ΔNDVI', color=LIMIT, fontsize=6.2, va='top', ha='left', linespacing=1.0)
    ax.text(1.5, 3.1, 'no imaging artefact', color=MUTED, fontsize=6.1, va='top', ha='left')
    # dot geometry in data units (axes is ~97 mm wide over 245 % and ~70 mm tall over 3.24 rows)
    bbox = ax.get_position(); fw, fh = ax.figure.get_size_inches()
    w_mm = bbox.width * fw * 25.4; h_mm = bbox.height * fh * 25.4
    dot_mm = 1.02
    dx = dot_mm / w_mm * (xmax - xmin); dy = dot_mm / h_mm * 3.90
    summary = []
    for tag, y0 in rows.items():
        z = q[q.tag == tag].reset_index(drop=True)
        off = swarm(z.share.values, dx, dy)
        for fam in ['rejected', 'sharpening', 'retained']:
            m = z.family == fam
            ax.scatter(z.share[m], y0 + off[m], s=4.9, c=FAM[fam], edgecolors='white', linewidths=.25, zorder=4 if fam == 'retained' else 3)
        r = z[z.family == 'retained']
        lo, hi = r.share.min(), r.share.max()
        # retained envelope and the unexplained remainder
        yb = y0 - .52
        ax.plot([lo, hi], [yb, yb], color=FAM['retained'], lw=1.6, solid_capstyle='butt', zorder=5)
        ax.annotate('', xy=(100, yb), xytext=(hi, yb),
                    arrowprops=dict(arrowstyle='<->,head_width=.12,head_length=.25', color=GROUND, lw=.75), zorder=5)
        ax.text((hi + 100) / 2, yb + .035, f'≥{np.floor(100 - hi):.0f}% not reproduced', color=GROUND, fontsize=6.0, ha='center', va='bottom')
        obs = z.obs.iloc[0]
        tr = ax.get_yaxis_transform()
        ax.text(-.012, y0 + .06, DATE_LABEL[tag], transform=tr, fontsize=7.0, color=INK, ha='right', va='bottom', path_effects=bold(INK, .2))
        ax.text(-.012, y0 + .04, f'observed ΔNDVI\n{signed(obs, 4)}\n= 100%', transform=tr, fontsize=6.1, color=SUB, ha='right', va='top', linespacing=1.05)
        summary.append(dict(tag=tag, observed=obs, retained_min=lo, retained_max=hi, sharpening_max=z[z.family == 'sharpening'].share.max(),
                            rejected_max=z[z.family == 'rejected'].share.max(), n_above_100=int((z.share > 100).sum()),
                            above_100_families=sorted(z[z.share > 100].family.unique().tolist())))
    ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xticks([-50, 0, 25, 50, 75, 100, 150])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'.replace('-', '−')))
    # two lines, so that it does not run into the axis label of (c)
    ax.set_xlabel('Share of the observed ΔNDVI reproduced by imaging alone' + '\n' + '(simulated ÷ observed, same units)')
    ax.grid(axis='x', color=GRID, lw=.45, zorder=0)
    handles = [Line2D([], [], marker='o', ls='', ms=3.8, mfc=FAM[f], mec='white', mew=.3, label=FAM_LABEL[f]) for f in ['retained', 'sharpening', 'rejected']]
    handles.append(Line2D([], [], color=FAM['retained'], lw=1.6, label='Retained range'))
    ax.figure.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .004), ncol=4, fontsize=6.1, handletextpad=.3, columnspacing=1.1, handlelength=1.3)
    return pd.DataFrame(summary)


def panel_b(ax):
    """Schematic of the forward model: blur carries panel darkness into the sampled ground."""
    x = np.linspace(-2.5, 7.5, 2001); dxm = x[1] - x[0]
    band = (x < 0) | (x > 5)
    truth = np.where(band, .05, .18)
    def blur(sig, tail=0.0, tail_sig=3.0):
        g = gaussian_filter1d(truth, sig / dxm, mode='nearest')
        if tail:
            g = (1 - tail) * g + tail * gaussian_filter1d(truth, tail_sig / dxm, mode='nearest')
        return g
    core = blur(.9); longt = blur(.9, .35, 2.6)
    ax.axvspan(-2.5, 0, color='#2B3A4A', alpha=.88, lw=0, zorder=0)
    ax.axvspan(5, 7.5, color='#2B3A4A', alpha=.88, lw=0, zorder=0)
    ax.axvspan(1.5, 3.5, color='#E3F2EC', lw=0, zorder=0)
    ax.text(-1.25, .155, 'row\nband', color='white', fontsize=5.8, ha='center', va='top', linespacing=1.0)
    ax.text(6.25, .155, 'row\nband', color='white', fontsize=5.8, ha='center', va='top', linespacing=1.0)
    ax.text(2.5, .1, 'sampled\ncentres', color=GROUND, fontsize=5.7, ha='center', va='bottom', linespacing=1.0)  # above the response labels
    ax.plot(x, truth, color=INK, lw=.9, drawstyle='steps-mid', zorder=3)
    ax.plot(x, core, color=FAM['retained'], lw=1.2, zorder=4)
    ax.plot(x, longt, color=FAM['rejected'], lw=1.1, ls=(0, (3, 1.6)), zorder=4)
    ax.text(2.5, .187, 'control spectrum\non array ground', fontsize=5.7, color=INK, va='bottom', ha='center', linespacing=.95)
    ax.annotate('measured-\nwidth\nresponse', xy=(.62, core[np.searchsorted(x, .62)]), xytext=(.18, .062), fontsize=5.6, color=FAM['retained'],
                arrowprops=dict(arrowstyle='-', color=FAM['retained'], lw=.5), ha='left', va='center', linespacing=.95)
    ax.annotate('long-tail\nresponse', xy=(4.45, longt[np.searchsorted(x, 4.45)]), xytext=(4.82, .062), fontsize=5.6, color='#6B7785',
                arrowprops=dict(arrowstyle='-', color='#9AA3AD', lw=.5), ha='right', va='center', linespacing=.95)
    ax.set_xlim(-2.5, 7.5); ax.set_ylim(0, .232)
    ax.set_xticks([0, 1.5, 3.5, 5]); ax.set_xticklabels(['0', '1.5', '3.5', '5'])
    ax.set_yticks([.05, .18]); ax.set_yticklabels(['panel', 'ground'])
    ax.set_xlabel('Distance across the gap (m)')
    ax.tick_params(axis='y', length=0, labelsize=6.0)
    for s in ['left']:
        ax.spines[s].set_visible(False)


def panel_c(ax, a):
    rows = []
    xs = {'1m': 0, '1p5m': 1, '2m': 2}
    for p, x0 in xs.items():
        for j, tag in enumerate(TAGS):
            z = a[(a.clearance == p) & (a.tag == tag)]
            r = z[z.family == 'retained'].share
            xx = x0 + (j - 1) * .22
            ax.plot([xx, xx], [r.min(), r.max()], color=DATE[tag], lw=2.6, solid_capstyle='round', alpha=.9, zorder=3)
            ax.plot(xx, r.median(), marker=DATE_MARK[tag], ms=3.2, mfc='white', mec=DATE[tag], mew=.9, zorder=4)
            rows.append(dict(clearance=p, tag=tag, retained_min=r.min(), retained_median=r.median(), retained_max=r.max()))
    ax.axhline(100, color=LIMIT, lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.axhspan(0, 100, color='#FBF3EA', lw=0, zorder=0)
    ax.axhline(0, color='#8A96A3', lw=.6, zorder=1)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['≥1 m', '≥1.5 m', '≥2 m'])
    ax.set_xlim(-.5, 2.5); ax.set_ylim(-60, 135)
    ax.set_yticks([-50, 0, 50, 100])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'.replace('-', '−')))
    ax.set_xlabel('Centre clearance from the row band')
    ax.set_ylabel('Share reproduced')
    ax.grid(axis='y', color=GRID, lw=.45, zorder=0)
    handles = [Line2D([], [], color=DATE[t], lw=2.4, label=DATE_SHORT[t]) for t in TAGS]
    ax.legend(handles=handles, loc='upper center', ncol=3, fontsize=5.9, handlelength=1.0, columnspacing=.7, borderaxespad=.1, bbox_to_anchor=(.5, 1.0))
    return pd.DataFrame(rows)


def main():
    a = scenarios()
    fig = figure(103.5)
    # panel rectangles were laid out for a 100 mm figure: 3.5 mm more height below for the two-line axis label
    # of (a); panels keep their size and place in mm
    pos = lambda r: [r[0], (r[1] * 100 + 3.5) / 103.5, r[2], r[3] * 100 / 103.5]
    axa = fig.add_axes(pos([.125, .175, .525, .755]))
    axb = fig.add_axes(pos([.735, .615, .25, .29]))
    axc = fig.add_axes(pos([.735, .175, .25, .30]))
    s = panel_a(axa, a)
    panel_b(axb)
    c = panel_c(axc, a)
    letter(axa, 'a', x=-.215, y=1.0, ha='left')
    axa.text(-.17, 1.0, 'Every imaging scenario, arrays aged ≥6 years, ≥1.5 m clearance', transform=axa.transAxes, fontsize=7.6, va='bottom', path_effects=bold(INK, .2))
    # letters and headings of (b, c) far enough left that the subtitle of (b) ends inside the page
    lx, tx = -.29, -.2
    letter(axb, 'b', x=lx, y=1.13, ha='left')
    axb.text(tx, 1.13, 'What a scenario simulates', transform=axb.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    v2fix(axb.text(tx, 1.035, 'schematic red reflectance, no ecological change', transform=axb.transAxes, fontsize=6.0, color=SUB, va='bottom'), 6.8)
    letter(axc, 'c', x=lx, y=1.05, ha='left')
    axc.text(tx, 1.05, 'Low at every clearance', transform=axc.transAxes, fontsize=7.4, va='bottom', path_effects=bold(INK, .2))
    s.to_csv(FIGDATA / 'fig7_summary.csv', index=False); c.to_csv(FIGDATA / 'fig7_clearance.csv', index=False)
    save(fig, 'Fig8')
    print(s.round(2).to_string()); print(c.round(1).to_string())


if __name__ == '__main__':
    main()
