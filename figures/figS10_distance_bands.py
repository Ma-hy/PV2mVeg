"""Fig. S10: ΔNDVI in centre-distance bands from the row band, 0-2 and >=6 years, three acquisitions.

Input (repository data): v12/cross_scene/distance.csv
Output: FigS10.png/.pdf; plotted values figdata/supp_diag/figS14_distance.csv
"""
from _diag import *


def main():
    d = rd(CSD / 'distance.csv')
    fig = figure(62)
    # axes high enough for the two rows of unit counts below them
    axs = [fig.add_axes([.072 + i * .318, .30, .262, .53]) for i in range(3)]
    # the axis includes the 4-6 m band, so that the 0-2-year estimates of that band (>=30 units, open) are shown
    xmax = 4.5
    bands = ['1–1.5', '1.5–2', '2–3', '3–4'] + ['4–6']
    rows = []
    for i, (ax, tag) in enumerate(zip(axs, TAGS)):
        for j, g in enumerate(['0–2', '≥6']):
            q = d[(d.tag == tag) & (d.group == g)].sort_values('d_lo').reset_index(drop=True)
            q = q[q.n >= 30]
            x = q.index.to_numpy() + (j - .5) * .16
            thin = (q.n < .5 * q.n.max()).to_numpy()
            ax.plot(x[~thin], q.estimate[~thin], color=AGE[g], lw=.7, alpha=.7, zorder=2)
            for xx, r, t in zip(x, q.itertuples(), thin):
                pr(ax, xx, r.estimate, r.lo, r.hi, AGE[g], filled=not t, ms=3.4)
                if xx > xmax:  # no count for a band outside the axis
                    continue
                ax.text(xx, -.045 - .08 * j, f'{int(r.n):,}', transform=ax.get_xaxis_transform(), fontsize=4.9, color=AGE[g], ha='center', va='top')
            rows.append(q.assign(thin=thin))
        tidy(ax); ax.set_xlim(-.5, xmax); ax.set_xticks(range(len(bands))); ax.set_xticklabels(bands); ax.tick_params(axis='x', pad=17)
        ax.set_xlabel('Centre distance from the row band (m)')
        if i == 0:
            ax.set_ylabel('ΔNDVI')
        head(ax, f'{"abc"[i]}   {DATE_LABEL[tag]}')
    lim = [min(a.get_ylim()[0] for a in axs), max(a.get_ylim()[1] for a in axs)]
    for a in axs:
        a.set_ylim(lim)
    keep('figS14_distance', pd.concat(rows))
    legend_row(fig, [H(AGE[g], f'{g} years') for g in ['0–2', '≥6']] + [H(MUTED, 'band with fewer than half of the units of its age group', filled=False)], .985)
    fig.text(.072, .03, 'Numbers below the axis: contributing units. ' + 'Bands with fewer than 30 units are not shown (4–6 m at ≥6 years).',
             fontsize=6.0, color=SUB)
    save(fig, 'FigS10')


if __name__ == '__main__':
    main()
