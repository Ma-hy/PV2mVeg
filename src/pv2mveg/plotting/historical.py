"""Text rules for the redrawn supplementary Figs. S1-S3.

Figs. S1-S3 are drawn by the figure code of earlier releases on their original inputs (figures/figS01-S03). Before
saving, the text on the figure is brought to the house style of the manuscript: ISO months and dates (a day and
month without a year is a 2025 acquisition), ΔNDVI for the baseline-adjusted contrast, panel labels "(a)" at the
panel-label size of style_v04, and no text smaller than the equivalent of 7.4 pt in a 165 mm figure (raised to the
equivalent of 7.5 pt). Black text and the 6.5 pt printed floor follow in the save pipeline (profile.save).
"""
import re

import matplotlib
from matplotlib.ticker import FixedFormatter, FixedLocator

from pv2mveg.plotting import profile as V2P
from pv2mveg.plotting.house import house

_ABBR = {m: i + 1 for i, m in enumerate(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])}
_A = '(' + '|'.join(_ABBR) + ')'


def fig_label(s):
    """House style for figure text, including the abbreviated forms used on the old figures."""
    s = s.replace(', Δ$_{adj}$', '')
    s = re.sub(r'\b' + _A + r'\.?\n(20\d\d)\b', lambda m: f'{m[2]}-{_ABBR[m[1]]:02d}', s)  # month and year on two lines
    s = re.sub(r'\b' + _A + r'\.? − ' + _A + r'\.? (20\d\d)\b', lambda m: f'{m[3]}-{_ABBR[m[1]]:02d} − {m[3]}-{_ABBR[m[2]]:02d}', s)
    s = re.sub(r'\b' + _A + r'\.? (2[0-9])\b(?![\d-])', lambda m: f'20{m[2]}-{_ABBR[m[1]]:02d}', s)
    s = fig_label_v2(s)
    return house(s)


def fig_label_v2(s):
    """Day-month dates in ISO form (without a year a 2025 acquisition); the dated test-combination labels of Fig. S2
    ("A06 (15 Jun)") as code and ISO date ("A06 2025-06-15")."""
    s = re.sub(r'\b(\d{1,2}) ' + _A + r'[a-z]*\.?(?: (20\d\d))?\b(?!-)', lambda m: f'{m[3] or 2025}-{_ABBR[m[2]]:02d}-{int(m[1]):02d}', s)
    s = re.sub(r'^([AB]\d{2,4}[ab]?) \((20\d\d-\d\d-\d\d)\)$', r'\1 \2', s)  # Fig. S2 (a, b): code and ISO date
    return s


def restyle(fig, log=None):
    """Apply fig_label to every text of the figure, including fixed tick labels."""
    for t in fig.findobj(matplotlib.text.Text):
        a = t.get_text()
        if a:
            b = fig_label(a)
            if b != a:
                t.set_text(b)
                if log is not None:
                    log.append((a, b))
    # Tick labels are regenerated at draw time (fixed, categorical or function formatters), so read
    # them after a draw and freeze the ticks of any axis whose labels change.
    fig.canvas.draw()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            old = [t.get_text() for t in axis.get_majorticklabels()]
            new = [fig_label(x) for x in old]
            if new != old:
                axis.set_major_locator(FixedLocator(axis.get_majorticklocs()))
                axis.set_major_formatter(FixedFormatter(new))
                if log is not None:
                    log.extend((a, b) for a, b in zip(old, new) if a != b)


def restyle_v2(fig, log, label_size, tight=False, rotate_first_two=False):
    """v2.0 sizes: panel labels "(a)" at `label_size`; texts smaller than the equivalent of 7.4 pt in a 165 mm figure
    raised to the equivalent of 7.5 pt (the tick-label and legend size of the v1.2 figures). `tight`: the figure is
    saved with bbox_inches='tight' (its composed width is that of the tight box). `rotate_first_two` (Fig. S2): the
    dated test-combination labels of (a, b) are longer in ISO form; at 50 degrees the first one stays on the page and
    the last ones clear the key below, and the long y labels of (a, b) at 8.2 pt fit the axes height."""
    _, composed = V2P.canvas(fig, 'tight' if tight else None)
    low, to = 7.4 * composed / V2P.BASE_MM, round(7.5 * composed / V2P.BASE_MM, 2)
    raised = 0
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text().strip()
        if not s:
            continue
        if re.fullmatch(r'\([a-z]\)', s):
            if t.get_fontsize() != label_size:
                log.append((f'{s} at {t.get_fontsize():g} pt', f'{s} at {label_size:g} pt'))
                t.set_fontsize(label_size)
        elif t.get_fontsize() < low - 1e-6:
            raised += 1
            t.set_fontsize(to)
    for ax in fig.axes:
        for axis in ('x', 'y'):
            ticks = getattr(ax, axis + 'axis').get_major_ticks()
            if ticks and ticks[0].label1.get_fontsize() < low - 1e-6:
                ax.tick_params(axis=axis, which='both', labelsize=to)
    if raised:
        log.append((f'{raised} texts below {low:.2f} pt', f'{to:g} pt (7.5 pt at 165 mm; composed width {composed:.0f} mm)'))
    if rotate_first_two:
        for ax in fig.axes[:2]:
            ax.tick_params(axis='x', labelrotation=50)
            ax.yaxis.label.set_fontsize(8.2)
        log.append(('Fig. S2 (a, b): tick labels at 40 degrees, y labels at the axis-label size', 'tick labels at 50 degrees, y labels at 8.2 pt'))


def finish(fig, log, label_size, tight=False, rotate_first_two=False):
    """Restyle the text of a redrawn figure (restyle, then restyle_v2); returns the list of changes."""
    restyle(fig, log)
    restyle_v2(fig, log, label_size, tight=tight, rotate_first_two=rotate_first_two)
    return log
