"""Save pipeline of the published (v2.0) figures.

Every figure is saved by `save()`, which

  * normalises the type sizes of the v1.2 figure family by role (`normalize`): base sizes in SIZES are font
    sizes of figures composed at BASE_MM = 165 mm (tick labels and legends 7.5 pt, axis labels 8.0, titles
    8.2, subtitles 7.0); figures composed at another width scale these by composed width / 165 mm;
  * draws text in pure black (`blacken`: neutral dark-grey text becomes #000000, with its fake-bold stroke);
  * enforces the printed floor (`apply_floor`): the manuscript inserts every figure at 450 pt = 158.75 mm
    (PRINT_MM) and no text is printed below FLOOR = 6.5 pt;
  * keeps every text readable on its ground (`legible`): a text whose WCAG contrast with the ground under it
    (measured on the drawn figure) is below 3:1 is recoloured: light text on a coloured box becomes black or
    white, coloured text on a light ground keeps its hue and is darkened to 4.5:1; symbols drawn as text keep
    their colour;
  * writes <name>.png at 500 dpi and <name>.pdf into the output folder (config: P.figures).

These steps change pixels and are part of the published figures. The text-overlap audit that the authors ran
alongside (overlap and canvas checks written to a report) does not change pixels and is not included.
"""
import os
import sys
import json
import hashlib
from pathlib import Path

from pv2mveg.config import P, REPO

V2 = True                 # the published figures are drawn under this profile only

DPI = 500
MIN_CONTRAST = 3.0        # smallest WCAG contrast ratio of a text with its ground (3:1, large text and graphics)
PRINT_MM = 158.75         # printed width: the Word files insert every figure at 450 pt (A4 text width 16 cm)
BASE_MM = 165.0           # composed width of the v1.2 figures; SIZES are font sizes at this width
FLOOR = 6.5               # smallest printed font size (pt): 6.76 pt in a 165 mm figure
BLACK = '#000000'
SIZES = dict(letter=9.5, title=8.2, sub=7.0, label=8.0, tick=7.5, legend=7.5,
             cbar=7.0, inset=6.8, note_min=6.8, note_max=7.2, note_lift=0.9)

# figures of the published set whose v1.2-family sizes are normalised by role (the author's table of the 19
# figures re-rendered for v2.0; Fig. 2, drawn for v2.0, takes its sizes from the style without normalisation)
NORMALISED = {'Fig1', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'Fig7', 'Fig8', 'FigS1', 'FigS2', 'FigS3', 'FigS4', 'FigS5',
              'FigS6', 'FigS7', 'FigS8', 'FigS9', 'FigS10', 'FigS11', 'FigS12'}


def v2pick(default, v2):
    """`v2` (the published v2.0 value); `default` is the value of the earlier release, kept for reference."""
    return v2


def v2fix(obj, size):
    """Exempt a text, a legend or an axes (its tick labels) from the role normalisation: use `size`."""
    obj._v2_size = float(size)
    return obj


def script():
    s = os.environ.get('PV2MVEG_FIG_SCRIPT') or (sys.argv[0] if sys.argv and sys.argv[0] else 'unknown')
    p = Path(s).resolve()
    try:
        return p.relative_to(REPO).as_posix()
    except ValueError:
        return p.name


# ------------------------------------------------------------------ colour helpers
def _rgba(c):
    from matplotlib.colors import to_rgba
    try:
        return to_rgba(c)
    except ValueError:
        return None


def _luminance(rgb):
    """WCAG relative luminance of sRGB colours in 0..1 (last axis R, G, B)."""
    import numpy as np
    c = np.asarray(rgb, dtype=float)
    lin = np.where(c <= .04045, c / 12.92, ((c + .055) / 1.055) ** 2.4)
    return lin @ np.array([.2126, .7152, .0722])


def _ratio(a, b):
    la, lb = float(_luminance(a[:3])), float(_luminance(b[:3]))
    return (max(la, lb) + .05) / (min(la, lb) + .05)


def _hex(v):
    return '#%02X%02X%02X' % tuple(int(round(min(1., max(0., u)) * 255)) for u in v[:3])


def _symbol(s):
    """A text without letters or digits (a legend swatch such as the square or the dot drawn as text)."""
    return not any(ch.isalnum() for ch in s)


def _grounds(fig, rows):
    """(row, ground colour, contrast ratio) of every drawn text, read from the drawn canvas: the ground is the median
    colour of the pixels in the text's box that differ from the text colour (glyphs and their anti-aliasing excluded; a
    label box, legend frame or halo under the text counts as ground). A box that is nearly all text colour, i.e. white
    text on white, gives the text colour itself and ratio 1."""
    import numpy as np
    try:
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3] / 255.
    except Exception:
        return []
    H, W = buf.shape[:2]
    out = []
    for x in rows:
        c = _rgba(x['t'].get_color())
        bb = x['bb']
        c0, c1 = max(0, int(np.floor(bb.x0))), min(W, int(np.ceil(bb.x1)))
        r0, r1 = max(0, int(np.floor(H - bb.y1))), min(H, int(np.ceil(H - bb.y0)))
        if c is None or c1 - c0 < 2 or r1 - r0 < 2:
            continue
        px = buf[r0:r1, c0:c1].reshape(-1, 3)
        far = np.abs(px - np.array(c[:3])).max(axis=1) > .25
        ground = np.median(px[far], axis=0) if far.mean() > .2 else np.array(c[:3])
        out.append((x, ground, _ratio(c, ground)))
    return out


def _recolour(t, new):
    """Set a text colour; a fake-bold stroke in the old colour follows (halos keep theirs)."""
    from matplotlib import patheffects as pe
    c = _rgba(t.get_color())
    effects = []
    for p in (t.get_path_effects() or []):
        if isinstance(p, pe.Stroke):
            fg = _rgba(p._gc.get('foreground', 'k'))
            if fg is not None and c is not None and max(abs(a - b) for a, b in zip(c[:3], fg[:3])) < .02:
                kw = dict(p._gc); kw['foreground'] = new
                p = type(p)(offset=p._offset, **kw)
        effects.append(p)
    t.set_color(new)
    if effects:
        t.set_path_effects(effects)


def legible(fig):
    """Recolour the texts whose contrast with their ground is below MIN_CONTRAST (see the module notes)."""
    import colorsys
    import matplotlib
    texts, r = _drawn(fig)
    rows = []
    for t in texts:
        bb = matplotlib.text.Text.get_window_extent(t, r)
        if bb.width > 0 and bb.height > 0:
            rows.append(dict(t=t, bb=bb, text=t.get_text().replace(chr(10), ' / ')[:70]))
    changed = []
    for x, g, ratio in _grounds(fig, rows):
        t = x['t']
        c = _rgba(t.get_color())
        if ratio >= MIN_CONTRAST or _symbol(t.get_text()) or c is None:
            continue
        if _luminance(c[:3]) > _luminance(g):        # light text on a darker ground (white on a coloured box)
            new = (0., 0., 0.) if _ratio((0., 0., 0.), g) >= _ratio((1., 1., 1.), g) else (1., 1., 1.)
        else:                                         # coloured text on a light ground: same hue, darker
            h, l, s = colorsys.rgb_to_hls(*c[:3])
            new = next((v for v in (colorsys.hls_to_rgb(h, l * (1 - k / 100), s) for k in range(1, 101))
                        if _ratio(v, g) >= 4.5), (0., 0., 0.))
        _recolour(t, (*new, c[3]))
        changed.append((x['text'], _hex(c), _hex(new), _hex(g), round(ratio, 2)))
    return changed


def _neutral(c):
    """(is neutral, mean level) of a colour: chroma below 0.12 on a 0-1 scale."""
    r = _rgba(c)
    if r is None:
        return False, 1.0
    mx, mn = max(r[:3]), min(r[:3])
    return (mx - mn) < .12, sum(r[:3]) / 3


def _strokes(t):
    from matplotlib import patheffects as pe
    return [p for p in (t.get_path_effects() or []) if isinstance(p, pe.Stroke)]


def _is_bold(t):
    """Fake bold: a thin stroke in the text colour (halos are white or dark and wide)."""
    c = _rgba(t.get_color())
    for p in _strokes(t):
        fg = _rgba(p._gc.get('foreground', 'k'))
        lw = p._gc.get('linewidth', 1.0)
        if c is not None and fg is not None and lw <= .6 and max(abs(a - b) for a, b in zip(c[:3], fg[:3])) < .02:
            return True
    return False


def blacken(fig):
    """Neutral dark-grey text -> black, its own-colour strokes too; light greys are only reported."""
    import matplotlib
    from matplotlib import patheffects as pe
    changed, light = [], []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible():
            continue
        c = _rgba(t.get_color())
        if c is None or c[3] == 0:
            continue
        neutral, level = _neutral(c)
        if not neutral or level < 1e-6:
            continue
        if level < .6:
            effects = []
            for p in (t.get_path_effects() or []):
                if isinstance(p, pe.Stroke):
                    fg = _rgba(p._gc.get('foreground', 'k'))
                    if fg is not None and max(abs(a - b) for a, b in zip(c[:3], fg[:3])) < .02:
                        kw = dict(p._gc); kw['foreground'] = BLACK
                        p = type(p)(offset=p._offset, **kw)
                effects.append(p)
            t.set_color(BLACK)
            if effects:
                t.set_path_effects(effects)
            changed.append(s[:60])
        elif level < .9:
            light.append(s[:60])
    return changed, light


# ------------------------------------------------------------------ role normalisation (v1.2 family)
def _axes_kinds(fig):
    kinds = {}

    def walk(ax, kind):
        if ax in kinds:
            return
        kinds[ax] = kind
        for ch in getattr(ax, 'child_axes', []):
            walk(ch, 'cbar' if hasattr(ch, '_colorbar') else 'inset')
    for ax in fig.axes:
        walk(ax, 'cbar' if hasattr(ax, '_colorbar') else 'main')
    return kinds


def note_size(s):
    """Annotations and direct labels: lifted by 0.9 pt into 6.8-7.2 pt, never made smaller."""
    return max(s, min(max(s + SIZES['note_lift'], SIZES['note_min']), SIZES['note_max']))


def normalize(fig):
    """One set of sizes for the v1.2 figures (font sizes at the composed width of 165 mm).

    tick labels 7.5 (insets 6.8, colour bars 7.0); axis labels 8.0; axes titles 8.2; legends keep the
    size the script gives them under the profile (7.5 by default; one below 6.9 is raised to 7.5);
    free texts: panel letters (>= 9 pt) unchanged, fake-bold headings of 7.05-8.5 pt -> 8.2, plain
    5.8-6.6 pt lines just above a panel (subtitles) -> 7.0, all other texts -> note_size().
    `v2fix(artist, size)` overrides any of these.
    """
    import matplotlib
    from matplotlib.legend import Legend
    S = SIZES
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    log, seen = [], set()

    def setsize(t, size, role):
        seen.add(t)
        if not t.get_text().strip():
            return
        size = getattr(t, '_v2_size', size)
        old = t.get_fontsize()
        if abs(size - old) > 1e-9:
            t.set_fontsize(size)
        log.append((role, t.get_text()[:50], round(old, 2), round(size, 2)))

    kinds = _axes_kinds(fig)
    for ax, kind in kinds.items():
        for name in ('x', 'y'):
            axis = getattr(ax, name + 'axis')
            size = getattr(axis, '_v2_size', getattr(ax, '_v2_size', {'main': S['tick'], 'inset': S['inset'], 'cbar': S['cbar']}[kind]))
            ticks = axis.get_major_ticks() + axis.get_minor_ticks()
            labels = [tk.label1 for tk in ticks] + [tk.label2 for tk in ticks]
            shown = [t for t in labels if t.get_visible() and t.get_text().strip()]
            if shown:
                log.append(('tick-' + kind, f'{name}: ' + ' '.join(t.get_text() for t in shown[:6])[:50], round(shown[0].get_fontsize(), 2), size))
            ax.tick_params(axis=name, which='both', labelsize=size)
            seen.update(labels)
            seen.add(axis.get_offset_text())
            axis.get_offset_text().set_fontsize(size)
            lab = axis.label
            if lab.get_text().strip():
                if kind == 'main':
                    setsize(lab, S['label'] if lab.get_fontsize() >= 7.0 else note_size(lab.get_fontsize()), 'label')
                else:
                    setsize(lab, S['cbar'] if kind == 'cbar' else S['inset'], 'label-' + kind)
            seen.add(lab)
        for tt in (ax.title, ax._left_title, ax._right_title):
            if kind == 'main':
                setsize(tt, S['title'], 'title')
            else:
                setsize(tt, S['cbar'] if kind == 'cbar' else note_size(tt.get_fontsize()), 'title-' + kind)
    for lg in fig.findobj(Legend):
        forced = getattr(lg, '_v2_size', None)
        for t in list(lg.get_texts()) + [lg.get_title()]:
            seen.add(t)
            if not t.get_text().strip():
                continue
            new = forced if forced is not None else (S['legend'] if t.get_fontsize() < 6.9 else t.get_fontsize())
            log.append(('legend', t.get_text()[:50], round(t.get_fontsize(), 2), new))
            if abs(new - t.get_fontsize()) > 1e-9:
                t.set_fontsize(new)
    for t in fig.findobj(matplotlib.text.Text):
        if t in seen or not t.get_text().strip() or not t.get_visible():
            continue
        s = t.get_fontsize()
        if hasattr(t, '_v2_size'):
            setsize(t, t._v2_size, 'fixed')
        elif s >= 9.0:
            setsize(t, s, 'letter')
        elif _is_bold(t) and 7.05 <= s <= 8.5:
            setsize(t, S['title'], 'heading')
        elif not _is_bold(t) and 5.8 <= s <= 6.6 and t.axes is not None and _above(t, t.axes, r):
            setsize(t, S['sub'], 'subtitle')
        else:
            setsize(t, note_size(s), 'note')
    return log


def _above(t, ax, r):
    import matplotlib
    bb = matplotlib.text.Text.get_window_extent(t, r)
    ab = ax.bbox
    return bb.y0 >= ab.y1 - 1.0 and bb.x1 > ab.x0 and bb.x0 < ab.x1


# ------------------------------------------------------------------ floor
def canvas(fig, bbox_inches=None, pad_inches=None):
    """Saved canvas in display pixels (figure dpi) and its width in mm."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    if bbox_inches == 'tight':
        import matplotlib as mpl
        pad = mpl.rcParams['savefig.pad_inches'] if pad_inches is None else pad_inches
        bb = fig.get_tightbbox(fig.canvas.get_renderer()).padded(pad)
        c = Bbox.from_extents(bb.x0 * fig.dpi, bb.y0 * fig.dpi, bb.x1 * fig.dpi, bb.y1 * fig.dpi)
    else:
        c = fig.bbox
    return c, c.width / fig.dpi * 25.4


def apply_floor(fig, composed_mm):
    import matplotlib
    import matplotlib.axes
    low = FLOOR * composed_mm / PRINT_MM
    raised = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text().strip() and t.get_fontsize() < low - 1e-6:
            raised.append((t.get_text()[:50], round(t.get_fontsize(), 2)))
            t.set_fontsize(round(low + .005, 2))
    for ax in fig.findobj(matplotlib.axes.Axes):
        for name in ('x', 'y'):
            axis = getattr(ax, name + 'axis')
            ticks = axis.get_major_ticks() + axis.get_minor_ticks()
            if ticks and ticks[0].label1.get_fontsize() < low - 1e-6:
                ax.tick_params(axis=name, which='both', labelsize=round(low + .005, 2))
    return raised


def _drawn(fig):
    """Draw the figure and return the texts actually drawn (visible, non-empty)."""
    import matplotlib
    drawn = []
    orig = matplotlib.text.Text.draw

    def draw(self, renderer):
        orig(self, renderer)
        if self.get_visible() and self.get_text().strip():
            drawn.append(self)
    matplotlib.text.Text.draw = draw
    try:
        fig.canvas.draw()
    finally:
        matplotlib.text.Text.draw = orig
    return list(dict.fromkeys(drawn)), fig.canvas.get_renderer()


# ------------------------------------------------------------------ save
def save(fig, name, apply_font, family, bbox_inches=None, pad_inches=None, normalize_roles=None, origin=None):
    """Normalise (v1.2 family), blacken, floor, recolour for contrast, then PNG (500 dpi) + PDF into P.figures."""
    import matplotlib
    import matplotlib.pyplot as plt
    folder = P.figures
    folder.mkdir(parents=True, exist_ok=True)
    apply_font(fig)
    if normalize_roles is None:
        normalize_roles = family == 'v12' and name in NORMALISED
    if normalize_roles:
        normalize(fig)
    greys, light = blacken(fig)
    cbox, composed = canvas(fig, bbox_inches, pad_inches)
    raised = apply_floor(fig, composed)
    cbox, composed = canvas(fig, bbox_inches, pad_inches)
    raised += apply_floor(fig, composed)
    readable = legible(fig)
    meta = {'Title': name, 'Script': script(), 'DPI': str(DPI), 'ComposedWidthMM': f'{composed:.2f}',
            'Software': f'Matplotlib version{matplotlib.__version__}, https://matplotlib.org/'}
    kw = dict(dpi=DPI, facecolor='white')
    if bbox_inches is not None:
        kw['bbox_inches'] = bbox_inches
        if pad_inches is not None:
            kw['pad_inches'] = pad_inches
    apply_font(fig)
    fig.savefig(folder / f'{name}.png', metadata=meta, **kw)
    apply_font(fig)
    fig.savefig(folder / f'{name}.pdf', metadata={'Title': name, 'Subject': f'{script()}; {DPI} dpi'}, **kw)
    plt.close(fig)
    rec = dict(name=name, script=script(), composed_mm=round(composed, 2), recoloured_black=len(greys),
               raised_to_floor=raised, recoloured_for_contrast=readable)
    P.checks.mkdir(parents=True, exist_ok=True)
    (P.checks / f'{name}.json').write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str), encoding='utf-8')
    print(f'SAVED {name}  composed {composed:.1f} mm, {len(greys)} recoloured black, {len(raised)} raised to floor, '
          f'{len(readable)} recoloured for contrast', flush=True)
    return folder / f'{name}.png', 0


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
