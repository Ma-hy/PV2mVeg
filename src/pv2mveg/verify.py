"""Comparison of the reproduced outputs with the published values (reference/).

numbers()      numbers.json, byte for byte and key by key
tables()       Tables 1-2 (Markdown and CSV) byte for byte, Table S1 cell by cell against the published supplement
figures()      every rendered PNG against the SHA-256 of the RGBA pixels of the published figure
figdata()      the plotted-value tables written by the figure scripts against those of the release
supplement()   the supplementary-table cells recomputed from the result files (pv2mveg.supplement)
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from pv2mveg.config import P, REPO
from pv2mveg.catalog import FIGURES

REF = REPO / 'reference'


def text_bytes(path):
    """File bytes with CRLF line ends read as LF (outputs written on Windows compare equal to the LF reference)."""
    return Path(path).read_bytes().replace(b'\r\n', b'\n')


def _flat(d, prefix=''):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flat(v, f'{prefix}{k}/'))
    else:
        out[prefix.rstrip('/')] = d
    return out


def numbers():
    ref = REF / 'numbers.json'; new = P.numbers / 'numbers.json'
    if not new.exists():
        return dict(status='missing', file=str(new))
    a, b = _flat(json.loads(ref.read_text(encoding='utf-8'))), _flat(json.loads(new.read_text(encoding='utf-8')))
    diff, worst = [], 0.0
    for k in sorted(set(a) | set(b)):
        x, y = a.get(k, '<missing>'), b.get(k, '<missing>')
        if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
            if not (x == y or (isinstance(x, float) and isinstance(y, float) and math.isnan(x) and math.isnan(y))):
                diff.append(k); worst = max(worst, abs(x - y))
        elif x != y:
            diff.append(k)
    return dict(status='identical' if text_bytes(ref) == text_bytes(new) else ('equal values' if not diff else 'differs'),
                keys=len(a), differing_keys=diff[:50], n_differing=len(diff), max_abs_diff=worst)


def _rows(md):
    import re
    return [[c.strip() for c in ln.strip().strip('|').split('|')] for ln in md.split('\n') if ln.startswith('|') and not re.match(r'^\|\s*:?-', ln)]


def supplement_blocks():
    """{'Table S1': (title, rows), ...} of the published supplement tables (rows include the header row)."""
    out, cur = {}, None
    for block in (REF / 'supplement_tables.md').read_text(encoding='utf-8').split('\n\n'):
        lines = block.strip().split('\n')
        if lines[0].startswith('**Table S'):
            cur = lines[0].split('**')[1].rstrip('.')
            out[cur] = (lines[0], _rows('\n'.join(lines[1:])))
        elif cur and lines[0].startswith('|'):
            out[cur] = (out[cur][0], out[cur][1] + _rows(block))
    return out


def tables():
    res = {}
    for n in ['table1.md', 'table1.csv', 'table2.md', 'table2.csv']:
        new = P.tables / n
        res[n] = 'missing' if not new.exists() else ('identical' if text_bytes(new) == text_bytes(REF / 'tables' / n) else 'differs')
    s1 = P.tables / 'tableS1.md'
    if s1.exists():
        res['tableS1'] = 'identical cells' if _rows(s1.read_text(encoding='utf-8')) == supplement_blocks()['Table S1'][1] else 'differs'
    else:
        res['tableS1'] = 'missing'
    return res


def pixel_sha(path):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    a = np.ascontiguousarray(np.asarray(Image.open(path).convert('RGBA')))
    return a.shape[1], a.shape[0], hashlib.sha256(a.tobytes()).hexdigest()


def figures(names=None):
    ref = json.loads((REF / 'figures.json').read_text(encoding='utf-8'))['figures']
    res = {}
    for name in names or FIGURES:
        p = P.figures / f'{name}.png'
        if not p.exists():
            res[name] = 'not rendered'
            continue
        w, h, s = pixel_sha(p)
        r = ref[name]
        res[name] = ('pixel-identical' if s == r['rgba_sha256'] else
                     f'size {w}x{h} (published {r["width"]}x{r["height"]})' if (w, h) != (r['width'], r['height']) else 'pixels differ')
    return res


def figdata(names=None):
    res = {}
    for name in names or FIGURES:
        for f in FIGURES[name][3]:
            new, ref = P.figdata / f, REF / 'figdata' / f
            if new.exists() and ref.exists():
                res[f] = 'identical' if text_bytes(new) == text_bytes(ref) else 'differs'
    return res


def supplement():
    from pv2mveg import supplement as S
    return S.check(supplement_blocks())


def report(names=None):
    from pv2mveg.config import SETTINGS
    from pv2mveg.plotting.fonts import FONT_FILE
    return dict(layout=P.layout, font=str(FONT_FILE) if FONT_FILE else 'fallback (STIX): pixel identity is not expected',
                output=str(P.out), numbers=numbers(), tables=tables(), figures=figures(names), figdata=figdata(names), supplement=supplement())
