"""Figure font: the TrueType file named in the configuration (PV2MVEG_FONT), else STIX.

The published figures used times+simsun.ttf, a combined Times New Roman / SimSun font of the first author's
computer. It is not redistributed (licence). Pixel-identical reproduction needs that file; with the fallback
(STIXGeneral, a Times-like face bundled with matplotlib) every figure is drawn with the same content and layout
rules, but glyph widths differ, so texts and the automatic size checks can shift slightly.
"""
from matplotlib import font_manager as fm

from pv2mveg.config import SETTINGS

FALLBACK = 'STIXGeneral'

FONT_FILE = SETTINGS.font if SETTINGS.font is not None and SETTINGS.font.exists() else None
if FONT_FILE is not None:
    fm.fontManager.addfont(str(FONT_FILE))
    FONT_NAME = fm.FontProperties(fname=str(FONT_FILE)).get_name()
else:
    FONT_NAME = FALLBACK


def properties(size):
    """FontProperties of the figure font at `size` (used by the scripts that set font properties per text)."""
    if FONT_FILE is not None:
        return fm.FontProperties(fname=str(FONT_FILE), size=size)
    return fm.FontProperties(family=FONT_NAME, size=size)
