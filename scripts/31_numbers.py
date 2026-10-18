"""Every number quoted in the text: <output>/numbers/numbers.json (pv2mveg.numbers).

Run after 30_tables.py. In the project layout, the Fig. 6 and Fig. S8 scripts must have run first (their side
outputs fig5_summary.json and figS27_coverage.csv are inputs); the repository layout ships these two summaries in
data/derived/summaries.
"""
import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg import numbers  # noqa: E402
from pv2mveg.config import P, TAGS  # noqa: E402

if __name__ == '__main__':
    N = numbers.write()
    print(json.dumps({t: {k: N[t][k] for k in ['units', 'km2', 'ndvi_≥6']} for t in TAGS}, indent=1, default=float))
    print('NUMBERS_DONE', P.numbers / 'numbers.json')
