"""Compare the reproduced outputs with the published values (reference/) and write <output>/verification.json.

Run after 30_tables.py, 31_numbers.py and 40_figures.py. Prints one line per item and exits with 1 if anything
that was produced differs (figures that were not rendered are listed, not counted as failures). Pixel identity of
the figures is expected with the figure font of the publication (PV2MVEG_FONT); with the fallback font only the
numbers, tables and plotted-value tables are expected to agree.
"""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg import verify  # noqa: E402
from pv2mveg.config import P  # noqa: E402

if __name__ == '__main__':
    names = sys.argv[1:] or None
    rep = verify.report(names)
    P.out.mkdir(parents=True, exist_ok=True)
    (P.out / 'verification.json').write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding='utf-8')
    n = rep['numbers']
    print(f'layout {rep["layout"]}; font {rep["font"]}')
    print(f'numbers.json: {n["status"]} ({n.get("keys", 0)} values, {n.get("n_differing", 0)} differing)')
    for k, v in rep['tables'].items():
        print(f'table {k}: {v}')
    for k, v in rep['figures'].items():
        print(f'{k:7s} {v}')
    for k, v in rep['figdata'].items():
        print(f'figdata {k}: {v}')
    s = rep['supplement']
    print(f'supplementary tables: {s["status"]} ({s["checked"]} cells or groups checked, {s["n_failed"]} failed)')
    for b in s['failed']:
        print('   ', b)
    bad = (n['status'] not in ('identical', 'missing')) or any(v not in ('identical', 'identical cells', 'missing') for v in rep['tables'].values()) \
        or any(v in ('pixels differ',) or v.startswith('size') for v in rep['figures'].values() if rep['font'].startswith('fallback') is False) \
        or any(v == 'differs' for v in rep['figdata'].values()) or s['n_failed']
    print('VERIFY', 'FAILED' if bad else 'PASSED', '->', P.out / 'verification.json')
    sys.exit(1 if bad else 0)
