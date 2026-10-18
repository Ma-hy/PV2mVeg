"""Main Tables 1-2, supplementary Table S1 and the held-out check (pv2mveg.tables).

Inputs: v12/cross_scene (cells, effects, optics), v12/buffer_primary/cells_08.csv,
        v12/cross_scene/secondary/1p5m/statistics/holdout.csv
Outputs: <output>/tables/table1.{csv,md}, table2.{csv,md}, tableS1.{csv,md}; <output>/figdata/heldout_primary.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg import tables  # noqa: E402

if __name__ == '__main__':
    out = tables.write_all()
    print(out['table1'].to_string(index=False))
    print(out['table2'].round(4).to_string(index=False))
    print(out['tableS1'].to_string(index=False))
    print(out['heldout'].round(4).to_string(index=False))
    print('TABLES_DONE')
