"""Three-clearance auxiliary spectra, using updated masks and frozen shifts."""
# Mirrored from experiments/bright_geometry/revision_aux_summary.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import json
import shutil
import ecology_revision as R
import revision_measure as M
import pandas as pd


def main():
    src=R.folder('1p5m')/'adjacency/spectral'
    dates=['20250528','20250628','20250721','20250824']
    veg=pd.concat([pd.read_csv(src/f'e1_veg_{d}.csv',dtype={'date':str}) for d in dates])
    veg.to_csv(src/'e1_veg.csv',index=False)
    for p in R.POLICIES:
        root=R.configure(p);dest=root/'adjacency/spectral';dest.mkdir(exist_ok=True)
        if p!='1p5m':
            for d in dates:
                for typ in ['prof','ctrl']:shutil.copy2(src/f'e1_{d}_{typ}.parquet',dest/f'e1_{d}_{typ}.parquet')
            shutil.copy2(src/'e1_veg.csv',dest/'e1_veg.csv')
        import e2_spec_analysis as E
        E.SP=str(dest);E.CLEARANCE_M=R.F.THRESHOLDS[R.POLICIES[p]];E.MIN_PIXELS=50
        E.ELIGIBLE_CELLS={d:set(M.tables('06' if d in dates[:2] else '08','1p5m').query('valid_same').cell) for d in dates}
        E.main()
        (dest/'complete.json').write_text(json.dumps(dict(complete=True,clearance_m=E.CLEARANCE_M,min_pixels=50,
            eligibility='Primary 2 m valid cells, with >=50 auxiliary footprints and local or matched available controls.',
            registration='Original per-block image-to-map registration held fixed; exact distance to revised raster-square boundary.',
            limitation='Independent acquisition, but map-guided registration, sparse footprints and map-date transfer are not independent boundary validation.'),indent=2))
        print('AUX_SUMMARY_COMPLETE',p,flush=True)
