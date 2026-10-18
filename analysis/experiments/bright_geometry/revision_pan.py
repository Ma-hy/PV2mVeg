"""PAN texture and across-gap profiles using each released clearance layer."""
# Mirrored from experiments/bright_geometry/revision_pan.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
from pathlib import Path
import ecology_revision as R
import revision_measure as M
import pandas as pd


def job(args):
    tag,clearance,row,scene,kind,cache=args
    root=R.configure(clearance)
    import e3_pan_texture as T
    import p1_extract as P
    cache=Path(cache)
    if cache.exists():return str(cache)
    if kind=='control' and clearance!='1p5m':
        original=R.folder('1p5m')/'panprof'/('cache_'+tag)/cache.name
        if original.exists():
            import shutil
            shutil.copy2(original,cache)
            return str(cache)
    # Selection is read from the released native support; the legacy PAN
    # routine also requires same-scene source coverage for all fine pixels.
    scen=R.A.scenes(tag)
    td=T.work(tag,[(int(row.cell),row.x,row.y,int(scene),kind)],scen)
    pd_=P.work(tag,[(int(row.cell),row.x,row.y,int(scene),kind,float(row.pitch_m))],scen)
    result=dict(texture=td[0],pixels=pd.concat(td[1]).to_dict('records') if td[1] else [],
                treated=pd_[0],control=pd_[1],profile=pd_[2],geometry=pd_[3])
    R.atomic_cache(cache,result)
    return str(cache)
