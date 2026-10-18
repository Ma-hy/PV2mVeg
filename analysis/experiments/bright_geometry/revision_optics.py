"""Recompute optical forward models on revised masks and all three supports.

Candidate kernel parameters remain frozen measurement assumptions. Their new
forward artefacts are recomputed; they are not interpreted as validated bounds.
"""
# Mirrored from experiments/bright_geometry/revision_optics.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
import json
import shutil
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import ecology_revision as R
import revision_measure as M
import pandas as pd


def seed(clearance):
    root=R.configure(clearance);dst=root/'adjacency';ew=dst/'ew_only';ew.mkdir(exist_ok=True)
    base=R.ROOT/'data/eco/adjacency'
    for name in ['k1_windows.csv','k1_kernels.npz','k2_kernels.csv','k2_np_kernels.npz']:
        if not (dst/name).exists():shutil.copy2(base/name,dst/name)
    # Provenance is retained: these initial values are replaced by newly
    # measured PAN summaries before the forward-model stage.
    for name in ['k3_pan_esf.csv','k3_holdout_summary.csv','a4_shadow_params.csv']:
        if not (ew/name).exists():shutil.copy2(base/'ew_only'/name,ew/name)
    return root


@lru_cache(maxsize=9)
def kernels(tag,clearance,kind):
    root=R.configure(clearance)
    import a2_cells as B
    import a4_ew_forward as C
    C.EW=str(root/'adjacency/ew_only')
    return B.load_kernels(tag) if kind=='a2' else C.kernels_for(tag)


def worker(args):
    tag,clearance,kind,row,path=args;path=Path(path)
    if path.exists():return str(path)
    root=R.configure(clearance)
    import a2_cells as B
    B.centre_dist=lambda pan,ro,co,h,w:R.F.exact_center_distance(pan,h,w,ro)
    import a4_ew_forward as C
    C.EW=str(root/'adjacency/ew_only')
    frame=pd.DataFrame([row]);sc=R.scenes(tag);kk=kernels(tag,clearance,kind)
    if kind=='a2':
        z=B.process(tag,frame,kk,sc)
        data=dict(zip(['cells','prof','pan','evar'],z))
    else:
        sh=pd.read_csv(root/'adjacency/ew_only/a4_shadow_params.csv',dtype={'tag':str})
        data={'cells':C.process(tag,frame,kk,sc,float(sh.loc[sh.tag==tag,'W'].iloc[0]))}
    R.atomic_cache(path,data)
    return str(path)


def run(tag,clearance,kind,workers=4,limit=None):
    root=seed(clearance);t=M.tables(tag,clearance);t=t[t.valid_same].copy()
    if limit:t=t.head(limit)
    out=root/'adjacency';dest=out if kind=='a2' else out/'ew_only'
    cache=dest/f'cache_{kind}_{tag}';cache.mkdir(exist_ok=True)
    jobs=[(tag,clearance,kind,row,cache/f'{int(row["cell"])}.json') for row in t.to_dict('records')]
    with ProcessPoolExecutor(workers) as pool:
        for i,_ in enumerate(pool.map(worker,jobs,chunksize=1),1):
            if i%100==0:print('OPTICS',kind,tag,clearance,i,len(jobs),flush=True)
    result={k:[] for k in (['cells','prof','pan','evar'] if kind=='a2' else ['cells'])}
    for *_,p in jobs:
        z=json.loads(p.read_text(encoding='utf-8'))
        for k in result:result[k].extend(z[k])
    for k,v in result.items():pd.DataFrame(v).to_parquet(dest/f'{kind}_{tag}_{k}.parquet',index=False)
    counts=pd.DataFrame(result['cells']).merge(t[['cell','ir_n']],on='cell',suffixes=('','_main'))
    ref='ir_n_main' if 'ir_n_main' in counts else 'ir_n'
    audit=dict(complete=limit is None,n=len(t),max_native_count_difference=float((counts.n_core-counts[ref]).abs().max()))
    (dest/f'complete_{kind}_{tag}.json').write_text(json.dumps(audit))
    print('OPTICS_COMPLETE',kind,tag,clearance,audit,flush=True)
