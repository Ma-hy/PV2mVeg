"""Resumable native-image measurements for the full ecological revision."""
# Mirrored from experiments/bright_geometry/revision_measure.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
import json
from functools import lru_cache
import ecology_revision as R
import numpy as np
import pandas as pd
from pathlib import Path
from rasterio.windows import from_bounds


def tables(tag, clearance):
    suf='_BJ3A1' if tag=='08' else ''
    return pd.read_csv(R.folder(clearance)/'chrono'/f'a1_{tag}_det_thr50{suf}_cells.csv')


@lru_cache(maxsize=12)
def canonical_tile(tag,r,c):
    g=R.F.B.grid_of(tag);h=min(256,g['height']-r);w=min(256,g['width']-c)
    return R.F.canonical_centres(tag,g,r,c,h,w)


def coords(tag,r,c,h,w):
    """Use the identical frozen projection knots even for arbitrary cell windows."""
    x=np.zeros((h,w));y=x.copy()
    for rr in range(r//256*256,(r+h-1)//256*256+1,256):
        for cc in range(c//256*256,(c+w-1)//256*256+1,256):
            a,b=canonical_tile(tag,rr,cc)
            rlo=max(r,rr);rhi=min(r+h,rr+a.shape[0]);clo=max(c,cc);chi=min(c+w,cc+a.shape[1])
            dst=np.s_[rlo-r:rhi-r,clo-c:chi-c];src=np.s_[rlo-rr:rhi-rr,clo-cc:chi-cc]
            x[dst]=a[src];y[dst]=b[src]
    return x,y


def cell_window(tag,row,pad=2):
    g=R.F.B.grid_of(tag)
    xx=row.x+np.array([-75,75,75,-75]);yy=row.y+np.array([-75,-75,75,75])
    x,y=R.F.D.R.from_utm(tag,xx,yy)
    c=max(0,int(np.floor((min(x)-g['x0'])/2))-pad)
    r=max(0,int(np.floor((g['y0']-max(y))/2))-pad)
    w=min(g['width'],int(np.ceil((max(x)-g['x0'])/2))+pad)-c
    h=min(g['height'],int(np.ceil((g['y0']-min(y))/2))+pad)-r
    ux,uy=coords(tag,r,c,h,w)
    inc=R.F.D.R.ts_cell_index(ux,uy)==row.cell
    return r,c,h,w,inc


def init(kind,clearance):
    global KIND,CLEARANCE
    KIND=kind;CLEARANCE=clearance;R.configure(clearance)


@lru_cache(maxsize=3)
def scene_data(tag):
    return R.A.scenes(tag),json.loads((R.ROOT/'data/eco/bj3cal/s2_bg/coefficients_bg.json').read_text())


@lru_cache(maxsize=6)
def grid10(tag,name):return R.F.D.R.read10(name,R.F.D.R.PER[tag]['config'])


def extra(tag,row):
    r,c,h,w,inc=cell_window(tag,row)
    rd=lambda key,fill=0:R.read(R.ROOT/'data/eco/bj3ndvi'/tag/(key+'_2m.tif'),r,c,h,w,fill=fill)
    red=rd('red',-32768).astype(float)/1e4;nir=rd('nir',-32768).astype(float)/1e4
    scene=rd('scene');qa=rd('qa');ok=inc&((qa&R.A.QA_EX)==0)&(red>-1)&(nir>-1)
    mask,known,pure,dist,gap,direction=R.geometry(tag,r,c,h,w)
    rid=R.read(R.ROOT/'data/eco/grid2m'/tag/'region_id_2m.tif',r,c,h,w)
    gr=R.F.B.grid_of(tag);ux,uy=coords(tag,r,c,h,w)
    pitch=R.F.D.R.sample10(grid10(tag,'pitch_m'),ux,uy,fill=np.nan)
    ori=R.F.D.R.sample10(grid10(tag,'orient'),ux,uy)
    base=ok&known&pure&(rid>0)&(ori==1)
    inter=base&np.isfinite(gap)&np.isfinite(pitch)&(gap<=pitch)
    lane=base&~(np.isfinite(gap)&np.isfinite(pitch)&(gap<=pitch))
    vals={'red':red,'nir':nir};sums=[]
    scenes,coef=scene_data(tag)
    blue=np.full((h,w),np.nan);green=blue.copy()
    x,y=gr['x0']+c*2,gr['y0']-r*2
    for sc in np.unique(scene[ok]):
        if sc<1:continue
        ss=scenes[int(sc)-1]
        if ss['sid'] not in coef:continue
        ds=R.dataset(ss['path']);win=from_bounds(x,y-h*2,x+w*2,y,ds.transform).round_offsets().round_lengths()
        arr=ds.read([1,2],window=win,boundless=True,fill_value=0).astype(float)
        sel=(scene==sc)&(arr>1).all(0);cf=coef[ss['sid']]
        blue[sel]=cf['blue'][0]+cf['blue'][1]*arr[0][sel]
        green[sel]=cf['green'][0]+cf['green'][1]*arr[1][sel]
    vals.update(blue=blue,green=green)
    def record(name,sel):
        rec=dict(cell=int(row.cell),subset=name,n=int(sel.sum()),mode=int(np.bincount(scene[sel],minlength=64).argmax()) if sel.any() else 0)
        for k in ['red','nir']:rec['s_'+k]=float(vals[k][sel].sum())
        bg=sel&np.isfinite(blue)&np.isfinite(green)
        rec['n_bg']=int(bg.sum())
        for k in vals:rec['bg_'+k]=float(vals[k][bg].sum())
        sums.append(rec)
    for clearance,policy in R.POLICIES.items():
        th=R.F.THRESHOLDS[policy];sel=inter&(dist>=th)
        layer=R.read(R.F.ecology_dir('clearance1p5')/tag/(R.LAYERS[clearance]+'.tif'),r,c,h,w,fill=255)==1
        assert np.array_equal(sel,layer&inc),(tag,int(row.cell),clearance,int((sel^(layer&inc)).sum()))
        record('ir_'+clearance,sel);record('lane_'+clearance,lane&(dist>=th))
    for lo,hi in [(1,1.5),(1.5,2),(2,3),(3,4),(4,6)]:
        m=inter&(dist>=lo)&(dist<hi);record(f'd{lo:g}_{hi:g}',m)
        for side in [1,3]:record(f'd{lo:g}_{hi:g}_side{side}',m&(direction==side))
        for g in np.arange(3,12,.5):record(f'd{lo:g}_{hi:g}_gap{g:g}',m&(gap>=g)&(gap<g+.5))
    # Exact-distance boundary perturbations, at the primary 1.5 m clearance.
    pad=R.F.HALO;sl=np.s_[pad:pad+h*4,pad:pad+w*4]
    pan=mask==1;nod=mask==255
    for variant in R.F.B.SCEN[:9]:
        if variant=='base':continue
        mp=R.F.B.perturb(pan,variant)
        d=R.F.exact_center_distance(mp,h,w,pad)
        pu=~mp[sl].reshape(h,4,w,4).any((1,3))
        gp=R.F.B.gaps(mp,nod,h,w,'ns')
        sel=ok&known&pu&(rid>0)&(ori==1)&(d>=1.5)&np.isfinite(gp)&np.isfinite(pitch)&(gp<=pitch)
        record('boundary_'+variant,sel)
    return sums


def worker(job):
    tag,row,cache=job;cache=Path(cache)
    if cache.exists():return str(cache)
    if KIND=='extra':out=extra(tag,row)
    else:raise ValueError(KIND)
    tmp=cache.with_suffix('.tmp.json');tmp.write_text(json.dumps(out),encoding='utf-8');os.replace(tmp,cache)
    return str(cache)
