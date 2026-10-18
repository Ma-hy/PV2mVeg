"""Versioned full-paper reanalysis; frozen matches and native-grid supports.

No historical dataset is overwritten. All derived files live under full_v4.
"""
# Mirrored from experiments/bright_geometry/ecology_revision.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
import sys
import json
import importlib.util
from pathlib import Path
from functools import lru_cache
import fixed_ecology as E
import fixed_core as F
import numpy as np
import rasterio
from rasterio.windows import Window, from_bounds
from scipy.ndimage import distance_transform_edt, maximum_filter, binary_dilation

ROOT = F.D.ROOT
OUT = F.ROOT / 'full_v4'
TAGS = ['06', '08', '2026_08']
POLICIES = {'1m':'clearance1p0', '1p5m':'clearance1p5', '2m':'exact'}
LAYERS = {'1m':'pure_interrow_2m', '1p5m':'clearance_1p5m_2m', '2m':'core_2m_2m'}
for p in ['eco/adjacency', 'eco/panprof', 'eco/chrono', 'eco/bj3ndvi']:
    sys.path.insert(0, str(CODE/p))
# cfg imports eco/common.py, whereas cc uses eco/ts/common.py. Keep them distinct.
saved_common = sys.modules.get('common')
spec = importlib.util.spec_from_file_location('eco_common_revision', CODE/'eco/common.py')
cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)
sys.modules['common'] = cm
import cfg
sys.modules['common'] = saved_common
import adjc as A
import pp


def folder(clearance='1p5m'):
    p = OUT/clearance
    p.mkdir(parents=True, exist_ok=True)
    return p


def load(tag, clearance='1p5m', threshold=50):
    """Keep the original control subsets and replace treatment measurements only."""
    d = E.A1.load(tag, 'det', threshold, 'all', 'BJ3A1' if tag=='08' else 'any')
    z = np.load(F.ecology_dir(POLICIES[clearance])/tag/'core_sums.npz')
    sums = z['sums'][2]; votes = z['scenes'][2]
    ref = d['tabs']['raw'].ir_scene.values.copy()
    same = (votes.argmax(1)==ref) & (votes.sum(1)>0) & (ref>0)
    for tab in d['tabs'].values():
        for field, arr in zip(['sr','sn','n'],sums):tab['ir_'+field] = arr
        tab['ir_scene_frac']=votes.max(1)/np.maximum(votes.sum(1),1)
        E.cc.add_index(tab,'ir')
    base = E.cc.scope_filter(d['cells'],tag)[d['tidx']] & same[d['tidx']] & (sums[2,d['tidx']]>=threshold)
    if tag=='08':
        sens = {s['index']:s['sensor'] for s in E.cc.bands_json(tag)['scenes']}
        base &= np.array([sens.get(int(s),'')=='BJ3A1' for s in ref[d['tidx']]])
    for v in d['valid']:d['valid'][v] = base & (d['sel'][v].sum(1)>=(2 if v=='same' else 1))
    # We do not expose stale scene-corrected treatment spectra as new estimates.
    d['valid']['fe_red'][:] = False;d['valid']['fe_nir'][:] = False
    return d


@lru_cache(maxsize=20)
def dataset(path):return rasterio.open(path)


def read(path,r,c,h,w,band=1,fill=0):
    return dataset(str(path)).read(band,window=Window(c,r,w,h),boundless=True,fill_value=fill)


def atomic_cache(path,payload):
    """Idempotent cache publication, including concurrent Windows readers."""
    import time
    path=Path(path);content=json.dumps(payload)
    tmp=path.with_suffix(f'.{os.getpid()}.tmp.json')
    tmp.write_text(content,encoding='utf-8')
    for attempt in range(4):
        try:
            if path.exists():
                assert path.read_text(encoding='utf-8')==content,('Conflicting cache payload',str(path))
                tmp.unlink();return
            os.replace(tmp,path);return
        except PermissionError:
            if attempt==3:raise
            time.sleep(.1*(attempt+1))


@lru_cache(maxsize=3)
def scenes(tag):return cfg.scenes(cfg.pcfg(tag))


@lru_cache(maxsize=3)
def geometry(tag,r,c,h,w):
    """Qualified row geometry on any native 2 m window; no stale distance grids."""
    g=F.B.grid_of(tag);pad=F.HALO;hb=pad*.5
    x,y=g['x0']+c*2,g['y0']-r*2
    bounds=(x-hb,y-h*2-hb,x+w*2+hb,y+hb);shape=(h*4+pad*2,w*4+pad*2)
    mask=F.read_bounds(dataset(str(F.ROOT/tag/'rowband_0p5m.tif')),bounds,255,shape)
    q=F.read_bounds(dataset(str(F.ROOT/tag/'certified_ew_0p5m.tif')),bounds,0,shape)==1
    if F.D.R.PER[tag].get('cloud'):
        cl=F.read_bounds(dataset(F.D.R.PER[tag]['cloud']),bounds,0,(shape[0]//4,shape[1]//4))
        cl=binary_dilation(np.isin(cl,[1,2]),np.ones((11,11),bool))
        q &= ~np.repeat(np.repeat(cl,4,0),4,1)
    mask=np.where(q,mask,255).astype('uint8');pan=mask==1;nod=mask==255
    sl=np.s_[pad:pad+h*4,pad:pad+w*4]
    blocks=lambda a:a[sl].reshape(h,4,w,4)
    known=~blocks(nod).any((1,3))
    known &= ~blocks(maximum_filter(nod,size=9))[:,1:3,:,1:3].any((1,3))
    pure=~blocks(pan).any((1,3))
    dist=F.exact_center_distance(pan,h,w,pad)
    gap=F.B.gaps(pan,nod,h,w,'ns')
    _,ind=distance_transform_edt(~pan,return_indices=True)
    rr=pad+np.arange(h)*4+2;cc=pad+np.arange(w)*4+2
    dr=ind[0][np.ix_(rr,cc)]-rr[:,None];dc=ind[1][np.ix_(rr,cc)]-cc[None,:]
    direction=np.where(abs(dr)>=abs(dc),np.where(dr<0,1,3),np.where(dc<0,4,2)).astype('uint8')
    return mask,known,pure,dist,gap,direction


def configure(clearance='1p5m'):
    """Adapt existing read-only measurement functions inside each worker."""
    root=folder(clearance);A.OUT=str(root/'adjacency');A.FIG=str(root/'figs');A.CHR=str(root/'chrono')
    pp.OUT=str(root/'panprof');pp.CHR_EW=A.CHR;pp.CHR=A.CHR
    pp.CLEARANCE_M=F.THRESHOLDS[POLICIES[clearance]]
    pp.DISCRETE_NATIVE_PHASE=True
    for p in [A.OUT,A.FIG,pp.OUT]:Path(p).mkdir(exist_ok=True)
    def pcfg(tag):
        c=cfg.pcfg(tag);c['mask']=str(F.ROOT/tag/'rowband_0p5m.tif');return c
    A.pcfg=pcfg;A.scenes=scenes;A.CENTER_DISTANCE_OFFSET=0.
    def adapted(path,r,c,h,w,band=1,fill=0):
        name=Path(path).name;tag=Path(path).parent.name
        if tag in TAGS and name in ['class_2m.tif','dist_ctr_2m.tif','dir_2m.tif','gap_ns_2m.tif']:
            layers=F.ecology_dir('clearance1p5')/tag
            sel=read(layers/(LAYERS[clearance]+'.tif'),r,c,h,w,fill=255)
            target=sel!=255
            a=read(path,r,c,h,w,band,fill).copy()
            if not target.any():return a
            mask,known,pure,dist,gap,direction=geometry(tag,r,c,h,w)
            if name=='class_2m.tif':
                a[target]=0
                a[target&known&~pure]=1
                a[target&known&pure]=2
                # Native primary support is authoritative. Lane/edge values
                # are measured separately by revision_measure with rid/pitch.
                # Do not infer lanes from missing geometry in this adapter.
                a[sel==1]=3
            else:a[target]={'dist_ctr_2m.tif':dist,'dir_2m.tif':direction,'gap_ns_2m.tif':gap}[name][target]
            return a
        return read(path,r,c,h,w,band,fill)
    A.read_win=adapted
    return root
