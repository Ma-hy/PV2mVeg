"""Common spatial support, without treating reprojected samples as independent."""
# Mirrored from experiments/bright_geometry/revision_crossperiod.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
import json
from pathlib import Path
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import numpy as np
from rasterio.windows import Window, from_bounds, transform
from rasterio.warp import reproject, transform_bounds, Resampling


def mapped(path,ref,r,c,h,w,method=Resampling.average,fill=0):
    ds=R.dataset(str(path));dt=transform(Window(c,r,w,h),ref.transform)
    if ds.crs==ref.crs and ds.transform==ref.transform:
        return R.read(path,r,c,h,w,fill=fill)
    bounds=(dt.c,dt.f-h*2,dt.c+w*2,dt.f)
    sb=transform_bounds(ref.crs,ds.crs,*bounds)
    win=from_bounds(*sb,ds.transform)
    x=int(np.floor(win.col_off))-3;y=int(np.floor(win.row_off))-3
    win=Window(x,y,int(np.ceil(win.width))+7,int(np.ceil(win.height))+7)
    arr=ds.read(1,window=win,boundless=True,fill_value=fill)
    dest=np.full((h,w),float(fill),np.float64)
    # No nodata renormalisation: an invalid contributing source must remain
    # visible to the conservative validity/support gate below.
    reproject(arr,dest,src_transform=transform(win,ds.transform),src_crs=ds.crs,
              dst_transform=dt,dst_crs=ref.crs,resampling=method,init_dest_nodata=False)
    return dest


def worker(args):
    a,b,row,path=args;path=Path(path)
    if path.exists():return str(path)
    r,c,h,w,inc=M.cell_window(a,row)
    ref=R.dataset(str(R.ROOT/'data/eco/bj3ndvi'/a/'red_2m.tif'))
    spectra={};scenes={};ok={};supports={}
    for tag in [a,b]:
        base=R.ROOT/'data/eco/bj3ndvi'/tag
        spectra[tag]={band:mapped(base/(band+'_2m.tif'),ref,r,c,h,w,fill=-32768)/1e4 for band in ['red','nir']}
        scenes[tag]=mapped(base/'scene_2m.tif',ref,r,c,h,w,Resampling.nearest,0).astype(int)
        for p in R.POLICIES:
            layer=R.F.ecology_dir('clearance1p5')/tag/(R.LAYERS[p]+'.tif')
            mn=mapped(layer,ref,r,c,h,w,Resampling.min,0)
            mx=mapped(layer,ref,r,c,h,w,Resampling.max,255)
            # Every source footprint contributing to a target footprint must
            # be selected; mixed 0/1 and 1/255 target footprints are excluded.
            supports[tag,p]=(mn==1)&(mx==1)
    out=[]
    for p in R.POLICIES:
        keep=inc&supports[a,p]&supports[b,p]
        for tag in [a,b]:
            rr=spectra[tag]['red'];nn=spectra[tag]['nir']
            good=keep&(rr>-1)&(nn>-1)
            assert good.sum()==keep.sum(),(a,b,row.cell,'invalid common reflectance')
            sc=scenes[tag][keep]
            out.append(dict(cell=int(row.cell),a=a,b=b,tag=tag,clearance=p,n=int(keep.sum()),
                            sr=float(rr[keep].sum()),sn=float(nn[keep].sum()),
                            scene=int(np.bincount(sc,minlength=64).argmax()) if keep.any() else 0))
    tmp=path.with_suffix('.tmp.json');tmp.write_text(json.dumps(out));os.replace(tmp,path)
    return str(path)
