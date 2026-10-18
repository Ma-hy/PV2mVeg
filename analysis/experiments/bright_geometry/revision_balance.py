"""Balance and held-out behaviour of the actually used same-scene subsets."""
# Mirrored from experiments/bright_geometry/revision_balance.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import numpy as np
import pandas as pd


def main():
    es=R.E.cc.es;C=R.E.cc.C;ts=R.ROOT/'data/eco/ts'
    Y=np.load(ts/'cell_y.npz');nd=Y['sum_ndvi'];red=Y['sum_red']
    cells=pd.read_parquet(ts/'ew/cells_pool.parquet');pool=np.flatnonzero(cells.pool.values)
    clim=C.climate();pz=np.array([clim.p_mjjas_z.get(y,np.nan) for y in es.YEARS])
    feature={}
    for year in range(2008,2027):
        fp,n=es.features(nd,red,cells,pool,year,pz,6)
        need=max(3,int(np.ceil(.6*(es.YEARS<=year-6).sum())))
        good=np.isfinite(fp).all(1)&(n>=need)
        feature[year]=(fp[good].mean(0),fp[good].std(0))
    for p in R.POLICIES:
        rows=[];held=[];counts=[]
        for tag in R.TAGS:
            d=R.load(tag,p);t=M.tables(tag,p);v=d['valid']['same'];t=t[v].copy()
            ids=d['tidx'][v];mat=d['M'][v];sel=d['sel']['same'][v];T=d['T'][v]
            yc=np.nanmean(np.where(sel[:,:,None],nd[mat],np.nan),axis=1)
            for kind,yy in [('same_scene',yc),('all_ten',np.nanmean(nd[mat],axis=1))]:
                hd,hr=es.holdout(nd[ids],yy,T)
                mean=np.nanmean(hd,axis=1)
                slope=np.where(np.isfinite(hd).all(1),((hd-np.nanmean(hd,axis=1)[:,None])*np.array([-1.5,-.5,.5,1.5])).sum(1)/5,np.nan)
                for quantity,x in [('heldout_mean',mean),('heldout_slope',slope)]:held.extend(dict(tag=tag,controls=kind,quantity=quantity,**a) for a in S.boot(x,t))
            for year in np.unique(T):
                use=T==year;mu,sd=feature[year]
                ft,_=es.features(nd,red,cells,ids[use],year,pz,6)
                fm=np.stack([es.features(nd,red,cells,mat[use,j],year,pz,6)[0] for j in range(mat.shape[1])],axis=1)
                for kind,mask in [('same_scene',sel[use]),('all_ten',np.ones_like(sel[use]))]:
                    cm=np.nanmean(np.where(mask[:,:,None],fm,np.nan),axis=1)
                    for j,name in enumerate(es.FEAT_NAMES):rows.append(dict(tag=tag,T=year,controls=kind,feature=name,n=int(use.sum()),smd_before=(np.nanmean(ft[:,j])-mu[j])/sd[j],smd_after=np.nanmean(ft[:,j]-cm[:,j])/sd[j]))
            dist=np.hypot(cells.x.values[ids,None]-cells.x.values[mat],cells.y.values[ids,None]-cells.y.values[mat])
            counts.append(dict(tag=tag,n_units=len(ids),median_same_controls=float(np.median(sel.sum(1))),median_same_distance=float(np.median(dist[sel])),median_all_distance=float(np.median(dist))))
        out=R.folder(p)/'statistics';pd.DataFrame(rows).to_csv(out/'balance.csv',index=False);pd.DataFrame(held).to_csv(out/'holdout.csv',index=False);pd.DataFrame(counts).to_csv(out/'controls.csv',index=False)
        print('BALANCE_COMPLETE',p,flush=True)
