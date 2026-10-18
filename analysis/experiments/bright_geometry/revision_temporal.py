"""Re-anchor seasonal and annual series; illumination scenarios by clearance.

These reconstructions impose an unchanging mixture and are diagnostics, not
independent evidence for recovery. Upstream reflectance, fractions and matches
remain frozen to isolate the changed ground anchor.
"""
# Mirrored from experiments/bright_geometry/revision_temporal.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'1')
import json
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import numpy as np
import pandas as pd


def temporal():
    ts=R.ROOT/'data/eco/ts';out=R.OUT/'temporal';out.mkdir(exist_ok=True)
    cells=pd.read_parquet(ts/'cells_pool.parquet')
    late=pd.read_parquet(ts/'ew/anchored/anchor_cells.parquet').set_index('cell').late_year
    yr=np.arange(2000,2027);sy=np.arange(2019,2027)
    annual=np.load(ts/'cell_y.npz')
    monthly={}
    for m in ['m06','m08']:
        for b in ['red','nir']:monthly[m,b]=np.stack([np.load(ts/'s2'/f'cell_s2_{y}.npz')[m+'_'+b] for y in sy],1)
    rows=[];season=[];paired=[];availability=[]
    for p in R.POLICIES:
        anchored={}
        for tag,mo,t0 in [('06','m06',2025),('08','m08',2025),('2026_08','m08',2026)]:
            d=R.load(tag,p);t=M.tables(tag,p);ids=t.cell.values;matches=d['M']
            f=t.pf_mean.values;T=t['T'].values;last=late.reindex(ids).fillna(0).values
            ok=t.valid_same.values&(f<=.6)&(np.maximum(T,last)+1<=2026)
            g={b:t['same_'+b+'_it'].values for b in ['red','nir']}
            age=sy[None,:]-T[:,None];post=(age>=1)&(sy[None,:]>=last[:,None]+1)
            spectra={b:g[b][:,None]+(monthly[mo,b][ids]-monthly[mo,b][ids,t0-2019][:,None])/(1-f[:,None]) for b in ['red','nir']}
            for k in ['ndvi','nirv','msavi']:
                fn=R.E.cc.C.INDEX_FUN[k]
                # Seasonal controls use the same frozen same-scene subset as
                # the anchor, unlike the historical all-neighbour series.
                cv=fn(monthly[mo,'red'],monthly[mo,'nir'])[matches]
                co=np.nanmean(np.where(d['sel']['same'][:,:,None],cv,np.nan),axis=1)
                val=fn(spectra['red'],spectra['nir'])-co
                val=np.where(ok[:,None]&post,val,np.nan)
                anchored[tag,k]=(t,val)
                for scope,yy in [('2025',sy==2025),('2019–2026',np.ones(8,bool))]:
                    for group,lo,hi in S.GROUPS[:-1]:
                        use=(age>=lo)&(age<=hi)&yy[None,:]
                        v=np.nanmean(np.where(use,val,np.nan),axis=1)
                        z=t.copy();z['age']=0
                        season.extend(dict(clearance=p,tag=tag,index=k,scope=scope,age_group=group,**a) for a in S.boot(v,z,[('all',0,99)]))
            availability.append(dict(clearance=p,tag=tag,n_anchor=int(ok.sum()),n_reference=int(t.valid_same.sum())))
            if tag=='2026_08':continue
            conv=pd.read_csv(ts/'ew/anchored'/f'conversion_{tag}.csv',dtype={'period':str}).query("group=='BJ3A1'").set_index('band')
            for fmul in [.9,1.,1.1]:
                ff=np.clip(f*fmul,0,.95)
                ar={b:(conv.loc[b,'a']+conv.loc[b,'b']*g[b])[:,None]+(annual['sum_'+b][ids]-annual['sum_'+b][ids,25][:,None])/(1-ff[:,None]) for b in ['red','nir']}
                age=yr[None,:]-T[:,None];post=(age>=1)&(yr[None,:]>=last[:,None]+1)
                for k in ['ndvi','nirv','msavi']:
                    fn=R.E.cc.C.INDEX_FUN[k];cv=fn(annual['sum_red'],annual['sum_nir'])[matches]
                    co=np.nanmean(np.where(d['sel']['same'][:,:,None],cv,np.nan),axis=1)
                    raw=fn(annual['sum_red'][ids],annual['sum_nir'][ids])-co
                    pre=np.nanmean(np.where(age<=-2,raw,np.nan),axis=1)
                    val=np.where(ok[:,None]&post,fn(ar['red'],ar['nir'])-co-pre[:,None],np.nan)
                    for event in range(1,14):
                        v=np.nanmean(np.where(age==event,val,np.nan),axis=1);z=t.copy();z['age']=0
                        rows.extend(dict(clearance=p,tag=tag,index=k,f_multiplier=fmul,event=event,**a) for a in S.boot(v,z,[('all',0,99)]))
        for a,b in [('06','08'),('06','2026_08')]:
            for k in ['ndvi','nirv','msavi']:
                ta,va=anchored[a,k];tb,vb=anchored[b,k]
                ia=pd.Series(np.arange(len(ta)),index=ta.cell);ib=pd.Series(np.arange(len(tb)),index=tb.cell)
                common=ta.cell[ta.cell.isin(tb.cell)];x=va[ia.loc[common]];y=vb[ib.loc[common]];z=ta.set_index('cell').loc[common].reset_index()
                use=np.isfinite(x)&np.isfinite(y)
                for scope,jj in [('2025',sy==2025),('2019–2026',np.ones(8,bool))]:
                    v=np.nanmean(np.where(use&jj[None,:],y-x,np.nan),axis=1)
                    paired.extend(dict(clearance=p,a=a,b=b,index=k,scope=scope,**a_) for a_ in S.boot(v,z))
    for name,x in [('annual',rows),('seasonal',season),('seasonal_paired',paired),('availability',availability)]:pd.DataFrame(x).to_csv(out/(name+'.csv'),index=False)
    (out/'protocol.json').write_text(json.dumps(dict(clearances=[1,1.5,2],minimum_pixels=50,
        fractions='Frozen original row-band fractions; vary by +/-10% for annual reconstruction.',
        conversion='Frozen outside-array BJ3A1-to-OLI coefficients for annual series.',
        controls='Frozen matched same-scene subsets of each anchor; no rematching.',
        inference='Conditional whole-station bootstrap. Reconstruction is model-dependent and not used as recovery evidence.'),indent=2))
    print('TEMPORAL_COMPLETE',flush=True)
