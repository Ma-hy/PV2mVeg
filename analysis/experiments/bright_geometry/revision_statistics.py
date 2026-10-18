"""Unit-level robustness, confounding, spatial summaries and paired contrasts."""
# Mirrored from experiments/bright_geometry/revision_statistics.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import json
import ecology_revision as R
import revision_measure as M
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

GROUPS=[('0–1',0,1),('0–2',0,2),('3–5',3,5),('≥6',6,99),('≥10',10,99),('all',0,99)]
J=R.E.J


def boot(x,t,groups=GROUPS):
    x=np.asarray(x,float);masks={n:np.asarray((t.age>=lo)&(t.age<=hi)) for n,lo,hi in groups}
    if not np.isfinite(x).any():return []
    return [dict(group=n,estimate=v['cell'],lo=v['cell_lo'],hi=v['cell_hi'],se=v['cell_se'],station_equal=v['stn'],n=v['n'],clusters=v['ncl'])
            for n,v in J.cboot(x,t.grp.values,masks).items()]


def base_stats(clearance):
    out=R.folder(clearance)/'statistics';out.mkdir(exist_ok=True)
    rows=[];reg=[];loo=[];stations=[];desc=[];cohorts=[];selection=[];bias=[];controls=[]
    tt={tag:M.tables(tag,clearance) for tag in R.TAGS}
    for tag,t in tt.items():
        v=t[t.valid_same].copy()
        for k in R.E.METRICS:
            for q in ['dadj','d','it','ic','raw','bias']:
                rows.extend(dict(tag=tag,index=k,quantity=q,**a) for a in boot(v[f'same_{k}_{q}'],v))
            for name,lo,hi in GROUPS:
                z=v[(v.age>=lo)&(v.age<=hi)];x=z[f'same_{k}_dadj']
                desc.append(dict(tag=tag,index=k,group=name,n=len(z),median=x.median(),q05=x.quantile(.05),q25=x.quantile(.25),q75=x.quantile(.75),q95=x.quantile(.95),positive=(x>0).mean(),treated=z[f'same_{k}_it'].mean(),control=z[f'same_{k}_ic'].mean()))
            z=v[v.age>=6];st=z.groupby('grp')[f'same_{k}_dadj'].agg(['sum','count'])
            for g,a in st.iterrows():loo.append(dict(tag=tag,index=k,omitted_cluster=g,estimate=(st['sum'].sum()-a['sum'])/(st['count'].sum()-a['count'])))
            for a,z in v.groupby('age'):
                cohorts.append(dict(tag=tag,index=k,age=a,cohort=int(z['T'].iloc[0]),n=len(z),estimate=z[f'same_{k}_dadj'].mean()))
            for b in [.2,.3,.4,.5,.6]:
                z=v[(v.pf_mean>=b)&(v.pf_mean<(b+.1 if b<.6 else 1.01))]
                bias.extend(dict(tag=tag,index=k,pf_lo=b,**a) for a in boot(z[f'same_{k}_bias'],z,[('all',0,99)]))
            if k not in ['ndvi','nirv','msavi']:continue
            z=v.dropna(subset=['pitch_m','gcr','ndvi_pre0010','elev','x','y',f'same_{k}_dadj']).copy()
            z['dep']=z[f'same_{k}_dadj'];z['x_km']=(z.x-z.x.mean())/1000;z['y_km']=(z.y-z.y.mean())/1000
            for name,f in [('unadjusted','dep ~ age'),('adjusted','dep ~ age + pitch_m + gcr + ndvi_pre0010 + elev + x_km + y_km')]:
                fit=smf.ols(f,z).fit(cov_type='cluster',cov_kwds={'groups':z.grp})
                for term in fit.params.index:
                    reg.append(dict(tag=tag,index=k,model=name,term=term,estimate=fit.params[term],lo=fit.conf_int().loc[term,0],hi=fit.conf_int().loc[term,1],p=fit.pvalues[term],n=int(fit.nobs),r2=fit.rsquared))
        for st,z in v.groupby('station'):
            stations.append(dict(tag=tag,station=st,n=len(z),x=z.x.mean(),y=z.y.mean(),ndvi=z.same_ndvi_dadj.mean(),nirv=z.same_nirv_dadj.mean(),msavi=z.same_msavi_dadj.mean(),age=z.age.mean()))
        for a,z in t.groupby('age'):selection.append(dict(tag=tag,age=a,parent=len(z),valid=int(z.valid_same.sum()),station_clusters=z.loc[z.valid_same,'grp'].nunique()))
        d=R.load(tag,clearance)
        for cutoff in [20,100,200]:
            dd=R.load(tag,clearance,cutoff);idx=dd['tidx'];z=t.copy()
            valid=dd['valid']['same'];z=z.loc[valid]
            for k in ['ndvi','nirv','msavi']:
                x,_=R.E.A1.per_cell(dd,'same',k)
                controls.extend(dict(tag=tag,index=k,scenario=f'minimum_{cutoff}_pixels',**a) for a in boot(x['dadj'][valid],z))
        # Control screening geometry is unchanged; update the treatment support
        # and use the same original matched neighbours, screened by each rule.
        for name,sub in [('screen_all',''),('screen_targeted','var_targeted'),('screen_dark','var_dark')]:
            path=R.ROOT/'data/eco/chrono/ew_ctrlclean'/sub/'celltab'/f'celltab_{tag}_raw.parquet'
            if not path.exists():raise FileNotFoundError(path)
            clean=pd.read_parquet(path);dd=dict(d);dd['tabs']=dict(d['tabs']);tab=d['tabs']['raw'].copy()
            for c in clean.columns:
                if c.startswith('ctrl_'):tab[c]=clean[c]
            R.E.cc.add_index(tab,'ctrl');dd['tabs']['raw']=tab
            dd['sel']=dict(d['sel']);dd['sel']['same']=d['sel']['same']&(tab.ctrl_n.values[d['M']]>=2000)
            valid=d['valid']['same']&(dd['sel']['same'].sum(1)>=2);dd['valid']=dict(d['valid']);dd['valid']['same']=valid
            for k in ['ndvi','nirv','msavi']:
                x,_=R.E.A1.per_cell(dd,'same',k);base,_=R.E.A1.per_cell(d,'same',k)
                for qty,a in [('effect',x['dadj']),('paired_change',x['dadj']-base['dadj'])]:
                    controls.extend(dict(tag=tag,index=k,scenario=name+'_'+qty,**b) for b in boot(a[valid],t.loc[valid]))
    pairs=[]
    for a,b in [('06','08'),('08','2026_08')]:
        x=tt[a];y=tt[b];z=x[x.valid_same].merge(y[y.valid_same],on='cell',suffixes=('_a','_b'))
        z['age']=z.age_a;z['grp']=z.grp_a
        for k in R.E.METRICS:
            pairs.extend(dict(period_a=a,period_b=b,index=k,support='same units; own-date pixels and same-scene controls',**v) for v in boot(z[f'same_{k}_dadj_b']-z[f'same_{k}_dadj_a'],z))
    for name,data in [('groups',rows),('regression',reg),('leave_station_out',loo),('station_effects',stations),('distribution',desc),('cohorts',cohorts),('selection',selection),('mixed_bias',bias),('robustness',controls),('period_pairs',pairs)]:
        pd.DataFrame(data).to_csv(out/(name+'.csv'),index=False)
    (out/'base_complete.json').write_text(json.dumps(dict(clearance=clearance,bootstrap=J.B,interval='Secondary conditional station bootstrap; original frozen controls travel with treatment. Main inference retains independent control-block resampling in ecology_comparison.csv.')))
    print('STATS_COMPLETE',clearance,flush=True)
