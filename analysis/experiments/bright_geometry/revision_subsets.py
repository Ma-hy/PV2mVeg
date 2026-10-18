"""Spectra, lanes, clearance strata and boundary scenarios on new supports."""
# Mirrored from experiments/bright_geometry/revision_subsets.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import json
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import k_bg_greenness as K
import numpy as np
import pandas as pd


def run(tag):
    out=R.OUT/'measurements'/tag
    assert json.loads((out/'complete.json').read_text())['complete']
    raw=pd.read_parquet(out/'subsets.parquet')
    d=R.load(tag,'1p5m');t=M.tables(tag,'1p5m').set_index('cell')
    pos=pd.Series(np.arange(len(d['tidx'])),index=d['tidx'])
    # Frozen original same-scene neighbours; only blue/green availability can
    # further restrict their use. Do not choose a new source scene or rematch.
    B=R.E.cc.read_bands(R.E.cc.prod_path(tag,'30m_bg'))
    ci,cn=K.cell_vals(B,R.E.cc.CTRL_BANDS)
    bgsel=d['sel']['same']&(cn[d['M']]>=2000)
    ctrl={k:R.E.cc.ctrl_mean(v,d['M'],bgsel)[0] for k,v in ci.items()}
    fixed_scope=R.E.cc.scope_filter(d['cells'],tag)
    stats=[];cells=[];changes=[];gaprows=[]
    for name,z in raw.groupby('subset',sort=False):
        z=z.set_index('cell').join(t,rsuffix='_meta').reset_index()
        p=pos.loc[z.cell].values
        ref=d['tabs']['raw'].ir_scene.values[z.cell]
        # Membership and source restrictions are fixed, while boundary/lane
        # eligibility is allowed to change under its own 50-pixel gate.
        scope=fixed_scope[z.cell]
        valid=scope&(z.n>=50)&(z['mode']==ref)&(ref>0)&(d['sel']['same'].sum(1)[p]>=2)
        if tag=='08':
            sensor={s['index']:s['sensor'] for s in R.E.cc.bands_json(tag)['scenes']}
            valid &= np.array([sensor.get(int(sc))=='BJ3A1' for sc in ref])
        indices=R.E.cc.index_from_sums(z.s_red.values,z.s_nir.values,z.n.values)
        isgap='_gap' in name
        keys=['ndvi'] if isgap else R.E.METRICS
        for k in keys:
            tr=indices[k];co=z['same_'+k+'_ic'].values
            # dadj = raw contrast - preconstruction difference.
            pre=z['same_'+k+'_d'].values-z['same_'+k+'_dadj'].values
            adj=tr-co-pre
            if isgap:
                for group,lo,hi in S.GROUPS:
                    keep=valid&(z.age>=lo)&(z.age<=hi)
                    gaprows.append(dict(tag=tag,subset=name,index=k,group=group,n=int(keep.sum()),estimate=float(np.mean(adj[keep])) if keep.any() else np.nan))
            else:
                stats.extend(dict(tag=tag,subset=name,index=k,quantity='dadj',**a) for a in S.boot(adj[valid],z.loc[valid]))
                if name.startswith('boundary_'):
                    keep=valid&z.valid_same
                    change=adj-z['same_'+k+'_dadj'].values
                    changes.extend(dict(tag=tag,subset=name,index=k,**a) for a in S.boot(change[keep],z.loc[keep]))
            if k=='ndvi':
                cell=z[['cell','age','grp','n']].copy();cell['subset']=name;cell['valid']=valid;cell['ndvi_dadj']=adj;cell['ndvi']=tr
                cells.append(cell)
        if name.startswith(('ir_','lane_')):
            bgok=valid&(z.n_bg>=50)&(bgsel.sum(1)[p]>=2)
            ix=K.indices(*[z['bg_'+b].values/np.maximum(z.n_bg,1) for b in K.BANDS])
            for k,tr in ix.items():
                diff=tr-ctrl[k][p]
                rows=S.boot(diff[bgok],z.loc[bgok])
                for a in rows:
                    lo,hi=next((lo,hi) for g,lo,hi in S.GROUPS if g==a['group'])
                    take=bgok&(z.age>=lo)&(z.age<=hi)
                    a['control']=float(np.mean(ctrl[k][p][take]));a['treated']=float(np.mean(tr[take]))
                    a['relative']=a['estimate']/a['control'] if k in K.BANDS else np.nan
                stats.extend(dict(tag=tag,subset=name,index=k,quantity='direct_spectral_difference',**a) for a in rows)
    pd.DataFrame(stats).to_csv(out/'subset_summary.csv',index=False)
    pd.DataFrame(changes).to_csv(out/'boundary_paired.csv',index=False)
    pd.DataFrame(gaprows).to_csv(out/'gap_stratified.csv',index=False)
    allcells=pd.concat(cells,ignore_index=True);allcells.to_parquet(out/'subset_cells.parquet',index=False)
    paired=[]
    for a,b in [('ir_1p5m','lane_1p5m'),('d1.5_2','d2_3'),('d2_3_side1','d2_3_side3')]:
        x=allcells[(allcells.subset==a)&allcells.valid];y=allcells[(allcells.subset==b)&allcells.valid]
        z=x.merge(y,on='cell',suffixes=('_a','_b'));z['age']=z.age_a;z['grp']=z.grp_a
        paired.extend(dict(tag=tag,a=a,b=b,quantity='b_minus_a',**v) for v in S.boot(z.ndvi_dadj_b-z.ndvi_dadj_a,z))
    pd.DataFrame(paired).to_csv(out/'within_unit.csv',index=False)
    (out/'summary_complete.json').write_text(json.dumps(dict(tag=tag,complete=True,bootstrap=S.J.B,cluster='whole station, controls frozen',scope='fixed parent and original same-scene controls')))
    print('SUBSETS_COMPLETE',tag,flush=True)
