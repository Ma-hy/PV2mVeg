"""One-year forward/backward test of the re-anchored reconstruction."""
# Mirrored from experiments/bright_geometry/revision_forward.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import ecology_revision as R
import revision_measure as M
import numpy as np
import pandas as pd


def main():
    ts=R.ROOT/'data/eco/ts';ls=np.load(ts/'cell_y.npz')
    s2={y:np.load(ts/'s2'/f'cell_s2_{y}.npz') for y in [2025,2026]}
    rows=[];records=[]
    for p in R.POLICIES:
        a=M.tables('08',p).query('valid_same');b=M.tables('2026_08',p).query('valid_same')
        z=a.merge(b,on='cell',suffixes=('_a','_b'));z=z[(z.age_a>=1)&(z.pf_mean_a<=.6)].copy();ids=z.cell.values
        for sensor in ['Sentinel-2','Landsat']:
            obs={}
            for side,tag in [('a','08'),('b','2026_08')]:
                cv=pd.read_csv(ts/'ew/anchored'/f'conversion_{tag}.csv').query("group=='BJ3A1'").set_index('band')
                for surface in ['it','ic']:
                    for band in ['red','nir']:
                        v=z[f'same_{band}_{surface}_{side}'].values
                        obs[side,surface,band]=v if sensor=='Sentinel-2' else cv.loc[band,'a']+cv.loc[band,'b']*v
            dm={band:(s2[2026]['m08_'+band][ids]-s2[2025]['m08_'+band][ids] if sensor=='Sentinel-2' else ls['sum_'+band][ids,26]-ls['sum_'+band][ids,25])/(1-z.pf_mean_a.values) for band in ['red','nir']}
            for k in ['ndvi','nirv','msavi']:
                fn=R.E.cc.C.INDEX_FUN[k]
                it={s:fn(obs[s,'it','red'],obs[s,'it','nir']) for s in ['a','b']}
                ic={s:fn(obs[s,'ic','red'],obs[s,'ic','nir']) for s in ['a','b']}
                observed=(it['b']-ic['b'])-(it['a']-ic['a'])
                forward=fn(obs['a','it','red']+dm['red'],obs['a','it','nir']+dm['nir'])-it['a']-(ic['b']-ic['a'])
                backward=it['b']-fn(obs['b','it','red']-dm['red'],obs['b','it','nir']-dm['nir'])-(ic['b']-ic['a'])
                # No change of the ground still allows the observed change
                # of controls; all comparisons use identical valid units.
                null=-(ic['b']-ic['a'])
                for direction,pred in [('forward',forward),('backward',backward)]:
                    for group,keep in [('all',np.ones(len(z),bool)),('≥6',z.age_a.values>=6)]:
                        good=keep&np.isfinite(observed)&np.isfinite(pred)&np.isfinite(null)
                        er=pred[good]-observed[good];en=null[good]-observed[good]
                        rows.append(dict(clearance=p,sensor=sensor,index=k,direction=direction,group=group,n=int(good.sum()),bias=float(np.mean(er)),rmse=float(np.sqrt(np.mean(er**2))),null_rmse=float(np.sqrt(np.mean(en**2))),skill=float(1-np.mean(er**2)/np.mean(en**2))))
    pd.DataFrame(rows).to_csv(R.OUT/'temporal/forward_backward.csv',index=False)
    print('FORWARD_COMPLETE',flush=True)
