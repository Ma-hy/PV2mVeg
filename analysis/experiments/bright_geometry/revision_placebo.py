"""Leave-one-control-out diagnostics for the revised eligible treatment units."""
# Mirrored from experiments/bright_geometry/revision_placebo.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import ecology_revision as R
import revision_measure as M
import numpy as np
import pandas as pd


def main():
    ys=np.arange(2000,2027);Y=R.E.cc.landsat_Y();rows=[];summary=[]
    for tag in R.TAGS:
        d=R.load(tag,'1p5m');t=M.tables(tag,'1p5m');ok=d['valid']['same'];tab=d['tabs']['raw'];mm=d['M'][ok];ss=d['sel']['same'][ok];tt=t[ok]
        for k in R.E.METRICS:
            val=tab['ctrl_'+k].values[mm]
            base=np.nanmean(np.where(ys[None,None,:]<=d['T'][ok,None,None]-2,Y[k][mm],np.nan),axis=2)
            for quantity,x in [('direct',val),('baseline_adjusted',val-base)]:
                total=np.where(ss,x,0).sum(1);n=ss.sum(1)
                placebo=x-(total[:,None]-x)/(n[:,None]-1)
                for i,z in enumerate(tt.itertuples()):
                    vals=placebo[i,ss[i]];assert abs(vals.mean())<1e-10
                    if k=='ndvi' and quantity=='baseline_adjusted':
                        for ctrl,v in zip(mm[i,ss[i]],vals):rows.append(dict(tag=tag,cell=z.cell,ctrl=int(ctrl),age=z.age,value=float(v)))
                for name,lo,hi in [('all',0,99),('0–2',0,2),('3–5',3,5),('≥6',6,99)]:
                    use=(tt.age.values>=lo)&(tt.age.values<=hi);p=placebo[use][ss[use]]
                    summary.append(dict(tag=tag,index=k,quantity=quantity,group=name,n_pairs=len(p),n_units=int(use.sum()),mean=float(p.mean()),sd=float(p.std()),q025=float(np.quantile(p,.025)),q975=float(np.quantile(p,.975))))
    out=R.folder()/'statistics';pd.DataFrame(rows).to_parquet(out/'placebo_pairs.parquet',index=False);pd.DataFrame(summary).to_csv(out/'placebo_summary.csv',index=False)
    print('PLACEBO_COMPLETE',flush=True)
