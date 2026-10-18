"""All-satellite inference adapter for the direction-repaired frozen cohort.

load/configure are the shared API for downstream diagnostics. No satellite
identity gate is used; original same-scene controls and scene identity persist.
"""
# Mirrored from experiments/bright_geometry/orientation_statistics.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(_k,'1')
import numpy as np
import pandas as pd
import ecology_revision as R
import revision_measure as M
from orientation_release import OLD_ROOT,NEW_ROOT,REPORT,TAGS,CLEARANCES,POLICIES,DIRS,LAYERS,old_ecology

E=R.E
GROUPS=E.GROUPS
METRICS=E.METRICS
OUT=NEW_ROOT/'analysis'
ORIGINAL_CONFIGURE=R.configure


def folder(clearance='1p5m'):
    assert clearance in CLEARANCES
    path=OUT/clearance;path.mkdir(parents=True,exist_ok=True);return path


def load(tag,clearance='1p5m',threshold=50,baseline=False):
    """Old data dictionary with replacement treatment sums and all-satellite validity."""
    assert clearance in CLEARANCES
    d=E.A1.load(tag,'det',threshold,'all','any')
    root=old_ecology(clearance,tag) if baseline else NEW_ROOT/DIRS[clearance]/tag
    z=np.load(root/'core_sums.npz');sums=z['sums'][2];votes=z['scenes'][2]
    ref=d['tabs']['raw'].ir_scene.to_numpy().copy()
    majority=np.where(votes.sum(1)>0,votes.argmax(1),0)
    same=(majority==ref)&(ref>0)&(votes.sum(1)>0)
    for tab in d['tabs'].values():
        for field,a in zip(['sr','sn','n'],sums):tab['ir_'+field]=a
        tab['ir_scene_frac']=votes.max(1)/np.maximum(votes.sum(1),1)
        E.cc.add_index(tab,'ir')
    idx=d['tidx'];parent=z['treated'][idx]
    base=E.cc.scope_filter(d['cells'],tag)[idx]&parent&same[idx]&(sums[2,idx]>=threshold)
    for v in d['valid']:d['valid'][v]=base&(d['sel'][v].sum(1)>=(2 if v=='same' else 1))
    d['valid']['fe_red'][:]=False;d['valid']['fe_nir'][:]=False
    d.update(release='fixed_ew_v1' if baseline else 'fixed_ew_v2_orientation',clearance=clearance,
             parent=parent,reference_scene=ref[idx],native_mode=majority[idx],native_sums=sums[:,idx],native_scene_votes=votes[idx])
    return d


def cell_table(tag,clearance='1p5m',baseline=False,write=True):
    d=load(tag,clearance,baseline=baseline);t=d['ct'].copy();idx=d['tidx']
    t['cell']=idx;t['T']=d['T'];t['age']=d['age'];t['grp']=np.where(t.station.to_numpy()>0,t.station,10**9+d['bt'])
    t['blk']=d['bt'];t['parent']=d['parent'];t['valid_same']=d['valid']['same'];t['all_satellite_valid']=d['valid']['same']
    raw=d['tabs']['raw']
    for k in ['sr','sn','n','scene','scene_frac','date']:t['ir_'+k]=raw['ir_'+k].to_numpy()[idx]
    t['native_mode']=d['native_mode'];t['n_same']=d['sel']['same'].sum(1);t['n_all']=d['sel']['all'].sum(1)
    scenes={s['index']:s for s in E.cc.bands_json(tag)['scenes']}
    t['sensor']=[scenes.get(int(s),{}).get('sensor','unknown') for s in d['reference_scene']]
    for k in METRICS:
        vals,valid=E.A1.per_cell(d,'same',k)
        for q,a in vals.items():t[f'same_{k}_{q}']=a
        t['finite_'+k]=valid
    if baseline:
        source=pd.read_csv(OLD_ROOT/'all_satellite_v2'/f'august_{clearance}_cells.csv') if tag=='08' else pd.read_csv(OLD_ROOT/'full_v4'/clearance/'chrono'/f'a1_{tag}_det_thr50_cells.csv')
        assert np.array_equal(t.cell,source.cell)
        flag='all_satellite_valid' if tag=='08' else 'valid_same'
        assert np.array_equal(t.valid_same,source[flag]),(tag,clearance,'baseline validity')
        for k in METRICS:assert np.allclose(t[f'same_{k}_dadj'],source[f'same_{k}_dadj'],atol=1e-12,rtol=0,equal_nan=True),(tag,clearance,k)
    if write:
        root=folder(clearance)/('baseline' if baseline else 'chrono');root.mkdir(exist_ok=True)
        path=root/f'a1_{tag}_det_thr50_cells.csv';tmp=path.with_suffix(f'.{os.getpid()}.tmp.csv')
        t.to_csv(tmp,index=False);os.replace(tmp,path)
        path=root/f'celltab_{tag}_raw.parquet';tmp=path.with_suffix(f'.{os.getpid()}.tmp.parquet')
        raw.to_parquet(tmp);os.replace(tmp,path)
    return t


def tables(tag,clearance='1p5m'):
    return pd.read_csv(folder(clearance)/'chrono'/f'a1_{tag}_det_thr50_cells.csv')


def configure(clearance='1p5m'):
    """Route existing local-window diagnostics into the repaired versioned release."""
    R.OUT=OUT;R.F.ROOT=NEW_ROOT;R.load=load;R.folder=folder;M.tables=tables
    # Every long-lived dataset/geometry cache belongs to a specific release.
    R.geometry.cache_clear();R.dataset.cache_clear()
    return ORIGINAL_CONFIGURE(clearance)
