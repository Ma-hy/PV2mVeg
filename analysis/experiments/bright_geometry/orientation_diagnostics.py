"""Release-isolated diagnostics for the verified local-orientation repair.

This adapter never changes frozen v1/v0.11 files. Treatment measurements are
recomputed on new masks; retained original same-scene control sets are frozen.
All qualified source satellites enter every diagnostic unless an actual band
or registration-coverage requirement prevents an observation.
"""
# Mirrored from experiments/bright_geometry/orientation_diagnostics.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '1')
import hashlib
import inspect
import json
from pathlib import Path
from functools import lru_cache
_SOURCE_HASH_AT_IMPORT=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
import numpy as np
import pandas as pd
import ecology_revision as R
import revision_measure as M
import revision_statistics as S
import revision_optics as O
import revision_pan as P
import revision_crossperiod as X

ROOT = PROJECT  # project data root (pv2mveg configuration)
NEW_ROOT = ROOT/'data/eco/releases/fixed_ew_v2_orientation'
OLD_ROOT = ROOT/'data/eco/releases/fixed_ew_v1'
CELLS = {t: f'a1_{t}_det_thr50_cells.csv' for t in R.TAGS}
_ORIGINALS = {'measure_init': M.init, 'optics_worker': O.worker,
              'pan_job': P.job, 'cross_worker': X.worker}
_INSTALLED = False
TRANSFORMS = []


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def rewrite(module, name, replacements, extra=None):
    """Small, exact, fail-closed edits to legacy functions, logged by SHA256."""
    original = getattr(module, name)
    source = inspect.getsource(original)
    edited = source
    for old, new in replacements:
        if edited.count(old) != 1:
            raise AssertionError((module.__name__, name, 'source drift', old))
        edited = edited.replace(old, new)
    # Keep the live module globals: EW/SP/CHR are deliberately changed between
    # clearance passes. A copied globals dictionary would pin a stale path.
    env = original.__globals__
    env.update(extra or {})
    exec(compile(edited, f'<orientation_adapter:{module.__name__}.{name}>', 'exec'), env)
    setattr(module, name, env[name])
    TRANSFORMS.append(dict(module=module.__name__, function=name,
                           original_sha256=hashlib.sha256(source.encode()).hexdigest(),
                           adapted_sha256=hashlib.sha256(edited.encode()).hexdigest()))


def configure(clearance='1p5m'):
    import orientation_statistics as N
    root = N.configure(clearance)
    # All legacy consumers now address a neutrally named, all-source table.
    import a2_cells as B
    import e3_pan_texture as T
    import e6_texture_analysis as E6
    R.pp.CELLS_CSV = dict(CELLS); B.CELLS_CSV = dict(CELLS); T.CELLS_CSV = dict(CELLS); E6.CELLS_CSV = dict(CELLS)
    return root


def init_measure(kind, clearance):
    install(); configure(clearance)
    _ORIGINALS['measure_init'](kind, clearance)


def optical_worker(args):
    install(); configure(args[1])
    tag, clearance, kind, row, path = args
    if kind == 'a4':
        path = Path(path)
        if path.exists(): return str(path)
        import a4_ew_forward as C
        import a2_cells as B
        root = R.folder(clearance)
        C.EW = str(root/'adjacency/ew_only')
        sc = R.scenes(tag)
        cam = B.camera_of(sc[int(row['ir_scene'])-1]['sid'], sc[int(row['ir_scene'])-1]['sensor'])
        params = pd.read_csv(root/'adjacency/ew_only/a4_shadow_params_camera.csv', dtype={'tag':str})
        pp = params[(params.tag == tag) & (params.camera == cam)]
        assert len(pp) == 1, (tag, cam, 'missing camera shadow parameters')
        width=float(pp.W.iloc[0]);kk=camera_kernels(tag,clearance)
        proof_path=root/'adjacency/ew_only'/f'a4_reuse_validation_{tag}.json'
        reuse=(cam in ['MS1','MS2'] and a4_reuse_ready(tag,clearance,cam) and unchanged_geometry(tag,int(row['cell'])))
        old=OLD_ROOT/'full_v4'/clearance/'adjacency/ew_only'/f'cache_a4_{tag}/{int(row["cell"])}.json'
        if reuse and old.exists():
            record=hybrid_a4(tag,clearance,row,kk,sc,width,old)
            data={'cells':[record]}
            R.atomic_cache(path.with_suffix('.reuse.json'),dict(kind='a4_base_weights_only',tag=tag,cell=int(row['cell']),camera=cam,
                clearance=clearance,source=str(old.relative_to(ROOT)),source_sha256=sha(old),validation_sha256=sha(proof_path),
                recomputed='All metadata, interior spectra, shadow wb/ws weights and every redabs weight; only PAN-parameter-independent w weights reused.'))
        else:data = {'cells': C.process(tag, pd.DataFrame([row]), kk, sc, width)}
        R.atomic_cache(path, data)
        return str(path)
    return _ORIGINALS['optics_worker'](args)


def pan_worker(args):
    install(); configure(args[1])
    tag,clearance,row,scene,kind,cache=args;cache=Path(cache)
    current_primary=R.folder('1p5m')/'panprof'/('cache_'+tag)/cache.name
    if kind=='control' and not cache.exists() and (clearance=='1p5m' or not current_primary.exists()):
        original=OLD_ROOT/'full_v4/1p5m/panprof'/('cache_'+tag)/cache.name
        if original.exists():
            signature=cache_signature(tag,row)
            if signature['identical']:
                R.atomic_cache(cache,json.loads(original.read_text(encoding='utf-8')))
                proof=cache.with_suffix('.reuse.json')
                R.atomic_cache(proof,dict(**signature,kind='control',scene=int(scene),source=str(original.relative_to(ROOT)),source_sha256=sha(original)))
                return str(cache)
    return _ORIGINALS['pan_job'](args)


def cross_worker(args):
    install(); configure('1p5m')
    return _ORIGINALS['cross_worker'](args)


def annual_conversion(tag, table, scene_column='ir_scene'):
    """Return per-unit, source-satellite OLI coefficients (never A1 for A4/N3)."""
    source = {s['index']:s['sensor'] for s in R.E.cc.bands_json(tag)['scenes']}
    sensor = table[scene_column].map(source)
    coef = pd.read_csv(R.ROOT/'data/eco/ts/ew/anchored'/f'conversion_{tag}.csv').set_index(['group','band'])
    suffix=scene_column.removeprefix('ir_scene')
    valid_column='valid_same'+suffix
    valid=table[valid_column].to_numpy(bool) if valid_column in table else np.ones(len(table),bool)
    values={}
    for band in ['red','nir']:
        for k in ['a','b']:
            array=np.full(len(table),np.nan)
            for i,(scene,sat) in enumerate(zip(table[scene_column],sensor)):
                if pd.isna(sat):
                    assert (pd.isna(scene) or int(scene)==0) and not valid[i], (tag,scene,valid[i],'unknown source for eligible unit')
                    continue
                assert (sat,band) in coef.index,(tag,sat,band,'missing source-specific OLI conversion')
                array[i]=coef.loc[(sat,band),k]
            values[band,k]=array
    class Conversion:
        def __init__(self):
            self.loc = values
    return Conversion()


def reference_band_available(tag, reference_scene):
    """A calibrated overlap fragment cannot qualify an uncalibrated reference."""
    scenes,coefficients=M.scene_data(tag)
    return np.array([int(i)>0 and scenes[int(i)-1]['sid'] in coefficients for i in reference_scene],bool)


def install():
    global _INSTALLED
    if _INSTALLED: return
    # N.configure must be called first: it owns new-mask geometry and main scope.
    import orientation_statistics as N
    N.configure('1p5m')
    # Retain its callable without mutating orientation_statistics.configure.
    R.configure = configure
    M.init = init_measure; O.worker = optical_worker; P.job = pan_worker; X.worker = cross_worker
    import revision_subsets as U
    rewrite(U, 'run', [("        if tag=='08':\n            sensor={s['index']:s['sensor'] for s in R.E.cc.bands_json(tag)['scenes']}\n            valid &= np.array([sensor.get(int(sc))=='BJ3A1' for sc in ref])\n", ''),
                       ("bgok=valid&(z.n_bg>=50)&(bgsel.sum(1)[p]>=2)",
                        "bgok=valid&(z.n_bg>=50)&(bgsel.sum(1)[p]>=2)&reference_band_available(tag,ref)")],
            {'reference_band_available':reference_band_available})
    import revision_temporal as T
    rewrite(T, 'temporal', [("conv=pd.read_csv(ts/'ew/anchored'/f'conversion_{tag}.csv',dtype={'period':str}).query(\"group=='BJ3A1'\").set_index('band')",
                            'conv=annual_conversion(tag,t)'),
                           ("conversion='Frozen outside-array BJ3A1-to-OLI coefficients for annual series.'",
                            "conversion='Frozen outside-array source-satellite-to-OLI coefficients, selected per unit scene; all qualified satellites.'")],
            {'annual_conversion': annual_conversion})
    rewrite(O,'run',[("    (dest/f'complete_{kind}_{tag}.json').write_text(json.dumps(audit))",
                       "    assert audit['max_native_count_difference']==0, (tag,clearance,kind,audit)\n    (dest/f'complete_{kind}_{tag}.json').write_text(json.dumps(audit))")])
    rewrite(S,'base_stats',[("for a,z in t.groupby('age'):selection.append",
                             "for a,z in t[t.parent].groupby('age'):selection.append")])
    import revision_forward as F
    rewrite(F, 'main', [("cv=pd.read_csv(ts/'ew/anchored'/f'conversion_{tag}.csv').query(\"group=='BJ3A1'\").set_index('band')",
                        "cv=annual_conversion(tag,z,'ir_scene_'+side)")], {'annual_conversion': annual_conversion})
    import revision_aux_summary as AS
    rewrite(AS,'main',[("M.tables('06' if d in dates[:2] else '08','1p5m')",
                        "M.tables('06' if d in dates[:2] else '08',p)"),
                       ("eligibility='Primary 2 m valid cells, with >=50 auxiliary footprints and local or matched available controls.'",
                        "eligibility='Corresponding-clearance 2 m valid cells from all qualified source satellites, with >=50 auxiliary footprints and local or matched available controls.'")])
    # Legacy PAN routines swallowed IO/registration exceptions and continued.
    # Empty returns for genuine support insufficiency remain valid; unexpected
    # exceptions must abort, never be converted to apparent satellite absence.
    import e3_pan_texture as E3
    rewrite(E3,'work',[("        try:\n            r, p = one(tag, *j, scen, tg, crs, tk, rok, cok, c['source'], rng)\n        except Exception as ex:  # noqa: BLE001\n            A.log(f'{tag} {j[0]} 失败 {ex}', 'log_e3.txt'); r = None",
                        "        r, p = one(tag, *j, scen, tg, crs, tk, rok, cok, c['source'], rng)")])
    import p1_extract as P1
    rewrite(P1,'work',[("        try:\n            r, p, g = one(tag, j, ctx)\n        except Exception as ex:  # noqa: BLE001\n            pp.log(f'{tag} {j[4]} {j[0]} 失败 {ex!r}')\n            continue",
                        "        r, p, g = one(tag, j, ctx)")])
    _INSTALLED = True
    configure('1p5m')


def cache_signature(tag, row):
    """Hash old/new native geometry with a conservative 150 m cell halo.

    This is larger than every measurement/optical convolution halo. Only exact
    equality of rowmask, certification and all three support layers permits
    reusing an old measurement. New source coefficients are immutable inputs.
    """
    r,c,h,w,_ = M.cell_window(tag,row,pad=75)
    g = R.F.B.grid_of(tag)
    bounds = (g['x0']+c*2,g['y0']-(r+h)*2,g['x0']+(c+w)*2,g['y0']-r*2)
    signatures = []
    for root in [OLD_ROOT, NEW_ROOT]:
        hs = hashlib.sha256()
        for name in ['rowband_0p5m.tif','certified_ew_0p5m.tif']:
            ds = R.dataset(str(root/tag/name))
            win = __import__('rasterio').windows.from_bounds(*bounds,ds.transform).round_offsets().round_lengths()
            arr = ds.read(1,window=win,boundless=True,fill_value=255 if name.startswith('rowband') else 0)
            hs.update(arr.tobytes())
        for name in R.LAYERS.values():
            hs.update(R.read(root/'ecology_1p5_v3'/tag/(name+'.tif'),r,c,h,w,fill=255).tobytes())
        signatures.append(hs.hexdigest())
    return dict(tag=tag,cell=int(row.cell),old=signatures[0],new=signatures[1],identical=signatures[0]==signatures[1])


@lru_cache(maxsize=9)
def camera_kernels(tag,clearance):
    import a2_cells as B
    root=configure(clearance);kk=B.load_kernels(tag)
    k2=pd.read_csv(root/'adjacency/k2_kernels.csv',dtype={'tag':str})
    k2=k2[(k2.tag==tag)&(k2.variant=='main')]
    ef=pd.read_csv(root/'adjacency/ew_only/k3_pan_esf_camera.csv',dtype={'tag':str})
    ef=ef[(ef.tag==tag)&(ef.age=='0–2')].set_index('camera')
    for _,r in k2.iterrows():
        if r.camera not in ef.index:continue
        sp=float(ef.loc[r.camera,'sigma_psf_upper']);th=[r.dx,r.dy,r.sx,r.sy,r.al,r.s2]
        kk['redabs',r.camera,r.band]=R.A.param_kernel([r.dx,r.dy,np.hypot(r.sx,sp),np.hypot(r.sy,sp),0,2] if r.band=='red' else th)
        for name,f in [('s08',.8),('s12',1.2)]:kk[name,r.camera,r.band]=R.A.param_kernel([r.dx,r.dy,r.sx*f,r.sy*f,0,2])
    return kk


@lru_cache(maxsize=3)
def unchanged_cells(tag):
    path=NEW_ROOT/'analysis/verification/cache_geometry_signatures.csv'
    d=pd.read_csv(path,dtype={'tag':str});return set(d[(d.tag==tag)&d.identical].cell.astype(int))


def unchanged_geometry(tag,cell):return cell in unchanged_cells(tag)


@lru_cache(maxsize=24)
def a4_reuse_ready(tag,clearance,camera):
    root=R.folder(clearance);dest=root/'adjacency/ew_only';path=dest/f'a4_reuse_validation_{tag}.json'
    if not path.exists():return False
    proof=json.loads(path.read_text())
    if proof.get('status')!='passed' or sum(x['camera']==camera for x in proof.get('cases',[]))<2:return False
    files={'k2_kernels.csv':root/'adjacency/k2_kernels.csv','k2_np_kernels.npz':root/'adjacency/k2_np_kernels.npz',
           'shadow_parameters':dest/'a4_shadow_params_camera.csv','esf_parameters':dest/'k3_pan_esf_camera.csv',
           'geometry_signatures':NEW_ROOT/'analysis/verification/cache_geometry_signatures.csv'}
    return all(proof['input_hashes'].get(name)==sha(p) for name,p in files.items())


@lru_cache(maxsize=1)
def partial_a4_function():
    """Recompute every camera-PAN-dependent field; omit only base w fields."""
    import a4_ew_forward as C
    source=inspect.getsource(C.process)
    start=source.index('                Wa = B.sample_w(Mall, kk, MF, MF, h, w)')
    end=source.index("                if var in ('main', 'wide', 'redabs', 'myoung', 'npkr'):",start)
    block=source[start:end]
    source=source[:start]+"                if var == 'redabs':\n"+''.join('    '+line if line.strip() else line for line in block.splitlines(keepends=True))+source[end:]
    source=source.replace('def process(','def process_changed_fields(',1)
    env=C.__dict__;exec(compile(source,'<orientation:a4_changed_fields>','exec'),env)
    return env['process_changed_fields']


def hybrid_a4(tag,clearance,row,kernels,scenes,width,old_path):
    old=json.loads(Path(old_path).read_text())['cells'][0]
    cam=old['camera'];assert cam in ['MS1','MS2']
    # These five families contain all shadow weights. Every wb/ws weight is
    # recalculated even when the rounded strip width happens to be unchanged.
    selected={k:v for k,v in kernels.items() if k[1]==cam and k[0] in ['main','wide','redabs','myoung','npkr']}
    current=partial_a4_function()(tag,pd.DataFrame([row]),selected,scenes,width)[0]
    current.update({k:v for k,v in old.items() if k.startswith('w_') and not k.startswith('w_redabs_')})
    return current
