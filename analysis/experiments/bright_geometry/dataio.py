"""Window-only access to frozen products; cache experimental inputs separately."""
from __future__ import annotations
# Mirrored from experiments/bright_geometry/dataio.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
from functools import lru_cache
from pathlib import Path
import sys

ROOT = PROJECT  # project data root (pv2mveg configuration)
import pv2mveg  # noqa: E402,F401  (PROJ/GDAL data of the rasterio wheel, no WMI)
sys.path.insert(0, str(CODE / 'eco' / 'rowband'))
import rb_common as R  # noqa: E402
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling, transform_bounds

OUT = ROOT / 'data' / 'experiments' / 'bright_geometry'
PROBS = {
    '06': ROOT / 'data/products_tmp_v5a/06/prob_0p5m.tif',
    '08': ROOT / 'data/products_tmp_v5a/08/prob_0p5m.tif',
    '2026_08': ROOT / 'data/eco/mask2026/work/tmp_final/2026_08/prob_0p5m.tif',
}


@lru_cache(maxsize=3)
def scenes(tag):
    return R.scenes(tag)


def warp_window(path, bounds, tf, crs, shape, bands=1, resampling=Resampling.nearest, fill=np.nan):
    bands = [bands] if isinstance(bands, int) else list(bands)
    out = np.full((len(bands), *shape), fill, np.float32)
    with rasterio.open(path) as ds:
        sb = transform_bounds(crs, ds.crs, *bounds, densify_pts=8)
        pad = 3 * abs(ds.transform.a)
        sb = (sb[0] - pad, sb[1] - pad, sb[2] + pad, sb[3] + pad)
        win = from_bounds(*sb, ds.transform).round_offsets().round_lengths()
        src = ds.read(bands, window=win, boundless=True, fill_value=ds.nodata or 0).astype(np.float32)
        for k in range(len(bands)):
            reproject(src[k], out[k], src_transform=ds.window_transform(win), src_crs=ds.crs,
                      src_nodata=ds.nodata, dst_transform=tf, dst_crs=crs,
                      dst_nodata=fill, resampling=resampling)
    return out[0] if len(bands) == 1 else out


def write_raster(path, arr, meta, nodata=255):
    from affine import Affine
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, 'w', driver='GTiff', height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype=arr.dtype, transform=Affine(*meta['transform']),
                       crs=f"EPSG:{meta['native_epsg']}", nodata=nodata, compress='deflate') as ds:
        ds.write(arr, 1)
