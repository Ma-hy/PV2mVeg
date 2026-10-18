"""Exact distances from native 2 m centres to occupied 0.5 m pixel squares."""
# Mirrored from experiments/bright_geometry/core_geometry.py of the research repository at d065677 (unused functions removed, roots configured, comments in English).
import sys as _sys, pathlib as _pathlib  # noqa: E401,E402
_SRC = str(_pathlib.Path(__file__).resolve().parents[3] / 'src')  # pv2mveg: configuration of the paths
_sys.path.insert(0, _SRC) if _SRC not in _sys.path else None
from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS  # noqa: E402,F401
import numpy as np
from scipy.ndimage import distance_transform_edt


def exact_center_distance(panel,h,w,halo=200):
    """Euclidean distance to the closed union of raster-pixel footprints.

    Both the 0.5 m pixel boundaries and native 2 m centres lie on the 0.5 m
    corner lattice. The closest point on any axis-aligned pixel-square edge
    to a centre lies on this lattice too. EDT on occupied corner nodes is
    therefore exact; no centre-to-edge scalar offset is needed.
    """
    if not panel.any():return np.full((h,w),np.inf)
    base=np.pad(np.asarray(panel,bool),((0,1),(0,1)))
    nodes=base.copy()
    nodes[1:]|=base[:-1];nodes[:,1:]|=base[:,:-1];nodes[1:,1:]|=base[:-1,:-1]
    distance=distance_transform_edt(~nodes,sampling=.5)
    return distance[np.ix_(halo+np.arange(h)*4+2,halo+np.arange(w)*4+2)]
