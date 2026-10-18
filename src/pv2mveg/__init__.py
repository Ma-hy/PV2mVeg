"""PV2mVeg: inter-row vegetation greenness of solar-park arrays from sub-metre satellite imagery.

Code and derived data behind the manuscript "Quantifying vegetation greenness between photovoltaic panel
rows using sub-metre satellite imagery: A case study of the Talatan solar park, China" (version 2.0).

Modules
  config        paths and settings (data layouts, output folder, figure font)
  catalog       the published figures: script, data tier and plotted-value tables
  plotting      figure style and the save pipeline of the published figures
  data          readers of the derived tables (shared by figures, numbers and tables)
  geometry      grids, projections and the Beijing-3 scene registry (restricted-data figures)
  numbers       every number quoted in the text (numbers.json)
  tables        main Tables 1-2 and the computed supplementary Table S1
  supplement    supplementary-table cells recomputed from the result files
  verify        comparison of the outputs with the published values (reference/)
  analysis      runner of the mirrored v1.2 analysis code (analysis/; restricted inputs)
  legacy        roots that the mirrored analysis code imports (project data, output folder, font)
"""
import os as _os
import sys as _sys
import importlib.util as _ilu

__version__ = '2.0'


def _gdal_proj_data():
    """Use the PROJ and GDAL data bundled with the rasterio wheel.

    PROJ_LIB/PROJ_DATA/GDAL_DATA left over from another installation (e.g. a conda environment) make rasterio and
    pyproj read an incompatible proj.db ("lacks DATABASE.LAYOUT.VERSION"); the analysis was run with the bundled data.
    """
    spec = _ilu.find_spec('rasterio')
    if spec is None or not spec.submodule_search_locations:
        return
    base = list(spec.submodule_search_locations)[0]
    for key, sub in (('PROJ_LIB', 'proj_data'), ('PROJ_DATA', 'proj_data'), ('GDAL_DATA', 'gdal_data')):
        p = _os.path.join(base, sub)
        if _os.path.isdir(p):
            _os.environ[key] = p


_gdal_proj_data()

if _sys.platform == 'win32':
    # Python 3.12 on Windows queries system information through WMI in the platform module; on the authors' machine
    # this occasionally failed with E_OUTOFMEMORY and ended the process (exit code 10). The fallback is equivalent.
    import platform as _platform

    def _no_wmi(*_a, **_k):
        raise OSError('WMI disabled')
    _platform._wmi_query = _no_wmi
