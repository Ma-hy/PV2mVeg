"""The figures of tier "repo" render from the repository data; their plotted-value tables equal those of the release,
and their pixels equal the published figures when the figure font of the publication is configured (PV2MVEG_FONT).

Rendering the 14 figures takes a few minutes; `pytest -m "not slow"` skips this test."""
import os

import pytest

from conftest import run, verify_in
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pv2mveg.catalog import FIGURES  # noqa: E402

REPO_FIGURES = [n for n, v in FIGURES.items() if v[1] == 'repo']


@pytest.mark.slow
def test_repository_figures(outputs):
    run('scripts/40_figures.py', outputs, *REPO_FIGURES)
    rep = verify_in(outputs, *REPO_FIGURES)
    for name in REPO_FIGURES:
        assert rep['figures'][name] != 'not rendered', name
    for f, status in rep['figdata'].items():
        assert status == 'identical', (f, status)
    if os.environ.get('PV2MVEG_FONT') and Path(os.environ['PV2MVEG_FONT']).exists():
        for name in REPO_FIGURES:
            assert rep['figures'][name] == 'pixel-identical', (name, rep['figures'][name])
