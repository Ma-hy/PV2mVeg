"""Shared fixtures: run the scripts of the repository in a fresh output folder with the repository data."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def run(script, out, *args, layout='repo'):
    env = dict(os.environ, PV2MVEG_LAYOUT=layout, PV2MVEG_OUTPUT=str(out), PYTHONIOENCODING='utf-8')
    if layout == 'repo':
        env['PV2MVEG_DATA_ROOT'] = str(REPO / 'data' / 'derived')
    r = subprocess.run([sys.executable, str(REPO / script), *args], cwd=str(REPO), env=env, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
    assert r.returncode == 0, f'{script} failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}'
    return r


@pytest.fixture(scope='session')
def outputs(tmp_path_factory):
    """Tables and numbers of the repository data, written once per test session."""
    out = tmp_path_factory.mktemp('output')
    run('scripts/30_tables.py', out)
    run('scripts/31_numbers.py', out)
    return out


def verify_in(out, *names):
    """pv2mveg.verify.report for an output folder (in a subprocess, so that the configuration is read afresh)."""
    code = ('import json, sys; sys.path.insert(0, r"%s"); from pv2mveg import verify; '
            'print(json.dumps(verify.report(%r or None)))') % (REPO / 'src', list(names))
    env = dict(os.environ, PV2MVEG_LAYOUT='repo', PV2MVEG_OUTPUT=str(out), PV2MVEG_DATA_ROOT=str(REPO / 'data' / 'derived'),
               PYTHONIOENCODING='utf-8')
    r = subprocess.run([sys.executable, '-c', code], cwd=str(REPO), env=env, capture_output=True, text=True, encoding='utf-8')
    assert r.returncode == 0, r.stderr[-3000:]
    return json.loads(r.stdout.strip().split('\n')[-1])
