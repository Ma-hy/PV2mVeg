"""Runner of the v1.2 analysis code (analysis/, mirrored from the research repository).

The analysis scripts run as in the archived run of the authors: one fresh Python process per script and, for the
supplementary diagnostics, one per stage. Several diagnostics resample with a module-level random generator
(analysis/eco/chrono/j_review.py), so their confidence intervals depend on the sequence of calls within a process;
one stage per process reproduces the archived values. The scripts need the project layout (restricted inputs) and
write to <output>/analysis/v12.

The child processes run as in the verification runs: UTF-8 mode on (PYTHONUTF8=1), because some of the analysis
code opens text files without naming an encoding, which would otherwise follow the locale; and single-threaded BLAS
(the scripts set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 1 only when they are unset; the runner
sets them, so that a different value in the caller's environment cannot change the order of floating-point sums).
"""
import os
import subprocess
import sys
import time

from pv2mveg.config import P, REPO

BG = REPO / 'analysis' / 'experiments' / 'bright_geometry'
OUT = P.out / 'analysis' / 'v12'


def run(script, *args):
    """Run analysis/experiments/bright_geometry/<script> with its arguments in a fresh process; stop on failure."""
    P.need_project('The analysis scripts (scripts/1x_*.py)')
    label = ' '.join([script, *args])
    t0 = time.time()
    print(f'== {label}', flush=True)
    env = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    r = subprocess.run([sys.executable, str(BG / script), *args], cwd=str(REPO), env=env)
    print(f'== {label}: exit {r.returncode} ({time.time() - t0:.0f} s)', flush=True)
    if r.returncode:
        sys.exit(r.returncode)
