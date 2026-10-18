"""Integrity of the repository: derived data match their manifest, no local paths or user names in tracked text."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_derived_data_match_manifest():
    man = json.loads((REPO / 'data/derived/MANIFEST.json').read_text(encoding='utf-8'))['files']
    for rel, rec in man.items():
        p = REPO / 'data/derived' / rel
        assert p.exists(), rel
        assert hashlib.sha256(p.read_bytes()).hexdigest() == rec['sha256'], rel


def test_no_identifying_columns_in_derived_data():
    import pandas as pd
    banned = {'x', 'y', 'lon', 'lat', 'row', 'col', 'T', 'T_hand', 'hand_year', 'det_year', 'station_label', 'lab_year', 'onset_year', 'array_id'}
    for p in (REPO / 'data/derived').rglob('*'):
        if p.suffix == '.csv':
            cols = set(pd.read_csv(p, nrows=0).columns)
        elif p.suffix == '.parquet':
            cols = set(pd.read_parquet(p).columns)
        else:
            continue
        if p.name == 'balance.csv':      # cohort-level table: T is the construction-year cohort of a group of units
            cols -= {'T'}
        assert not cols & banned, (p.relative_to(REPO), cols & banned)


def tracked_text_files():
    try:
        out = subprocess.run(['git', 'ls-files'], cwd=str(REPO), capture_output=True, text=True, check=True).stdout.split()
    except Exception:
        out = [str(p.relative_to(REPO)) for p in REPO.rglob('*') if p.is_file() and '.git' not in p.parts]
    return [REPO / f for f in out if Path(f).suffix in ('.py', '.md', '.toml', '.txt', '.json', '.yml', '.cfg') and 'derived' not in Path(f).parts]


def test_no_local_paths():
    pat = re.compile(r'[A-Za-z]:\\\\?(?:RUN|Users)|TalatanAna\\|C:/RUN|AppData|\bDELL\b')
    bad = []
    for p in tracked_text_files():
        if p.resolve() == Path(__file__).resolve():      # the pattern above
            continue
        for i, line in enumerate(p.read_text(encoding='utf-8', errors='replace').split('\n'), 1):
            if pat.search(line):
                bad.append(f'{p.relative_to(REPO)}:{i}: {line.strip()[:120]}')
    assert not bad, '\n'.join(bad)
