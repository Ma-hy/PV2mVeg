"""The mirrored analysis code (analysis/): every module parses, takes its roots from pv2mveg.legacy, and has English
comments and docstrings; the analysis wrappers stop with a clear message without the restricted inputs.

Running the analysis itself needs the project data (scripts/10_units.py ... 13_secondary.py, docs/VERIFICATION.md)."""
import ast
import io
import os
import re
import subprocess
import sys
import tokenize
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / 'analysis'
CJK = re.compile(r'[　-〿㐀-鿿＀-￯]')


def modules():
    return sorted(CODE.rglob('*.py'))


def docstring_lines(tree):
    spans = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                spans.add((first.value.lineno, first.value.end_lineno))
    return spans


def test_modules_parse_and_take_their_roots_from_the_configuration():
    mods = modules()
    assert len(mods) == 48
    for p in mods:
        src = p.read_text(encoding='utf-8')
        ast.parse(src)
        assert 'from pv2mveg.legacy import PROJECT, CODE, V12_OUT, FONT_FILE, SELECTED_STATIONS' in src, p


def test_comments_and_docstrings_in_english():
    for p in modules():
        src = p.read_text(encoding='utf-8')
        docs = docstring_lines(ast.parse(src))
        for t in tokenize.generate_tokens(io.StringIO(src).readline):
            if CJK.search(t.string):
                assert t.type != tokenize.COMMENT, (p.relative_to(REPO), t.start)
                assert not (t.type == tokenize.STRING and (t.start[0], t.end[0]) in docs), (p.relative_to(REPO), t.start)


def test_analysis_needs_the_project_layout(tmp_path):
    env = dict(os.environ, PV2MVEG_LAYOUT='repo', PV2MVEG_OUTPUT=str(tmp_path), PYTHONIOENCODING='utf-8')
    r = subprocess.run([sys.executable, str(REPO / 'scripts' / '10_units.py')], cwd=str(REPO), env=env,
                       capture_output=True, text=True, encoding='utf-8', errors='replace')
    assert r.returncode != 0
    assert 'needs restricted inputs' in r.stdout + r.stderr
    assert not (tmp_path / 'analysis').exists()
