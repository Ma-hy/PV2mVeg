"""Numbers of the text, main tables and supplementary-table cells from the repository data equal the published ones."""
from conftest import verify_in


def test_numbers_json_identical(outputs):
    rep = verify_in(outputs)
    assert rep['numbers']['status'] == 'identical', rep['numbers']


def test_main_tables_identical(outputs):
    rep = verify_in(outputs)
    for name, status in rep['tables'].items():
        assert status in ('identical', 'identical cells'), (name, status)


def test_supplementary_table_cells(outputs):
    rep = verify_in(outputs)
    s = rep['supplement']
    assert s['n_failed'] == 0 and s['checked'] > 250, s
