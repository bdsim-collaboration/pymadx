import os
import pymadx
import pytest

_skip_cpymad_tests = False
try:
    import cpymad
    import cpymad.madx
except ImportError:
    _skip_cpymad_tests = True

if _skip_cpymad_tests:
    pytestmark = pytest.mark.xfail(run=True, reason="requires cpymad")


def _fn(filename):
    return os.path.join(os.path.dirname(__file__), filename)


def test_Tfs_from_cpymad_table():
    m = cpymad.madx.Madx(stdout=None)
    m.call(file=_fn("test_input/h6.madx"))
    tt = m.table['twiss']
    tfs = pymadx.Data.Tfs(tt)
    print(tfs)
