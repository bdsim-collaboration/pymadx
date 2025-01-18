import os
import pymadx


def _fn(filename):
    return os.path.join(os.path.dirname(__file__), "test_input", filename)


def test_data_tfs_load():
    d = pymadx.Data.Tfs(_fn("h6-simple.tfs"))
    assert len(d) == 226


def test_data_tfs_survey_load():
    d = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    assert len(d) == 226


def test_data_tfs_load_gz():
    d = pymadx.Data.Tfs(_fn("h6-positive-120gev-fm.tfs.gz"))
    assert len(d) == 226