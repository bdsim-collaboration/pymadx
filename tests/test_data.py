import os
import numpy as np
import pymadx


abs_difference = 1e-6

def within_tolerance(val, ref, tolerance=abs_difference):
    v = np.abs(val - ref) <= tolerance
    return np.all(v)


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


def test_roto_translation():
    r = pymadx.Data.RotoTranslation2D(np.pi/6, np.array([3, 4]))
    w, l = 3, 5
    x0, y0 = 10, 1
    points = np.array([[x0, y0+0.5*w], [x0, y0-0.5*w], [x0-l, y0-0.5*w], [x0-l, y0+0.5*w]])
    rpoints = r * points
    ref = np.array([[10.41025404, 11.16506351],
                    [11.91025404,  8.5669873],
                    [7.58012702,   6.0669873],
                    [6.08012702,   8.66506351]])
    assert within_tolerance(rpoints, ref)


def test_roto_translation_rot_origin():
    x0, y0 = 10, 1
    r = pymadx.Data.RotoTranslation2D(np.pi/6, np.array([3, 4]), rotationOrigin=np.array([x0, y0]))
    w, l = 3, 5
    points = np.array([[x0, y0+0.5*w], [x0, y0-0.5*w], [x0-l, y0-0.5*w], [x0-l, y0+0.5*w]])
    rpoints = r * points
    ref = np.array([[12.25,        6.29903811],
                    [13.75,        3.70096189],
                    [9.41987298,   1.20096189],
                    [7.91987298,   3.79903811]])
    assert within_tolerance(rpoints, ref)

def test_roto_tranlsate_from_survey():
    d = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    r = d.GetRotoTranslationFromElementZX("MBXHC.X0410132")
    localPos = np.array([[0, 0]])
    globalPos = r * localPos
    ref = np.array([d['MBXHC.X0410132']['Z'], d['MBXHC.X0410132']['X']])
    assert within_tolerance(globalPos, ref)


def test_roto_tranlsate_from_survey_inverse():
    d = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    r = d.GetRotoTranslationFromElementZX("MBXHC.X0410132")
    rinv = r.Inverse()
    ele = d['MBXHC.X0410132']
    globalPos = np.array([[ele['Z'], ele['X']]])
    localPos = rinv * globalPos
    ref = np.array([[0, 0]])
    assert within_tolerance(localPos, ref)

test_roto_translation()
test_roto_translation_rot_origin()
test_roto_tranlsate_from_survey()
test_roto_tranlsate_from_survey_inverse()