import os
import pymadx


def _fn(filename):
    return os.path.join(os.path.dirname(__file__), "test_input", filename)
