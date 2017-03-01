import unittest
import pandas as pd

from extract_coverage_caveat import extract_qc_call


class TestMethods(unittest.TestCase):
    def test_extract_coverage_caveat_ok1(self):
        data = [[110, 3, 548, 0.04],
                [111, 7, 548, 0.21],
                [112, 7, 548, 0.25],
                [113, 10, 548, 0.25],
                [120, 7, 548, 0.25]]
        test_table = pd.DataFrame(data)

        self.assertEquals(extract_qc_call(test_table, 0.95, 100, 0.95, 50), "OK")

    def test_extract_coverage_caveat_ok2(self):
        data = [[40, 3, 548, 0.04],
                [101, 7, 548, 0.21],
                [112, 7, 548, 0.25],
                [113, 10, 548, 0.25],
                [120, 7, 548, 0.25]]
        test_table = pd.DataFrame(data)

        self.assertEquals(extract_qc_call(test_table, 0.95, 100, 0.95, 50), "OK")

    def test_extract_coverage_caveat_warn1(self):
        data = [[40, 3, 548, 0.04],
                [80, 7, 548, 0.21],
                [112, 7, 548, 0.25],
                [113, 10, 548, 0.25],
                [120, 7, 548, 0.25]]
        test_table = pd.DataFrame(data)

        self.assertEquals(extract_qc_call(test_table, 0.95, 100, 0.95, 50), "WARN")

    def test_extract_coverage_caveat_fail1(self):
        data = [[40, 3, 548, 0.04],
                [45, 7, 548, 0.21],
                [112, 7, 548, 0.25],
                [113, 10, 548, 0.25],
                [120, 7, 548, 0.25]]
        test_table = pd.DataFrame(data)

        self.assertEquals(extract_qc_call(test_table, 0.95, 100, 0.95, 50), "FAIL")
