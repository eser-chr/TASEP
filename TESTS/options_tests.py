import unittest
from lib.parameters import Parameters

class TestOptions(unittest.TestCase):

    def test_valid_types(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.are_valid_types())

    def test_invalid_types(self):
        options = Parameters(10.5, "1000", "invalid", 0.05, 0.01, 0.2, 0.2)
        self.assertFalse(options.are_valid_types())

    def test_mathematically_valid_ranges(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.mathematically_valid_ranges())

    def test_mathematically_invalid_ranges(self):
        options = Parameters(-10, 1000, 0.1, -0.05, 0.01, 0, 0)
        self.assertFalse(options.mathematically_valid_ranges())

    def test_q_and_kq_invalid(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0)
        self.assertFalse(options.is_q_and_kq_valid())

    def test_q_and_kq_valid_a(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 1, 0)
        self.assertTrue(options.is_q_and_kq_valid())

    def test_q_and_kq_valid_b(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 2, 10)
        self.assertTrue(options.is_q_and_kq_valid())

    def test_no_debug_check(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.no_debug_check())

    def test_debug_check(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.debug_check())

    def test_check_with_debug(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.check(debug=True))

    def test_check_without_debug(self):
        options = Parameters(10, 1000, 0.1, 0.05, 0.01, 0.2, 0.2)
        self.assertTrue(options.check())
