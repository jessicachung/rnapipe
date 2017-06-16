import unittest
from unittest.mock import patch, mock_open
from rnapipe.samples import *

class TestParseCSV(unittest.TestCase):
    '''Unit tests for parsing CSV files'''
    def test_remove_comments(self):
        data = "#\n#\nsample-A,c1\nsample-B,c2 # comment\n"
        expect = ["sample-A,c1", "sample-B,c2"]
        with patch('builtins.open', mock_open(read_data=data)) as m:
            contents = read_csv("foo")
        self.assertEqual(contents, expect)

    def test_samples_csv(self):
        data = "sample-A,c1,x\nsample-B,c2,x\n"
        expect = [["sample-A", "c1", "x"], ["sample-B", "c2", "x"]]
        with patch('builtins.open', mock_open(read_data=data)) as m:
            contents = parse_samples_csv("foo")
        self.assertEqual(contents, expect) 

    def test_fail_bad_csv_1(self):
        data = "sample-A,c1,\nsample-B,c2\n"
        with self.assertRaises(SystemExit):
            with patch('builtins.open', mock_open(read_data=data)) as m:
                contents = parse_samples_csv("foo")

    def test_fail_bad_csv_2(self):
        data = "sample-A,c1\nsample-B\n"
        with self.assertRaises(SystemExit):
            with patch('builtins.open', mock_open(read_data=data)) as m:
                contents = parse_samples_csv("foo")


