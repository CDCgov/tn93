#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#

"""
Title:test_tn93.py
Description:
Usage:
Date Created: 2022-08-15 11:43
Last Modified: Mon 15 Aug 2022 12:08:11 PM EDT
Author: Reagan Kelly (ylb9@cdc.gov)
"""

import unittest
from unittest.mock import Mock
import tn93


class TestTN93(unittest.TestCase):
    def setUp(self):
        self.unambig_seq1 = "AACTCATG"
        self.unambig_seq2 = "ACATACCT"
        self.ambig_seq1 = ""
        self.ambig_seq2 = ""

    def test_skip_unambig(self):
        pass

    def test_gapmm_unambig(self):
        pass

    def test_average_unambig(self):
        pass

    def test_resolve_unambig(self):
        pass

    def test_skip_ambig(self):
        pass

    def test_gapmm_ambig(self):
        pass

    def test_average_ambig(self):
        pass

    def test_resolve_ambig(self):
        pass

    def test_mode_skip(self):
        pass

    def test_mode_gapmm(self):
        pass

    def test_mode_average(self):
        pass

    def test_mode_resolve(self):
        pass

    def test_mode_unknown(self):
        pass

    def test_get_nucleotide_frequency(self):
        pass

    def test_get_distance(self):
        pass


if __name__ == "__main__":
    unittest.main()
