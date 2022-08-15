#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#

"""
Title:test_tn93.py
Description:
Usage:
Date Created: 2022-08-15 11:43
Last Modified: Mon 15 Aug 2022 04:30:23 PM EDT
Author: Reagan Kelly (ylb9@cdc.gov)
"""

import unittest
import json
from unittest.mock import Mock
import tn93
from Bio import SeqIO


class TestTN93(unittest.TestCase):
    def setUp(self):
        self.test_seqs = {
            x.id: x for x in SeqIO.parse("test/test_sequences.fasta", format="fasta")
        }
        self.tn93 = tn93.TN93()
        with open("test/expected_counts.json", "r") as json_file:
            self.counts = json.load(json_file)
        with open("test/expected_distances.json", "r") as json_file:
            self.distances = json.load(json_file)

    def test_skip_unambig(self):
        pairwise_counts = self.tn93.get_counts_skip(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["unambig_skip"])

    def test_gapmm_unambig(self):
        pairwise_counts = self.tn93.get_counts_gapmm(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["unambig_gapmm"])

    def test_average_unambig(self):
        pairwise_counts = self.tn93.get_counts_resolve(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["unambig_average"])

    def test_resolve_unambig(self):
        pairwise_counts = self.tn93.get_counts_average(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["unambig_resolve"])

    def test_skip_ambig(self):
        pairwise_counts = self.tn93.get_counts_skip(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["ambig_skip"])
        pass

    def test_gapmm_ambig(self):
        pairwise_counts = self.tn93.get_counts_gapmm(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["ambig_gapmm"])
        pass

    def test_average_ambig(self):
        pairwise_counts = self.tn93.get_counts_average(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["ambig_average"])
        pass

    def test_resolve_ambig(self):
        pairwise_counts = self.tn93.get_counts_resolve(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        self.assertEqual(pairwise_counts, self.counts["ambig_resolve"])
        pass

    def test_mode_skip(self):
        self.tn93.get_counts_skip = Mock()
        mode_value = "SKIP"
        self.tn93.get_counts(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"], mode_value
        )
        self.tn93.get_counts_skip.assert_called_once()

    def test_mode_gapmm(self):
        self.tn93.get_counts_gapmm = Mock()
        mode_value = "GAPMM"
        self.tn93.get_counts(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"], mode_value
        )
        self.tn93.get_counts_gapmm.assert_called_once()
        pass

    def test_mode_average(self):
        self.tn93.get_counts_average = Mock()
        mode_value = "AVERAGE"
        self.tn93.get_counts(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"], mode_value
        )
        self.tn93.get_counts_average.assert_called_once()
        pass

    def test_mode_resolve(self):
        self.tn93.get_counts_resolve = Mock()
        mode_value = "RESOLVE"
        self.tn93.get_counts(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"], mode_value
        )
        self.tn93.get_counts_resolve.assert_called_once()
        pass

    def test_mode_unknown(self):
        mode_value = "UNKNOWN"
        with self.assertRaises(SystemExit) as cm:
            self.tn93.get_counts(
                self.test_seqs["unambig_1"], self.test_seqs["unambig_2"], mode_value
            )
            self.assertEqual(cm.exception_code, 1)

    def test_get_nucleotide_frequency(self):
        nuc_freq = self.tn93.get_nucleotide_frequency(self.counts["unambig_resolve"])
        self.assertEqual(nuc_freq, [740.5, 262.0, 377.0, 438.5])

    def test_get_distance_ambig_skip(self):
        pairwise = self.tn93.get_counts_skip(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["ambig_skip"])

    def test_get_distance_ambig_gapmm(self):
        pairwise = self.tn93.get_counts_gapmm(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["ambig_gapmm"])

    def test_get_distance_ambig_average(self):
        pairwise = self.tn93.get_counts_average(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["ambig_average"])

    def test_get_distance_ambig_resolve(self):
        pairwise = self.tn93.get_counts_resolve(
            self.test_seqs["ambig_1"], self.test_seqs["ambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["ambig_resolve"])

    def test_get_distance_unambig_skip(self):
        pairwise = self.tn93.get_counts_skip(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["unambig_skip"])

    def test_get_distance_unambig_gapmm(self):
        pairwise = self.tn93.get_counts_gapmm(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["unambig_gapmm"])

    def test_get_distance_unambig_average(self):
        pairwise = self.tn93.get_counts_average(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["unambig_average"])

    def test_get_distance_unambig_resolve(self):
        pairwise = self.tn93.get_counts_resolve(
            self.test_seqs["unambig_1"], self.test_seqs["unambig_2"]
        )
        nuc_freq = self.tn93.get_nucleotide_frequency(pairwise)
        distance = self.tn93.get_distance(pairwise, nuc_freq)
        self.assertEqual(distance, self.distances["unambig_resolve"])


if __name__ == "__main__":
    unittest.main()
