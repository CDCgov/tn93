#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#

"""
Title: argerparser.py
Description: argument parser for tn93
Usage: imported by tn93
Date Created: 2022-09-21 12:15
Last Modified: Thu 22 Sep 2022 05:20:20 PM EDT
Author: Reagan Kelly (ylb9@cdc.gov)
"""

import argparse


def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        action="store",
        required=True,
        help="Path to the input fasta file",
    )
    parser.add_argument(
        "-m",
        "--match_mode",
        action="store",
        required=True,
        help="""
        How to handle ambiguities. This can be one of four options:
        average - Averages the possible nucleotide values for each ambiguity in a sequence;
        resolve - Tries to resolve ambiguities;
        skip - Ignores gaps and ambiguities;
        gapmm - Treats gaps in only one sequence as 'N's;
        """,
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        required=True,
        help="The name of the output file to create",
    )
    parser.add_argument(
        "-g",
        "--max_ambig_fraction",
        action="store",
        required=False,
        default=1.0,
        help="Sequences that have proportions of ambiguities lower than this value will be resolved, otherwise they will be averaged (Default: 1.0)",
    )
    parser.add_argument(
        "-c",
        "--show_counts",
        action="store_true",
        default=False,
        help="Show counts and other debugging information (Default: False)",
    )
    parser.add_argument(
        "-n",
        "--ignore_terminal_gaps",
        action="store_true",
        default=False,
        help="Should gaps at the beginning and end of a sequence be ignored (GAPMM only)? (Default: False)",
    )
    parser.add_argument(
        "-j",
        "--json_output",
        action="store_true",
        default=False,
        help="Should the output be in JSON format? (Default: False)",
    )
    return parser
