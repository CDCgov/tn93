#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#

"""
Title: tn93.py
Description: Implementation of Takamura-Nei distance calculation for pair of HIV sequences
Usage: Used by other software
Date Created: 2022-08-09 18:11
Last Modified: Tue 09 Aug 2022 06:19:56 PM EDT
Author: Reagan Kelly (ylb9@cdc.gov)
"""

import sys
import math
import logging


def main(args):
    seq1 = args[0]
    seq2 = args[1]
    match_mode = args[2]
    tn93_distance(seq1, seq2, match_mode)


def tn93_distance(seq1, seq2, match_mode):
    length = math.min(seq1.len, seq2.len)


if __name__ == "__main__":
    main(sys.argv[1:])
