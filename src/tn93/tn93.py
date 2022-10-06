#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#

"""
Title: tn93.py
Description: Implementation of Tamura-Nei distance calculation for pair of HIV sequences
Usage: Used by other software
Date Created: 2022-08-09 18:11
Last Modified: Thu 06 Oct 2022 02:29:25 PM EDT
Author: Reagan Kelly (ylb9@cdc.gov)
"""

import copy
import string
import random
import csv
import json
import logging
import math
import re
import sys
from itertools import chain
from Bio import SeqIO
import numpy as np
import argparser


def main(args):
    parser = argparser.setup_parser()
    args = parser.parse_args(args)
    tn93 = TN93(
        verbose=args.verbose,
        ignore_gaps=args.ignore_terminal_gaps,
        json_output=args.json_output,
        max_ambig_fraction=args.max_ambig_fraction,
    )
    fasta_file = args.input_file
    fasta_sequences = [x for x in SeqIO.parse(fasta_file, format="fasta")]
    if (
        len(
            [
                x
                for x in fasta_sequences
                if len(str(x.seq)) != len(str(fasta_sequences[0].seq))
            ]
        )
        > 0
    ):
        logging.error(
            "Sequence lengths not identical - all sequences should be aligned to a common reference and should be the same length"
        )
        sys.exit(1)
    match_mode = args.match_mode
    final_distance = []
    for i in range(len(fasta_sequences) - 1):
        for j in range(i + 1, len(fasta_sequences)):
            dist_value = tn93.tn93_distance(
                fasta_sequences[i], fasta_sequences[j], match_mode
            )
            if not args.json_output:
                dist_holder = dist_value.split(",")
            final_distance += [dist_value if args.json_output else dist_holder]
    with open(args.output, "w") as output_file:
        if args.json_output:
            json.dump(final_distance, output_file)
        else:
            writer = csv.writer(output_file)
            writer.writerow(["ID1", "ID2", "Distance"])
            writer.writerows(final_distance)


class TN93(object):
    def __init__(
        self,
        verbose=False,
        ignore_gaps=False,
        json_output=False,
        max_ambig_fraction=1.0,
    ):
        (
            self.map_character,
            self.resolutions,
            self.resolutionsCount,
        ) = self.get_constants()
        self.verbose = verbose
        self.ignore_gaps = ignore_gaps
        self.json_output = json_output
        self.max_ambig_fraction = max_ambig_fraction
        self.first_nongap = False
        self.last_non_gap = False

    def tn93_distance(self, seq1, seq2, match_mode):
        try:  # If sequences are passed as SeqRecord objects
            pairwise_counts = self.get_counts(str(seq1.seq), str(seq2.seq), match_mode)
            seq_ids = [seq1.id, seq2.id]
        except AttributeError:  # If sequences are passed as strings
            pairwise_counts = self.get_counts(seq1, seq2, match_mode)
            seq_ids = [1, 2]
        nucleotide_frequency = self.get_nucleotide_frequency(pairwise_counts)
        dist = self.get_distance(pairwise_counts, nucleotide_frequency)
        if self.verbose > 0:
            logging.error(dist)
        if self.json_output:
            output = {"ID1": seq_ids[0], "ID2": seq_ids[1], "Distance": dist}
        else:
            output = f"{seq_ids[0]},{seq_ids[1]},{dist}"
        return output

    def get_distance(self, pairwise_counts, nucleotide_frequency):
        if self.verbose > 0:
            logging.error(pairwise_counts[0])
            logging.error(pairwise_counts[1])
            logging.error(pairwise_counts[2])
            logging.error(pairwise_counts[3])
            logging.error(nucleotide_frequency)
            logging.error(sum(nucleotide_frequency) / 2)
        dist = 0
        total_non_gap = 2 / sum(nucleotide_frequency)
        AG_counts = float(pairwise_counts[0][2] + pairwise_counts[2][0])
        AG = AG_counts * total_non_gap
        CT_counts = float(pairwise_counts[1][3] + pairwise_counts[3][1])
        CT = CT_counts * total_non_gap
        matching = (
            pairwise_counts[0][0]
            + pairwise_counts[1][1]
            + pairwise_counts[2][2]
            + pairwise_counts[3][3]
        ) * total_non_gap
        tv = 1 - (AG + CT + matching)
        if self.verbose > 0:
            logging.error(f"Initially: AG={AG} CT={CT} tv={tv}")
        if 0 in nucleotide_frequency:
            AG = 1 - 2 * (AG + CT) - tv
            CT = 1 - 2 * tv
            if AG > 0 and CT > 0:
                dist = -0.5 * math.log(AG) - 0.25 * math.log(CT)
            else:
                dist = 1.0
        else:
            auxd = 1 / sum(nucleotide_frequency)
            nucF = [0, 0, 0, 0]
            for j in range(4):
                nucF[j] = float(nucleotide_frequency[j]) * auxd
            fR = nucF[0] + nucF[2]
            fY = nucF[1] + nucF[3]
            K1 = 2 * nucF[0] * nucF[2] / fR
            K2 = 2 * nucF[1] * nucF[3] / fY
            K3 = 2 * (
                fR * fY - nucF[0] * nucF[2] * fY / fR - nucF[1] * nucF[3] * fR / fY
            )
            AG = 1 - AG / K1 - 0.5 * tv / fR
            CT = 1 - CT / K2 - 0.5 * tv / fY
            tv = 1 - 0.5 * tv / fY / fR
            dist = self.round(
                -K1 * math.log(AG) - K2 * math.log(CT) - K3 * math.log(tv)
            )

        if self.verbose > 0:
            logging.error(
                f"Finally: auxd={auxd} fR={fR} fY={fY} K1={K1} K2={K2} K3={K3} AG={AG} CT={CT} tv={tv}"
            )
            logging.error(f"dist={dist}")

        return self.round(dist) if dist > -0.0 else 0.0  # If it's really -0.0 return 0

    def get_nucleotide_frequency(self, pairwise_counts):
        nucleotide_frequency = [0, 0, 0, 0]
        for j in range(4):
            for k in range(4):
                nucleotide_frequency[j] += pairwise_counts[j][k]
                nucleotide_frequency[k] += pairwise_counts[j][k]
        return nucleotide_frequency

    def find_terminal_gaps(self, seq1, seq2):
        s1 = re.match(r"-*", seq1).end()
        s2 = re.match(r"-*", seq2).end()
        e1 = re.search(r"-*$", seq1).start()
        e2 = re.search(r"-*$", seq2).start()
        self.first_nongap = max(s1, s2)
        self.last_nongap = min(e1, e2) - 1
        if self.verbose > 0:
            logging.error(
                f"{self.first_nongap}-{seq1[self.first_nongap]},{seq2[self.first_nongap]} {self.last_nongap}-{seq1[self.last_nongap]},{seq2[self.last_nongap]}"
            )

    def ambig_fraction_too_high(self, seq):
        ambig_count = len([x for x in seq if x not in ["A", "C", "G", "T", "U", "-"]])
        if (ambig_count / len(seq)) > self.max_ambig_fraction:
            return True
        return False

    def get_counts(self, seq1, seq2, match_mode):
        self.find_terminal_gaps(seq1, seq2)
        if match_mode.upper() == "RESOLVE":
            if self.ambig_fraction_too_high(seq1) or self.ambig_fraction_too_high(seq2):
                pairwise_counts = self.get_counts_average(seq1, seq2)
            else:
                pairwise_counts = self.get_counts_resolve(seq1, seq2)
        elif match_mode.upper() == "AVERAGE":
            pairwise_counts = self.get_counts_average(seq1, seq2)
        elif match_mode.upper() == "GAPMM":
            pairwise_counts = self.get_counts_gapmm(seq1, seq2)
        elif match_mode.upper() == "SKIP":
            pairwise_counts = self.get_counts_skip(seq1, seq2)
        else:
            logging.error(f"Match mode {match_mode} is not recognized")
            sys.exit(1)
        return pairwise_counts

    def round(self, number):
        val = np.format_float_positional(
            number, precision=6, unique=False, fractional=False, trim="k"
        )
        return float(val)

    def get_counts_resolve(self, seq1, seq2):
        all_pairwise_count_arrays = {}
        length = min(len(seq1), len(seq2))
        pairwise_counts = [
            # A, C, G, T
            [0, 0, 0, 0],  # A
            [0, 0, 0, 0],  # C
            [0, 0, 0, 0],  # G
            [0, 0, 0, 0],  # T
        ]
        for p in range(0, length):
            nuc1 = self.map_character[ord(seq1[p])]
            nuc2 = self.map_character[ord(seq2[p])]
            if nuc1 < 4 and nuc2 < 4:  # If neither nucleotide is ambiguous
                if self.verbose and nuc1 == 3 and nuc2 == 1:
                    logging.error(f"{p} is TC")
                pairwise_counts[nuc1][nuc2] += 1
            else:
                if nuc1 == 17 and nuc2 == 17:
                    continue
                elif nuc1 < 4:  # nuc1 is resolved, nuc2 is ambiguous
                    if self.resolutionsCount[nuc2] > 0:
                        if self.resolutions[nuc2][
                            nuc1
                        ]:  # Resolve nuc2 to nuc1 if possible
                            pairwise_counts[nuc1][nuc1] += 1
                            continue
                        for j in range(4):
                            if self.resolutions[nuc2][j]:
                                pairwise_counts[nuc1][j] += self.resolutionsCount[nuc2]
                elif nuc2 < 4:  # nuc2 is resolved, nuc1 is ambiguous
                    if self.resolutionsCount[nuc1] > 0:
                        if self.resolutions[nuc1][
                            nuc2
                        ]:  # Resolve nuc1 to nuc2 if possible
                            pairwise_counts[nuc2][nuc2] += 1
                            continue
                        for j in range(4):
                            if self.resolutions[nuc1][j]:
                                pairwise_counts[j][nuc2] += self.resolutionsCount[nuc1]
                else:  # Both nuc1 and nuc2 are ambiguous
                    norm = self.resolutionsCount[nuc1] * self.resolutionsCount[nuc2]
                    if norm > 0.0:
                        matched_count = 0
                        positive_match = [False, False, False, False]
                        for j in range(4):
                            if self.resolutions[nuc1][j] and self.resolutions[nuc2][j]:
                                matched_count += 1
                                positive_match[j] = True
                        if matched_count > 0:
                            norm2 = 1 / matched_count
                            for j in range(4):
                                if positive_match[j]:
                                    pairwise_counts[j][j] += norm2
                            continue
                        for j in range(4):
                            if self.resolutions[nuc1][j]:
                                for k in range(4):
                                    if self.resolutions[nuc2][k]:
                                        pairwise_counts[j][k] += norm
            all_pairwise_count_arrays[str(p)] = copy.deepcopy(pairwise_counts)
        if self.verbose > 1:
            count_holder = []
            for k in all_pairwise_count_arrays.keys():
                count_line = list(chain.from_iterable(all_pairwise_count_arrays[k]))
                final_line = [int(k)] + count_line
                count_holder += [final_line]
            filename_slug = self.make_filename()
            with open(f"{filename_slug}_python_resolve_counts.csv", "w") as cf:
                writer = csv.writer(cf)
                writer.writerows(count_holder)
        return pairwise_counts

    def get_counts_average(self, seq1, seq2):
        all_pairwise_count_arrays = {}
        length = min(len(seq1), len(seq2))
        pairwise_counts = [
            # A, C, G, T
            [0, 0, 0, 0],  # A
            [0, 0, 0, 0],  # C
            [0, 0, 0, 0],  # G
            [0, 0, 0, 0],  # T
        ]
        for p in range(0, length):
            nuc1 = self.map_character[ord(seq1[p])]
            nuc2 = self.map_character[ord(seq2[p])]
            if nuc1 < 4 and nuc2 < 4:  # If neither nucleotide is ambiguous
                pairwise_counts[nuc1][nuc2] += 1
            else:
                if nuc1 == 17 or nuc2 == 17:  # One or both sequences have a gap
                    continue
                elif nuc1 < 4:  # nuc1 is resolved, nuc2 is ambiguous
                    if self.resolutionsCount[nuc2] > 0:
                        for j in range(4):
                            if self.resolutions[nuc2][j]:
                                pairwise_counts[nuc1][j] += self.resolutionsCount[nuc2]
                elif nuc2 < 4:  # nuc2 is resolved, nuc1 is ambiguous
                    if self.resolutionsCount[nuc1] > 0:
                        for j in range(4):
                            if self.resolutions[nuc1][j]:
                                pairwise_counts[j][nuc2] += self.resolutionsCount[nuc1]
                else:  # nuc1 and nuc2 are ambiguous
                    norm = self.resolutionsCount[nuc1] * self.resolutionsCount[nuc2]
                    if norm > 0.0:
                        for j in range(4):
                            if self.resolutions[nuc1][j]:
                                for k in range(4):
                                    if self.resolutions[nuc2][k]:
                                        pairwise_counts[j][k] += norm
            all_pairwise_count_arrays[str(p)] = copy.deepcopy(pairwise_counts)
        if self.verbose > 1:
            count_holder = []
            for k in all_pairwise_count_arrays.keys():
                count_line = list(chain.from_iterable(all_pairwise_count_arrays[k]))
                final_line = [int(k)] + count_line
                count_holder += [final_line]
            filename_slug = self.make_filename()
            with open(f"{filename_slug}_python_avg_counts.csv", "w") as cf:
                writer = csv.writer(cf)
                writer.writerows(count_holder)
        return pairwise_counts

    def get_counts_gapmm(self, seq1, seq2):
        all_pairwise_count_arrays = {}
        length = min(len(seq1), len(seq2))
        pairwise_counts = [
            # A, C, G, T
            [0, 0, 0, 0],  # A
            [0, 0, 0, 0],  # C
            [0, 0, 0, 0],  # G
            [0, 0, 0, 0],  # T
        ]
        start = self.first_nongap if self.ignore_gaps else 0
        end = self.last_nongap if self.ignore_gaps else length
        for p in range(start, end):
            nuc1 = self.map_character[ord(seq1[p])]
            nuc2 = self.map_character[ord(seq2[p])]
            if nuc1 < 4 and nuc2 < 4:  # If neither nucleotide is ambiguous
                pairwise_counts[nuc1][nuc2] += 1
            else:
                if nuc1 == 17 or nuc2 == 17:  # One or both sequences have a gap
                    if nuc1 == 17 and nuc2 == 17:  # Both sequences have a gap
                        continue
                    else:
                        if nuc1 == 17:
                            nuc1 = 15
                        else:
                            nuc2 = 15
                if nuc1 < 4:  # nuc1 is resolved, nuc2 ambiguous
                    if self.resolutionsCount[nuc2] > 0:
                        for j in range(4):
                            if self.resolutions[nuc2][j]:
                                pairwise_counts[nuc1][j] += self.resolutionsCount[nuc2]
                else:
                    if nuc2 < 4:  # nuc2 is resolved, nuc1 is ambiguous
                        if self.resolutionsCount[nuc1] > 0:
                            for j in range(4):
                                if self.resolutions[nuc1][j]:
                                    pairwise_counts[j][nuc2] += self.resolutionsCount[
                                        nuc1
                                    ]
                    else:  # Both nuc1 and nuc2 are ambiguous
                        norm = self.resolutionsCount[nuc1] * self.resolutionsCount[nuc2]
                        if norm > 0.0:
                            for j in range(4):
                                if self.resolutions[nuc1][j]:
                                    for k in range(4):
                                        if self.resolutions[nuc2][k]:
                                            pairwise_counts[j][k] += norm
            all_pairwise_count_arrays[str(p)] = copy.deepcopy(pairwise_counts)
        if self.verbose > 1:
            count_holder = []
            for k in all_pairwise_count_arrays.keys():
                count_line = list(chain.from_iterable(all_pairwise_count_arrays[k]))
                final_line = [int(k)] + count_line
                count_holder += [final_line]
            filename_slug = self.make_filename()
            with open(f"{filename_slug}_python_gapmm_counts.csv", "w") as cf:
                writer = csv.writer(cf)
                writer.writerows(count_holder)
        return pairwise_counts

    def get_counts_skip(self, seq1, seq2):
        all_pairwise_count_arrays = (
            {}
        )  # This is for debugging the difference between impl
        length = min(len(seq1), len(seq2))
        pairwise_counts = [
            # A, C, G, T
            [0, 0, 0, 0],  # A
            [0, 0, 0, 0],  # C
            [0, 0, 0, 0],  # G
            [0, 0, 0, 0],  # T
        ]
        for p in range(0, length):
            nuc1 = self.map_character[ord(seq1[p])]
            nuc2 = self.map_character[ord(seq2[p])]
            if nuc1 < 4 and nuc2 < 4:  # If neither nucleotide is ambiguous
                pairwise_counts[nuc1][nuc2] += 1
            all_pairwise_count_arrays[str(p)] = copy.deepcopy(pairwise_counts)
        if self.verbose > 1:
            count_holder = []
            for k in all_pairwise_count_arrays.keys():
                count_holder += [p] + [
                    list(chain.from_iterable(all_pairwise_count_arrays[k]))
                ]
            filename_slug = self.make_filename()
            with open(f"{filename_slug}_python_skip_counts.csv", "w") as cf:
                writer = csv.writer(cf)
                writer.writerows(count_holder)

        return pairwise_counts

    def make_filename(self, N=8):
        return "".join(random.choices(string.ascii_uppercase + string.digits, k=N))

    def get_constants(self):
        map_character = [16 for x in range(256)]
        map_character[45] = 17  # GAP
        map_character[65] = 0  # A
        map_character[66] = 11  # B
        map_character[67] = 1  # C
        map_character[68] = 12  # D
        map_character[71] = 2  # G
        map_character[72] = 13  # H
        map_character[75] = 9  # K
        map_character[77] = 10  # M
        map_character[78] = 15  # N
        map_character[82] = 5  # R
        map_character[83] = 7  # S
        map_character[84] = 3  # T
        map_character[85] = 4  # U
        map_character[86] = 14  # V
        map_character[87] = 8  # W
        map_character[89] = 6  # Y
        map_character[97] = 0  # a
        map_character[98] = 11  # b
        map_character[99] = 1  # c
        map_character[100] = 12  # d
        map_character[103] = 2  # g
        map_character[104] = 13  # h
        map_character[107] = 9  # k
        map_character[109] = 10  # m
        map_character[110] = 15  # n
        map_character[114] = 5  # r
        map_character[115] = 7  # s
        map_character[116] = 3  # t
        map_character[117] = 4  # u
        map_character[118] = 14  # v
        map_character[119] = 8  # w
        map_character[121] = 6  # y

        resolutions = [
            # A,C,G,T
            [1, 0, 0, 0],  # A             -> A (0) (Adenine)
            [0, 1, 0, 0],  # C             -> C (1) (Cytosine)
            [0, 0, 1, 0],  # G             -> G (2) (Guanine)
            [0, 0, 0, 1],  # T             -> T (3) (Thymine)
            [0, 0, 0, 1],  # T             -> U (4) (Uracil)
            [1, 0, 1, 0],  # A | G         -> R (5) (Either Purine)
            [0, 1, 0, 1],  # C | T         -> Y (6) (Either Pyrimidine)
            [0, 1, 1, 0],  # C | G         -> S (7)
            [1, 0, 0, 1],  # A | T         -> W (8)
            [0, 0, 1, 1],  # G | T         -> K (9)
            [1, 1, 0, 0],  # A | C         -> M (10)
            [0, 1, 1, 1],  # C | G | T     -> B (11) (Not Adenine)
            [1, 0, 1, 1],  # A | G | T     -> D (12) (Not Cytosine)
            [1, 1, 0, 1],  # A | C | T     -> H (13) (Not Guanine)
            [1, 1, 1, 0],  # A | C | G     -> V (14) (Not Thymine)
            [1, 1, 1, 1],  # A | C | G | T -> N (15)
            [1, 1, 1, 1],  # A | C | G | T -> ? (16)
            [0, 0, 0, 0],  # GAP
        ]

        resolutionsCount = [
            1.0,  # A
            1.0,  # C
            1.0,  # G
            1.0,  # T
            1.0,  # U
            1.0 / 2.0,  # R
            1.0 / 2.0,  # Y
            1.0 / 2.0,  # S
            1.0 / 2.0,  # W
            1.0 / 2.0,  # K
            1.0 / 2.0,  # M
            1.0 / 3.0,  # B
            1.0 / 3.0,  # D
            1.0 / 3.0,  # H
            1.0 / 3.0,  # V
            1.0 / 4.0,  # N
            1.0 / 4.0,  # ?
            0.0,  # GAP
        ]
        return map_character, resolutions, resolutionsCount


if __name__ == "__main__":
    main(sys.argv[1:])
