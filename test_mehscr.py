"""Fast regression tests for the calculation and context helpers."""

import random

import numpy as np
import pandas as pd

import MeHscr as meh


class FakeFasta:
    def __init__(self, sequence):
        self.sequence = sequence

    def get_reference_length(self, chrom):
        return len(self.sequence)

    def fetch(self, chrom, start, end):
        return self.sequence[start:end]


def test_patterns_and_hamming_distance():
    patterns = meh._all_patterns(4)
    distance = meh.PattoDis(pd.DataFrame(patterns), 1)
    assert patterns.shape == (16, 4)
    assert distance.shape == (16, 16)
    assert distance[0, 15] == 4
    assert np.array_equal(distance, distance.T)


def test_uniform_patterns_pairwise_score():
    patterns = meh._all_patterns(4)
    distance = meh.PattoDis(pd.DataFrame(patterns), 1)
    result = meh.MeHperwindow(
        pd.DataFrame(patterns), 10, 30, "chr1", distance, 4, False
    )
    assert result.loc[0, "MeH"] == 22.62742
    assert list(result.columns) == meh.RESULT_COLUMNS


def test_optional_pattern_counts_use_legacy_bit_order():
    patterns = pd.DataFrame([[0, 0, 0], [1, 0, 0], [1, 0, 0]])
    distance = meh.PattoDis(pd.DataFrame(meh._all_patterns(3)), 1)
    _, composition = meh.MeHperwindow(
        patterns, 1, 2, "chr1", distance, 3, True, MeH=1
    )
    assert composition.loc[0, "p01"] == 1
    assert composition.loc[0, "p02"] == 2


def test_imputation_can_be_reproduced():
    window = np.array([[0, 0, 0], [1, 1, 1], [0, 1, 0], [1, 0, 1], [0, np.nan, 1]])
    first = meh.impute(window, 3, random.Random(7))
    second = meh.impute(window, 3, random.Random(7))
    assert np.array_equal(first, second)
    assert not np.isnan(first).any()


def test_context_matching_and_reference_boundaries():
    fasta = FakeFasta("AACGACTGACATG")
    assert meh._context_matches(fasta, "chr1", 2, "CG", "f")
    assert meh._context_matches(fasta, "chr1", 5, "CHG", "f")
    assert meh._context_matches(fasta, "chr1", 9, "CHH", "f")
    assert not meh._context_matches(fasta, "chr1", len(fasta.sequence), "CHG", "r")


def test_read_thresholds():
    complete = meh._all_patterns(3)
    assert meh.enough_reads(complete, 3, complete=True)
    assert not meh.enough_reads(complete[:7], 3, complete=True)


def test_progress_message():
    message = meh._progress_message("CG", "sample1", 25, "Chr1", 1234, 56)
    assert "[#####---------------]" in message
    assert "25%" in message
    assert "Chr1:1,234" in message
    assert "windows=56" in message
