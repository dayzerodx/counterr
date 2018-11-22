"""
This file contains dependencies and static data structures
accessed throughout the command line tool.
"""
from functools import partial
import pandas as pd
import pysam
import numpy as np
import os
from time import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gzip
import argparse
import csv
import seaborn as sns
import sys
import math
import warnings
np.random.seed(42)
warnings.filterwarnings("ignore", category=UserWarning) # Intended to ignore pesky Matplotlib kwarg warning.
if sys.version_info.major > 2:
    xrange = range

N_bases = 4
base2int = {"A" : 0, "C" : 1, "G" : 2, "T" : 3, "-": 4}
int2base = {0: "A", 1: "C", 2: "G", 3: "T"}
bases = {"A", "C", "G", "T"}
bases_ext = {"A", "C", "G", "T", "-"}
bases_list = ["A", "C", "G", "T"]
bases_ext_list = ["A", "C", "G", "T", "-"]
base2complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
np_int_type = np.int32


def reverse_complement(seq):
    """
    Given a sequence, return it's reverse complements.
    If the character is not in "bases", do nothing.
    """
    len_seq = len(seq)
    seq_out = []
    for i in xrange(len_seq-1, -1, -1):
        c = seq[i]
        if c in bases:
            seq_out.append(base2complement[c])
        else:
            seq_out.append(c)

    return "".join(seq_out)

def kmer_to_idx(kmer):
    """
    Turns a string of nucleotides to 64 bit integers with the coversion {00: A, 01: C, 10: G, 11: T}.
    """
    assert len(kmer) <= 30
    idx = 0
    for c in kmer:
        idx = idx << 2
        idx += base2int[c]
    return idx

def idx_to_kmer(K, idx):
    """
    Inverse operation of the string conversion kmer_to_idx.
    """
    letters = []
    for i in xrange(K):
        letters.append(int2base[idx & 3])
        idx = idx >> 2
    return ''.join(letters[::-1])

def merge_dict(dict1, dict2):
    rtrn_dict = {}
    for key in dict1.keys():
        rtrn_dict[key] = dict1[key] + dict2[key]
    return rtrn_dict

def match_rate2phredQ(rate):
    rate = np.asarray(rate)
    with np.errstate(divide='ignore'):
        return -10 * np.log10((1-rate))

def count_dict2bar_chart(count_dict):
    length = []
    count = []
    for key in count_dict.keys():
        length.append(key)
        count.append(count_dict[key])

    return length, count

def base_counts2np_arr(base_counts):
    arr = np.zeros(4, dtype=float)
    arr[0] = base_counts["A"]
    arr[1] = base_counts["C"]
    arr[2] = base_counts["G"]
    arr[3] = base_counts["T"]

    return arr

def compute_mean_med_std(Q):
    """
    Given an array return its mean, med, std
    """
    return np.mean(Q), np.median(Q), np.std(Q)

def has_hp(kmer):
    """
    Return True if the kmer has a homopolymer.
    """
    idx = 0
    while idx < (len(kmer)-1):
        if kmer[idx] == kmer[idx+1]:
            n_repeat = 1
            while (idx < (len(kmer)-1)) and (kmer[idx] == kmer[idx+1]):
                n_repeat += 1
                idx += 1
            if n_repeat >= 3:
                return True
        else:
            idx += 1
    return False