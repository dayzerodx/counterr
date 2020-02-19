from .util import *


def make_read_filter(mapq_threshold, len_min_read, len_min_aln, bitflag):
    """
    Returns a read_filter to filter out the reads that do not meet the
    defined standards.

    Parameters
    ----------
    - mapq_threshold: Minimum mapQ threshold to pass
    - len_min_read  : Minimum length of the read to pass
    - len_min_aln   : Minimum length of the read that must align to pass
    - read          : Pysam.AlignmentSegment Object, the read that will be
                      checked against the filters
    - bitflag       : Used to determine which filters to use.
                      See http://samtools.github.io/hts-specs/SAMv1.pdf

    1    0x1   template having multiple segments in sequencing
    2    0x2   each segment properly aligned according to the aligner
    4    0x4   segment unmapped
    8    0x8   next segment in the template unmapped
    16   0x10  SEQ being reverse complemented
    32   0x20  SEQ of the next segment in the template being reverse complemented
    64   0x40  the first segment in the template
    128  0x80  the last segment in the template
    256  0x100 secondary alignment
    512  0x200 not passing filters, such as platform/vendor quality controls
    1024 0x400 PCR or optical duplicate
    2048 0x800 supplementary alignment

    If bitflag is set, then ~(bitflag & ---) will evaluate to False. In order
    for the read to pass, then the corresponding read property must be False
    so that ~property is evaluated to True.
    """
    def read_filter(read):
        cigar = read.cigarstring

        if ((read.mapping_quality > mapq_threshold)
                and ((read.query_length > len_min_read))
                and ((read.query_alignment_length > len_min_aln))
                and ((cigar.count("S") <= 2) and (cigar.count("H") <= 2))
                and (~read.is_unmapped or ~(bitflag & 4))
                and (~read.is_reverse or ~(bitflag & 16))
                and (~read.is_secondary or ~(bitflag & 256))
                and (~read.is_qcfail or ~(bitflag & 512))
                and (~read.is_duplicate or ~(bitflag & 1024))
                and (~read.is_supplementary or ~(bitflag & 2048))):

            return True
        return False

    return read_filter


def mapQ_stats_per_read(reads, verbose=False, comment=""):
    """
    For each read compute mean/median/std
    """
    if verbose:
        print("Computing per-read Q-score statistics " + comment)
        start = time()

    pile = []
    for name in reads.keys():
        for read in reads[name]:  # Iterate through the reads
            Qs = read.query_qualities
            if Qs is not None:
                pile.append(compute_mean_med_std(Qs))
            else:
                pile.append((0, 0, 0))

    means = [x[0] for x in pile]
    meds  = [x[1] for x in pile]
    stds  = [x[2] for x in pile]

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return means, meds, stds


def mapQ_stats_aligned_readsegment(reads, verbose):
    """
    Compute mean/std/len of each read grouped by aligned vs. unaligned regions
    """
    if verbose:
        print("Computing Q-score statistics by aligned/unaligned region")
        start = time()

    means_in = []
    stds_in = []
    lens_in = []
    means_out = []
    stds_out = []
    lens_out = []
    for name in reads.keys():
        for read in reads[name]:  # Iterate through the reads
            # Inside aligned region
            Q_in = read.query_alignment_qualities
            # Q_in is None implies that there are no alignment qualities, so we
            # append dummy values
            if Q_in is None:
                means_in.append(0)
                stds_in.append(0)
            else:
                means_in.append(np.mean(Q_in))
                stds_in.append(np.std(Q_in))
            lens_in.append(read.query_alignment_length)
            # Outside aligned region
            Q = read.query_qualities
            # There may not be qualities for the read, so creating an empty
            # array if there aren't any.
            if Q is None:
                Q_out = np.array([])
            else:
                Q_out = np.concatenate([Q[:read.query_alignment_start],
                                        Q[read.query_alignment_end:]])
            if (read.query_length != read.query_alignment_length) and (Q_out.size > 0): # Sometimes the returned Q_out is an emptry array.
                means_out.append(np.mean(Q_out))
                stds_out.append(np.std(Q_out))
                lens_out.append(read.query_length - read.query_alignment_length)
            else:
                means_out.append(0)
                stds_out.append(0)
                lens_out.append(0)

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return means_in, stds_in, lens_in, means_out, stds_out, lens_out


def reconstruct_all_alignments(contigs, reads, len_max_indel=20, correct_orientation=True, len_trim_contig_edge=0, verbose=False):
    """
    See reconstruct_alignment function for more detail regarding the output.

    Returns
    -------
    alns: A list of reconstructed alignments that passed the filter
    lens: The original query length.
    """
    if verbose:
        print("Reconstructing all alignments that pass the filter")
        start = time()

    alns = []
    lens = []
    for (name, contig) in contigs:
        if name in reads.keys():
            for read in reads[name]:  # Iterate through the reads
                alns.append(reconstruct_alignment(contig, read, len_max_indel, correct_orientation, len_trim_contig_edge))
                lens.append(read.query_length)

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return alns, lens


def reconstruct_alignment(reference, read, len_max_indel=20, correct_orientation=True, len_trim_contig_edge=0):
    """
    Given reference string and PySam readsegment object read, return a reconstructed alignment of the read and the reference. Example:

    Ref  : AGCT--GTCA--AAACCC
    Read : AGCATTG-CACCAAAGCC
    CIGAR: ===XIIMD==II===X==
    PhQ  : 333373303335333332

    Args
    ----
    - correct_orientation: Reverse complement reads whose reverse complement was mapped to the assembly and recorded.
    - len_trim_contig_edge: If this is non-zero, then trim reconstructed alingments if mapped at the edges of the reference contigs by len_trim_contig_edge amount.
    """
    # Place holder for returned values
    ref = []
    rd = []
    cig = []
    Qs = []

    # Get the query sequence
    query = read.query_sequence

    # Get the Q array
    Q_arr = read.query_alignment_qualities  # Full qualities
    if Q_arr is None:
        Q_arr = np.zeros(read.query_alignment_length)

    # Index variables
    pos_ref = read.reference_start # Refers to the full reference sequence
    pos_rd = read.query_alignment_start # Start of the read. Used to index both query sequence AND Phred Q-scores
    pos_Q = 0

    # Iterate through each cigar tuple and reconstruct the alignment
    for tup in read.cigar:
        code, len_cig = tup
        if code == 0: # match
            ref.append(reference[pos_ref:pos_ref+len_cig])
            rd.append(query[pos_rd:pos_rd+len_cig])
            cig.append("M" * len_cig)
            Qs.append(Q_arr[pos_Q:pos_Q+len_cig])
            pos_ref += len_cig
            pos_rd += len_cig
            pos_Q += len_cig
        elif code == 1: # Insertion
            if len_cig > len_max_indel:
                pass
            else:
                ref.append("-" * len_cig)
                rd.append(query[pos_rd:pos_rd+len_cig])
                cig.append("I" * len_cig)
                Qs.append(Q_arr[pos_Q:pos_Q+len_cig])
            pos_rd += len_cig
            pos_Q += len_cig
        elif code == 2: # Deletion
            if len_cig > len_max_indel:
                pass
            else:
                rd.append("-" * len_cig)
                ref.append(reference[pos_ref:pos_ref+len_cig])
                cig.append("D" * len_cig)
                Qs.append(np.zeros(len_cig, dtype=int))
            pos_ref += len_cig
        # elif (code != 4) and (code != 5): # Skip the initial cuts
            # pass
    ref = "".join(ref)
    rd = "".join(rd)
    cig = "".join(cig)
    Qs = np.concatenate(Qs)

    # --- Trim the edges as necessary
    # This is not exact but should be sufficiently accurate
    len_trim, len_trim2 = None, None
    if (read.reference_start < len_trim_contig_edge):
        len_trim = len_trim_contig_edge - read.reference_start
        ref = ref[len_trim:]
        rd  = rd[len_trim:]
        cig = cig[len_trim:]
        Qs  = Qs[len_trim:]
    if (read.reference_end >= (len(reference) - len_trim_contig_edge)):
        len_trim2 = (read.reference_end + 1) - (len(reference) - len_trim_contig_edge)
        ref = ref[:-len_trim2]
        rd  = rd[:-len_trim2]
        cig = cig[:-len_trim2]
        Qs  = Qs[:-len_trim2]

    if read.is_reverse and correct_orientation:
        ref = reverse_complement(ref)
        rd = reverse_complement(rd)
        cig = cig[::-1]
        Qs = Qs[::-1]

    # Modify the cigar to include "=" and "X" instead of "M".
    cig = modify_cigar(ref, rd, cig)

    return ref, rd, cig, Qs


def modify_cigar(ref, rd, cig):
    """
    Instead of "M" for "match", use "=" and "X".
    """
    cig_modified = []
    for i, c in enumerate(cig):
        if c == "M":
            if ref[i] == rd[i]:
                cig_modified.append("=")
            else:
                cig_modified.append("X")
        else:
            cig_modified.append(c)
    return "".join(cig_modified)


def find_hp(rd, h=3):
    """
    Find a homopolymer string corresponding to a single string.

    Example:
    rd = "ACGTTT-CG"
    print(rd)
    print("".join(find_hp(rd, h=3)))
    """
    len_rd = len(rd)
    hp = ["-"] * len_rd

    idx = 0
    while idx <= (len_rd-1):
        ltr = rd[idx] # Get the current letter
        if ltr in bases:
            idx_end = idx
            n_repeat = 0
            while rd[idx_end] in [ltr, "-"]: # Update idx_end until the character does not match
                if rd[idx_end] == ltr:
                    n_repeat += 1
                idx_end += 1
                if idx_end > (len_rd-1):
                    break

            if n_repeat >= h: # If number of repeats greater than min homopolymer length
                # If the last "match" was "-", then revert back
                while rd[idx_end-1] == "-":
                    idx_end -= 1
                hp[idx] = "S"  # Start of the homopolymer
                for m in range(idx+1, idx_end-1, 1):
                    hp[m] = "H"
                hp[idx_end-1] = "E"

            # Update the new start position
            idx = idx_end
        else:
            idx += 1

    return hp


def find_hp_region(ref, rd, cig, len_min_hp=3):
    """
    Given a reconstructed alignment, returns the corresponding homopolymer string.

    Args
    ----
    - len_min_hp: minimum repeated character in order for a region to be considered a homopolymer.
    """
    # Find homopolymer region solely based on the reference.
    hp = find_hp(ref, len_min_hp)

    len_ref = len(ref)
    idx = 0
    while idx < (len_ref-1):
        hc = hp[idx] # Current hp character
        ltr = ref[idx] # Current letter
        if hc == "S": # If start of the homopolymer region, check left to see if there is an insertion of the same character
            idx_hc = idx-1
            while (idx_hc >= 0) and (cig[idx_hc] in ["I", "="]) and (rd[idx_hc] == ltr):
                hp[idx_hc] = "S"
                hp[idx_hc+1] = "H"
                idx_hc -= 1
            idx += 1
        elif hc == "E":
            idx_hc = idx+1
            while (idx_hc <= (len_ref-1)) and (cig[idx_hc] in ["I", "="]) and (rd[idx_hc] == ltr):
                hp[idx_hc-1] = "H"
                hp[idx_hc] = "E"
                idx_hc += 1
            idx = idx_hc
        else:
            idx += 1

    return "".join(hp)


def get_hist_len_hp(contigs, len_min_hp, len_max_hp, len_trim_contig_edge, verbose=False):
    """
    Given contigs list, return the counts of homopolymers of various lengths for each of the letter in ["A", "C", "G", "T"]
    """
    if verbose:
        print("Computing homopolymer length histogram in the assembly")
        start = time()

    hist_len_hp = np.zeros((N_bases, len_max_hp))

    for (name, ref) in contigs: # For every contig in the reference
        len_ref = len(ref)
        idx = len_trim_contig_edge
        while idx < (len_ref-len_trim_contig_edge-1): # Could be a tighter bound
            ltr = ref[idx] # Get the current letter
            if ltr in bases:
                idx_end = idx
                while ref[idx_end] == ltr: # Update idx_end until the character does not match
                    idx_end += 1
                    if idx_end > (len_ref-len_trim_contig_edge-1):
                        break
                n_repeat = idx_end - idx # Get the number of repeat characters

                if (n_repeat >= len_min_hp) and (n_repeat < len_max_hp): # If number of repeats greater than min homopolymer length
                    hist_len_hp[base2int[ltr], n_repeat] += 1

                # Update the new start position
                idx = idx_end
            else:
                idx += 1

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return hist_len_hp


def phredQ_vs_error(alns, max_Q=100, verbose=False):
    """
    Given the reconstructed alignments, for each recorded Q value compute the various error rates. Note that deleted characters can't be associated with a phred score.
    """
    if verbose:
        print("Computing computed Q vs. empirical Q")
        start = time()

    # Counter dictionary for the Q-score
    # - First level: Integer key corresponding to the recorded Q-socre
    # - Second level: Counter for various types of errors.
    # Note that Q-score corresponding to deletin is zero.
    Q_dict = {}
    for i in range(max_Q):
        Q_dict[i] = {"X": 0, "I": 0, "=": 0, "D": 0}

    # Count errors corresponding to each q
    for aln in alns:
        ref, _, cig, Qs = aln
        for i, q in enumerate(Qs):
            if ref[i] in bases_ext:
                Q_dict[q][cig[i]] += 1

    # Make a total tally in numpy arrays and return
    qs = np.arange(1, max_Q)
    len_qs = qs.size
    nums_sub = np.zeros(len_qs, dtype=np_int_type)
    nums_ins = np.zeros(len_qs, dtype=np_int_type)
    nums_match = np.zeros(len_qs, dtype=np_int_type)
    nums_err = np.zeros(len_qs, dtype=np_int_type)
    nums_tot = np.zeros(len_qs, dtype=np_int_type)

    for i, q in enumerate(qs):
        counts = Q_dict[q]
        nums_match[i] = (counts["="])
        nums_sub[i] = (counts["X"])
        nums_ins[i] = (counts["I"])
        nums_err[i] = (counts["X"] + counts["I"])
        nums_tot[i] = (counts["X"]+counts["I"]+counts["="])

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return qs, nums_match, nums_sub, nums_ins, nums_err, nums_tot


def count_errors_per_read(alns, verbose, skip_nonACGT=True):
    """
    Given reconstructed alignments, count errors of various types on per read basis. By default, if the read maps to non-ACGT region, then skip them when counting errors. Return alignment lengths.
    """
    if verbose:
        print("Counting errors per read")
        start = time()

    nums_match = np.zeros(len(alns))
    nums_sub = np.zeros(len(alns))
    nums_ins = np.zeros(len(alns))
    nums_del = np.zeros(len(alns))
    nums_skip = np.zeros(len(alns))

    for k, aln in enumerate(alns):
        ref, rd, cig, _ = aln

        for i in range(len(ref)):
            if ref[i] in bases_ext:
                c = cig[i]
                if c == "=":
                    nums_match[k] += 1
                elif c == "X":
                    nums_sub[k] += 1
                elif c == "D":
                    nums_del[k] += 1
                elif c == "I":
                    nums_ins[k] += 1
            else:
                nums_skip[k] += 1
    # Alignment length
    lens_aligned = nums_match + nums_sub + nums_ins + nums_del

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned


def count_errors_hp(alns_scrubbed, len_min_hp=3, len_max_hp=20, verbose=False):
    """
    Given a set of scrubbed reconstructed alignments, compute truth vs. observed homopolymer length histogram.

    Args
    ----
    - alns_scrubbed: a list of reconstructed alignments
    - len_min_hp: minimum length for a string to be considered a homopolymer.

    Returns
    -------
    - Homopolymer length distribution
    """
    if verbose:
        print("Computing homopolymer true vs. observed length distribution")
        start = time()

    # For tallying homopolymer lengths distributions.
    hist_len_hp = np.zeros((N_bases, len_max_hp, len_max_hp))

    for aln in alns_scrubbed:
        ref, rd, cig, _, hp = aln

        # Iterate through each homopolymer
        len_ref = len(ref)
        idx = 0
        while idx < len_ref-1:
            c_hp = None
            if hp[idx] == "S": # If you are at a start of the homopolymer region
                num_del = 0
                num_ins = 0

                # First boundary
                if cig[idx] == "D":
                    num_del += 1
                    c_hp = ref[idx] # Hmopolymer character
                elif cig[idx] == "I":
                    num_ins += 1
                    c_hp = rd[idx]
                elif cig[idx] == "=": # If it's a match
                    c_hp = ref[idx]
                elif cig[idx] == "X": # If it's a mis-match
                    c_hp = ref[idx]
                else:
                    assert False
                len_hp = 1
                assert c_hp is not None

                while hp[idx+len_hp] != "E":
                    if cig[idx+len_hp] == "D":
                        num_del += 1
                    elif cig[idx+len_hp] == "I":
                        num_ins += 1
                    len_hp += 1
                assert (hp[idx+len_hp] == "E")

                # Last boundary
                if cig[idx+len_hp] == "D":
                    num_del += 1
                elif cig[idx+len_hp] == "I":
                    num_ins += 1
                len_hp += 1 # This includes the last basepair in the homopolymer

                # Record the length distribution
                len_hp_orig = len_hp - num_ins
                len_observed = len_hp - num_del
                if (len_hp_orig >= len_min_hp) and (len_hp_orig < len_max_hp) and (len_observed < len_max_hp):
                    hist_len_hp[base2int[c_hp], len_hp_orig, len_observed] += 1

                # Increment index
                idx += len_hp
            else:
                idx += 1

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return hist_len_hp


def cutout_hp_region(ref, rd, cig, hp):
    len_hp = len(hp)
    ref_ex = []
    rd_ex = []
    cig_ex = []
    idx_start = 0
    idx_end = 0
    while idx_end < (len_hp-1):
        if hp[idx_end] == "S": # Start of a homopolymer
            ref_ex.append(ref[idx_start:idx_end])
            rd_ex.append(rd[idx_start:idx_end])
            cig_ex.append(cig[idx_start:idx_end])
            while hp[idx_end] != "E":
                idx_end += 1
            idx_end += 1 # The first character after homopolymer region
            idx_start = idx_end
        else:
            idx_end += 1
    ref = "".join(ref_ex)
    rd = "".join(rd_ex)
    cig = "".join(cig_ex)

    return ref, rd, cig


def count_errors(alns_scrubbed, len_min_hp=3, exclude_hp = True, len_max_indel=20, verbose=False):
    """
    Given a set of scrubbed reconstructed alignments, count errors of various types. This can be computed including or excluding homopolymer region.

    Args
    ----
    - alns_scrubbed: a list of reconstructed alignments
    - exclude_hp: If True, exclude homopolymer region.
    - len_min_hp: Minimum length for something to be considered polymer region.

    Return
    ------
    - sub_matrix: Substitution matrix
    - hist_len_ins/del: indel length distribution
    """
    if verbose:
        if exclude_hp:
            print("Counting errors, excluding homopolymer regions")
        else:
            print("Counting errors, including homopolymer regions")
        start = time()

    # Tally mismatches outside indel regions
    sub_matrix = np.zeros((5, 5), dtype=int)

    # Length distribution of insertions and deletions
    hist_len_del = np.zeros(len_max_indel)
    hist_len_ins = np.zeros(len_max_indel)

    for aln in alns_scrubbed:
        ref, rd, cig, _, hp = aln

        if exclude_hp:
            ref, rd, cig = cutout_hp_region(ref, rd, cig, hp)

        # ----- start counting errors
        idx = 0
        len_ref = len(ref)
        while idx <= (len_ref-1):
            c = cig[idx]
            if c in ["X", "="]: # Match
               sub_matrix[base2int[ref[idx]], base2int[rd[idx]]] += 1
               idx += 1
            elif c == "D": # Deletion
                num_del = 0
                while cig[idx] == "D":
                    sub_matrix[base2int[ref[idx]], -1] += 1
                    idx += 1
                    num_del += 1
                    if idx > (len_ref-1):
                        break
                if num_del < len_max_indel:
                    hist_len_del[num_del] += 1
            elif c == "I":
                num_ins = 0
                while cig[idx] == "I":
                    sub_matrix[-1, base2int[rd[idx]]] += 1
                    idx += 1
                    num_ins += 1
                    if idx > (len_ref-1):
                        break
                if num_ins < len_max_indel:
                    hist_len_ins[num_ins] += 1
            else:
                assert False, "Something went wrong"

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return sub_matrix, hist_len_del, hist_len_ins


def context_sub_table(alns_scrubbed, len_context_sub=7, verbose=False):
    """
    Compile a table of errors where each row corresponds to a substitution table for the middle character of a g-mer (where g = len_context_sub and odd).
    """
    if verbose:
        print("Computing %d-mer substitution matrix" % len_context_sub)
        start = time()

    assert (len_context_sub % 2) == 1, "len_context_sub must be odd"
    sub_matrix = np.zeros((4**len_context_sub, 5))

    len_side_mer = len_context_sub // 2
    acceptable_cigar = {"D", "=", "X"}
    for aln in alns_scrubbed:
        try: # Wrap the counting text inside try-except to avoid stopping the program due to a pesky indexing error.
            ref, rd, cig, _, _ = aln
            len_ref = len(ref)
            # Find start and end points
            counter = 0
            idx_start = 0
            while counter < len_context_sub:
                if ref[idx_start] != "-":
                    counter += 1
                idx_start += 1
            counter = 0
            idx_end = len_ref-1
            while counter < len_context_sub:
                if ref[idx_end] != "-":
                    counter += 1
                    if counter == len_context_sub:
                        break
                idx_end -= 1

            for idx in xrange(idx_start, idx_end, 1):
                c = cig[idx]
                if c in acceptable_cigar: # If the center character of the read is not an insertion
                    # --- Find the g-mer centered around idx.
                    # Finding left index
                    len_left = 0
                    s = idx-1
                    while True:
                        if ref[s] != "-":
                            len_left += 1
                        if (len_left == len_side_mer) or (s == 0):
                            break
                        s -= 1

                    # Finding the right index
                    len_right = 0
                    e = idx+1
                    while True:
                        if ref[e] != "-":
                            len_right += 1
                        if (len_right == len_side_mer) or (e == (len_ref-1)):
                            break
                        e += 1

                    context = ref[s:e+1].replace("-", "")
                    if len(context) == len_context_sub:
                        # Record what the error is!
                        context_idx = kmer_to_idx(context)
                        sub_matrix[context_idx][base2int[rd[idx]]] += 1
        except:
            continue

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return sub_matrix


def context_ins_table(alns_scrubbed, len_context_ins=8, len_max_ins=20, verbose=False):
    """
    Compile a table of errors where each row corresponds to insertions occuring in the middle of f-mer (f = len_context_ins and even).
    """
    if verbose:
        print("Computing %d-mer insertion error table" % len_context_ins)
        start = time()

    assert (len_context_ins % 2) == 0, "len_context_ins must be even"

    # 1) tally of inserted characters
    # 2) tally of observed insertion length.
    hist_len_ins = np.zeros((4**len_context_ins, len_max_ins), dtype=int)
    hist_char_ins = np.zeros((4**len_context_ins, 4), dtype=int)

    acceptable_cigar = {"D", "=", "X"}
    len_side_mer = len_context_ins // 2
    for aln in alns_scrubbed:
        try: # Wrap the counting text inside try-except to avoid stopping the program due to a pesky indexing error.
            ref, rd, cig, _, _ = aln
            len_ref = len(ref)
            # Find start and end points
            counter = 0
            idx_start = 0
            while counter < len_context_ins:
                if ref[idx_start] != "-":
                    counter += 1
                idx_start += 1
            counter = 0
            idx_end = len_ref-1
            while counter < len_context_ins:
                if ref[idx_end] != "-":
                    counter += 1
                    if counter == len_context_ins:
                        break
                idx_end -= 1

            # --- Record errors
            # idx indicates the left character next to insertion. For instance, G in ACG----TCC.
            for idx in xrange(idx_start, idx_end, 1):
                c = cig[idx] # Cigar character
                if c in acceptable_cigar: # If the current character is not an insertion.
                    # Finding the left context
                    len_left = 0
                    s = idx
                    while True:
                        if ref[s] != "-":
                            len_left += 1
                        if (len_left == len_side_mer) or (s == 0):
                            break
                        s -= 1

                    # Finding the right context
                    len_right = 0
                    e = idx+1
                    while True:
                        if ref[e] != "-":
                            len_right += 1
                        if (len_right == len_side_mer) or (e == (len_ref-1)):
                            break
                        e += 1

                    # Total context
                    context = ref[s:e+1].replace("-", "")
                    if len(context) == len_context_ins:
                        # Find the inserted string
                        right_mer = ref[idx+1:e+1]
                        rd_right_mer = rd[idx+1:e+1]
                        len_ins = 0
                        while right_mer[len_ins] == "-":
                            len_ins += 1

                        # Record what the error is!
                        context_idx = kmer_to_idx(context)
                        if len_ins < len_max_ins:
                            hist_len_ins[context_idx, len_ins] += 1
                            for ci in rd_right_mer[:len_ins]:
                                hist_char_ins[context_idx, base2int[ci]] += 1
        except:
            continue

    # The total length of insertions must be equal to the total number of characters inserted.
    assert (hist_len_ins * np.arange(hist_len_ins.shape[1])).sum() == hist_char_ins.sum()

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return hist_len_ins, hist_char_ins

def compute_indel_rates(sub_matrix, dist_len_ins, dist_len_del):
    """
    Given substitution matrix that includes indel and indel length distributions, return indel rates.
    """
    # Total number of transition opportunitites.
    N_trans = np.sum(sub_matrix[:4, :4])

    # Total number of indels
    N_ins = np.sum(dist_len_ins)
    N_del = np.sum(dist_len_del)

    # Indel rates
    R_ins = N_ins/float(N_trans)
    R_del = N_del/float(N_trans)

    return R_ins, R_del


def append_hp(aln, h):
    """
    Given a reconstructed alignment, add homopolymer string.
    """
    # Unpack
    ref, rd, cig, Qs = aln

    # Compute the homopolymer string
    hp = find_hp_region(ref, rd, cig, len_min_hp=h)

    return (ref, rd, cig, Qs, hp)


def split_recontructed_alignment(aln):
    """
    Break happens at non-ACGT sites
    """
    aln_split = [] # Alignment split into list of alignments
    ref, rd, cig, Qs = aln
    len_ref = len(ref)
    idx_start = 0
    idx_end = 0
    # Keep skipping until non-ACGT letter encountered. Tuck away portions [idx_start:idx_end] and repeat.
    while idx_end < (len_ref-1):
        if (ref[idx_end] in bases_ext) and (rd[idx_end] in bases_ext):
            idx_end += 1
        else:
            aln_split.append((ref[idx_start:idx_end], rd[idx_start:idx_end], cig[idx_start:idx_end], Qs[idx_start:idx_end]))
            # Skip over non-ACGT letters
            while (idx_end <= (len_ref-1)) and ((ref[idx_end] not in bases_ext) or (rd[idx_end] not in bases_ext)):
                idx_end += 1
            idx_start = idx_end
    # Treat the last case separately
    aln_split.append((ref[idx_start:idx_end], rd[idx_start:idx_end], cig[idx_start:idx_end], Qs[idx_start:idx_end]))
    return aln_split


def scrub_reconstructed_alignments(alns, h=3, verbose=False):
    """
    1) Break reconstructed reads in regions where there are non-ACGT basepairs.
    2) Add homopolymer strings
    """
    if verbose:
        print("Splitting alignments at non-ACGT sites and homopolymer marker string")
        start = time()

    alns_scrubbed = []

    for aln in alns:
        # Break
        aln_split = split_recontructed_alignment(aln)

        # Add homopolymer string for each and save
        for a in aln_split:
            a = append_hp(a, h)
            alns_scrubbed.append(a)

    if verbose:
        end = time()
        print("Time taken: %.2f seconds" % (end-start))

    return alns_scrubbed