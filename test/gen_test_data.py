from counterr.util import *
from counterr.counters import modify_cigar

def gen_rand_seq(L):
    return "".join(np.random.choice(bases_list, size=L))

def errorful_copy(seq, error_rate):
    """
    Given a sequence of DNA letters seq (1D numpy array), return a copy with errors.
    """
    N_errors = int(error_rate * seq.size)
    idx_change = np.random.randint(low=0, high=seq.size, size=N_errors)
    seq_w_errors = np.copy(seq)
    seq_w_errors[idx_change] = gen_rand_seq(N_errors)

    return seq_w_errors

class MockPysamRead(object):
    """
    Mock pysam read class used to test reconstruct_alignment function.
    """
    def __init__(self, contig, L_read, idx_start, R_ins, R_del, R_x):
        """
        Return a mock pysam read object, making an errorful copy of the contig.

        Parameters
        ----------
        - contig    : Generating sequence
        - L_read    : Length of the desired read
        - idx_start : Start position of the alignment on the contig
        - R_ins     : Insertion rate
        - R_del     : Deletion rate
        - R_x       : Mis-match rate

        Properties
        ----------
        - reference_start          : Same as idx_start
        - query_alignment_start    : Start position of the alignment on the read
        - query                    : Read sequence
        - cigar                    : Cigar tuple
        - ref, rd, cig, Qs         : Sequences corresponding to the alignment
        - query_alignment_qualities: Quality scores corresponding to the aligned portion of the read
        """

        # Constant params
        self.reference_start = idx_start 
        self.query_alignment_start = 100 # The first hundred letters are random junk.
        self.is_reverse = False # Assume that mapping is in the correct alignment

        # Initialize
        self.query = []
        self.query_alignment_qualities = []
        self.cigar = []
        # Answer
        self.ref = []
        self.rd = []
        self.cig = []
        self.Qs = []        

        # Generate the sequence
        idx = self.reference_start # index on the contig
        len_read = 0 # Length of the read generated so far
        while len_read < L_read:
            c = np.random.choice(["M", "I", "D"], p=[1-R_ins-R_del, R_ins, R_del])
            if c == "M": # If it's a match
                len_match = np.random.randint(low=1, high=10, size=1)[0]
                # Cigar tuple
                self.cigar.append((0, len_match))
                # Query update
                ref_add = contig[idx:idx+len_match]
                rd_add = "".join(errorful_copy(np.asarray(list(ref_add)), R_x))
                Q_add = np.random.randint(low=0, high=50, size=len_match)
                self.query.append(rd_add)
                self.query_alignment_qualities.append(Q_add)
                # Answer
                self.ref.append(ref_add)
                self.rd.append(rd_add)
                self.cig.append("M" * len_match)
                self.Qs.append(Q_add)
                # Update contig index/read length
                idx += len_match
                len_read += len_match
            elif c == "I":
                len_ins = np.random.randint(low=1, high=5, size=1)[0]
                # Cigar tuple
                self.cigar.append((1, len_ins))
                # Query update
                rd_add = gen_rand_seq(len_ins)
                Q_add = np.random.randint(low=0, high=50, size=len_ins)
                self.query.append(rd_add)
                self.query_alignment_qualities.append(Q_add)
                # Answer
                self.ref.append("-" * len_ins)
                self.rd.append(rd_add)
                self.cig.append("I" * len_ins)
                self.Qs.append(Q_add)
                # Update contig index/read length
                len_read += len_ins             
            elif c == "D":
                len_del = np.random.randint(low=1, high=5, size=1)[0]
                # Cigar tuple
                self.cigar.append((2, len_del))
                # Query update
                ref_add = contig[idx:idx+len_del]
                # Answer
                self.ref.append(ref_add)
                self.rd.append("-" * len_del)
                self.cig.append("D" * len_del)
                self.Qs.append(np.zeros(len_del))
                # Update contig index/read length
                idx += len_del
            else:
                assert False
        # Concatenate the results
        self.query_sequence = "".join([gen_rand_seq(self.query_alignment_start), "".join(self.query), gen_rand_seq(self.query_alignment_start)])
        self.query_alignment_qualities = np.concatenate(self.query_alignment_qualities)
        # Answer
        self.ref = "".join(self.ref)
        self.rd  = "".join(self.rd)
        self.cig = "".join(self.cig)
        self.Qs  = np.concatenate(self.Qs)
        # Modify the cigar
        self.cig = modify_cigar(self.ref, self.rd, self.cig)
        # Boundary conditions
        self.reference_end = idx

        return

def gen_mock_hp(L, len_max_hp, return_hist_len_hp = False):
    """
    Generate mock example that has the correct hp string.
    """
    # Initialization
    idx = 0
    ref = []
    rd  = []
    cig = []
    hp_truth  = []
    ltr_previous = "A"

    # Recorded observed vs. homopolymer length
    hist_len_hp = np.zeros((N_bases, len_max_hp, len_max_hp))

    while idx < L:
        case = np.random.choice(["hp", "not_hp"], p=[0.1, 0.9])
        bases_list = ["A", "C", "G", "T"]
        bases_list.remove(ltr_previous)
        ltr = np.random.choice(bases_list)
        ltr_previous = ltr

        if case == "hp":
            len_hp = np.random.randint(low=3, high=len_max_hp, size=1)[0]
            len_hp_obs = np.random.randint(low=max(0, len_hp-3), high=min(len_hp+3, len_max_hp), size=1)[0]
            hist_len_hp[base2int[ltr], len_hp, len_hp_obs] += 1 # Record the true vs. obs length
            len_hp_region = max(len_hp, len_hp_obs) # Larger of the two number
            ref_add = [ltr] * len_hp_region
            rd_add = [ltr] * len_hp_region
            cig_add = ["="] * len_hp_region
            idx_edits = np.random.choice(range(len_hp_region), size=np.abs(len_hp-len_hp_obs), replace=False)
            if len_hp > len_hp_obs:
                for i in idx_edits: 
                    rd_add[i] = "-"
                    cig_add[i] = "D"
            elif len_hp < len_hp_obs:
                # Insertion would normally be all collected at one end but thi is more general
                for i in idx_edits:
                    ref_add[i] = "-"
                    cig_add[i] = "I"
            else:
                pass
            ref.append("".join(ref_add))
            rd.append("".join(rd_add))
            cig.append("".join(cig_add))                
            hp_truth.append("".join(["S", "H" * (len_hp_region-2), "E"]))
            idx += len_hp_region
        else:
            ref.append(ltr)
            rd.append(ltr)
            cig.append("=")
            hp_truth.append("-")
            idx += 1
    ref = "".join(ref)
    rd = "".join(rd)
    cig = "".join(cig)
    hp_truth = "".join(hp_truth)

    if return_hist_len_hp:
        return ref, rd, cig, hp_truth, hist_len_hp
    else:
        return ref, rd, cig, hp_truth

def gen_mock_errors(L, len_max_indel):
    # Initialization
    idx = 0
    ref = []
    rd  = []
    cig = []
    sub_matrix = np.zeros((5, 5))
    hist_len_del = np.zeros(len_max_indel)
    hist_len_ins = np.zeros(len_max_indel)
    R_ins = 0.05
    R_del = 0.05
    R_x = 0.2
    R_m = 1 - R_ins - R_del - R_x
    case_previous = "="

    while idx < L:
        cases = ["I", "D", "=", "X"]
        cases.remove(case_previous)
        case = np.random.choice(cases)
        case_previous = case
        L_error = np.random.randint(low=1, high=len_max_indel, size=1)[0]
        if case == "=":
            add_str = "".join(np.random.choice(["A", "C", "G", "T"], size=L_error))
            ref.append(add_str)
            rd.append(add_str)
            cig.append("=" * L_error)
            for c in add_str:
                sub_matrix[base2int[c], base2int[c]] += 1
        elif case == "X":
            ref_add = "".join(np.random.choice(["A", "C", "G", "T"], size=L_error))
            ref.append(ref_add)           
            rd_add = []
            for c in ref_add:
                bases = ["A", "C", "G", "T"]
                bases.remove(c)
                rd_add.append(np.random.choice(bases))
            rd_add = "".join(rd_add)
            rd.append(rd_add)
            cig.append("X" * L_error)
            for i in range(L_error):
                sub_matrix[base2int[ref_add[i]], base2int[rd_add[i]]] += 1
        elif case == "I":
            ref_add = "-" * L_error
            rd_add = "".join(np.random.choice(["A", "C", "G", "T"], size=L_error))
            ref.append(ref_add)           
            rd.append(rd_add)
            cig.append("I" * L_error)
            for i in range(L_error):
                sub_matrix[-1, base2int[rd_add[i]]] += 1
            hist_len_ins[L_error] += 1
        elif case == "D":
            ref_add = "".join(np.random.choice(["A", "C", "G", "T"], size=L_error))
            rd_add = "-" * L_error            
            ref.append(ref_add)
            rd.append(rd_add)
            cig.append("D" * L_error)
            for i in range(L_error):
                sub_matrix[base2int[ref_add[i]], -1] += 1
            hist_len_del[L_error] += 1
        idx += L_error

    ref = "".join(ref)
    rd = "".join(rd)
    cig = "".join(cig)

    return ref, rd, cig, sub_matrix, hist_len_del, hist_len_ins


def gen_mock_context_sub(L, len_context_sub=5):
    # Initialization
    cig = []
    sub_matrix = np.zeros((4**len_context_sub, 5))
    R_ins = 0.1
    R_del = 0.1
    R_x = 0.2

    # --- Procedure
    # Generate a long reference
    ref = np.random.choice(bases_list, size=L)

    # Generate read introducing "X" and "X" but not "I". Update sub_mtrix
    idx_start = len_context_sub
    idx_end = L - len_context_sub

    rd = errorful_copy(ref, R_x)
    N_del = int(L * R_del)
    idx_del = np.random.choice(range(idx_start, idx_end), replace=False, size=N_del)
    rd[idx_del] = "-"
    for i in range(L):
        if ref[i] == rd[i]:
            cig.append("=")
        elif rd[i] == "-":
            cig.append("D")
        elif ref[i] != rd[i]:
            cig.append("X")

    # Update submatrx
    len_side = len_context_sub // 2
    for i in range(idx_start, idx_end, 1):
        idx = kmer_to_idx(ref[i-len_side:i+len_side+1])
        sub_matrix[idx, base2int[rd[i]]] += 1

    # At various places, append insertion. This should not change the sub_matrix
    ref = list(ref)
    rd = list(rd)
    N_ins = int(L * R_del)
    idx_ins = np.sort(np.random.choice(range(L), replace=False, size=N_ins))[::-1]
    for i in idx_ins:
        L_ins = np.random.randint(low=1, high=5, size=1)[0]
        ref = ref[:i] + ["-"*L_ins] + ref[i:]
        rd = rd[:i] + list(gen_rand_seq(L_ins)) + rd[i:]        
        cig = cig[:i] + ["I"*L_ins] + cig[i:]

    ref = "".join(ref)
    rd = "".join(rd)
    cig = "".join(cig)

    return ref, rd, cig, sub_matrix


if __name__ == "__main__":
    # Set numpy random seed
    np.random.seed(42)
    L_contig = int(1e6)

    # Generate reference sequence
    contig = gen_rand_seq(L_contig)

    # Generate mock read
    L_read = int(5e5) 
    idx_start = int(1e5)
    R_ins = 1e-1
    R_del = 1e-1
    R_x = 5e-1
    read = MockPysamRead(contig, L_read, idx_start, R_ins, R_del, R_x)

    # Print the answer
    print(read.ref[:100])
    print(read.rd[:100])
    print(read.cig[:100])
    print(read.Qs[:100])