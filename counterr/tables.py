from .util import *

def context_sub_matrix_row2str(idx, context_sub_matrix , kmers, error_rates, N_obs):
    N_obs = float(N_obs)
    return "\t".join([kmers[idx], "%2.2f" % (error_rates[idx] * 100), "%d" % N_obs, "%2.2f" % (context_sub_matrix[idx, 0]/N_obs * 100), "%2.2f" % (context_sub_matrix[idx, 1]/N_obs * 100), "%.2f" % (context_sub_matrix[idx, 2]/N_obs * 100), "%2.2f" % (context_sub_matrix[idx, 3]/N_obs * 100), "%2.2f" % (context_sub_matrix[idx, 4]/N_obs * 100)]) + "\n"    

def context_ins_row2str(idx, idx_N1, kmers, error_rates, N_obs, context_hist_char_ins, context_hist_len_ins, len_max_ins):
    N_obs = float(N_obs)
    N_obs_char = float(np.sum(context_hist_char_ins[idx])) + 1e-6
    return "\t".join([kmers[idx], kmers[idx][idx_N1: idx_N1+2], 
                    "%2.2f" % (error_rates[idx] * 100), 
                    "%d" % N_obs, 
                    "%2.2f" % (context_hist_char_ins[idx, 0]/N_obs_char * 100) , 
                    "%2.2f" % (context_hist_char_ins[idx, 1]/N_obs_char * 100) , 
                    "%2.2f" % (context_hist_char_ins[idx, 2]/N_obs_char * 100) , 
                    "%2.2f" % (context_hist_char_ins[idx, 3]/N_obs_char * 100) , 
                    "\t".join(["%2.2f" % (context_hist_len_ins[idx, i]/N_obs * 100) for i in xrange(len_max_ins)])])+ "\n"

def tabulate_context_sub_matrix(context_sub_matrix, len_context_sub, output_dir, ex_hp=False, N_top_bottom=None, suffix=""):
    """
    Output in .tsv the context dependent substitution matrix ordered by error rate.

    kmer, error rate, N_obs, A, C, G, T, -

    The last five numbers are in counts. If ex_hp = True, then exclude k-mer that contain any homopolymers.
    """
    assert (len_context_sub % 2) == 1, "Subtitution context must be odd."
    
    N_kmer = N_bases ** len_context_sub
    assert N_kmer == context_sub_matrix.shape[0]

    idx_center = len_context_sub // 2 
    kmers = np.array([idx_to_kmer(len_context_sub, idx) for idx in xrange(N_kmer)])
    Ns_obs = np.array([context_sub_matrix[idx].sum() for idx in xrange(N_kmer)])
    error_rates = np.array([1. - (context_sub_matrix[idx, base2int[kmers[idx][idx_center]]] / float(Ns_obs[idx]+1e-6)) for idx in xrange(N_kmer)])

    # Sort all arrays according to error rates
    idx_sort = np.argsort(error_rates)[::-1]
    error_rates = error_rates[idx_sort]
    kmers = kmers[idx_sort]
    context_sub_matrix = context_sub_matrix[idx_sort, :]
    Ns_obs = Ns_obs[idx_sort]

    if ex_hp:
        iselect = ~np.array([has_hp(kmer) for kmer in kmers])
        error_rates = error_rates[iselect]
        kmers = kmers[iselect]
        context_sub_matrix = context_sub_matrix[iselect, :]
        Ns_obs = Ns_obs[iselect]

    f = open(os.path.join(output_dir, "".join(["context_sub", suffix, ".tsv"])), "w")
    f.write("Context\t%Error\tCounts\t%A\t%C\t%G\t%T\t%-\n")
    if N_top_bottom is None:
        for idx in xrange(Ns_obs.size):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                f.write(context_sub_matrix_row2str(idx, context_sub_matrix, kmers, error_rates, N_obs))
    else:
        strings_to_write_top = []
        # Bottom N
        strings_to_write_bottom = []
        counter = 0
        for idx in xrange(-1, -Ns_obs.size, -1):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                strings_to_write_bottom.append(context_sub_matrix_row2str(idx, context_sub_matrix, kmers, error_rates, N_obs))
                counter += 1
                if counter == N_top_bottom:
                    break        
        # Top N
        counter = 0
        for idx in xrange(Ns_obs.size):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                strings_to_write_top.append(context_sub_matrix_row2str(idx, context_sub_matrix, kmers, error_rates, N_obs))
                counter += 1                
                if counter == N_top_bottom:
                    break
        strings_to_write = strings_to_write_top + strings_to_write_bottom[::-1]
        for x in strings_to_write:
            f.write(x)
    f.close()

    return


def tabulate_context_ins_hist(context_hist_len_ins, context_hist_char_ins, len_context_ins, len_max_ins, output_dir, ex_hp=False, N_top_bottom=None, suffix=""):
    """
    Output in .tsv the context dependent insertion distribution ordered by error rate

    kmer, NN, error rate, N_obs, A, C, G, T, I0, I1, I2, ... Ilen_max_ins

    NN -- the two letters where the insertions are being recorded.
    """
    assert (len_context_ins % 2) == 0, "Insertion context must be even."
    
    N_kmer = N_bases ** len_context_ins
    assert N_kmer == context_hist_char_ins.shape[0]

    lengths = np.arange(len_max_ins)
    assert len_max_ins == context_hist_len_ins.shape[1]

    idx_N1 = (len_context_ins // 2) - 1 # The left character before insertion
    kmers = np.array([idx_to_kmer(len_context_ins, idx) for idx in xrange(N_kmer)])
    Ns_obs = np.array([context_hist_len_ins[idx].sum() for idx in xrange(N_kmer)])
    error_rates = np.array([1. - (context_hist_len_ins[idx, 0] / float(Ns_obs[idx]+1e-6)) for idx in xrange(N_kmer)])

    # Sort all arrays according to error rates
    idx_sort = np.argsort(error_rates)[::-1]
    error_rates = error_rates[idx_sort]
    kmers = kmers[idx_sort]
    context_hist_len_ins = context_hist_len_ins[idx_sort, :]
    context_hist_char_ins = context_hist_char_ins[idx_sort, :]    
    Ns_obs = Ns_obs[idx_sort]

    if ex_hp:
        iselect = ~np.array([has_hp(kmer) for kmer in kmers])
        error_rates = error_rates[iselect]
        kmers = kmers[iselect]
        context_hist_len_ins = context_hist_len_ins[iselect, :]
        context_hist_char_ins = context_hist_char_ins[iselect, :]    
        Ns_obs = Ns_obs[iselect]

    f = open(os.path.join(output_dir, "".join(["context_ins", suffix, ".tsv"])), "w")
    f.write("Context\tNN\t%Error\tCounts\t%A\t%C\t%G\t%T\t" + "\t".join(["%%I%d" % l for l in lengths]) + "\n")
    if N_top_bottom is None:
        for idx in xrange(Ns_obs.size):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                f.write(context_ins_row2str(idx, idx_N1, kmers, error_rates, N_obs, context_hist_char_ins, context_hist_len_ins, len_max_ins))
    else:
        strings_to_write_top = []
        # Bottom N
        strings_to_write_bottom = []
        counter = 0
        for idx in xrange(-1, -Ns_obs.size, -1):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                strings_to_write_bottom.append(context_ins_row2str(idx, idx_N1, kmers, error_rates, N_obs, context_hist_char_ins, context_hist_len_ins, len_max_ins))
                counter += 1
                if counter == N_top_bottom:
                    break
        # Top N
        counter = 0
        for idx in xrange(Ns_obs.size):
            N_obs = float(Ns_obs[idx])
            if N_obs > 0:
                strings_to_write_top.append(context_ins_row2str(idx, idx_N1, kmers, error_rates, N_obs, context_hist_char_ins, context_hist_len_ins, len_max_ins))
                counter += 1                
                if counter == N_top_bottom:
                    break
        strings_to_write = strings_to_write_top + strings_to_write_bottom[::-1]
        for x in strings_to_write:
            f.write(x)
    f.close()
    return

