from .util import *

def load_files(asm, bam, bai, len_min_contig, read_filter, lim, verbose, cram=False):
    """
    Load the assembly contigs and their corresopnding reads, applying variosu filters.

    Parameters
    ----------
    - lim: If -1, load and store all reads. If not, randomly select lim number of samples.

    Return
    ------
    - contigs   : A list of tuples (contig_name, contig_str)
    - reads_pass: A dictionary with key = contig name and val = list of reads that passed the filter. 
    - reads_fail: Same as above except for those did not pass the filter.
    """
    # ---- Load the reference and eliminate contigs that do not meet the cut.
    fasta = pysam.FastaFile(asm)
    contigs = [(name, fasta.fetch(name)) for name in fasta.references]
    for i in xrange(len(contigs)-1, -1, -1):
        name, contig = contigs[i]
        if len(contig) < len_min_contig:
            del contigs[i]
    contig_names = [x[0] for x in contigs]
    fasta.close()

    # ---- Load the alignment file, save only those that mapped to the contigs saved above. 
    if cram:
        alignment = pysam.AlignmentFile(bam, index_filename = bai, mode ="rc")
    else:
        alignment = pysam.AlignmentFile(bam, index_filename = bai, mode ="rb")
    reads_pass = {} 
    reads_fail = {} 

    N_pass = 0
    N_fail = 0
    if lim == -1:
        for (name, _) in contigs: # For every contig in the reference
            reads_pass[name] = []
            reads_fail[name] = []
            for read in alignment.fetch(name):
                if read_filter(read):
                    reads_pass[name].append(read)
                    N_pass += 1
                else:
                    reads_fail[name].append(read)
                    N_fail += 1
    else:
        # Get index statistics
        index_stats = alignment.get_index_statistics()
        nums_mapped = [index.mapped for index in index_stats if index.contig in contig_names]
        names = [index.contig for index in index_stats if index.contig in contig_names]        
        num_total = sum(nums_mapped)

        if lim > num_total:
            lim = num_total

        # Select a fixed number of reads to choose for each contig
        nums_pick = np.random.multinomial(lim, np.asarray(nums_mapped)/float(num_total), size=1)[0]

        for (name, _) in contigs:
            idx = names.index(name)
            num_pick = nums_pick[idx]
            num_tot = nums_mapped[idx]
            if (num_tot == 0) or (num_pick == 0):
                continue
            assert num_pick <= num_tot
            reads_pass[name] = []
            reads_fail[name] = []
            N_reads = 0
            indices = np.random.choice(range(num_tot), size=num_pick, replace=False)
            for idx, read in enumerate(alignment.fetch(name)):
                if idx in indices:
                    if read_filter(read):
                        reads_pass[name].append(read)
                        N_pass += 1
                    else:
                        reads_fail[name].append(read)
                        N_fail += 1
                    N_reads += 1
                    if N_reads == num_pick:
                        break
    alignment.close()

    if verbose:
        print("Number of reads that passed/failed the filter: %d/%d" % (N_pass, N_fail))

    if N_pass == 0:
        print("No read passed the filter.\nTo allow more reads, please adjust the filter parameters: mapq_thres, len_min_read, len_min_aln, len_min_contig")
        sys.exit()        

    return contigs, reads_pass, reads_fail

def save_per_read_stats_fail(item_to_save, output_dir):
    """
    Parameters
    ----------
    - item_to_save: A list of arrays corresponding to the reads that did not pass the read filter
    """
    recarr = np.rec.fromarrays(item_to_save, dtype=[('mean', np.float32), ('med', np.float32), ('std', np.float32)])
    np.save(os.path.join(output_dir, "per_read_stats_fail"), recarr)

    return

def save_per_read_stats_pass(item_to_save, output_dir):
    """
    Parameters
    ----------
    - item_to_save: A list of arrays corresponding to the reads that did pass the read filter
    """
    recarr = np.rec.fromarrays(item_to_save, dtype=[('mean', np.float32), ('med', np.float32), ('std', np.float32), ('mean_in' , np.float32), ('std_in',  np.float32), ('len_in', np.int64), ('mean_out',  np.float32), ('std_out',  np.float32), ('len_out', np.int64), ('num_match', np.int64), ('num_sub', np.int64), ('num_ins', np.int64), ('num_del', np.int64), ('num_skip', np.int64), ('len_aligned', np.int64), ('len', np.int64)])
    np.save(os.path.join(output_dir, "per_read_stats_pass"), recarr)

    return    

def save_phredQ_stats(item_to_save, output_dir):
    recarr = np.rec.fromarrays(item_to_save, dtype=[('q', np.int64), ('num_match', np.int64), ('num_sub', np.int64), ('num_ins', np.int64), ('num_err', np.int64), ('num_tot', np.int64)])
    np.save(os.path.join(output_dir, "phredQ_stats"), recarr)
    return

def save_hist_len_hp_asm(hist_len_hp_asm, len_max_hp, output_dir):
    lengths = np.arange(len_max_hp)
    item_to_save = [lengths] + [hist_len_hp_asm[i, :] for i in xrange(N_bases)]
    # Reformat for clarity
    recarr = np.rec.fromarrays(item_to_save, dtype=[('length', np.int64), ('A', np.int64), ('C', np.int64), ('G', np.int64), ('T', np.int64)])
    np.save(os.path.join(output_dir, "hist_len_hp_asm"), recarr)

    return

def save_len_hp(hist_len_hp, output_dir):
    """
    Output format:
    - axis=0: Letters A, C, G, T
    - axis=1: Truth length
    - axis=2: Observed length
    """
    np.save(os.path.join(output_dir, "hist_len_hp"), hist_len_hp)

    return

def save_context_independent_errors(item_to_save, output_dir, suffix=""):
    sub_matrix, hist_len_del, hist_len_ins, R_ins, R_del = item_to_save
    np.savez(os.path.join(output_dir, "".join(["context_independent_errors", suffix])), sub_matrix=sub_matrix, hist_len_del=hist_len_del, hist_len_ins=hist_len_ins, R_ins=R_ins, R_del=R_del)
    return

def save_context_dependent_errors(item_to_save, output_dir):
    context_sub_matrix, context_hist_len_ins, context_hist_char_ins = item_to_save
    np.savez(os.path.join(output_dir, "".join(["context_dependent_errors"])), sub_matrix=context_sub_matrix, hist_len_ins=context_hist_len_ins, hist_char_ins=context_hist_char_ins) 
    return
