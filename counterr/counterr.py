from .IO import *
from .counters import *
from .plots import *
from .tables import *

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # --- Required arguments
    parser.add_argument("-bam", required=True, help="the input bam file")
    parser.add_argument("-genome", required=True, help="the input fasta file")
    parser.add_argument("-output_dir", required=True, help="the output directory for figures and stats")
    # --- Optional arguments
    parser.add_argument("-no_figures", help="pass this flag to not generate figures", action="store_true")
    parser.add_argument("-bai", help="the input bai filename if non-conventionally named", type=str, default=None)
    parser.add_argument("-cram", help="treat bam and bai inputs as cram and crai", action="store_true")
    parser.add_argument("-verbose", help="pass this flag to follow progress in the terminal", action="store_true")
    parser.add_argument("-len_min_contig", help="minimum contig length", type=int, default=1500)    
    parser.add_argument("-mapq_thres", help="minimum mapq threshold", type=int, default=40)
    parser.add_argument("-len_min_read", help="minimum read length", type=int, default=1500)
    parser.add_argument("-len_min_aln", help="minimum length aligned", type=int, default=1000)
    parser.add_argument("-bitflag", help="bit flag for read filter (as specified in SAM doc) ", type=int, default=3845)
    parser.add_argument("-len_min_hp", help="minimum homopolymer length", type=int, default=3)
    parser.add_argument("-len_max_hp", help="maximum homopolymer length", type=int, default=20)
    parser.add_argument("-len_context_sub", help="length of the k-mer context for context dependent substitution table", type=int, default=5)
    parser.add_argument("-len_context_ins", help="length of the k-mer context for context dependent insertion table", type=int, default=6)
    parser.add_argument("-len_max_indel", help="maximum length of indels to consider", type=int, default=20)
    parser.add_argument("-len_trim_contig_edge", help="length of the contig edge to trim before computing various statistics", type=int, default=1)
    parser.add_argument("-use_recorded", help="pass this flag to NOT perform reverse complementing of the reverse complement mapped reads", action="store_true")
    parser.add_argument("-lim", help="pass this flag to run the program with 'lim' randomly selected reads (both pass and fail)", type=int, default=-1)
    parser.add_argument("-num_pts_max", help="maximum number of points to be plotted for any scatter plots", type=int, default=100000)
    parser.add_argument("-report_name", help="the name of the output PDF report if the user wishes to use a non-default name", type=str, default="report.pdf")
    parser.add_argument("-only_png", help="save all figures in png format", action="store_true")
    parser.add_argument("-illumina", help="use this option to make figures look nicer with Illumina data", action="store_true")
    args = parser.parse_args()

    # Create variables of the same name
    bam = args.bam
    bai = args.bai
    asm = args.genome
    output_dir = args.output_dir
    verbose = args.verbose
    mapq_thres = args.mapq_thres
    len_min_contig = args.len_min_contig
    len_min_read = args.len_min_read
    len_min_hp = args.len_min_hp
    len_max_hp = args.len_max_hp
    len_context_ins = args.len_context_ins
    len_context_sub = args.len_context_sub
    len_max_indel = args.len_max_indel
    len_min_aln = args.len_min_aln
    len_trim_contig_edge = args.len_trim_contig_edge
    lim = args.lim
    num_pts_max = args.num_pts_max
    bitflag = args.bitflag
    report_name = args.report_name
    only_png = args.only_png
    cram = args.cram
    illumina = args.illumina

    if not args.no_figures:
        generate_figures = True
    else:
        generate_figures = False
    if args.use_recorded:
        correct_orientation = False
    else:
        correct_orientation = True

    assert bai is None or os.path.exists(bai), "The index file (.bai) was not found."

    # Check to make sure output_dir exists. If not, make one
    if not os.path.exists(output_dir):
        print("The output directory does not exist. Creating one.")
        os.mkdir(output_dir)
    else:
        if len(os.listdir(output_dir)) != 0:
            print("The output directory is not empty. Files will be overwritten.")
    # Make figures directory
    output_dir_figures = os.path.join(output_dir, "figures")
    if not os.path.exists(output_dir_figures):
        os.mkdir(output_dir_figures)
    # Make output stats directory
    output_dir_stats = os.path.join(output_dir, "stats")
    if not os.path.exists(output_dir_stats):
        os.mkdir(output_dir_stats)
    if generate_figures:
        # Blank PDF to collect all figures into a single file
        if not report_name.endswith(".pdf"):
            report_name = ".".join([report_name, "pdf"])
        report = PdfPages(os.path.join(output_dir, report_name))
    if verbose:
        start = time() # Measure the total time

    # ---- Read filter used throughout the program
    read_filter = make_read_filter(mapq_thres, len_min_read, len_min_aln, bitflag)


    # --- Load the data
    contigs, reads_pass, reads_fail = load_files(asm, bam, bai, len_min_contig, read_filter, lim, verbose, cram)


    # ---- Compute Q-score related statistics.
    # Mean/median/std Q-scores per read by pass/fail
    means_pass, meds_pass, stds_pass = mapQ_stats_per_read(reads_pass, verbose=verbose, comment="for pass reads")
    means_fail, meds_fail, stds_fail = mapQ_stats_per_read(reads_fail, verbose=verbose, comment="for fail reads")
    if generate_figures:
        plot_per_read_Q_stats(means_pass, meds_pass, stds_pass, means_fail, meds_fail, stds_fail, output_dir_figures, report=report, num_pts_max=num_pts_max, illumina=illumina)

    # Mean/std Q-score per read by aligned/unaligned region
    means_in, stds_in, lens_in, means_out, stds_out, lens_out = mapQ_stats_aligned_readsegment(reads_pass, verbose=verbose)
    if generate_figures:
        plot_per_read_Q_stats_aligned(means_in, stds_in, lens_in, means_out, stds_out, lens_out, output_dir_figures, report=report, num_pts_max=num_pts_max, illumina=illumina)


    # ---- Reconstruct all alignments that pass the read_filter defined above
    alns, lens = reconstruct_all_alignments(contigs, reads_pass, len_max_indel, correct_orientation, len_trim_contig_edge, verbose=verbose)


    # ---- Compute per-read error rates
    # Count the number of various errors in each read
    nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned = count_errors_per_read(alns, verbose=verbose)
    lens -= nums_skip # Exclude regions where reference is non-ACGT.
    if generate_figures:
        plot_per_read_error_stats(lens, lens_aligned, nums_match, nums_sub, nums_ins, nums_del, output_dir_figures, report=report, num_pts_max=num_pts_max)


    # ---- Save per-read statsitics
    save_per_read_stats_pass([means_pass, meds_pass, stds_pass, means_in, stds_in, lens_in, means_out, stds_out, lens_out, nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned, lens], output_dir_stats)
    save_per_read_stats_fail([means_fail, meds_fail, stds_fail], output_dir_stats)
    del reads_pass, reads_fail, means_pass, meds_pass, stds_pass, means_in, stds_in, lens_in, means_out, stds_out, lens_out, nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned, lens, means_fail, meds_fail, stds_fail

    # ---- PhredQ: Empirical vs. Computed Q
    qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q = phredQ_vs_error(alns, verbose=verbose)
    save_phredQ_stats([qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q], output_dir_stats)
    if generate_figures:
        plot_phredQ_stats(qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q, output_dir_figures, report=report, only_png=only_png)
    del qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q


    # ---- Homopolymer length histogram based only on the assembly
    hist_len_hp_asm = get_hist_len_hp(contigs, len_min_hp, len_max_hp, len_trim_contig_edge, verbose=verbose)
    save_hist_len_hp_asm(hist_len_hp_asm, len_max_hp, output_dir_stats)
    if generate_figures:
        plot_hist_len_hp_asm(hist_len_hp_asm, output_dir_figures, report=report, only_png=only_png)
    del hist_len_hp_asm

    # ---- Scrub the reconstructed alignments
    alns_scrubbed = scrub_reconstructed_alignments(alns, len_min_hp, verbose=verbose)
    del alns


    # ----- Compute errors in homopolymer regions
    hist_len_hp = count_errors_hp(alns_scrubbed, len_min_hp, len_max_hp, verbose=verbose)
    save_len_hp(hist_len_hp, output_dir_stats)
    if generate_figures:
        plot_dist_len_hp(hist_len_hp, output_dir_figures, report=report, only_png=only_png)
    del hist_len_hp


    # --- Context independent errors
    # Exluding hp regions
    sub_matrix, hist_len_del, hist_len_ins = count_errors(alns_scrubbed, len_min_hp=len_min_hp, exclude_hp = True, len_max_indel=len_max_indel, verbose=verbose)
    # Including all regions
    sub_matrix_all, hist_len_del_all, hist_len_ins_all = count_errors(alns_scrubbed, len_min_hp=len_min_hp, exclude_hp = False, len_max_indel=len_max_indel, verbose=verbose)
    # Only homopolymer region
    sub_matrix_hp = sub_matrix_all - sub_matrix
    hist_len_del_hp = hist_len_del_all - hist_len_del
    hist_len_ins_hp = hist_len_ins_all - hist_len_ins
    # Compute rate of indels
    R_ins, R_del = compute_indel_rates(sub_matrix, hist_len_del, hist_len_ins) # ex hp
    R_ins_all, R_del_all = compute_indel_rates(sub_matrix_all, hist_len_del_all, hist_len_ins_all) # all
    R_ins_hp, R_del_hp = compute_indel_rates(sub_matrix_hp, hist_len_del_hp, hist_len_ins_hp) # hp
    # Save
    save_context_independent_errors([sub_matrix, hist_len_del, hist_len_ins, R_ins, R_del], output_dir_stats, suffix="_ex_hp")
    save_context_independent_errors([sub_matrix_all, hist_len_del_all, hist_len_ins_all, R_ins_all, R_del_all], output_dir_stats, suffix="_all")
    save_context_independent_errors([sub_matrix_hp, hist_len_del_hp, hist_len_ins_hp, R_ins_hp, R_del_hp], output_dir_stats, suffix="_hp")

    if generate_figures:
        for (sub, R_i, R_d, hist_ins, hist_del, suffix) in [
            (sub_matrix, R_ins, R_del, hist_len_del, hist_len_ins, "_ex_hp"),
            (sub_matrix_all, R_ins_all, R_del_all, hist_len_del_all, hist_len_ins_all, "_all"),
            (sub_matrix_hp, R_ins_hp, R_del_hp, hist_len_del_hp, hist_len_ins_hp, "_hp")]:
            plot_sub_heatmap(sub, output_dir_figures, vmin=0., vmax=8, fname="sub_matrix%s.png" % suffix, report=report, only_png=only_png)
            plot_dist_indel(R_i, R_d, hist_del, hist_ins, output_dir_figures, fname="dist_indel%s.png" %suffix, report=report, only_png=only_png)

    # ---- Context dependent errors
    context_sub_matrix = context_sub_table(alns_scrubbed, len_context_sub=len_context_sub, verbose=verbose)
    context_hist_len_ins, context_hist_char_ins = context_ins_table(alns_scrubbed, len_context_ins=len_context_ins, len_max_ins=len_max_indel, verbose=verbose)
    save_context_dependent_errors([context_sub_matrix, context_hist_len_ins, context_hist_char_ins], output_dir_stats)
    tabulate_context_sub_matrix(context_sub_matrix, len_context_sub, output_dir)
    tabulate_context_ins_hist(context_hist_len_ins, context_hist_char_ins, len_context_ins, len_max_indel, output_dir)

    if generate_figures:
        # Close the report
        report.close()
    if verbose:    
        end = time()
        print("Total time taken: %.1f seconds" % (end-start))


if __name__ == "__main__":
    main()