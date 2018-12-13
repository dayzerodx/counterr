from .util import *
ft_size = 15

def plot_per_read_Q_stats(means_pass, meds_pass, stds_pass, means_fail, meds_fail, stds_fail, output_dir, report=None, num_pts_max=50000, illumina=False):
    # -- Q mean/med distribution
    fig_name = os.path.join(output_dir, "per_read_Q_mean_med_pass_vs_fail.png")
    if illumina:
        bins = np.arange(-0.5, 50, 1.)
    else:
        bins = np.arange(-0.5, 30.5, 1.)        
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    # pass
    ax1.hist(means_pass, bins=bins, histtype="step", color="black", label="pass/mean", lw=1.5, normed=True)
    ax1.hist(meds_pass, bins=bins, histtype="step", color="red", label="pass/med", lw=1.5, normed=True)
    ax1.set_title("Pass reads", fontsize=ft_size)
    ax1.set_xlabel("Q-score stat", fontsize=ft_size)
    if illumina:
        ax1.set_xlim([0, 50])
    else:
        ax1.set_xlim([0, 30])        
    ax1.legend(loc="upper right", fontsize=ft_size)
    # fail
    ax2.hist(means_fail, bins=bins, histtype="step", color="black", label="fail/mean", lw=1.5, normed=True)
    ax2.hist(meds_fail, bins=bins, histtype="step", color="red", label="fail/med", lw=1.5, normed=True)
    ax2.set_title("Fail reads", fontsize=ft_size)
    ax2.set_xlabel("Q-score stat", fontsize=ft_size)    
    if illumina:
        ax2.set_xlim([0, 50])
    else:
        ax2.set_xlim([0, 30])
    ax2.legend(loc="upper right", fontsize=ft_size)
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # -- Q mean/med distribution -- II
    fig_name = os.path.join(output_dir, "per_read_Q_pass_fail_mean_vs_med.png")
    if illumina:
        bins = np.arange(-0.5, 50, 1.)
    else:
        bins = np.arange(-0.5, 30.5, 1.)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    # pass
    ax1.hist(means_pass, bins=bins, histtype="step", color="black", label="pass/mean", lw=1.5, normed=True)
    ax1.hist(means_fail, bins=bins, histtype="step", color="red", label="fail/mean", lw=1.5, normed=True)
    ax1.set_title("Mean", fontsize=ft_size)
    ax1.set_xlabel("Q-score stat", fontsize=ft_size)        
    if illumina:
        ax1.set_xlim([0, 50])
    else:
        ax1.set_xlim([0, 30])                
    ax1.legend(loc="upper right", fontsize=ft_size)
    # fail
    ax2.hist(meds_pass, bins=bins, histtype="step", color="black", label="pass/med", lw=1.5, normed=True)
    ax2.hist(meds_fail, bins=bins, histtype="step", color="red", label="fail/med", lw=1.5, normed=True)
    ax2.set_title("Med", fontsize=ft_size)
    ax2.set_xlabel("Q-score stat", fontsize=ft_size)    
    if illumina:
        ax2.set_xlim([0, 50])
    else:
        ax2.set_xlim([0, 30])     
    ax2.legend(loc="upper right", fontsize=ft_size)
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # -- Q med vs std.
    means_pass, stds_pass = uniform_downsample([means_pass, stds_pass], num_pts_max)
    means_fail, stds_fail = uniform_downsample([means_fail, stds_fail], num_pts_max)

    fig_name = os.path.join(output_dir, "per_read_Q_pass_fail_mean_vs_std.png")
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 7))
    # pass
    ax.scatter(means_pass, stds_pass, color="black", label="pass", s=2, edgecolors="none")
    ax.scatter(means_fail, stds_fail, color="red", label="fail", s=2, edgecolors="none")
    ax.set_xlabel("Mean", fontsize=ft_size)
    ax.set_ylabel("Std", fontsize=ft_size)
    lgnd = plt.legend(loc="upper right", fontsize=ft_size, numpoints=1)
    #change the marker size manually
    lgnd.legendHandles[0]._sizes = [30]
    lgnd.legendHandles[1]._sizes = [30]    
    if illumina:
        ax.set_xlim([0, 50])
        ax.set_ylim([0, 25])
    else:
        ax.set_xlim([0, 30])
        ax.set_ylim([0, 15])
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    return 

def plot_per_read_Q_stats_aligned(means_in, stds_in, lens_in, means_out, stds_out, lens_out, output_dir, report=None, num_pts_max=50000, illumina=False):
    # -- Q med vs std.
    means_out, stds_out = uniform_downsample([means_out, stds_out], num_pts_max)
    means_in, stds_in = uniform_downsample([means_in, stds_in], num_pts_max)    
    fig_name = os.path.join(output_dir, "per_read_in_or_out_align_mean_vs_std.png")
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 7))
    # pass
    ax.scatter(means_out, stds_out, color="red", label="Unaligned", s=2, edgecolors="none")
    ax.scatter(means_in, stds_in, color="black", label="Aligned", s=2, edgecolors="none")
    ax.set_xlabel("Mean", fontsize=ft_size)
    ax.set_ylabel("Std", fontsize=ft_size)
    # ax.set_title("Mean", fontsize = 10)
    lgnd = plt.legend(loc="upper right", fontsize=ft_size, numpoints=1)
    #change the marker size manually
    lgnd.legendHandles[0]._sizes = [30]
    lgnd.legendHandles[1]._sizes = [30]    
    if illumina:
        ax.set_xlim([0, 50])
        ax.set_ylim([0, 25])
    else:
        ax.set_xlim([0, 30])
        ax.set_ylim([0, 15])
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # -- Length distribution of aligned and unaligned regions
    plt.close()    
    fig_name = os.path.join(output_dir, "per_read_dist_len_in_or_out_align.png")
    fig, ax = plt.subplots(1, figsize=(10, 7))
    len_max_in, len_max_out = 0, 0
    if len(lens_in) > 0:
        len_max_in = max(lens_in)
    if len(lens_out) > 0:
        len_max_out = max(lens_out)
    len_max = max(len_max_in, len_max_out)*1.05
    if len_max > 10: # Otherwise what's the point?
        bins = np.arange(0, len_max, len_max/100.)        
        ax.hist(lens_in, bins=bins, color="black", label="Aligned", histtype="step", normed=True)
        ax.hist(lens_out, bins=bins, color="red", label="Unaligned", histtype="step", normed=True)
        ax.set_xlabel("Length", fontsize=ft_size)
        ax.set_xlim([0, len_max])
        ax.legend(loc="upper right", fontsize=ft_size)    
        plt.savefig(fig_name, dpi=200, bbox_inches="tight")
        if report is not None:
            plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
            report.savefig(fig, dpi=200)

        plt.close()

def plot_per_read_error_stats(lens, lens_aligned, nums_match, nums_sub, nums_ins, nums_del, output_dir, report=None, num_pts_max=50000):
    len_max_lim = np.percentile(lens, 99.9) * 1.1
    w = (np.percentile(lens, 90) - np.percentile(lens, 10))/50.
    lens = lens.astype(float) + 1e-6
    lens_aligned = lens_aligned.astype(float) + 1e-6

    # -- Length distribution
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 5))
    ax.hist(lens, bins=np.arange(0, len_max_lim, w), histtype="step", lw=1, label="Total", color="black")
    ax.hist(lens_aligned, bins=np.arange(0, len_max_lim, w), histtype="step", lw=1, label="Aligned", color="red")
    fig_name = os.path.join(output_dir, "per_read_dist_len.png")
    ax.legend(loc="upper right", fontsize=ft_size)
    ax.set_xlabel("Length", fontsize=ft_size)
    ax.set_ylabel("Counts", fontsize=ft_size)
    ax.set_xlim([0, len_max_lim])
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # ---- accuracy and error histogram
    plt.close()
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))
    bins = np.arange(0, 1.005, 0.005)
    ax1.hist(nums_match/lens_aligned, color="green", histtype="step", label="Match", lw=1, bins=bins)
    ax1.hist(nums_sub/lens_aligned, color="red", histtype="step", label="Sub.", lw=1, bins=bins)
    ax1.hist(nums_del/lens_aligned, color="orange", histtype="step", label="Deletion", lw=1, bins=bins)
    ax1.hist(nums_ins/lens_aligned, color="blue", histtype="step", label="Insertion", lw=1, bins=bins)
    ax1.legend(loc="upper right", fontsize=ft_size)
    ax1.set_xlabel("Proportion", fontsize=ft_size)
    ax1.set_ylabel("Counts", fontsize=ft_size)
    ax1.set_xlim([0, 1])
    fig_name = os.path.join(output_dir, "per_read_hist_errors.png")
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()    

    lens, lens_aligned, nums_match, nums_sub, nums_ins, nums_del = uniform_downsample([lens, lens_aligned, nums_match, nums_sub, nums_ins, nums_del], num_pts_max)

    # -- Length vs. aligned length
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 7))
    ax.scatter(lens, lens_aligned, c="black", edgecolors="none", s=1)
    ax.plot([0, len_max_lim], [0, len_max_lim], ls="--", lw=0.5, c="red")
    fig_name = os.path.join(output_dir, "per_read_len_vs_len_aligned.png")
    # ax.legend(loc="upper right", fontsize=ft_size)
    ax.set_xlabel("Length", fontsize=ft_size)
    ax.set_ylabel("Length aligned", fontsize=ft_size)
    ax.set_xlim([0, len_max_lim])
    ax.set_ylim([0, len_max_lim])
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    lens_ratio = lens_aligned/lens.astype(float)
    # -- Length vs. aligned length/length
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 7))
    ax.scatter(lens, lens_ratio, c="black", edgecolors="none", s=1)
    ax.plot([0, len_max_lim], [1, 1], ls="--", lw=0.5, c="red")
    fig_name = os.path.join(output_dir, "per_read_len_vs_len_aligned_div_len.png")
    # ax.legend(loc="upper right", fontsize=ft_size)
    ax.set_xlabel("Length", fontsize=ft_size)
    ax.set_ylabel("Length aligned / length", fontsize=ft_size)
    ax.set_xlim([0, len_max_lim])
    ax.set_ylim([0, np.max(lens_ratio)])
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # -- Length vs. error
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 7))
    ax.scatter(lens, nums_match/lens_aligned, c="green", edgecolors="none", s=1.5, label="Match")
    ax.scatter(lens, nums_del/lens_aligned, c="orange", edgecolors="none", s=1.5, label="Deletion")
    ax.scatter(lens, nums_ins/lens_aligned, c="blue", edgecolors="none", s=1.5, label="Insertion")
    ax.scatter(lens, nums_sub/lens_aligned, c="red", edgecolors="none", s=1.5, label="Sub.")
    ax.plot([0, len_max_lim], [0, len_max_lim], ls="--", lw=0.5, c="red")
    fig_name = os.path.join(output_dir, "per_read_len_vs_error_div_len_aligned.png")
    lgnd = ax.legend(loc="center right", fontsize=ft_size, scatterpoints=1)
    lgnd.legendHandles[0]._sizes = [20]
    lgnd.legendHandles[1]._sizes = [20]
    lgnd.legendHandles[2]._sizes = [20]
    lgnd.legendHandles[3]._sizes = [20]
    ax.set_xlabel("Length", fontsize=ft_size)
    ax.set_ylabel("(# occurrence) / length aligned", fontsize=ft_size)
    ax.set_xlim([0, len_max_lim])
    ax.set_ylim([0, 1.])
    plt.savefig(fig_name, dpi=400, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)
    plt.close()

    return

def plot_phredQ_stats(qs, nums_match, nums_sub, nums_ins, nums_err, nums_tot, output_dir, report=None, only_png=False):
    # Plot only that has at least one observation
    iselect = nums_tot > 0 
    qs = qs[iselect]
    if qs.size == 0: # If it appears that phred-Q scores are not available, then terminate early.
        return

    q_max = np.max(qs) + 3
    nums_match = nums_match[iselect]
    nums_sub = nums_sub[iselect]
    nums_ins = nums_ins[iselect]
    nums_err = nums_err[iselect]
    nums_tot = nums_tot[iselect].astype(float)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # -- Q-score vs. error rates
    ax1.scatter(qs, nums_match / nums_tot, c="green", s=10, edgecolors="none", label="Match")
    ax1.scatter(qs, nums_sub / nums_tot, c="red", s=10, edgecolors="none", label="Sub.")
    ax1.scatter(qs, nums_ins / nums_tot, c="blue", s=10, edgecolors="none", label="Insertion")
    ax1.scatter(qs, nums_err / nums_tot, c="black", s=10, edgecolors="none", label="Error")
    ax1.set_xlabel("Computed Phred Q", fontsize=15)
    ax1.set_ylim([0, 1])
    ax1.set_xlim([0, q_max])
    ax1.set_ylabel("Rate", fontsize=15)
    ax1.legend(loc="upper right", fontsize=10)
    # Predicted Q-score vs. empirical
    ax2.scatter(qs, match_rate2phredQ(nums_match/ nums_tot), s=10, c="black")
    ax2.plot([0, q_max], [0, q_max], ls="--", c="red", lw=1.)
    ax2.set_xlabel("Computed Phred Q", fontsize=15)
    ax2.set_ylabel("Empirical Phred Q", fontsize=15)
    ax2.set_ylim([0, q_max])
    ax2.set_xlim([0, q_max])
    fig_name = os.path.join(output_dir, "phredQ_vs_error.pdf")
    fig_name = change_ext_to_png(fig_name, only_png)    
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    # Histgoram of q-score distribution
    plt.close()
    fig, ax = plt.subplots(1, figsize=(10, 5))
    # -- Q-score vs. error rates
    ax.bar(qs, nums_tot)
    ax.set_xlabel("Phred Q", fontsize=15)
    ax.set_ylabel("Counts", fontsize=15)
    ax.set_xlim([0, q_max])
    fig_name = os.path.join(output_dir, "phredQ.pdf")    
    fig_name = change_ext_to_png(fig_name, only_png)
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    return

def plot_hist_len_hp_asm(len_hist_hp, output_dir, report=None, only_png=False):
    len_max = len_hist_hp.shape[1]
    plt.close()
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for idx, ltr in enumerate(bases_list):
        ax.plot(np.arange(len_max), len_hist_hp[idx], lw=1., marker="o", label=ltr)
    ax.set_xlabel("Homopolymer length")
    ax.legend(loc="upper right", fontsize=ft_size)
    ax.set_xlim([0, len_max])
    ax.set_ylim([0, np.max(len_hist_hp) * 1.05])
    ax.set_xticks(np.arange(len_max))
    ax.set_ylabel("Counts")
    fig_name = os.path.join(output_dir, "asm_hist_len_hp.pdf")
    fig_name = change_ext_to_png(fig_name, only_png)
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    return

def plot_sub_heatmap(sub_matrix, output_dir, vmin=0., vmax=10., fname="sub_table.pdf", report=None, only_png=False):
    """
    Note that units are in percentage.
    """
    # Heatmap of substitution matrix
    mask = np.zeros((5, 5), dtype=bool)
    mask[-1, -1] = True
    plt.close()
    fig, ax = plt.subplots(1, figsize=(5, 5))
    matrix = sub_matrix/sub_matrix.sum(axis=1).reshape(5, 1).astype(float) * 100
    ax = sns.heatmap(matrix, vmax=vmax, vmin=vmin, fmt="2.2f", annot=True, xticklabels=bases_ext_list, yticklabels=bases_ext_list, cmap="Reds", mask=mask)
    ax.set_ylabel("Truth", fontsize=ft_size)
    ax.set_xlabel("Observed", fontsize=ft_size)
    fig_name = os.path.join(output_dir, fname)
    fig_name = change_ext_to_png(fig_name, only_png)
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    return

def plot_dist_len_hp(dist_len_hp, output_dir, len_min=3, len_max=11, report=None, only_png=False):
    assert (len_max - len_min) == 8, "The range must be 8."
    # Lengths distributions
    plt.close()
    fig, ax_list = plt.subplots(4, 2, figsize=(17, 17))
    length = np.arange(dist_len_hp.shape[2])
    for i, len_hp in enumerate(range(len_min, len_max)):
        idx_col = i // 4
        idx_row = i % 4
        for ltr in bases_list: # For each letter homopolymer        
            counts = dist_len_hp[base2int[ltr], len_hp]
            count_total = counts.sum()
            counts_normed = counts / (count_total + 1e-6)
            ax_list[idx_row, idx_col].plot(length, counts_normed, marker="o", label="%s: %d" % (ltr, count_total), ls="--", lw=1)
        ax_list[idx_row, idx_col].axvline(x=len_hp, ls="--", c="black", lw=1, label="True length: %d" % len_hp)
        ax_list[idx_row, idx_col].legend(loc="upper right", fontsize=13)
        ax_list[idx_row, idx_col].set_xlim([0, 15])
        ax_list[idx_row, idx_col].set_ylim([0, 1])
        ax_list[idx_row, idx_col].set_xticks(range(15))
        ax_list[idx_row, idx_col].set_xlabel("Observed homopolymer length", fontsize=15)
    fig_name = os.path.join(output_dir, "dist_len_hp.pdf")
    fig_name = change_ext_to_png(fig_name, only_png)    
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

def plot_dist_indel(R_ins, R_del, dist_len_del, dist_len_ins, output_dir, fname="dist_indel.pdf", report=None, only_png=False):
    plt.close()
    # Normalize to unit-norm
    dist_len_del = dist_len_del[:]/float(dist_len_del.sum())
    dist_len_ins = dist_len_ins[:]/float(dist_len_ins.sum())

    length = np.arange(1, dist_len_del.size)
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.plot(length, dist_len_ins[1:], marker="o", ls="--", lw=1, color="black", label="Ins: %.3f%%" % (R_ins * 100))
    ax.plot(length, dist_len_del[1:], marker="o", ls="--", lw=1, color="red", label="del: %.3f%%" % (R_del * 100))
    ax.legend(loc="upper right", fontsize=ft_size)
    ax.set_xticks(range(1, 10))
    ax.set_xlim([1, 10])
    ax.set_xlabel("Length of indel", fontsize=ft_size)
    ax.set_ylabel("Probability", fontsize=ft_size)
    max_count = max(np.max(dist_len_del), np.max(dist_len_ins))
    ax.set_ylim([0, 1.0])
    ax.set_title("Ins/del: rate per match", fontsize=ft_size)
    fig_name = os.path.join(output_dir, fname)
    fig_name = change_ext_to_png(fig_name, only_png)    
    plt.savefig(fig_name, dpi=200, bbox_inches="tight")
    if report is not None:
        plt.suptitle(fig_name.split("/")[-1].split(".")[0], fontsize=20)
        report.savefig(fig, dpi=200)

    plt.close()

    return None
