"""
Given a set of output directories produced using the same set of options, the script outputs aggregate results.

Notes:
- The script assumes that all statistics have been computed.
- Certain figures are not included in reports.pdf.
"""
from counterr.util import *
from counterr.tables import *
from counterr.plots import *
from counterr.counters import compute_indel_rates
from counterr.IO import *

# ---- User must specify the output directories (out_dir/stats/) to be aggregated here.
output_dirs = []

# ---- User must specify the output driectory for the aggregate results.
aggregate_dir = ""

output_dir_figures = os.path.join(aggregate_dir, "figures")
output_dir_stats = os.path.join(aggregate_dir, "stats")
if not os.path.exists(aggregate_dir):
    os.mkdir(aggregate_dir)
if not os.path.exists(output_dir_figures):
    os.mkdir(output_dir_figures)
if not os.path.exists(output_dir_stats):
    os.mkdir(output_dir_stats)

# ---- PDF report output
report = PdfPages(os.path.join(aggregate_dir, "report.pdf"))


# ---- Placeholder for aggregate results
reads_pass = []
reads_fail = []
hist_len_hp_asm = None
hist_len_hp = None
phredQ_stats = None
sub_matrix_all = None
hist_len_del_all = None
hist_len_ins_all = None
sub_matrix_ex_hp = None
hist_len_del_ex_hp = None
hist_len_ins_ex_hp = None
sub_matrix_hp = None
hist_len_del_hp = None
hist_len_ins_hp = None
context_sub_matrix = None
context_hist_len_ins = None
context_hist_char_ins = None

for path in output_dirs:
    # per-read stats for pass reads
    data = np.load(os.path.join(path, "per_read_stats_pass.npy"))
    reads_pass.append(data)

    # per-read stats for fail reads
    data = np.load(os.path.join(path, "per_read_stats_fail.npy"))
    reads_fail.append(data)

    # Homopolymer length distribution (assembly)
    data = np.load(os.path.join(path, "hist_len_hp_asm.npy"))
    if hist_len_hp_asm is None:
        hist_len_hp_asm = data
    else:
        for ltr in bases_list:
            hist_len_hp_asm[ltr] += data[ltr]

    # Homopolymer true vs. obs length distribution
    data = np.load(os.path.join(path, "hist_len_hp.npy"))
    if hist_len_hp is None:
        hist_len_hp = data
    else:
        hist_len_hp += data

    # Phred Q-score stats
    data = np.load(os.path.join(path, "phredQ_stats.npy"))
    if phredQ_stats is None:
        phredQ_stats = data
    else:
        for col in ["num_match", "num_sub", "num_ins", "num_err", "num_tot"]:
            phredQ_stats[col] += data[col]

    # Context independent errors -- all
    data = np.load(os.path.join(path, "context_independent_errors_all.npz"))
    if sub_matrix_all is None:
        sub_matrix_all = data["sub_matrix"]
        hist_len_del_all = data["hist_len_del"]
        hist_len_ins_all = data["hist_len_ins"]
    else:
        sub_matrix_all += data["sub_matrix"]
        hist_len_del_all += data["hist_len_del"]
        hist_len_ins_all += data["hist_len_ins"]

    # Context independent errors -- ex_hp
    data = np.load(os.path.join(path, "context_independent_errors_ex_hp.npz"))
    if sub_matrix_ex_hp is None:
        sub_matrix_ex_hp = data["sub_matrix"]
        hist_len_del_ex_hp = data["hist_len_del"]
        hist_len_ins_ex_hp = data["hist_len_ins"]
    else:
        sub_matrix_ex_hp += data["sub_matrix"]
        hist_len_del_ex_hp += data["hist_len_del"]
        hist_len_ins_ex_hp += data["hist_len_ins"]

    # Context independent errors -- hp
    data = np.load(os.path.join(path, "context_independent_errors_hp.npz"))
    if sub_matrix_hp is None:
        sub_matrix_hp = data["sub_matrix"]
        hist_len_del_hp = data["hist_len_del"]
        hist_len_ins_hp = data["hist_len_ins"]
    else:
        sub_matrix_hp += data["sub_matrix"]
        hist_len_del_hp += data["hist_len_del"]
        hist_len_ins_hp += data["hist_len_ins"]

    # Context dependent errors
    data = np.load(os.path.join(path, "context_dependent_errors.npz"))
    if context_sub_matrix is None:
        context_sub_matrix = data["sub_matrix"]
        context_hist_len_ins = data["hist_len_ins"] 
        context_hist_char_ins = data["hist_char_ins"]
    else:
        context_sub_matrix += data["sub_matrix"]
        context_hist_len_ins += data["hist_len_ins"] 
        context_hist_char_ins += data["hist_char_ins"]


# ---- Unpack/compute stats based on the aggregated results. suffix "_" is for unpacked.
# Pass reads
reads_pass = np.hstack(reads_pass)
means_pass, meds_pass, stds_pass, means_in, stds_in, lens_in, means_out, stds_out, lens_out, nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned, lens = reads_pass["mean"], reads_pass["med"], reads_pass["std"], reads_pass['mean_in'], reads_pass['std_in'], reads_pass['len_in'], reads_pass['mean_out'],reads_pass['std_out'], reads_pass['len_out'], reads_pass['num_match'], reads_pass['num_sub'], reads_pass['num_ins'], reads_pass['num_del'], reads_pass['num_skip'], reads_pass['len_aligned'], reads_pass['len']

# Fail reads
reads_fail = np.hstack(reads_fail)
means_fail, meds_fail, stds_fail = reads_fail["mean"], reads_fail["med"], reads_fail["std"]

# Homopolymer length in assembly
hist_len_hp_asm_ = np.zeros((4, hist_len_hp_asm["length"].size))
for ltr in bases_list:
    hist_len_hp_asm_[base2int[ltr]] = hist_len_hp_asm[ltr]

# Q-statas
qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q = phredQ_stats["q"], phredQ_stats["num_match"], phredQ_stats["num_sub"], phredQ_stats["num_ins"], phredQ_stats["num_err"], phredQ_stats["num_tot"]

# Compute indel rates
R_ins_all, R_del_all = compute_indel_rates(sub_matrix_all, hist_len_del_all, hist_len_ins_all)
R_ins_ex_hp, R_del_ex_hp = compute_indel_rates(sub_matrix_ex_hp, hist_len_del_ex_hp, hist_len_ins_ex_hp)
R_ins_hp, R_del_hp = compute_indel_rates(sub_matrix_hp, hist_len_del_hp, hist_len_ins_hp)

# len_context
len_context_sub = int(math.log(context_sub_matrix.shape[0], N_bases))
len_context_ins = int(math.log(context_hist_len_ins.shape[0], N_bases))
len_max_indel = int(context_hist_len_ins.shape[1])


# ---- Plots and figures based on the aggregated results
plot_per_read_Q_stats(means_pass, meds_pass, stds_pass, means_fail, meds_fail, stds_fail, output_dir_figures, report=None)
plot_per_read_Q_stats_aligned(means_in, stds_in, lens_in, means_out, stds_out, lens_out, output_dir_figures, report=None)
plot_per_read_error_stats(lens, lens_aligned, nums_match, nums_sub, nums_ins, nums_del, output_dir_figures, report=None)
plot_phredQ_stats(qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q, output_dir_figures, report=report)
plot_hist_len_hp_asm(hist_len_hp_asm_, output_dir_figures, report=report)
plot_dist_len_hp(hist_len_hp, output_dir_figures, report=report)
plot_sub_heatmap(sub_matrix_all, output_dir_figures, vmin=0., vmax=8, fname="sub_matrix_all.png", report=report)
plot_dist_indel(R_ins_all, R_del_all, hist_len_del_all, hist_len_ins_all, output_dir_figures, fname="dist_indel_all.png", report=report)
plot_sub_heatmap(sub_matrix_ex_hp, output_dir_figures, vmin=0., vmax=8, fname="sub_matrix_ex_hp.png", report=report)
plot_dist_indel(R_ins_ex_hp, R_del_ex_hp, hist_len_del_ex_hp, hist_len_ins_ex_hp, output_dir_figures, fname="dist_indel_ex_hp.png", report=report)
plot_sub_heatmap(sub_matrix_hp, output_dir_figures, vmin=0., vmax=8, fname="sub_matrix_hp.png", report=report)
plot_dist_indel(R_ins_hp, R_del_hp, hist_len_del_hp, hist_len_ins_hp, output_dir_figures, fname="dist_indel_hp.png", report=report)
tabulate_context_sub_matrix(context_sub_matrix, len_context_sub, aggregate_dir, ex_hp=True, N_top_bottom=5, suffix="_5")
tabulate_context_ins_hist(context_hist_len_ins, context_hist_char_ins, len_context_ins, len_max_indel, aggregate_dir, ex_hp=True, N_top_bottom=5, suffix="_5")
tabulate_context_sub_matrix(context_sub_matrix, len_context_sub, aggregate_dir)
tabulate_context_ins_hist(context_hist_len_ins, context_hist_char_ins, len_context_ins, len_max_indel, aggregate_dir)


# ---- Save the aggregated results
save_per_read_stats_pass([means_pass, meds_pass, stds_pass, means_in, stds_in, lens_in, means_out, stds_out, lens_out, nums_match, nums_sub, nums_ins, nums_del, nums_skip, lens_aligned, lens], output_dir_stats)
save_per_read_stats_fail([means_fail, meds_fail, stds_fail], output_dir_stats)
save_phredQ_stats([qs, nums_match_q, nums_sub_q, nums_ins_q, nums_err_q, nums_tot_q], output_dir_stats)
save_len_hp(hist_len_hp, output_dir_stats)
save_hist_len_hp_asm(hist_len_hp_asm_, hist_len_hp_asm_.shape[1], output_dir_stats)
save_context_independent_errors([sub_matrix_all, hist_len_del_all, hist_len_ins_all, R_ins_all, R_del_all], output_dir_stats, suffix="_all")
save_context_independent_errors([sub_matrix_ex_hp, hist_len_del_ex_hp, hist_len_ins_ex_hp, R_ins_ex_hp, R_del_ex_hp], output_dir_stats, suffix="_ex_hp")
save_context_independent_errors([sub_matrix_hp, hist_len_del_hp, hist_len_ins_hp, R_ins_hp, R_del_hp], output_dir_stats, suffix="_hp")
save_context_dependent_errors([context_sub_matrix, context_hist_len_ins, context_hist_char_ins], output_dir_stats)


# ---- Close the report
report.close()