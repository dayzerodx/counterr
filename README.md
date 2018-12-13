# Counterr
Counterr is a light-weight command line tool that computes errors in sequencing data by comparing the reads to a reference genome. Counterr takes as input an alignment of the reads (bam + bai) and the corresponding reference (fasta) and outputs summaries (figures/tables) of errors in the reads. The tool computes both context independent and context dependent error statistics of the reads. The latter can be particularly helpful for uncovering systematic errors latent in the sequencing technology (due to hardware, basecaller, etc.).

The tool was developed with Oxford Nanopore Technology reads in mind but is general and applicable to other sequencing platforms such as Illumina and PacBio.

## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Quick Examples](#quick-examples)
4. [Full Usage](#full-usage)
5. [Outputs](#outputs)
6. [Program Description](#program-description)
7. [Caveats](#caveats)
8. [Miscellaneous](#miscellaneous)
9. [Known Issues and Work-arounds](#known-issues-and-work-arounds)
10. [Contact Information](#contact-information)
11. [License](#license)

## Requirements
- Python 2.7/Python 3.4 or later
- pysam>=0.14.1
- seaborn>=0.9.0
- matplotlib>=2.2.3
- numpy>=1.15.4
- pandas>=0.23.4

The user must manually install these packages before installing counterr as below. 

Notes:
- For installing pysam, please refer to its [documentation](https://pysam.readthedocs.io/en/latest/installation.html). 
- matplotlib, numpy, and pandas are prerequisites for seaborn. If you use pip, those packages will be automatically installed.
- If you use conda, you may want to create a new virtual environment in order to avoid package version collisions.

## Installation
```
git clone git@github.com:dayzerodx/counterr.git
cd counterr
python setup.py install
counterr -h
```

To use the latest version of counterr, you can either remove the counterr directory and obtain the latest repo as above or move to the root directory and use the following commands:
```
git pull
rm -r dist/ build/ counterr.egg-info/
python setup.py install
```

## Quick Examples
**Basic use case with terminal outputs**
```
counterr -bam file.bam -genome assembly.fa -output_dir output -verbose
```
**Use only 100 reads randomly selected from the bam file (useful for testing)**
```
counterr -bam file.bam -genome assembly.fa -output_dir output -lim 100
```
**Apply custom read filtering conditions**
```
counterr -bam file.bam -genome assembly.fa -output_dir output -mapq_thres 20 -len_min_aln 2000 -len_min_read 2000
```

## Full Usage
```
usage: counterr [-h] -bam BAM -genome GENOME -output_dir OUTPUT_DIR
                 [-no_figures] [-bai BAI] [-verbose]
                 [-len_min_contig LEN_MIN_CONTIG] [-mapq_thres MAPQ_THRES]
                 [-len_min_read LEN_MIN_READ] [-len_min_aln LEN_MIN_ALN]
                 [-bitflag BITFLAG] [-len_min_hp LEN_MIN_HP]
                 [-len_max_hp LEN_MAX_HP] [-len_context_sub LEN_CONTEXT_SUB]
                 [-len_context_ins LEN_CONTEXT_INS]
                 [-len_max_indel LEN_MAX_INDEL]
                 [-len_trim_contig_edge LEN_TRIM_CONTIG_EDGE] [-use_recorded]
                 [-lim LIM]

optional arguments:
  -h, --help            show this help message and exit
  -bam BAM              the input bam file (default: None)
  -genome GENOME        the input fasta file (default: None)
  -output_dir OUTPUT_DIR
                        the output directory for figures and stats (default:
                        None)
  -no_figures           pass this flag to not generate figures (default:
                        False)
  -bai BAI              the input bai filename if non-conventionally named
                        (default: )
  -verbose              pass this flag to follow progress in the terminal
                        (default: False)
  -len_min_contig LEN_MIN_CONTIG
                        minimum contig length (default: 1500)
  -mapq_thres MAPQ_THRES
                        minimum mapq threshold (default: 40)
  -len_min_read LEN_MIN_READ
                        minimum read length (default: 1500)
  -len_min_aln LEN_MIN_ALN
                        minimum length aligned (default: 1000)
  -bitflag BITFLAG      bit flag for read filter (as specified in SAM doc)
                        (default: 3845)
  -len_min_hp LEN_MIN_HP
                        minimum homopolymer length (default: 3)
  -len_max_hp LEN_MAX_HP
                        maximum homopolymer length (default: 20)
  -len_context_sub LEN_CONTEXT_SUB
                        length of the k-mer context for context dependent
                        substitution table (default: 5)
  -len_context_ins LEN_CONTEXT_INS
                        length of the k-mer context for context dependent
                        insertion table (default: 6)
  -len_max_indel LEN_MAX_INDEL
                        maximum length of indels to consider (default: 20)
  -len_trim_contig_edge LEN_TRIM_CONTIG_EDGE
                        length of the contig edge to trim before computing
                        various statistics (default: 1)
  -use_recorded         pass this flag to NOT perform reverse complementing of
                        the reverse complement mapped reads (default: False)
  -lim LIM              pass this flag to run the program with 'lim' randomly
                        selected reads (both pass and fail) (default: -1)
  -num_pts_max NUM_PTS_MAX
                        maximum number of points to be plotted for any scatter
                        plots (default: 50000)
  -report_name REPORT_NAME
                        the name of the output PDF report if the user wishes
                        to use a non-default name. (default: report.pdf)                        
```

## Outputs
Example outputs are available in "example_outdir" directory. 

### Figures
All generated figures are collected into a single PDF file (default name is “reports.pdf”). The individual figures are saved under the "figures" sub-directory.
- per_read_Q_mean_med_pass_vs_fail.png: Unit-normalized density histogram of per-read mean/median Phred-Q scores grouped by read filter status. 
- per_read_Q_pass_fail_mean_vs_med.png: Unit-normalized density histogram of per-read mean/median Phred-Q scores grouped by the statistics.
- per_read_Q_pass_fail_mean_vs_std.png: Scatter plot of per-read Phred-Q score mean vs. std grouped by read filter status.
- per_read_in_or_out_align_mean_vs_std.png: Scatter plot of per-read Phred-Q score mean vs. std grouped by aligned vs. unaligned regions.
- per_read_dist_len.png: Histograms of the original read length and the length of the aligned portion of the read.
- per_read_dist_len_in_or_out_align.png: Unit-normalized density histogram of read length grouped by aligned vs. unaligned regions.
- per_read_len_vs_len_aligned.png: Scatter plot of recorded read length vs. **alignment length** defined as the sum of matches, mis-matches, insertions, and deletions in the alignment. Any regions with non-ACGT letters in the assembly are skipped and not counted in either of the lengths.
- per_read_len_vs_len_aligned_div_len.png: Scatter plot of recorded read length vs. alignment length divided by the former.
- per_read_len_vs_error_div_len_aligned.png: Scatter plot of original read length vs. per-read error rates grouped by error type.
- per_read_hist_errors.png: Histogram of per-read error rate grouped by error type.
- phredQ.pdf: Histogram of Phred-Q scores over all reads.
- phredQ_vs_error.pdf: (left) Computed Phred-Q score vs. error rates. (right) Computed Phred-Q score vs. empirical Phred-Q score. The former is what is recorded in the fastq files of the reads and the latter is based on the alignments.
- asm_hist_len_hp.pdf: Histogram of homopolymer length in the reference genome grouped by DNA base.
- dist_len_hp.pdf: Histogram of observed homopolymer length grouped by the true length and DNA base.
- dist_indel_**.pdf: Distribution of insertions and deletions and their rates (computed as the ratio of the number of insertions/deletions to the number of matches or mismatches).
- sub_matrix_**.pdf: Substitution matrices including insertions and deletions.

The suffixes for the last two figures take one of the following values: 
- all: No region in the assembly is excluded when computing various statistics.
- hp: Only homopolymer regions in the assembly are considered when computing various statistics.
- ex_hp: Excludes homopolymer regions in the assembly when computing various statistics.

### Tables
- context_ins.tsv: Context-dependent insertion table with each row corresponding to a particular k-mer context (e.g., ACTTCA). Ordered by error rate from top to bottom. If a k-mer is not observed, then it is excluded from the table.
- context_sub.tsv: Context-dependent substitution table with each row corresponding to a particular k-mer context (e.g., ACGCA). Ordered by error rate from top to bottom. If a k-mer is not observed, then it is excluded from the table.

### Data
Various computed statistics are saved in the "stats" subdirectory. These files are necessary for aggregating results from multiple runs (see "Miscellaneous" section). At the moment, there is no complete documentation for the files.

## Program Description
Given input data, counterr progresses in a linear fashion, executing the following steps:
- Load all contigs and their corresponding PySam read objects and store them in a Python list and dictionary, respectively. Group reads into "pass" and "fail" piles based on the read filter chosen by the user or following the default settings.
- Compute *per-read* Phred quality score statistics.
- Create an "alignment reconstruction” for each "pass" read and the contig to which it aligns. The alignment reconstruction is a lightweight string representation of the alignment. Creating the reconstruction entails inserting “-“’s in appropriate places within the read and reference strings according to the alignment CIGAR string. Here is an example of an alignment reconstruction:
```
Ref  : AGCT--GTCA--AAACCC
Read : AGCATTG-CACCAAAGCC
CIGAR: ===XIIMD==II===X==
PhQ  : 333373303335333332
```
The reconstruction structure makes it easier to write and debug routines that count various context independent and dependent errors. By default, for reads that mapped in the reverse orientation to the reference, the read sequence recorded in the bam is reverse complemented in the alignment reconstruction. This step is necessary since for reads that map to the reverse strand of the reference, the bam file typically records the reverse complement of the read and not the read in its original orientation. The user may disable this default setting by passing -use_recorded flag but this will likely result in an incorrect error profile characterization.

- Based on alignment reconstructions, count errors (i.e., match, mismatch, insertion, and deletion) in each read, skipping over non-ACGT letter regions in the reference.
- Compute *aggregate* Phred quality score statistics.
- Compute the histogram of homopolymer length grouped by DNA base in the reference.
- "Scrub" the alignment reconstructions: 1) Split alignments at non-ACGT sites and 2) add a "homopolymer string" to each reconstructed alignment. The homopolymer string is used to mark homopolymer regions. After scrubbing, the example above becomes
```
Ref  : AGCT--GTCA--AAACCC
Read : AGCATTG-CACCAAAGCC
CIGAR: ===XIIMD==II===X==
PhQ  : 333373303335333332
HP   : ------------SHE---
```
The start of a homopolymer region is denoted by "S", the end "E", and the in-between homopolymer region "H". Note that the length of a homopolymer region can be greater than or equal to the length of the homopolymer in the reference. The latter situation occurs with an insertion in the read, for example:
```
Ref  : TCGAAA-ACG
Read : TCTAAAAACG
CIGAR: ==X===I===
HP   : ---SHHHE--
```
- Compute true vs. observed homopolymer length distribution.
- Compute context independent error statistics.
- Compute context dependent error statistics.

## Caveats
**High quality reference must be used.**

In order for the computed statistics to be accurate estimates of the true error profiles, the reference assembly must have a high (e.g., >99.999%) accuracy. If such reference is not readily available, users can create an assembly based on high (>100x) depth sequencing data followed by several rounds of polishing (e.g., using high quality Illumina reads). If a short-reads assembly is used as the reference, then the assembly may contain many contigs whose edges are substantially less accurate than the rest of the reference. For that reason, counterr has an optional argument (-len_trim_contig_edge) that allows the user to exclude from consideration regions close to contig edges when computing various statistics. The program also splits contigs where non-ACGT letters (e.g., "S", "H", "N") are present in the assembly.


**Computational requirements**

counterr is not yet highly optimized for speed or memory. In particular, the tool does not support the use of multiple threads. For large alignment/assembly inputs, the user may run the tool on a cluster with a sufficiently large RAM to avoid an out-of-memory error. As an example, on a MacBook Pro 2015 (2.7 GHz Intel Core i5), counterr took approximately 21 minutes to process a 137 MB bam file and its maximum RAM footprint was 1.2 GB. The run time and memory footprint scale linearly with the input bam file size. If you do not have a large memory cluster available, down-sampling the bam file using either samtools ("view -s") or "-lim" option is recommended.


## Miscellaneous
### Aggregating results
To aggregate results from various runs, please refer to the script "/misc/aggregator.py".


## Known Issues and Work-arounds
The following issues will be fixed in the near future. Until then, work-arounds are provided for some of them.

  The list is currently empty.
## Contact Information
Please use the GitHub issues page for the repository if you have questions. Inquiries about the tool can also be made to jae@dayzerodiagnostics.com.

## License
GNU Public License, V3