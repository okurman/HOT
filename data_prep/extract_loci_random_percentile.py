import glob
import shutil
import subprocess
import os
import sys
import tempfile
from os.path import join, basename, exists
from random import shuffle

import numpy as np
from pybedtools import BedTool

import warnings
warnings.filterwarnings('ignore')

DATA_PATH = "/panfs/pan1/devdcode/sanjar/overbinders"
PERC_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/perc_bins/")
# LOG_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/log_bins_random/subsample_to_half/")

LOG_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/log_bins_random/")
PEAKS_DIR = join(DATA_PATH, "chipseq_files/peaks/")

from basic import load_metadata

a = load_metadata()


def run_locus_extraction(frac="30", cell_line="HepG2"):

    global LOG_BINS_DIR
    LOG_BINS_DIR += "%s/" % frac

    if not os.path.exists(LOG_BINS_DIR):
        os.mkdir(LOG_BINS_DIR)

    if not os.path.exists(join(LOG_BINS_DIR, cell_line)):
        os.mkdir(join(LOG_BINS_DIR, cell_line))

    for i in range(10):

        print("Concat")
        concat_all_rep_loci(i, cell_line, frac)
        print("extract_bound_loci")
        extract_bound_loci(i, cell_line)
        print("add_stats_columns")
        add_stats_columns(i, cell_line)
        print("extract_to_bins")
        extract_to_log_bins(i, cell_line)
        print("extract to HOTs")
        extract_to_HOTs(i, cell_line)

    clear_directory(cell_line)


def concat_all_rep_loci(i, cell_line="HepG2", frac="20"):

    cell_line2tf2file_id = load_metadata()

    tfs = list(cell_line2tf2file_id[cell_line].keys())

    shuffle(tfs)
    frac = int(frac)/100
    rand_tfs = tfs[:int(len(tfs) * frac)]

    all_files = [cell_line2tf2file_id[cell_line][tf] for tf in rand_tfs]

    beds = [BedTool(join(PEAKS_DIR, "%s_hg19_peak_8bp.bed" % _)) for _ in all_files]

    concat_bed = beds[0]
    for bed in beds[1:]:
        concat_bed = concat_bed.cat(bed, postmerge=False)

    save_file = join(LOG_BINS_DIR, cell_line, "%d_all_tfbs.bed" % i)
    print("Saving to:", save_file)
    concat_bed.sort().groupby(g=[1, 2, 3], c=[4, 5], o="collapse").saveas(save_file)


def extract_bound_loci(i, cell_line="HepG2"):

    all_rep_loci_file = join(LOG_BINS_DIR, cell_line, "%d_all_tfbs.bed" % i)

    windows_file = os.path.join(DATA_PATH, "chipseq_files/homer/windows_400_hg19.bed")
    windows_grp_file = os.path.join(LOG_BINS_DIR, cell_line, "%d_windows_400_groupby.bed" % i)

    print(windows_grp_file)

    BedTool(windows_file).\
        intersect(all_rep_loci_file, wo=True).\
        groupby(g=[1, 2, 3], c=[5, 6, 7, 8], o="collapse").\
        saveas(windows_grp_file)


def add_stats_columns(i, cell_line="HepG2"):

    windows_grp_file = join(LOG_BINS_DIR, cell_line, "%d_windows_400_groupby.bed" % i)
    save_file = join(LOG_BINS_DIR, cell_line, "%d_400_loci.bed" % i)

    print(save_file)

    with open(save_file, "w") as of:

        header_line = "#chr\tstart\tstop\tcov\tuq_tfbs\tuq_tfs\ttf_starts\ttf_stops\tencode_ids\tsignal_values\n"
        of.write(header_line)

        out_fmt = "%s\t%d\t%d\t%f\t%d\t%d\t%s\t%s\t%s\t%s\n"

        for r in BedTool(windows_grp_file):

            total_range = set([_ for _ in range(r.start, r.stop)])

            tf_ranges = []
            tf_starts = [int(_) for _ in r.fields[3].split(",")]
            tf_stops = [int(_) for _ in r.fields[4].split(",")]
            for (_start, _stop) in zip(tf_starts, tf_stops):
                tf_ranges += [_ for _ in range(_start, _stop)]

            overlap_cnt = len(total_range.intersection(tf_ranges))
            cov = np.round(overlap_cnt / len(total_range), 5)

            tfs = len(tf_starts)
            tfs_uq = len(set(r.fields[5].split(",")))

            out_line = out_fmt % (r.chrom,
                                  r.start,
                                  r.stop,
                                  cov,
                                  tfs,
                                  tfs_uq,
                                  r.fields[3],
                                  r.fields[4],
                                  r.fields[5],
                                  r.fields[6])

            of.write(out_line)


def extract_to_log_bins(i, cell_line="HepG2"):

    locus_file = os.path.join(LOG_BINS_DIR, cell_line, "%d_400_loci.bed" % i)
    lines = open(locus_file).readlines()

    max_tfs = max([int(l.split("\t")[5]) for l in lines[1:]])
    # max_tfs = 20

    bin_edges = [1, 2, 3, 4] + list(np.logspace(np.log10(5), np.log10(max_tfs), 11, dtype=int))
    bin_edges = np.asarray(bin_edges)

    fnames = ["%d_400_loci.bin.%d.bed" % (i, j) for j in range(bin_edges.shape[0])]

    lines_list = [[] for _ in range(np.size(bin_edges, 0))]

    # skip the header
    for l in lines[1:]:
        parts = l.split()
        tfs = int(parts[5])
        ind = np.argmax((bin_edges - tfs) > 0) - 1
        lines_list[ind].append(l)

    for _lines, _fname in zip(lines_list, fnames):
        if not _lines:
            continue
        with open(os.path.join(LOG_BINS_DIR, cell_line, _fname), "w") as of:
            for line in _lines:
                of.write(line)


def extract_to_HOTs(i, cell_line="HepG2"):

    blacklist_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/blacklisted_regions/hg19-blacklist.v2.bed"

    files = glob.glob(join(LOG_BINS_DIR, cell_line, "%d_400_loci.bin.*.bed" % i))
    max_ind = max([int(basename(_).split(".")[2]) for _ in files])

    files = [join(LOG_BINS_DIR, cell_line, "%d_400_loci.bin.%d.bed" % (i, j)) for j in range(max_ind-3, max_ind+1)]
    beds = [BedTool(f) for f in files if exists(f)]
    bed = beds[0]
    for _bed in beds[1:]:
        bed = bed.cat(_bed, postmerge=False)

    bed.sort().intersect(blacklist_file, v=True, wa=True).saveas(join(LOG_BINS_DIR, cell_line, "%d_HOTs.bed.gz" % i))


def clear_directory(cell_line="HepG2"):

    _d = join(LOG_BINS_DIR, cell_line)

    files = [f for f in os.listdir(_d) if not f.endswith("HOTs.bed.gz")]
    [os.remove(join(_d, f)) for f in files]


if __name__ == "__main__":

    run_locus_extraction(frac=sys.argv[1], cell_line=sys.argv[2])