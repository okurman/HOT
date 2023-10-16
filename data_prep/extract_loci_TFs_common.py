import glob
import shutil
import subprocess
import os
import tempfile
from os.path import join, basename
from random import shuffle

import numpy as np
from pybedtools import BedTool

import warnings
warnings.filterwarnings('ignore')

DATA_PATH = "/panfs/pan1/devdcode/sanjar/overbinders"
PERC_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/perc_bins/")
LOG_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/log_bins_random/common_TFs/")
PEAKS_DIR = join(DATA_PATH, "chipseq_files/peaks/")

from matplotlib_venn import venn3
import matplotlib.pyplot as plt

from basic import load_metadata
a = load_metadata()

common_tfs = set(a["H1"].keys()).intersection(set(a["HepG2"].keys())).intersection(set(a["K562"].keys()))


def run_locus_extraction(cell_line="HepG2"):

    print("Concat")
    concat_all_rep_loci(cell_line)
    print("extract_bound_loci")
    extract_bound_loci(cell_line)
    print("add_stats_columns")
    add_stats_columns(cell_line)
    print("extract_to_bins")
    extract_to_log_bins(cell_line)
    print("extract to HOTs")
    extract_to_HOTs(cell_line)

    cl_dir = join(LOG_BINS_DIR, cell_line)
    for f in os.listdir(cl_dir):

        if f.endswith("HOTs.bed.gz"):
            continue
        os.remove(join(cl_dir, f))


def concat_all_rep_loci(cell_line="HepG2"):

    cell_line2tf2file_id = load_metadata()

    all_files = [cell_line2tf2file_id[cell_line][tf] for tf in common_tfs]
    print(len(common_tfs), len(all_files))

    beds = [BedTool(join(PEAKS_DIR, "%s_hg19_peak_8bp.bed" % _)) for _ in all_files]

    concat_bed = beds[0]
    for bed in beds[1:]:
        concat_bed = concat_bed.cat(bed, postmerge=False)

    save_file = join(LOG_BINS_DIR, cell_line, "all_tfbs.bed")
    concat_bed.sort().groupby(g=[1, 2, 3], c=[4, 5], o="collapse").saveas(save_file)


def extract_bound_loci(cell_line="HepG2"):

    all_rep_loci_file = join(LOG_BINS_DIR, cell_line, "all_tfbs.bed")

    windows_file = os.path.join(DATA_PATH, "chipseq_files/homer/windows_400_hg19.bed")
    windows_grp_file = os.path.join(LOG_BINS_DIR, cell_line, "%windows_400_groupby.bed")

    BedTool(windows_file).\
        intersect(all_rep_loci_file, wo=True).\
        groupby(g=[1, 2, 3], c=[5, 6, 7, 8], o="collapse").\
        saveas(windows_grp_file)


def add_stats_columns(cell_line="HepG2"):

    windows_grp_file = join(LOG_BINS_DIR, cell_line, "%windows_400_groupby.bed")
    save_file = join(LOG_BINS_DIR, cell_line, "400_loci.bed")

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


def extract_to_log_bins(cell_line="HepG2"):

    locus_file = os.path.join(LOG_BINS_DIR, cell_line, "400_loci.bed")
    lines = open(locus_file).readlines()

    max_tfs = max([int(l.split("\t")[5]) for l in lines[1:]])

    bin_edges = [1, 2, 3, 4] + list(np.logspace(np.log10(5), np.log10(max_tfs), 11, dtype=int))
    bin_edges = np.asarray(bin_edges)

    fnames = ["400_loci.bin.%d.bed" % j for j in range(bin_edges.shape[0])]

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


def extract_to_HOTs(cell_line="HepG2"):

    blacklist_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/blacklisted_regions/hg19-blacklist.v2.bed"

    files = glob.glob(join(LOG_BINS_DIR, cell_line, "400_loci.bin.*.bed"))
    max_ind = max([int(basename(_).split(".")[2]) for _ in files])

    files = [join(LOG_BINS_DIR, cell_line, "400_loci.bin.%d.bed" % j) for j in range(max_ind-3, max_ind+1)]
    beds = [BedTool(f) for f in files]
    bed = beds[0].cat(beds[1], postmerge=False).cat(beds[2], postmerge=False).cat(beds[3], postmerge=False)

    bed.sort().intersect(blacklist_file, v=True, wa=True).saveas(join(LOG_BINS_DIR, cell_line, "HOTs.bed.gz"))


def load_overlap_data(matched=True):

    if matched:

        loci_dir = "/panfs/pan1/devdcode/sanjar/overbinders/definitions/peak_8bp_v2/log_bins_random/common_TFs/"
        hepg2 = BedTool(join(loci_dir, "HepG2/HOTs.bed.gz"))
        k562 = BedTool(join(loci_dir, "K562/HOTs.bed.gz"))
        h1 = BedTool(join(loci_dir, "H1/HOTs.bed.gz"))

    else:

        loci_dir = "/panfs/pan1/devdcode/sanjar/overbinders/definitions/peak_8bp_v2/HOTs/"
        hepg2 = BedTool(join(loci_dir, "HepG2_HOTs.bed"))
        k562 = BedTool(join(loci_dir, "K562_HOTs.bed"))
        h1 = BedTool(join(loci_dir, "H1_HOTs.bed"))

    data = [
            hepg2.intersect(k562.cat(h1, postmerge=False), v=True, f=0.5).count(),
            k562.intersect(hepg2.cat(h1, postmerge=False), v=True, f=0.5).count(),
            hepg2.intersect(k562, wa=True, f=0.5).cat(hepg2.intersect(k562, wb=True, f=0.5), postmerge=False).intersect(h1, v=True, f=0.5).count(),
            h1.intersect(k562.cat(hepg2, postmerge=False), v=True, f=0.5).count(),
            hepg2.intersect(h1, wa=True, f=0.5).cat(hepg2.intersect(h1, wb=True, f=0.5), postmerge=False ).intersect(k562, v=True, f=0.5).count(),
            k562.intersect(h1, wa=True, f=0.5).cat(k562.intersect(h1, wb=True, f=0.5), postmerge=False).intersect(hepg2, v=True, f=0.5).count(),
            hepg2.intersect(k562, wa=True, f=0.5).cat(hepg2.intersect(k562, wb=True, f=0.5), postmerge=False).intersect(h1, wa=True, f=0.5).count()
        ]

    total_count = hepg2.count() + k562.count() + h1.count()

    data_perc = [100 * (_/total_count) for _ in data]

    return data, data_perc


def get_overlaps_H1():

    loci_dir = "/panfs/pan1/devdcode/sanjar/overbinders/definitions/peak_8bp_v2/log_bins_random/common_TFs/"
    hepg2 = BedTool(join(loci_dir, "HepG2/HOTs.bed.gz"))
    k562 = BedTool(join(loci_dir, "K562/HOTs.bed.gz"))
    h1 = BedTool(join(loci_dir, "H1/HOTs.bed.gz"))

    common_TFs_ovp = h1.intersect(hepg2.cat(k562), wa=True).merge(d=-1).count()/h1.count()
    print("common TFs: %f" % common_TFs_ovp)

    loci_dir = "/panfs/pan1/devdcode/sanjar/overbinders/definitions/peak_8bp_v2/HOTs/"
    hepg2 = BedTool(join(loci_dir, "HepG2_HOTs.bed"))
    k562 = BedTool(join(loci_dir, "K562_HOTs.bed"))
    h1 = BedTool(join(loci_dir, "H1_HOTs.bed"))

    all_TFs_ovp = h1.intersect(hepg2.cat(k562), wa=True).merge(d=-1).count() / h1.count()
    print("All TFs: %f" % all_TFs_ovp)

    return common_TFs_ovp, all_TFs_ovp



def plot_venn_diagram(data, data_perc, common=True):

    PLOTS_DIR = "/panfs/pan1/devdcode/sanjar/overbinders/plots/HOTs/"

    plt.figure(figsize=(4, 4))
    venn3(subsets=data, set_labels=('HepG2', 'K562', 'H1'))
    save_file = join(PLOTS_DIR, "VENN_HOTs_%s.pdf" % ("common" if common else "all"))
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(4, 4))
    venn3(subsets=data_perc, set_labels=('HepG2', 'K562', 'H1'), subset_label_formatter=lambda x: "%.1f%%" % np.round(x,1))
    save_file = join(PLOTS_DIR, "VENN_HOTs_%s_perc.pdf" % ("common" if common else "all"))
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()

    # save_file = join(PLOTS_DIR, "VENN_overbinders_3CL_nolabels.pdf")
    # for _id in ['100', '110', '010', '101', '111', '011', '001']:
    #     v.get_label_by_id(_id).set_text("")



