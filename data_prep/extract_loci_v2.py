import glob
import shutil
import subprocess
import os
import tempfile
from os.path import join

import numpy as np
from pybedtools import BedTool

import warnings
warnings.filterwarnings('ignore')

EXONS_FILE = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/exon_regions.bed"
PROMOTERS_FILE = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed"
INTRONS_FILE = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/genes/gencode/knownGene_introns.liftOver_hg19.bed"
CODING_PROMOTERS_FILE = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed.coding"
NONCODING_PROMOTERS_FILE = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed.noncoding"

DATA_PATH = "/net/intdev/devdcode/sanjar/overbinders"
PERC_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/perc_bins/")
LOG_BINS_DIR = os.path.join(DATA_PATH, "definitions/peak_8bp_v2/log_bins/")
PEAKS_DIR = join(DATA_PATH, "chipseq_files/peaks/")
LOCI_DIR = join(DATA_PATH, "chipseq_files/loci_v2/")

import subprocess as sp
# from overbinders.data_prep.basic import load_metadata
# from basic import load_metadata
# from data_prep.basic import load_metadata


"""
Order of execution:

- peakwise_annotation.annotate_peaks()
- concat_all_rep_loci()
- extract_bound_loci()
- add_stats_columns()
- extract_to_bins():

"""


def run_locus_extraction(cell_line="HepG2"):

    print("Concat")
    concat_all_rep_loci(cell_line)
    print("extract_bound_loci")
    extract_bound_loci(cell_line)
    print("add_stats_columns")
    add_stats_columns(cell_line)
    print("extract_to_bins")
    extract_to_bins(cell_line)
    print("overlap_promoters")
    overlap_promoters(cell_line)


def concat_all_rep_loci(cell_line="HepG2"):

    cell_line2tf2file_id = load_metadata()
    all_files = cell_line2tf2file_id[cell_line].values()

    tmp_file = tempfile.NamedTemporaryFile()

    for i, d in enumerate(all_files):
        print(i, len(all_files))
        rep_file = join(PEAKS_DIR, "%s_hg19_peak_8bp.bed" % d)
        cmd = "cat %s >> %s" % (rep_file, tmp_file.name)
        sp.call(cmd, shell=True)

    save_file = join(LOCI_DIR, "%s_all_tfbs.bed" % cell_line)
    print("Saving to:", save_file)
    BedTool(tmp_file.name).sort().groupby(g=[1, 2, 3], c=[4, 5], o="collapse").saveas(save_file)


def extract_bound_loci(cell_line="HepG2"):

    all_rep_loci_file = join(LOCI_DIR, "%s_all_tfbs.bed" % cell_line)

    windows_file = os.path.join(DATA_PATH, "chipseq_files/homer/windows_400_hg19.bed")
    windows_grp_file = os.path.join(LOCI_DIR, "%s_windows_400_groupby.bed" % cell_line)

    print(windows_grp_file)

    BedTool(windows_file).\
        intersect(all_rep_loci_file, wo=True).\
        groupby(g=[1, 2, 3], c=[5, 6, 7, 8], o="collapse").\
        saveas(windows_grp_file)


def add_stats_columns(cell_line="HepG2"):

    windows_grp_file = join(LOCI_DIR, "%s_windows_400_groupby.bed" % cell_line)
    save_file = join(LOCI_DIR, "%s_400_loci.bed" % cell_line)

    print(windows_grp_file)
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


def extract_to_bins_perc(cell_line="HepG2"):

    locus_file = os.path.join(LOCI_DIR, "%s_400_loci.bed" % cell_line)

    tmp_file = join(LOCI_DIR, "tmp.txt")
    cmd = "sed -n '1d;p' %s | sort -k6,6n > %s"
    subprocess.call(cmd % (locus_file, tmp_file), shell=True)
    shutil.move(tmp_file, locus_file)

    fnames = ["%s_400_loci.%d.bed" % (cell_line, i) for i in range(10)]
    lines = open(locus_file).readlines()

    block_size = len(lines)//10
    blocks = [lines[i: i+block_size] for i in range(0, len(lines), block_size)]

    for block, fname in zip(blocks, fnames):
        print(fname)
        with open(join(PERC_BINS_DIR, fname), "w") as of:
            [of.write(l) for l in block]


def extract_to_log_bins(cell_line="HepG2"):

    locus_file = os.path.join(LOCI_DIR, "%s_400_loci.bed" % cell_line)
    lines = open(locus_file).readlines()
    if cell_line == "H1":
        max_tfs = int(lines[-1].split("\t")[5]) + 1
    else:
        max_tfs = int(lines[-1].split("\t")[5]) + 10

    bin_edges = [1, 2, 3, 4] + list(np.logspace(np.log10(5), np.log10(max_tfs), 11, dtype=int))
    bin_edges = np.asarray(bin_edges)

    fnames = ["%s_400_loci.%d.bed" % (cell_line, i) for i in range(bin_edges.shape[0])]

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
        print(_fname)
        with open(os.path.join(LOG_BINS_DIR, _fname), "w") as of:
            for line in _lines:
                of.write(line)


def overlap_promoters(cell_line, bins_dir=LOG_BINS_DIR):

    # fnames = ["%s_400_loci.%d.bed" % (cell_line, i) for i in range(14)]
    fnames = glob.glob(join(bins_dir, "%s_400_loci.*.bed" % cell_line))

    tmp_file = tempfile.NamedTemporaryFile()

    for loci_file in fnames:
        # loci_file = join(LOG_BINS_DIR, fname)
        loci_bed = BedTool(loci_file)
        noprom_file = loci_file + ".noprom"
        prom_file = loci_file + ".prom"
        loci_bed.sort().intersect(PROMOTERS_FILE, v=True).sort().saveas(tmp_file.name)
        cmd = "uniq %s > %s" % (tmp_file.name, noprom_file)
        sp.call(cmd, shell=True)
        loci_bed.sort().intersect(PROMOTERS_FILE, wa=True).sort().saveas(tmp_file.name)
        cmd = "uniq %s > %s" % (tmp_file.name, prom_file)
        sp.call(cmd, shell=True)
        print(prom_file)
        print(noprom_file)
        print("\n")


def extract_to_HOTs(cell_line="HepG2"):

    d = "/panfs/pan1/devdcode/sanjar/overbinders/definitions/peak_8bp_v2/HOTs/"

    locus_file = os.path.join(LOCI_DIR, "%s_400_loci.bed" % cell_line)
    HOTs_file = join(d, "%s_HOTs.bed" % cell_line)
    proms_file = join(d, "%s_HOTs.proms.bed" % cell_line)
    enhs_file = join(d, "%s_HOTs.noproms.bed" % cell_line)

    # lines = open(locus_file).readlines()
    # outf = open(save_file, "w")
    # for l in lines:
    #     parts = l.split()
    #     tfs = int(parts[4])
    #     if tfs >= 50:
    #         outf.write(l)

    BedTool(HOTs_file).sort().intersect(PROMOTERS_FILE, v=True).sort().saveas(enhs_file)

    BedTool(HOTs_file).\
        sort().\
        intersect(PROMOTERS_FILE, wo=True).\
        sort().\
        groupby(g=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], c=[15], o="collapse").\
        saveas(proms_file)

