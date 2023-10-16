import math
import shutil
import sys
import os
from os.path import join

import numpy as np

from pybedtools import BedTool

import gzip

exons = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/exon_regions.bed"
promoters = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/from_wei/output_ENSG_promoter.UCSC.alternativetss"

from basic import load_metadata

import pybedtools
pybedtools.set_bedtools_path("/home/hudaiber/progs/bedtools2/bin/")
import subprocess as sp
import glob
from collections import defaultdict
import matplotlib.pyplot as plt

DATA_PATH = "/panfs/pan1/devdcode/sanjar/overbinders"
HOMER_DIR = os.path.join(DATA_PATH, "chipseq_files/homer/homer_dirs/")
PEAKS_DIR = os.path.join(DATA_PATH, "chipseq_files/peaks/")

WORK_DIR = os.path.join(DATA_PATH, "definitions/knownMotif")

"""
Homer detects and maps back the locations of all of the known motifs. 
For every ChIP-seq peak:
 1 - location of the top (pvalue) known motif is assigned as TFBS
 2 - if no known motif was bound, the peak center is assigned.  

"""


def concat_known_motifs():

    m = load_metadata()

    ids = list(m["HepG2"].values()) + list(m["K562"].values())

    for i, _id in enumerate(ids):

        print(i+1, len(ids), _id)

        _dir = os.path.join(HOMER_DIR, _id, "knownResults")

        if not os.path.exists(_dir):
            continue

        tmp_file = os.path.join(_dir, "tmp.bed")
        save_file = os.path.join(_dir, "all_known.peaks.bed")

        if os.path.exists(save_file):
            continue

        with open(tmp_file, "w") as of:
            for f in glob.glob(os.path.join(_dir, "known*.peaks.bed")):
                _motif = os.path.basename(f).split(".")[0]
                for line in open(f):
                    out_line = line.rstrip() + "\t%s\n" % _motif
                    of.write(out_line)

        BedTool(tmp_file).sort().saveas(save_file)


def annotate_peaks(data_id):

    peak_file = os.path.join(PEAKS_DIR, "%s.bed" % data_id)
    motifs_file = os.path.join(HOMER_DIR, data_id, "knownResults/all_known.peaks.bed")

    if not os.path.exists(motifs_file):
        raise OSError("concat_known_motifs() should have created: %s" % motifs_file)

    motif_annotated_file = os.path.join(HOMER_DIR, data_id, "annotations.bed")
    peak_tfbs_file = os.path.join(HOMER_DIR, data_id, "peak_tfbs.bed")

    if os.stat(motifs_file).st_size > 0:

        tmp_file = BedTool(peak_file).\
                    intersect(motifs_file, F=1, wo=True).\
                    groupby(g=[1, 2, 3, 7], c=[12, 13, 14, 15, 17], o="collapse")

        with open(motif_annotated_file, "w") as of:
            for r in tmp_file:
                parts = r.fields
                starts = [_ for _ in parts[4].split(",")]
                stops = [_ for _ in parts[5].split(",")]
                names = [_ for _ in parts[6].split(",")]
                scores = [float(_) for _ in parts[7].split(",")]
                codes = [_ for _ in parts[8].split(",")]

                ind = np.argsort(scores)[-1]

                out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (parts[0],
                                                             parts[1],
                                                             parts[2],
                                                             starts[ind],
                                                             stops[ind],
                                                             names[ind],
                                                             codes[ind],
                                                                 parts[3])

                of.write(out_line)

    tmp_lines = []

    with open(motif_annotated_file, "a+") as of:
        cnt = 0
        for r in BedTool(peak_file).\
                intersect(motif_annotated_file, v=True):

            peak = int(r.fields[-1])
            out_line = "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" % (r.chrom,
                                                             r.start,
                                                             r.stop,
                                                             r.start + peak - 4,
                                                             r.start + peak + 4,
                                                             "None",
                                                             "None",
                                                             r.fields[6])
            of.write(out_line)
            cnt += 1

            out_line = "\t".join(r.fields) + "\n"
            tmp_lines.append(out_line)

    with open(motif_annotated_file) as inf, open(peak_tfbs_file, "w") as outf:
        for line in inf:
            parts = line.strip().split()
            out_line = "\t".join([parts[0], parts[3], parts[4], data_id, parts[-1]]) + "\n"
            outf.write(out_line)



def annotate_peaks_8bp(src_file):

    data_id = os.path.basename(src_file).split("_")[0]

    annotation_file = src_file.replace(".bed", "_peak_8bp.bed")
    if os.path.exists(annotation_file):
        return

    bed_str = []
    for r in BedTool(src_file):
        peak = int(r.fields[-1])
        out_line = "%s\t%d\t%d\t%s\t%s" % (r.chrom,
                                           r.start + peak - 4,
                                           r.start + peak + 4,
                                           data_id,
                                           r.fields[6])
        bed_str.append(out_line)

    BedTool("\n".join(bed_str), from_string=True).sort().saveas(annotation_file)


def _annotate_peaks():

    source_files = glob.glob(join(PEAKS_DIR, "*_hg19.bed"))

    for i, f in enumerate(source_files):
        print(i, len(source_files), f)
        annotate_peaks_8bp(f)


def locus_analysis(locus_id="chr7-141438000-141438400", cell_line = "HepG2"):

    # locus_id = "chr1-162217200-162217600"
    # locus_id = "chr10-114847200-114847600"
    # locus_id = "chr12-2861600-2862000"

    a = load_metadata()
    work_dir = DATA_PATH + "/chipseq_files/homer/locus_analysis/"

    l_peaks_ovlp_file = work_dir + "%s.peaks.bed" % locus_id

    l_bed = BedTool(locus_id.replace("-", "\t"), from_string=True)
    overlap_strig_pool = []
    for tf, _id in a[cell_line].items():

        # peak_file = HOMER_DIR + "%s/annotations.bed" % _id
        peak_file = HOMER_DIR + "%s/annotations_8bp.bed" % _id

        if not os.path.exists(peak_file):
            continue

        l_ovrlp = l_bed.intersect(peak_file, wo=True)
        if l_ovrlp.count() == 0:
            continue
        
        for r in l_ovrlp:
            fields = r.fields
            out_line = fields[:3] + fields[4:-2] + [tf]
            out_line = "\t".join(out_line)
            overlap_strig_pool.append(out_line)

    print("%d peaks overlapped. Saving to: %s" % (len(overlap_strig_pool), l_peaks_ovlp_file))
    BedTool("\n".join(overlap_strig_pool), from_string=True).saveas(l_peaks_ovlp_file)

    plot_locus(locus_id)


def millify(n):

    millnames = ['', ' K', ' M', ' B']
    n = float(n)
    millidx = max(0, min(len(millnames) - 1,
                         int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))))

    return '{:.0f}{}'.format(n / 10 ** (3 * millidx), millnames[millidx])


def plot_locus(locus_id):

    work_dir = DATA_PATH + "/chipseq_files/homer/locus_analysis/"
    l_peaks_ovlp_file = work_dir + "%s.peaks.bed" % locus_id
    plot_file = work_dir + "%s.peaks.pdf" % locus_id

    [chrom, l_start, l_stop] = locus_id.split("-")
    l_start = int(l_start)
    l_stop = int(l_stop)

    lines = ["\t".join([chrom] + r.fields[3:]) for r in BedTool(l_peaks_ovlp_file)]
    ovlps = BedTool("\n".join(lines), from_string=True).sort()

    plt.figure(figsize=(10, 30))

    cnt = 3
    for r in ovlps:

        f = r.fields
        p_start, p_stop = int(f[1]), int(f[2])
        m_start, m_stop = int(f[3]), int(f[4])
        tf = f[-1]

        plt.plot([p_start, p_stop], [cnt, cnt], color="grey")
        plt.text(p_start, cnt + 0.1, tf, fontsize=7)
        plt.plot([m_start, m_stop], [cnt, cnt], color="red")

        cnt += 2

        # if cnt > 30:
        #     break

    ax = plt.gca()
    plt.plot([l_start, l_stop], [1, 1], color="blue")
    plt.text(l_start + (l_stop - l_start)/2, 1.1, "400bp")
    ax.axvspan(l_start, l_stop, alpha=0.1, color='blue')

    x_ticks = ax.get_xticks()
    x_tick_labels = [millify(_) for _ in x_ticks]
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel(chrom)
    ax.set_yticks([])
    plt.tight_layout()
    print("Saving to:", plot_file)
    plt.savefig(plot_file)
    plt.close()
