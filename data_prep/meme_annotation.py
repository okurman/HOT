import glob
import gzip
import shutil
import subprocess
import sys
import os
import tempfile
from collections import defaultdict, OrderedDict, Counter
from os.path import join

import networkx
import numpy as np
import configparser

###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
from lib import phyloP
from lib import phastCons

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

PROJECT_DIR = "/panfs/pan1/devdcode/sanjar/overbinders/"
RESULTS_DIR = join(PROJECT_DIR, "chipseq_files/meme/results/")
PEAKS_DIR = join(PROJECT_DIR, "chipseq_files/peaks/")

from overbinders.data_prep.basic import load_metadata


def fimo2bed(fimo_file):

    bed_str = []
    for l in open(fimo_file):
        if l.startswith("motif_id") or l.startswith("#") or not l.strip():
            continue
        parts = l.split("\t")
        if float(parts[8]) > 0.01:
            continue
        out_line = "%s\t%s\t%s\t%s\t%s\n" % (parts[2], parts[3], parts[4], parts[7], parts[8])
        bed_str.append(out_line)

    return BedTool("\n".join(bed_str), from_string=True).sort()


def annotate_single_dir(d):

    # print("Starting for:", d)
    peak_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
    peaks = BedTool(peak_file)

    pool = []

    for n in range(3):
        _d = join(RESULTS_DIR, d, str(n + 1))
        for root, dirs, files in os.walk(_d):
            for f in files:
                if f == "fimo.tsv":
                    fimo_file = join(root, f)
                    bed = fimo2bed(fimo_file)
                    if not bed.count():
                        continue
                    bed_file = fimo_file.replace(".tsv", ".bed.gz")
                    bed.saveas(bed_file)
                    ann_cnt = peaks.intersect(bed, F=1.0, wo=True).groupby(c=[4], o="collapse").count()
                    pool.append([ann_cnt, bed_file])

    if not pool:
        return

    pool = sorted(pool, key=lambda x: x[0], reverse=True)

    res_motif_file = pool[0][1]
    annotation_lines = [str(l) for l in gzip.open(res_motif_file, mode="rt").readlines()]
    annotation_lines = [str(l).replace("\n", "\t%s\n" % d.split("_")[0]) for l in annotation_lines]
    BedTool("\n".join(annotation_lines), from_string=True).saveas(join(RESULTS_DIR, d, "annotations.bed"))

    with open(join(RESULTS_DIR, d, "annotations_info.txt"), "w") as of:
        of.write("%s\n" % res_motif_file)
    print("Extraction completed!")


def expand_motifs_runs():

    dirs = [d for d in os.listdir(RESULTS_DIR) if d.endswith("_controls")]

    for cnt, d in enumerate(dirs):
        if cnt < 7: continue
        print(cnt, len(dirs), d)
        annotate_single_dir(d)


def limit_annotations():

    dirs = [d for d in os.listdir(RESULTS_DIR) if d.endswith("_controls")]

    for cnt, d in enumerate(dirs):

        print(cnt, len(dirs), d)

        peak_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
        old_file = join(RESULTS_DIR, d, "annotations.bed")
        new_file = join(RESULTS_DIR, d, "annotations_peaks.bed")
        try:
            tmp_file = tempfile.NamedTemporaryFile()
            BedTool(old_file).intersect(peak_file, f=1.0, wa=True).saveas(tmp_file.name)
            cmd = "sort -k1,1 -k2,2n %s | uniq > %s " % (tmp_file.name, new_file)
            subprocess.call(cmd, shell=True)
        except:
            print("skipping")


def tmp():

    dirs = [d for d in os.listdir(RESULTS_DIR) if d.endswith("_controls")]

    for cnt, d in enumerate(dirs):
        print(cnt, len(dirs), d)
        f = join(RESULTS_DIR, d, "annotations_info.txt")
        if not os.path.exists(f):
            for f in [join(RESULTS_DIR, d, "annotations.bed"), join(RESULTS_DIR, d, "annotations_peaks.bed")]:
                if os.path.exists(f):
                    os.remove(f)


if __name__ == "__main__":

    annotate_single_dir(sys.argv[1])