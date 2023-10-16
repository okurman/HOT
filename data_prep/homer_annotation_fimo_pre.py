#!/home/hudaiber/miniconda3/bin/python

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


import pybedtools
pybedtools.helpers.set_bedtools_path("/home/hudaiber/progs/bedtools2/bin/")

from pybedtools import BedTool
from lib import phyloP
from lib import phastCons

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

PROJECT_DIR = "/panfs/pan1/devdcode/sanjar/overbinders/"
RESULTS_DIR = join(PROJECT_DIR, "chipseq_files/homer/results/")
PEAKS_DIR = join(PROJECT_DIR, "chipseq_files/peaks/")
from overbinders.data_prep.basic import load_metadata


def fimo_to_beds(fimo_file, save_dir):

    lines = [l for l in gzip.open(fimo_file, "rt")]
    lines = [l for l in lines if not (l.startswith("motif_id") or l.startswith("#") or not l.strip())]

    motif2lines = defaultdict(list)

    for l in lines:
        parts = l.strip().split()
        pvalue = parts[7]
        qvalue = parts[8]

        if float(qvalue) > 0.05:
            continue

        motif_id = parts[1]
        chrom = parts[2].split(":")[0]
        peak_start = int(parts[2].split(":")[1].split("-")[0])
        motif_start = peak_start + int(parts[3])
        motif_stop = peak_start + int(parts[4])

        out_line = "\t".join([str(_) for _ in [chrom, motif_start, motif_stop, pvalue, qvalue]])
        motif2lines[motif_id].append(out_line)

    for motif, lines in motif2lines.items():
        save_file = join(save_dir, "motif_%s.bed.gz" % motif)
        motif_bed = BedTool("".join(lines), from_string=True)
        motif_bed.sort().saveas(save_file)


def get_top_motif_from_fimo(fimo_file, peaks_bed, suffix):

    lines = [l for l in gzip.open(fimo_file, "rt")]
    lines = [l for l in lines if not (l.startswith("motif_id") or l.startswith("#") or not l.strip())]

    motif2lines = defaultdict(list)

    for l in lines:
        parts = l.strip().split()
        pvalue = parts[7]
        qvalue = parts[8]

        if float(qvalue) > 0.05:
            continue

        motif_id = parts[1]
        chrom = parts[2].split(":")[0]
        peak_start = int(parts[2].split(":")[1].split("-")[0])
        motif_start = peak_start + int(parts[3])
        motif_stop = peak_start + int(parts[4])

        out_line = "\t".join([str(_) for _ in [chrom, motif_start, motif_stop, pvalue, qvalue, suffix]])
        motif2lines[motif_id].append(out_line)

    # for k, v in motif2lines.items():
    #     print(k, len(v))

    if not motif2lines:
        return [None, None, None]

    pool = []

    for motif, lines in motif2lines.items():

        motif_bed = BedTool("\n".join(set(lines)), from_string=True).sort()
        cnt = peaks_bed.intersect(motif_bed, wo=True).groupby(g=[1,2,3], c=[4], o=["collapse"]).count()
        pool.append([cnt, motif_bed, motif])

    top = sorted(pool, key=lambda x: x[0], reverse=True)[0]
    
    return top


def expand_single_dir(d):
    
    print("Starting for:", d)
    peak_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
    peaks = BedTool(peak_file)

    pool = []

    for n in range(3):

        _d = join(RESULTS_DIR, d, str(n + 1))

        fimo_target_dir = join(_d, "fimo_target")
        # fimo_ctr_dir = join(_d, "fimo_control")
        # fimo_to_beds(join(fimo_target_dir, "fimo.tsv.gz"), fimo_target_dir)
        # fimo_to_beds(join(fimo_ctr_dir, "fimo.tsv.gz"), fimo_ctr_dir)

        fimo_file = join(fimo_target_dir, "fimo.tsv.gz")
        if not os.path.exists(fimo_file):
            continue
        [cnt, bed, motif] = get_top_motif_from_fimo(fimo_file, peaks, d.replace("_hg19_controls", ""))
        if not cnt and not bed:
            continue
        pool.append([cnt, bed, motif, n+1])

    pool = sorted(pool, key=lambda x: x[0], reverse=True)[0]

    pool[1].saveas(join(RESULTS_DIR, d, "annotations_fimo.bed"))

    with open(join(RESULTS_DIR, d, "annotations_fimo_info.txt"), "w") as of:
        of.write("%d\n" % pool[3])
        of.write("%s\n" % pool[2])
        of.write("%f\n" % (pool[0]/peaks.count()))

    print("Extraction completed!")


def extract_motif_concats():

    a = load_metadata()

    for cell_line in ["K562", "HepG2"]:
        print(cell_line)
        tmp_file = tempfile.NamedTemporaryFile()
        for tf_id in a[cell_line].values():
            ann_file = join(RESULTS_DIR, "%s_hg19_controls/annotations_fimo.bed" % tf_id)
            cmd = "cat %s >> %s" % (ann_file, tmp_file.name)
            subprocess.call(cmd, shell=True)

        concat_file = join(PROJECT_DIR, "chipseq_files/%s_motif_concat_fimo.bed.gz" % cell_line)
        print(concat_file)
        cmd = "sort -k1,1 -k2,2n %s | gzip > %s" % (tmp_file.name, concat_file)
        subprocess.call(cmd, shell=True)



if __name__ == "__main__":

    expand_single_dir(sys.argv[1])