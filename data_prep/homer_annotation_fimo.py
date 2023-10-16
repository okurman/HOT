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


def get_top_motif_from_fimo(fimo_file, suffix):

    lines = [l for l in gzip.open(fimo_file, "rt")]
    lines = [l for l in lines if not (l.startswith("motif_id") or l.startswith("#") or not l.strip())]

    motif2lines = defaultdict(list)

    for l in lines:
        parts = l.strip().split()
        pvalue = parts[7]
        qvalue = parts[8]

        if float(qvalue) > 0.05:
            continue
        if float(pvalue) > 10**-4:
            continue

        motif_id = parts[1].split(",")[0]
        chrom = parts[2].split(":")[0]
        start = int(parts[3])
        stop = int(parts[4])

        out_line = "\t".join([str(_) for _ in [chrom, start, stop, pvalue, qvalue, suffix]])
        motif2lines[motif_id].append(out_line)

    return motif2lines


def run_fimo(d, fasta_file):

    cmd_fmt = "cat %s | /home/hudaiber/meme/libexec/meme-5.0.2/chen2meme > %s"
    m1 = join(d, "homerResults/motif1.motif")
    m2 = join(d, "homerResults/motif1.motif.meme")
    subprocess.call(cmd_fmt % (m1, m2), shell=True)

    fimo_dir = join(d, "fimo_motif1")
    if not os.path.exists(fimo_dir):
        os.mkdir(fimo_dir)
    cmd_fmt = "/home/hudaiber/meme/bin/fimo --max-strand --max-stored-scores 1000000 --parse-genomic-coord --oc %s %s %s >%s 2>%s"
    cmd = cmd_fmt % (fimo_dir,
                     m2,
                     fasta_file,
                     join(fimo_dir, "std.out"),
                     join(fimo_dir, "std.err"))

    subprocess.call(cmd, shell=True)
    os.remove(join(fimo_dir, "cisml.xml"))
    os.remove(join(fimo_dir, "fimo.gff"))
    os.remove(join(fimo_dir, "fimo.xml"))
    os.remove(join(fimo_dir, "fimo.html"))
    subprocess.call("gzip -f %s" % join(fimo_dir, "fimo.tsv"), shell=True)


def expand_single_dir(d):

    print("Starting for:", d)

    peak_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
    peak_bed = BedTool(peak_file)
    fasta_file = "/dev/shm/%s.fa" % d
    print("Extracting fasta")
    peak_bed.getfasta(fi="/panfs/pan1/patternquest/data/genomes/hg19/hg19.fa").save_seqs(fasta_file)

    _d = join(RESULTS_DIR, d, "4")

    fimo_target_dir = join(_d, "fimo_motif1")
    fimo_file = join(fimo_target_dir, "fimo.tsv.gz")
    if not os.path.exists(fimo_file):
        print("Running FIMO")
        run_fimo(_d, fasta_file)

    motif2lines = get_top_motif_from_fimo(fimo_file, d.replace("_hg19_controls", ""))

    with open(join(_d, "homerResults/motif1.motif")) as inf:
        l = inf.readline()
        motif_name = l.split()[1].split(",")[0]
    bed_str = "\n".join(motif2lines[motif_name])
    motif_bed = BedTool(bed_str, from_string=True)

    # best match if reverse complement is also a match
    motif_bed = motif_bed.groupby(g=[1, 2, 3, 6], c=[4, 5], o="min")

    bed_str = "\n".join(["\t".join(r[:3]+r[4:]+[r[3]]) for r in motif_bed])
    BedTool(bed_str, from_string=True).saveas(join(RESULTS_DIR, d, "annotations_fimo_c4.bed"))

    print("Motifs:", motif_bed.count())
    print("Peaks:", peak_bed.count())
    print("Coveraga:", motif_bed.count()/peak_bed.count()*100)
    print("Extraction completed!")


def extract_motif_concats():

    a = load_metadata()

    for cell_line in ["K562", "HepG2"]:
        print(cell_line)
        tmp_file = tempfile.NamedTemporaryFile()
        for tf, tf_id in a[cell_line].items():
            if "POL" in tf:
                continue
            ann_file = join(RESULTS_DIR, "%s_hg19_controls/annotations_fimo_c4.bed" % tf_id)
            if not os.path.exists(ann_file):
                continue
            cmd = "cat %s >> %s" % (ann_file, tmp_file.name)
            subprocess.call(cmd, shell=True)

        concat_file = join(PROJECT_DIR, "chipseq_files/%s_motif_concat_fimo.bed.gz" % cell_line)
        print(concat_file)
        cmd = "sort -k1,1 -k2,2n %s | gzip > %s" % (tmp_file.name, concat_file)
        subprocess.call(cmd, shell=True)


def check_annotations():
    
    # unannotated TFs

    a = load_metadata()

    for cl in ["HepG2", "K562"]:
        for tf, id in a[cl].items():
            ann_file = join(RESULTS_DIR, "%s_hg19_controls/annotations_fimo_c4.bed" % id)
            if not os.path.exists(ann_file):
                print(cl, tf, id)


if __name__ == "__main__":

    expand_single_dir(sys.argv[1]+"_hg19_controls")