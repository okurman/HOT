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
RESULTS_DIR = join(PROJECT_DIR, "chipseq_files/homer/results/")
PEAKS_DIR = join(PROJECT_DIR, "chipseq_files/peaks/")
from overbinders.data_prep.basic import load_metadata


def expand_single_dir(d):

    print("Starting for:", d)
    peak_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
    peaks = BedTool(peak_file)

    pool = []

    for n in range(3):
        _d = join(RESULTS_DIR, d, str(n + 1))

        all_motifs_file = join(_d, "homerMotifs.all.motifs.bed")

        motif2lines = defaultdict(list)
        for l in open(all_motifs_file):
            _motif = l.split()[3]
            motif2lines[_motif].append(l)

        if not motif2lines:
            continue

        for _motif, _lines in motif2lines.items():

            save_file = join(_d, "motif_%s.bed.gz" % _motif)
            print(save_file)
            try:
                motif_bed = BedTool("".join(_lines), from_string=True)
                motif_bed.sort().saveas(save_file)
                motifs_cnt = peaks.intersect(motif_bed, F=1.0, wa=True).groupby(c=[4], o="collapse").count()
                pool.append([motifs_cnt, save_file])
            except:
                print("Skipping")
                continue

    pool = sorted(pool, key=lambda x: x[0], reverse=True)

    res_motif_file = pool[0][1]

    annotation_lines = [str(l) for l in gzip.open(res_motif_file, mode="rt").readlines()]
    annotation_lines = [str(l).replace("\n", "\t%s\n" % d.split("_")[0]) for l in annotation_lines]
    BedTool("\n".join(annotation_lines), from_string=True).saveas(join(RESULTS_DIR, d, "annotations.bed"))

    with open(join(RESULTS_DIR, d, "annotations_info.txt"), "w") as of:
        of.write("%s\n" % res_motif_file)

    print("Extraction completed!")


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


def merge_annotations():

    homer_dir = join(PROJECT_DIR, "chipseq_files/homer/results/")
    meme_dir = join(PROJECT_DIR, "chipseq_files/meme/results/")

    dirs = [d for d in os.listdir(RESULTS_DIR) if d.endswith("_controls")]

    pool = []
    for cnt, d in enumerate(dirs):
        print(cnt, len(dirs), d)

        peaks_file = join(PEAKS_DIR, "%s.bed" % d.replace("_controls", ""))
        peaks = BedTool(peaks_file)
        peaks_cnt = peaks.groupby(c=[4], o="collapse").count()

        homer_file = join(homer_dir, d, "annotations_peaks.bed")
        homer_cnt = peaks.intersect(homer_file, F=1.0, wa=True).groupby(c=[4], o="collapse").count()

        res_file = homer_file
        res_frac = homer_cnt / peaks_cnt

        meme_file = join(meme_dir, d, "annotations_peaks.bed")
        if os.path.exists(meme_file):
            meme_cnt = peaks.intersect(meme_file, F=1.0, wa=True).groupby(c=[4], o="collapse").count()
            if meme_cnt > homer_cnt:
                res_file = meme_file
                res_frac = meme_cnt / peaks_cnt

        pool.append([d.replace("_controls", ""), res_frac, res_file])

    report_file = join(PROJECT_DIR, "chipseq_files/motif_annotations_info.txt")
    with open(report_file, "w") as outf:
        for p in pool:
            out_line = "\t".join([str(_) for _ in p]) + "\n"
            outf.write(out_line)


def extract_motif_concats():

    a = load_metadata()

    report_file = join(PROJECT_DIR, "chipseq_files/motif_annotations_info.txt")
    id2ann_file = {l.split()[0]:l.split()[2].rstrip() for l in open(report_file)}

    for cell_line in ["K562", "HepG2"]:
        print(cell_line)
        concat_file = join(PROJECT_DIR, "chipseq_files/%s_motif_concat.bed" % cell_line)
        print(concat_file)
        with open(concat_file, "w") as outf:
            for tf_id in a[cell_line].values():
                ann_file = id2ann_file[tf_id+"_hg19"]
                for l in open(ann_file):
                    parts = l.split()
                    # parts = parts[:3] + [tf_id]
                    logodds = "-" if "meme" in ann_file else parts[4]
                    parts = parts[:3] + [logodds, parts[-1]]
                    out_line = "\t".join(parts) + "\n"
                    outf.write(out_line)
                outf.flush()


def extract_motif_concats_v2():

    a = load_metadata()

    for cell_line in ["K562", "HepG2"]:
        print(cell_line)
        concat_file = join(PROJECT_DIR, "chipseq_files/%s_motif_concat.bed.gz" % cell_line)
        print(concat_file)
        tmp_file = tempfile.NamedTemporaryFile()
        with open(tmp_file.name, "w") as outf:
            for tf_id in a[cell_line].values():
                ann_file = join(RESULTS_DIR, "%s_hg19_controls/annotations_fimo.bed" % tf_id)
                [outf.write(l) for l in open(ann_file)]
                outf.flush()
        BedTool(tmp_file.name).sort().saveas(concat_file)


def extract_from_fimo(d):

    _dir = join(RESULTS_DIR, d)

    info_parts = open(join(_dir, "annotations_info.txt")).readline().split("/")
    motif_name = info_parts[-1].split("_")[-1].replace(".bed.gz", "").rstrip()
    fimo_file = join(_dir, info_parts[-2], "fimo_target/fimo.tsv.gz")

    fimo_lines = [l for l in gzip.open(fimo_file, "rt")]
    fimo_lines = [l for l in fimo_lines if not (l.startswith("motif_id") or l.startswith("#") or not l.strip())]
    fimo_lines = [l for l in fimo_lines if motif_name in l]

    bed_lines = []
    for l in fimo_lines:
        parts = l.split("\t")
        # if float(parts[8]) > 0.05:
        #     continue
        out_line = "%s\t%s\t%s\t%s\t%s\n" % (parts[2], parts[3], parts[4], parts[7], parts[8])
        bed_lines.append(out_line)

    BedTool("".join(bed_lines), from_string=True).saveas(join(_dir, "annotations_pvalue.bed"))




if __name__ == "__main__":

    expand_single_dir(sys.argv[1])