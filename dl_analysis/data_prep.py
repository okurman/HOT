import subprocess
import tempfile
import warnings
warnings.filterwarnings('ignore')

import sys
import os
from os.path import join, basename
import configparser
import pandas
import numpy as np

###############################################################

config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)
code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool

BINS_DIR = ROOT_DIR+"/definitions/peak_8bp_v2/log_bins/"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom" % (x, i)) for i in range(14)]

PROJECT_DIR = ROOT_DIR+"overbinders/"
WORK_DIR = join(PROJECT_DIR, "DL_analysis/datasets")

all_dhs_file = "..."
genome_file = "hg19_22.genome"

val_chroms = ["chr7"]
test_chroms = ["chr8", "chr9"]
train_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

from lib.data_io.tools import get_chrom2seq, seq2one_hot
import h5py


def create_dataset_multitask_3_classes(cl="K562", chrom2seq=None, seq_len=400):

    if not chrom2seq:
        chrom2seq = get_chrom2seq()
        # return chrom2seq

    save_file = join(WORK_DIR, "%s_3class_%d.hdf5" % (cl, seq_len))
    print("Output file to:", save_file)

    loci_files = get_loci_files(cl)

    cat_files = [tempfile.NamedTemporaryFile(), tempfile.NamedTemporaryFile(), tempfile.NamedTemporaryFile()]
    cat_ranges = [loci_files[:5], loci_files[5:10], loci_files[10:]]

    for cat_file, cat_range in zip(cat_files, cat_ranges):
        cmd = "cat %s | sort -k1,1 -k2,2n > %s" % (" ".join(cat_range), cat_file.name)
        subprocess.call(cmd, shell=True)

    cat_loci = [BedTool(f.name) for f in cat_files]
    all_pos_loci = BedTool([r for loci in cat_loci for r in loci])

    all_pos_file = tempfile.NamedTemporaryFile()
    all_pos_loci.saveas(all_pos_file.name)

    control_loci_file = tempfile.NamedTemporaryFile()
    cmd = "bedtools shuffle -i %s -g %s -incl %s -excl %s > %s" % (all_pos_file.name,
                                                                   genome_file,
                                                                   all_dhs_file,
                                                                   all_pos_file.name,
                                                                   control_loci_file.name)
    subprocess.call(cmd, shell=True)

    control_loci = BedTool(control_loci_file.name).subtract(all_pos_loci)

    train_loci = [[], [], [], []]
    validation_loci = [[], [], [], []]
    test_loci = [[], [], [], []]

    for i, loci in enumerate(cat_loci):

        for r in loci:
            if r.chrom in train_chroms:
                train_loci[i].append(r)
            elif r.chrom in test_chroms:
                test_loci[i].append(r)
            elif r.chrom == val_chroms[0]:
                validation_loci[i].append(r)

    for r in control_loci:
        if r.chrom in train_chroms:
            train_loci[3].append(r)
        elif r.chrom in test_chroms:
            test_loci[3].append(r)
        elif r.chrom == val_chroms[0]:
            validation_loci[3].append(r)

    with h5py.File(save_file, "w") as of:

        sub_cat_loci = [train_loci, validation_loci, test_loci]
        prefixes = ["train", "validation", "test"]
        labels = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]]

        for _cat_loci, prefix in zip(sub_cat_loci, prefixes):
            print(prefix)
            total_regions = sum([len(_) for _ in _cat_loci])

            X = of.create_dataset("%s_data" % prefix, (total_regions, seq_len, 4), dtype=np.int8)
            Y = of.create_dataset("%s_labels" % prefix, (total_regions, 3), dtype=np.int8)

            cnt = 0
            for loci, label in zip(_cat_loci, labels):
                for r in loci:
                    flank_len = (seq_len - r.length)//2
                    _seq = chrom2seq[r.chrom][r.start-flank_len: r.stop+flank_len]

                    if not len(_seq) == seq_len:
                        print("Skipping:", r, "sequence length not equal")
                        continue

                    _mat = seq2one_hot(_seq)
                    X[cnt, :, :] = _mat
                    Y[cnt, :] = label

                    cnt += 1


def create_dataset_multitask(cl="K562", chrom2seq=None, seq_len=400):

    if not chrom2seq:
        chrom2seq = get_chrom2seq()
        # return chrom2seq

    save_file = join(WORK_DIR, "%s_14class_%d_tmp.hdf5" % (cl, seq_len))
    print("Output file to:", save_file)
    loci_files = get_loci_files(cl)

    ######### generate control regions
    all_pos_loci = BedTool([r for loci_file in loci_files for r in BedTool(loci_file)])
    all_pos_file = tempfile.NamedTemporaryFile()
    all_pos_loci.saveas(all_pos_file.name)

    control_loci_file = tempfile.NamedTemporaryFile()
    cmd = "bedtools shuffle -i %s -g %s -incl %s -excl %s > %s" % (all_pos_file.name,
                                                                   genome_file,
                                                                   all_dhs_file,
                                                                   all_pos_file.name,
                                                                   control_loci_file.name)
    subprocess.call(cmd, shell=True)
    control_loci = BedTool(control_loci_file.name).subtract(all_pos_loci)
    ####################################################################################

    print("Pos loci:", all_pos_loci.count())
    print("Control loci:", control_loci.count())

    train_loci = [[] for _ in range(len(loci_files)+1)]
    validation_loci = [[] for _ in range(len(loci_files)+1)]
    test_loci = [[] for _ in range(len(loci_files)+1)]

    for i, loci_file in enumerate(loci_files):

        loci = BedTool(loci_file)

        for r in loci:
            if r.chrom in train_chroms:
                train_loci[i].append(r)
            elif r.chrom in test_chroms:
                test_loci[i].append(r)
            elif r.chrom == val_chroms[0]:
                validation_loci[i].append(r)

    for r in control_loci:
        if r.chrom in train_chroms:
            train_loci[-1].append(r)
        elif r.chrom in test_chroms:
            test_loci[-1].append(r)
        elif r.chrom == val_chroms[0]:
            validation_loci[-1].append(r)

    with h5py.File(save_file, "w") as of:

        sub_cat_loci = [train_loci, validation_loci, test_loci]
        prefixes = ["train", "validation", "test"]

        for _cat_loci, prefix in zip(sub_cat_loci, prefixes):
            print(prefix)

            total_regions = sum([len(_) for _ in _cat_loci])
            X = of.create_dataset("%s_data" % prefix, (total_regions, seq_len, 4), dtype=np.int8)
            Y = of.create_dataset("%s_labels" % prefix, (total_regions, len(loci_files)), dtype=np.int8)

            cnt = 0
            for label_ind, loci in enumerate(_cat_loci):
                label = np.zeros(len(loci_files), dtype=np.int8)
                if label_ind < len(loci_files):
                    label[label_ind] = 1

                print(label_ind, len(loci), label)

                for r in loci:
                    flank_len = (seq_len - r.length) // 2
                    _seq = chrom2seq[r.chrom][r.start - flank_len: r.stop + flank_len]

                    if not len(_seq) == seq_len:
                        print("Skipping:", r, "sequence length not equal")
                        continue

                    _mat = seq2one_hot(_seq)
                    X[cnt, :, :] = _mat
                    Y[cnt, :] = label
                    cnt += 1

                print(cnt)



if __name__ == "__main__":

    # create_dataset_multitask_3_classes(cl=sys.argv[1], seq_len=int(sys.argv[2]))
    create_dataset_multitask(cl=sys.argv[1], seq_len=int(sys.argv[2]))
    # create_dataset_multitask(cl="HepG2")
