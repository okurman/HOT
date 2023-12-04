#!/usr/bin/env python

import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import warnings
warnings.filterwarnings('ignore')
import random
from os.path import join
import numpy as np
from pybedtools import BedTool
import h5py
from pathlib import Path

DATA_PATH = Path(os.environ["HOT_DATA"])
BINS_DIR = DATA_PATH/"log_bins"
HOTS_DIR = DATA_PATH/"HOTs"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

FASTA_FILE = DATA_PATH/"src_files/hg19_files/hg19.fa"
VAL_CHROMS = ["chr6", "chr7"]
TEST_CHROMS = ["chr8", "chr9"]
TRAIN_CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

from data_prep.basic import get_chrom2seq, seq2one_hot, download_hg19_fasta

SAVE_DIR = DATA_PATH / "classification/datasets"
SAVE_DIR.mkdir(exist_ok=True, parents=True)

(SAVE_DIR/"fasta_datasets").mkdir(exist_ok=True)
(SAVE_DIR/"features_datasets").mkdir(exist_ok=True)
(SAVE_DIR/"one_hot_datasets").mkdir(exist_ok=True)

random.seed(1024)
import multiprocessing as mp


def create_dataset_binary(cl="HepG2", chrom2seq=None, seq_len=400, control_type="dhs"):

    assert control_type in ["dhs", "proms", "re"]

    if not chrom2seq:
        print("Loading hg19")
        chrom2seq = get_chrom2seq(FASTA_FILE)

    hots_file = join(SAVE_DIR, f"control_regions/{cl}_hots_{seq_len}.bed.gz")
    controls_file = join(SAVE_DIR, f"control_regions/{cl}_{control_type}_controls_{seq_len}.bed.gz")
    save_file = join(SAVE_DIR, f"one_hot_datasets/{cl}_binary_{control_type}_{seq_len}.hdf5")

    hots = BedTool(hots_file)
    controls = BedTool(controls_file)

    train_pool = [BedTool([r for r in hots if r.chrom in TRAIN_CHROMS]),
                  BedTool([r for r in controls if r.chrom in TRAIN_CHROMS])]
    val_pool =   [BedTool([r for r in hots if r.chrom in VAL_CHROMS]),
                  BedTool([r for r in controls if r.chrom in VAL_CHROMS])]
    test_pool =  [BedTool([r for r in hots if r.chrom in TEST_CHROMS]),
                  BedTool([r for r in controls if r.chrom in TEST_CHROMS])]

    Y_train = np.concatenate((np.ones(train_pool[0].count()), np.zeros(train_pool[1].count())))
    Y_val = np.concatenate((np.ones(val_pool[0].count()), np.zeros(val_pool[1].count())))
    Y_test = np.concatenate((np.ones(test_pool[0].count()), np.zeros(test_pool[1].count())))

    train_bed = train_pool[0].cat(train_pool[1], postmerge=False)
    val_bed = val_pool[0].cat(val_pool[1], postmerge=False)
    test_bed = test_pool[0].cat(test_pool[1], postmerge=False)

    print(f"\t\t\tSaving to: {save_file}")
    with h5py.File(save_file, "w") as of:

        prefixes = ["train", "validation", "test"]

        for prefix, bed, Y in zip(prefixes,
                                [train_bed, val_bed, test_bed],
                                [Y_train, Y_val, Y_test]):

            X = np.zeros((len(Y), seq_len, 4))

            for i, r in enumerate(bed):
                seq = chrom2seq[r.chrom][r.start: r.stop]
                one_hot = seq2one_hot(seq)
                try:
                    X[i, :, :] = one_hot
                except:
                    print("skipping", i, r, len(seq))
                    continue

            of.create_dataset(name="%s_data" % prefix, data=X)
            of.create_dataset(name="%s_labels" % prefix, data=Y)


def create_dataset_binary_features(cl="HepG2", chrom2seq=None, control_type="dhs"):

    seq_len = 400

    if not chrom2seq:
        print("Loading hg19")
        chrom2seq = get_chrom2seq(FASTA_FILE)

    hots_file = join(SAVE_DIR, f"control_regions/{cl}_hots_{seq_len}.bed.gz")
    controls_file = join(SAVE_DIR, f"control_regions/{cl}_{control_type}_controls_{seq_len}.bed.gz")
    save_file = join(SAVE_DIR, f"features_datasets/{cl}_binary_{control_type}_{seq_len}_features.hdf5")

    hots = BedTool(hots_file)
    controls = BedTool(controls_file)

    cpg_map = get_cpg_map(hots)
    cpg_map.update(get_cpg_map(controls))

    train_pool = [BedTool([r for r in hots if r.chrom in TRAIN_CHROMS]),
                  BedTool([r for r in controls if r.chrom in TRAIN_CHROMS])]
    val_pool = [BedTool([r for r in hots if r.chrom in VAL_CHROMS]),
                BedTool([r for r in controls if r.chrom in VAL_CHROMS])]
    test_pool = [BedTool([r for r in hots if r.chrom in TEST_CHROMS]),
                 BedTool([r for r in controls if r.chrom in TEST_CHROMS])]

    Y_train = np.concatenate((np.ones(train_pool[0].count()), np.zeros(train_pool[1].count())))
    Y_val = np.concatenate((np.ones(val_pool[0].count()), np.zeros(val_pool[1].count())))
    Y_test = np.concatenate((np.ones(test_pool[0].count()), np.zeros(test_pool[1].count())))

    train_bed = train_pool[0].cat(train_pool[1], postmerge=False)
    val_bed = val_pool[0].cat(val_pool[1], postmerge=False)
    test_bed = test_pool[0].cat(test_pool[1], postmerge=False)

    print(f"\t\t\tSaving to: {save_file}")
    with h5py.File(save_file, "w") as of:

        prefixes = ["train", "validation", "test"]

        for prefix, bed, Y in zip(prefixes,
                                [train_bed, val_bed, test_bed],
                                [Y_train, Y_val, Y_test]):

            X = np.zeros((len(Y), 4))

            for i, r in enumerate(bed):
                seq = chrom2seq[r.chrom][r.start: r.stop]
                X[i, :] = get_features(r, seq, cpg_map)

            of.create_dataset(name="%s_data" % prefix, data=X, compression="gzip")
            of.create_dataset(name="%s_labels" % prefix, data=Y, compression="gzip")


def get_cpg_map(loci):

    cpg = BedTool(join(DATA_PATH, "src_files/hg19_files/cpgIslandExt.bed"))

    loci_cpg = loci.intersect(cpg, wo=True)
    loci_non_cpg = loci.intersect(cpg, v=True)

    cpg_map = dict()
    for r in loci_cpg:
        cpg_map[f"{r.chrom}-{r.start}-{r.stop}"] = int(r.fields[-1]) / r.length
    for r in loci_non_cpg:
        cpg_map[f"{r.chrom}-{r.start}-{r.stop}"] = 0

    return cpg_map


def get_features(r, seq, cpg_map):

    seq = str(seq).upper()
    r_id = f"{r.chrom}-{r.start}-{r.stop}"

    gc = (seq.count("G") + seq.count("C"))/len(seq)
    cpg = seq.count("CG")
    gpc = seq.count("GC")
    cpg_cov = cpg_map[r_id]

    return np.asarray([gc, cpg, gpc, cpg_cov])


def create_dataset_binary_fasta(cl="HepG2", seq_len=400, control_type="dhs"):

    hots_file = join(SAVE_DIR, f"control_regions/{cl}_hots_{seq_len}.bed.gz")
    controls_file = join(SAVE_DIR, f"control_regions/{cl}_{control_type}_controls_{seq_len}.bed.gz")
    save_dir = SAVE_DIR/f"fasta_datasets/{cl}_binary_{control_type}_{seq_len}_lsgkm"
    save_dir.mkdir(exist_ok=True, parents=True)

    hots = BedTool(hots_file)
    controls = BedTool(controls_file)

    hots_train = BedTool([r for r in hots if r.chrom in TRAIN_CHROMS + VAL_CHROMS])
    hots_test = BedTool([r for r in hots if r.chrom in TEST_CHROMS])

    controls_train = BedTool([r for r in controls if r.chrom in TRAIN_CHROMS + VAL_CHROMS])
    controls_test = BedTool([r for r in controls if r.chrom in TEST_CHROMS])

    fnames = ["hots_train.fa", "hots_test.fa", "controls_train.fa", "controls_test.fa"]
    for fname, bed in zip(fnames, [hots_train, hots_test, controls_train, controls_test]):
        save_file = join(save_dir, fname)
        print(save_file)
        bed.sequence(fi=str(FASTA_FILE), fo=save_file)


if __name__ == "__main__":

    if not FASTA_FILE.exists():
        download_hg19_fasta(FASTA_FILE)
        quit()

    cl = sys.argv[1]
    with mp.Pool() as pool:
        results = []
        for ctr in ["proms", "re", "dhs"]:

            for seq_len in [400, 1000]:
                _kwargs = {"cl": cl, "seq_len": seq_len, "control_type": ctr}
                results.append(pool.apply_async(create_dataset_binary, (), _kwargs))
                results.append(pool.apply_async(create_dataset_binary_fasta, (), _kwargs))

            _kwargs = {"cl": cl, "control_type": ctr}
            results.append(pool.apply_async(create_dataset_binary_features, (), _kwargs))

        for r in results:
            r.wait()
        pool.close()
        pool.join()

    ####################################################################################################
    #### To run individually, in a non-parallel way uncomment the following part:  #####################
    ####################################################################################################

    # chrom2seq = get_chrom2seq(FASTA_FILE)
    # for ctr in ["proms", "re", "dhs"]:
    #     print(f"\t{ctr}")
    #
    #     print("\t\tone_hot dataset")
    #     for seq_len in [400, 1000]:
    #         print(f"\t\t\t{seq_len}")
    #         create_dataset_binary(cl=cl, seq_len=seq_len, control_type=ctr, chrom2seq=chrom2seq)
    #
    #     print("\t\tfasta dataset")
    #     for seq_len in [400, 1000]:
    #         print(f"\t\t\t{seq_len}")
    #         create_dataset_binary_fasta(cl=cl, seq_len=seq_len, control_type=ctr)
    #
    #     # if seq_len == 400:
    #     print("\t\tsequence features dataset")
    #     create_dataset_binary_features(cl=cl, chrom2seq=chrom2seq, control_type=ctr)
