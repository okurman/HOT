#!/usr/bin/env python

import argparse
import warnings
warnings.filterwarnings('ignore')
import random
from pathlib import Path
import subprocess as sp

DATA_PATH = Path("../data")
DATASETS_DIR = DATA_PATH/"classification/datasets/fasta_datasets"
MODELS_DIR = DATA_PATH/"classification/models/lsgkm_models"
MODELS_DIR.mkdir(exist_ok=True, parents=True)

random.seed(1024)


def run_training(args):

    cl = args.cl
    ctr = args.ctr
    seq_len = args.len

    data_dir = DATASETS_DIR / f"{cl}_binary_{ctr}_{seq_len}_lsgkm"
    save_dir = MODELS_DIR / f"{cl}_{ctr}_{seq_len}"
    save_dir.mkdir(exist_ok=True)

    pos_file = data_dir/"hots_train.fa"
    neg_file = data_dir/"controls_train.fa"
    model_prefix = save_dir/f"model.kernel_{args.kernel}"

    cmd_fmt = "gkmtrain -t {} -T {} {} {} {}"
    command_str = cmd_fmt.format(args.kernel, args.threads, pos_file, neg_file, model_prefix)
    sp.call(command_str, shell=True)

    pos_file = data_dir / "hots_test.fa"
    pos_scores_file = save_dir / "hots_test_scores.txt"
    neg_file = data_dir / "controls_test.fa"
    neg_scores_file = save_dir / "controls_test_scores.txt"
    model_file = str(model_prefix) + ".svm.txt"

    for seq_file, score_file in zip([pos_file, neg_file], [pos_scores_file, neg_scores_file]):
        command_str = "gkmpredict -T {} {} {} {}".format(args.threads, seq_file, model_file, score_file)
        sp.call(command_str, shell=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-cl',
                        choices=["HepG2", "K562"],
                        type=str,
                        default="HepG2")
    parser.add_argument('-ctr',
                        choices=["proms", "dhs", "re"],
                        type=str,
                        default="dhs")
    parser.add_argument('-len',
                        choices=[400, 1000],
                        type=int,
                        default=400)
    parser.add_argument('-kernel',
                        choices=range(6),
                        type=int,
                        default=4)
    parser.add_argument('-threads',
                        choices=[1, 4, 16],
                        type=int,
                        default=4)
    args = parser.parse_args()

    print("argparse arguments:")
    for arg in vars(args):
        print(f"{arg}:\t{getattr(args, arg)}")
    print("\n")

    run_training(args)

