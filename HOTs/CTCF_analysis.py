import warnings
from scipy.stats import wilcoxon, mannwhitneyu, kruskal, spearmanr

warnings.filterwarnings('ignore')

import sys
import os
from os.path import join

import numpy as np
import configparser

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
import pandas
import matplotlib.pyplot as plt

import seaborn as sns
import pandas as pd
from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

TF_CLASS_FILE = "ROOT_DIR/chipseq_files/human_transcription_factors.txt"


def get_TF_distributions():

    # nearest_tss_loci = BedTool(join(LOCI_DIR, "HepG2_HOTs.noprom.nearest_tss"))

    tf_ids = ["ENCFF648HFC", "ENCFF049KLQ"]
    tf_names = ["CTCF", "Cohesin"]

    data = []

    for tf_name, tf_id in zip(tf_names, tf_ids):

        hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.bed"))
        proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed"))
        enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.noproms.bed"))

        for name, loci in zip(["HOT loci", "HOT enh", "HOT prom"], [hots, enhs, proms]):

            ctcf = BedTool([r for r in loci if tf_id in r.fields[8]])
            no_ctcf = BedTool([r for r in loci if tf_id not in r.fields[8]])

            ctcf_tfs = [int(r.fields[5]) for r in ctcf]
            no_ctcf_tfs = [int(r.fields[5]) for r in no_ctcf]

            data += [[name, tf_name, "bound", _] for _ in ctcf_tfs]
            data += [[name, tf_name, "not bound", _] for _ in no_ctcf_tfs]

            print("%s  %s  bound_tfs: %d, non-bound_tfs: %d" % (tf_name, name, np.mean(ctcf_tfs), np.mean(no_ctcf_tfs)))

    df = pd.DataFrame(data=data, columns=["locus", "tf", "status", "tfs"])
    return df


def plot_TF_distributions(df):

    plt.figure(figsize=(4, 2.5))
    ax = plt.gca()

    df = df[df["locus"] == "HOT loci"]

    g = sns.boxplot(x="tf", y="tfs", hue="status", data=df, showfliers=False, palette="Set1", ax=ax)

    g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=2, frameon=False)
    g.set_xlabel("")
    g.set_ylabel("bound TFs")

    plt.grid(which="both", axis="y", alpha=0.5, linestyle="--")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "CTCF_bound_tfs.pdf")
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def get_bound_TFs_stats(df):

    df = df[df["locus"] == "HOT loci"]

    for tf in ["CTCF", "Cohesin"]:
        print(tf)
        _df = df[df["tf"] == tf]
        bound_tfs = _df[_df["status"]=="bound"].tfs
        not_bound_tfs = _df[_df["status"]=="not bound"].tfs

        print("Bound TFs: %f  Not bound TFs: %f  Ratio: %f p-value %f "
              % (bound_tfs.mean(),
                 not_bound_tfs.mean(),
                 bound_tfs.mean()/not_bound_tfs.mean(),
                 mannwhitneyu(bound_tfs, not_bound_tfs).pvalue))
        print("")


def get_nearest_TSS_distances():

    tss = BedTool(join(LOCI_DIR, "HepG2_HOTs.noprom.nearest_tss"))

    tf_ids = ["ENCFF648HFC", "ENCFF049KLQ"]
    tf_names = ["CTCF", "Cohesin"]

    data = []

    for tf_name, tf_id in zip(tf_names, tf_ids):

        enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.noproms.bed"))

        bound = BedTool([r for r in enhs if tf_id in r.fields[8]]).intersect(tss, wb=True)
        not_bound = BedTool([r for r in enhs if tf_id not in r.fields[8]]).intersect(tss, wb=True)

        bound_dist = set([int(r.fields[-1]) for r in bound])
        not_bound_dist = set([int(r.fields[-1]) for r in not_bound])

        data += [[tf_name, "bound", _] for _ in bound_dist]
        data += [[tf_name, "not bound", _] for _ in not_bound_dist]

        print("%s  bound_tfs: %d, non-bound_tfs: %d" % (tf_name, np.mean(list(bound_dist)), np.mean(list(not_bound_dist))))

    df = pd.DataFrame(data=data, columns=["tf", "status", "dist"])
    return df


def plot_nearest_TSS_distances(df):

    plt.figure(figsize=(4, 2.5))
    ax = plt.gca()

    g = sns.boxplot(x="tf", y="dist", hue="status", data=df, showfliers=False, palette="Set1", ax=ax)

    g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=2, frameon=False)
    g.set_xlabel("")
    g.set_ylabel("nearest TSS (kbps)")
    g.set_yticklabels(["", "0", "20", "40", "60", "80"])

    plt.grid(which="both", axis="y", alpha=0.5, linestyle="--")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "CTCF_bound_nearest_TSS.pdf")
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def get_nearest_TSS_stats(df):

    for tf in ["CTCF", "Cohesin"]:
        print(tf)
        _df = df[df["tf"] == tf]
        bound_tfs = _df[_df["status"]=="bound"].dist
        not_bound_tfs = _df[_df["status"]=="not bound"].dist

        pvalue = mannwhitneyu(bound_tfs, not_bound_tfs).pvalue

        bound_tfs = bound_tfs.mean()
        not_bound_tfs = not_bound_tfs.mean()

        print("Bound TFs: %f  Not bound TFs: %f  Ratio: %f Diff: %d Perc: %f p-value %f "
              % (bound_tfs.mean(),
                 not_bound_tfs.mean(),
                 bound_tfs/not_bound_tfs,
                 not_bound_tfs - bound_tfs,
                 (not_bound_tfs - bound_tfs)/not_bound_tfs,
                 pvalue))
        print("")
