
import warnings
from collections import Counter, defaultdict
import seaborn as sns
import pandas
from scipy.stats import mannwhitneyu, kruskal, kstest

warnings.filterwarnings('ignore')

import sys
import os
from os.path import join, exists

import numpy as np
import configparser

###############################################################

config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool

import matplotlib.pyplot as plt

INTRONS_FILE ="/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/knownGene.introns.merged.bed.gz"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"
# PRIM_PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.primary_TSS.liftOver_hg19.bed"

PROJECT_DIR = "ROOT_DIR/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"
LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"


def get_overlaps_subsample_h1():

    print("merged")
    h1_hots = BedTool(join(LOCI_DIR, "H1_HOTs.bed")).merge(d=-1)
    h1_count = h1_hots.count()

    ovp_list = []

    d = "ROOT_DIR/definitions/peak_8bp_v2/log_bins_random/H1/"
    for i in range(10):
        hepg2_hots = BedTool(join(d, "HepG2", "%d_HOTs.bed.gz" % i))
        k562_hots = BedTool(join(d, "K562", "%d_HOTs.bed.gz" % i))
        cl_hots = hepg2_hots.cat(k562_hots).sort().merge()
        ovlp_fr = h1_hots.intersect(cl_hots, wa=True).merge(d=-1).count() / h1_count
        ovp_list.append(ovlp_fr)

    return ovp_list


def get_overlaps_H1():

    loci_dir = "ROOT_DIR/definitions/peak_8bp_v2/log_bins_random/common_TFs/"
    hepg2 = BedTool(join(loci_dir, "HepG2/HOTs.bed.gz"))
    k562 = BedTool(join(loci_dir, "K562/HOTs.bed.gz"))
    h1 = BedTool(join(loci_dir, "H1/HOTs.bed.gz"))

    common_TFs_ovp = h1.intersect(hepg2.cat(k562), wa=True).merge(d=-1).count()/h1.count()

    loci_dir = "ROOT_DIR/definitions/peak_8bp_v2/HOTs/"
    hepg2 = BedTool(join(loci_dir, "HepG2_HOTs.bed"))
    k562 = BedTool(join(loci_dir, "K562_HOTs.bed"))
    h1 = BedTool(join(loci_dir, "H1_HOTs.bed"))

    all_TFs_ovp = h1.intersect(hepg2.cat(k562), wa=True).merge(d=-1).count() / h1.count()

    return common_TFs_ovp, all_TFs_ovp


def load_data():

    data = []

    common_TFs_ovp, all_TFs_ovp = get_overlaps_H1()

    data += [
        ["complete", all_TFs_ovp],
        ["common DAPs", common_TFs_ovp]
    ]

    print("H1 size")
    ovp_values = get_overlaps_subsample_h1()
    data += [["subsampled to H1", _] for _ in ovp_values]

    df = pandas.DataFrame(data=data, columns=["sample", "value"])
    df["value"] *= 100
    return df


def get_overlaps_subsample_perc(perc=20):

    ovp_list = []

    d = "ROOT_DIR/definitions/peak_8bp_v2/log_bins_random/%s/" % str(perc)
    for i in range(10):
        f1 = join(d, "HepG2", "%d_HOTs.bed.gz" % i)
        f2 = join(d, "K562", "%d_HOTs.bed.gz" % i)
        f3 = join(d, "H1", "%d_HOTs.bed.gz" % i)

        hepg2_hots = BedTool(f1).merge(d=-1)
        k562_hots = BedTool(f2).merge(d=-1)
        h1_hots = BedTool(f3).merge(d=-1)

        cl_hots = hepg2_hots.cat(k562_hots).sort().merge()
        ovlp_fr = h1_hots.intersect(cl_hots, wa=True).merge(d=-1).count() / h1_hots.count()

        ovp_list.append(ovlp_fr)

    return ovp_list


def load_data_perc():

    data = []

    for perc in [20, 30, 40, 50, 60, 70]:
        print(perc)
        ovp_list = get_overlaps_subsample_perc(perc)
        data += [[perc, _] for _ in ovp_list]

    df = pandas.DataFrame(data=data, columns=["sample", "value"])
    df["value"] *= 100

    return df


def plot_subsample_h1(df):

    # df = load_data()

    # data = [["H1 size", _] for _ in ovp_list]
    # data += [["All DAPs", 0.84]]
    # df = pandas.DataFrame(data=data, columns=["type", "value"])
    # df["value"] *= 100

    plt.figure(figsize=(3.5, 3))
    ax = plt.gca()
    # sns.barplot(x="type", y="value", data=df, ci="sd", ax=ax, order=["All DAPs", "H1 size"])
    sns.barplot(x="sample", y="value", data=df, ci="sd", ax=ax, order=['complete', 'common DAPs', 'subsampled to H1'])
    ax.set_ylim([0, 100])
    ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax.set_xticklabels(['all\nDAPs', 'common\nDAPs', 'subsampled\nto H1'])
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.4)
    plt.tight_layout()
    save_file = join(PLOTS_DIR, "H1_complete_common_subsampled.pdf")
    plt.savefig(save_file, bbox_inches='tight')
    print(save_file)
    plt.close()


def plot_subsample_overlaps(df):

    plt.figure(figsize=(3, 3))
    ax = plt.gca()
    sns.boxplot(x="sample", y="value", data=df, ax=ax, order=["H1 size", "20%", "50%"])
    ax.set_ylim([0, 100])
    ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax.set_ylabel("")
    ax.set_xlabel("")
    plt.tight_layout()
    ax.grid(axis="y", alpha=0.7, linestyle='--')
    save_file = join(PLOTS_DIR, "H1_subsampled_overlaps_3vers.pdf")
    plt.savefig(save_file, bbox_inches='tight')
    print(save_file)
    plt.close()





