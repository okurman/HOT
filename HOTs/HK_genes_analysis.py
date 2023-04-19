import gzip
import warnings
from collections import defaultdict, Counter
from random import shuffle

import h5py
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
from lib import phastCons, phyloP
import seaborn as sns
from upsetplot import from_memberships
from upsetplot import plot as upset_plot

from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"
LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"
TSS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/TSS.bed.gz"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"
HK_GENES_FILE = "/panfs/pan1/devdcode/sanjar/housekeeping_genes/Housekeeping_GenesHuman.csv"

SE_DIR = "/panfs/pan1/devdcode/sanjar/superenhancers/"


def extract_nearest_genes(regions, window=False, distance=25000):

    """ONLY FOR ENHANCERS"""

    nearest_tss = regions.closest(PROMOTERS_FILE).groupby(c=7, o="collapse")

    if window:
        _window = regions.window(TSS_FILE, w=distance).groupby(c=7, o="collapse")
        nearest_tss = nearest_tss.cat(_window, postmerge=False).sort().groupby(c=4, o="distinct")

    return nearest_tss


def extract_HK_fractions():

    hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed")).merge(d=-1)
    hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed")).merge(d=-1)
    hots = hot_proms.cat(hot_enhs).merge(d=-1)

    f = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac.bed"
    re = BedTool(f).intersect(hots, v=True)

    f = "/panfs/pan1/devdcode/sanjar/superenhancers/HepG2.bed"
    se = BedTool(f).sort().merge()

    hk_genes = set([l.split(";")[0] for l in open(HK_GENES_FILE).readlines()[1:]])

    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True, wa=True)

    names = ["HOT promoters", "HOT enhancers", "Regular promoters",  "Regular enhancers", "Super-enhancers"]
    regions_list = [hot_proms, hot_enhs, rp, re, se]

    nearest_genes_list = [hot_proms.intersect(PROMOTERS_FILE, wo=True).groupby(c=7, o="distinct"),
                          extract_nearest_genes(hot_enhs),
                          rp,
                          extract_nearest_genes(re),
                          extract_nearest_genes(se)]

    data = []

    for name, regions, nearest_genes in zip(names, regions_list, nearest_genes_list):
        print(name)
        assert(regions.count() == nearest_genes.count())

        nearest_tss_pool = set([_.split(".")[0] for r in nearest_genes for _ in r.fields[3].split(",")])
        found_hk_genes = set(hk_genes).intersection(nearest_tss_pool)
        hk_count = 0

        for r in nearest_genes:
            target_genes = r.fields[3].split(",")
            target_genes = set([g.split(".")[0] for g in target_genes])
            _ovlp = target_genes.intersection(hk_genes)
            hk_count += 1 if _ovlp else 0
            found_hk_genes.update(_ovlp)

        hk_fr = len(found_hk_genes)/len(hk_genes)
        regions_fr = hk_count/regions.count()
        data.append([name, hk_fr, regions_fr])

    return data


def plot_fractions(data):

    # data = extract_HK_fractions()

    data = np.asarray(data)
    names = data[:, 0]
    values = 100*np.asarray(data[:, 1:], dtype=float)

    fig, axs = plt.subplots(2, 1, figsize=(3, 5))

    x_range = np.arange(len(names))
    bars = axs[0].bar(x_range, values[:, 0])
    axs[0].bar_label(bars, label=values[:, 0], fmt="%d%%", size=9)
    axs[0].set_ylabel("% of regulated \n HK genes")
    axs[0].set_ylim([0, 60])
    axs[0].set_xticklabels([])
    axs[0].grid(axis="y")

    bars = axs[1].bar(x_range, values[:, 1], color="green")
    axs[1].bar_label(bars, label=values[:, 1], fmt="%d%%", size=9)
    axs[1].set_ylabel("% of loci\n regulating HK genes")
    axs[1].set_ylim([0, 20])
    axs[1].grid(axis="y")

    axs[1].set_xticks(x_range)
    axs[1].set_xticklabels(names, rotation=45, horizontalalignment="right")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "barplot_houskeeping_genes.pdf")
    print(save_file)
    plt.savefig(save_file)


def plot_fractions_merged(data):

    # data = extract_HK_fractions()

    data = np.asarray(data)
    names = data[:, 0]
    hk_fracs = 100 * np.asarray(data[:, 1], dtype=float)
    loci_fracs = 100 * np.asarray(data[:, 2], dtype=float)

    data_list = []

    for n, hk_frac, l_frac in zip(names, hk_fracs, loci_fracs):
        data_list.append([n, "% of HK genes", hk_frac])
        data_list.append([n, "% of regions", l_frac])

    df = pandas.DataFrame(data=data_list, columns=["cat", "frac", "value"])
    _order = ['HOT promoters', 'Regular promoters', 'HOT enhancers', 'Regular enhancers', 'Super-enhancers']
    x_labels = ["HOT prom", "reg prom", "HOT enh", "reg enh", "super-enh"]

    plt.figure(figsize=(4, 2.6))
    ax = plt.gca()
    sns.barplot(x="cat", y="value", hue="frac", data=df, ax=ax, palette="tab10", order=_order)
    plt.legend(frameon=False)
    ax.grid(axis="y", alpha=0.5)
    ax.set_ylim([0, 95])
    for bars in ax.containers:
        ax.bar_label(bars, fmt='%.0f%%', fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("%")
    ax.set_xticklabels(x_labels, rotation=45, ha="right")

    plt.tight_layout()
    # plt.show()

    save_file = join(PLOTS_DIR, "barplot_houskeeping_genes_merged.pdf")
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def tau_data():

    tau_file = "ROOT_DIR/gene_expressions/tau/Tau_transcript_V8.csv"

    enst2tau = {}

    for l in open(tau_file):
        parts = l.split(",")
        try:
            enst2tau[parts[0]] = float(parts[1])
        except:
            continue

    re_file = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac.bed"
    se_file = "/panfs/pan1/devdcode/sanjar/superenhancers/HepG2.bed"

    hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed")).merge(d=-1)
    hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed")).merge(d=-1)

    re = BedTool(re_file).intersect(hot_enhs, v=True, wa=True)
    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True, wa=True)
    se = BedTool(se_file).sort().merge()

    names = ["HOT promoters", "HOT enhancers", "Regular promoters",  "Regular enhancers", "Super-enhancers"]

    nearest_genes_list = [hot_proms.intersect(PROMOTERS_FILE, wo=True).groupby(c=7, o="distinct"),
                          extract_nearest_genes(hot_enhs),
                          rp,
                          extract_nearest_genes(re),
                          extract_nearest_genes(se)]

    data = []

    for name, nearest_genes in zip(names, nearest_genes_list):
        print(name)
        nearest_tss_pool = set([_.split(".")[0] for r in nearest_genes for _ in r.fields[3].split(",")])

        for enst in nearest_tss_pool:
            if enst in enst2tau:
                data.append([name, enst2tau[enst]])

    df = pandas.DataFrame(data=data, columns=["category", "tau"])
    return df


def plot_tau_single_cl(df):

    # df = tau_data_single_cl()

    _order = ['HOT promoters', 'Regular promoters', 'HOT enhancers',  'Regular enhancers', 'Super-enhancers']

    plt.figure(figsize=(3.8, 2.7))
    ax = plt.gca()
    g = sns.violinplot(data=df, x="category", order=_order, y="tau", hue="category", dodge=False, legend=False)
    g.set_ylim([0, 1.25])
    g.legend_.remove()
    g.set_xlabel("")

    # ax.set_xticklabels(["HOT\npromoters", "HOT\nenhancers", "regular\nenhancers", "regular\npromoters"])
    ax.set_xticklabels(["HOT prom", "reg prom", "HOT enh", "reg enh", "super-enh"], rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.5)
    ax.set_ylim([0, 1.4])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    plt.tight_layout()
    save_file = join(PLOTS_DIR, "HepG2_tau_plot.pdf")
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')



if __name__ == "__main__":

    pass






