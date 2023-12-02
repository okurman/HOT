import gzip
from os.path import join
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from pybedtools import BedTool
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

import os
DATA_PATH = Path(os.environ["HOT_DATA"])
TSS_FILE = DATA_PATH / "src_files/hg19_files/TSS.bed.gz"
PROMOTERS_FILE = DATA_PATH / "src_files/hg19_files/promoters.merged.bed.gz"
HK_GENES_FILE = DATA_PATH / "src_files/Housekeeping_GenesHuman.csv"
TAU_FILE = DATA_PATH / "src_files/Tau_gene_V8.csv.gz"


def extract_nearest_genes(regions, window=False, distance=25000):

    """ONLY FOR ENHANCERS"""

    nearest_tss = regions.closest(str(PROMOTERS_FILE)).groupby(c=7, o="collapse")

    if window:
        _window = regions.window(TSS_FILE, w=distance).groupby(c=7, o="collapse")
        nearest_tss = nearest_tss.cat(_window, postmerge=False).sort().groupby(c=4, o="distinct")

    return nearest_tss


def extract_HK_fractions():

    hot_proms = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.proms.bed.gz")).merge(d=-1)
    hot_enhs = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.noproms.bed.gz")).merge(d=-1)
    hots = hot_proms.cat(hot_enhs).merge(d=-1)

    f = join(DATA_PATH, "src_files/HepG2_enhancers_DHS_H3K27ac.bed.gz")
    re = BedTool(f).intersect(hots, v=True)

    f = join(DATA_PATH, "src_files/HepG2_superenhancers.bed.gz")
    se = BedTool(f).sort().merge()

    hk_genes = set([l.split(";")[0] for l in open(HK_GENES_FILE).readlines()[1:]])

    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True, wa=True)

    names = ["HOT promoters", "HOT enhancers", "Regular promoters",  "Regular enhancers", "Super-enhancers"]
    regions_list = [hot_proms, hot_enhs, rp, re, se]

    nearest_genes_list = [hot_proms.intersect(str(PROMOTERS_FILE), wo=True).groupby(c=7, o="distinct"),
                          extract_nearest_genes(hot_enhs),
                          rp,
                          extract_nearest_genes(re),
                          extract_nearest_genes(se)]

    data = []

    for name, regions, nearest_genes in zip(names, regions_list, nearest_genes_list):

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


def plot_fractions_merged(save_file):

    data = extract_HK_fractions()

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
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def load_enst2ensg():

    enst2ensg = {}

    f = DATA_PATH/"src_files/hg19_files/ensGene.txt.gz"
    for l in gzip.open(f, "rt"):
        parts = l.split("\t")
        enst2ensg[parts[1]] = parts[12]

    return enst2ensg


def tau_data():

    ensg2tau = {}

    for l in gzip.open(TAU_FILE, "rt"):

        if l.startswith("gene_id"):
            continue

        parts = l.split(",")
        expr = parts[1]
        if expr == "NA":
            continue
        ensg2tau[parts[0]] = float(expr)

    hots = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.bed.gz")).merge(d=-1).sort()
    hot_enhs = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.noproms.bed.gz")).merge(d=-1).sort()

    re = BedTool(join(DATA_PATH, "src_files/HepG2_enhancers_DHS_H3K27ac.bed.gz")).intersect(hots, v=True, wa=True)
    se = BedTool(join(DATA_PATH, "src_files/HepG2_superenhancers.bed.gz")).sort().merge()
    rp = BedTool(PROMOTERS_FILE).intersect(hots, v=True, wa=True)

    names = ["HOT promoters", "HOT enhancers", "Regular promoters",  "Regular enhancers", "Super-enhancers"]

    nearest_genes_list = [hots.intersect(str(PROMOTERS_FILE), wo=True).groupby(c=7, o="distinct"),
                          extract_nearest_genes(hot_enhs),
                          rp,
                          extract_nearest_genes(re),
                          extract_nearest_genes(se)]

    data = []

    enst2ensg = load_enst2ensg()

    for name, nearest_genes in zip(names, nearest_genes_list):
        nearest_tss_pool = set([_.split(".")[0] for r in nearest_genes for _ in r.fields[3].split(",")])

        for enst in nearest_tss_pool:

            if enst not in enst2ensg:
                continue

            ensg = enst2ensg[enst]
            if ensg in ensg2tau:
                data.append([name, ensg2tau[ensg]])

    df = pandas.DataFrame(data=data, columns=["category", "tau"])
    return df


def plot_tau_single_cl(save_file):

    df = tau_data()

    _order = ['HOT promoters', 'Regular promoters', 'HOT enhancers',  'Regular enhancers', 'Super-enhancers']

    plt.figure(figsize=(3.8, 2.7))
    ax = plt.gca()
    g = sns.violinplot(data=df, x="category", order=_order, y="tau", hue="category", dodge=False, legend=False)
    g.set_ylim([0, 1.25])
    g.legend_.remove()
    g.set_xlabel("")

    ax.set_xticklabels(["HOT prom", "reg prom", "HOT enh", "reg enh", "super-enh"], rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.5)
    ax.set_ylim([0, 1.4])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    plt.tight_layout()

    plt.savefig(save_file)



if __name__ == "__main__":

    pass






