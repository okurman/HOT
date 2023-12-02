
import warnings
warnings.filterwarnings('ignore')
import seaborn as sns
import pandas
from os.path import join
from pybedtools import BedTool
import matplotlib.pyplot as plt
from pathlib import Path

import os
DATA_PATH = Path(os.environ["HOT_DATA"])
LOCI_DIR = DATA_PATH / "HOTs"


def get_overlaps_subsample_h1():

    h1_hots = BedTool(join(LOCI_DIR, "H1_HOTs.bed.gz")).merge(d=-1)
    h1_count = h1_hots.count()

    ovp_list = []

    _d = DATA_PATH / "src_files/HOTs_with_subsampled_to_H1"
    for i in range(10):
        hepg2_hots = BedTool(str(_d / "HepG2" / ("%d_HOTs.bed.gz" % i)))
        k562_hots = BedTool(str(_d / "K562" / ("%d_HOTs.bed.gz" % i)))
        cl_hots = hepg2_hots.cat(k562_hots, postmerge=False).sort().merge()
        ovlp_fr = h1_hots.intersect(cl_hots, wa=True).merge(d=-1).count() / h1_count
        ovp_list.append(ovlp_fr)

    return ovp_list


def get_overlaps_H1():

    loci_dir = DATA_PATH / "src_files/HOTs_with_common_tfs/"
    hepg2 = BedTool(str(join(loci_dir, "HepG2_HOTs.bed.gz")))
    k562 = BedTool(str(join(loci_dir, "K562_HOTs.bed.gz")))
    h1 = BedTool(str(join(loci_dir, "H1_HOTs.bed.gz")))

    common_TFs_ovp = h1.intersect(hepg2.cat(k562, postmerge=False), wa=True).merge(d=-1).count()/h1.count()

    hepg2 = BedTool(str(join(LOCI_DIR, "HepG2_HOTs.bed.gz")))
    k562 = BedTool(str(join(LOCI_DIR, "K562_HOTs.bed.gz")))
    h1 = BedTool(str(join(LOCI_DIR, "H1_HOTs.bed.gz")))

    all_TFs_ovp = h1.intersect(hepg2.cat(k562, postmerge=False), wa=True).merge(d=-1).count() / h1.count()

    return common_TFs_ovp, all_TFs_ovp


def load_data():

    data = []

    common_TFs_ovp, all_TFs_ovp = get_overlaps_H1()

    data += [
        ["complete", all_TFs_ovp],
        ["common DAPs", common_TFs_ovp]
    ]

    ovp_values = get_overlaps_subsample_h1()
    data += [["subsampled to H1", _] for _ in ovp_values]

    df = pandas.DataFrame(data=data, columns=["sample", "value"])
    df["value"] *= 100
    return df


def plot_subsample_h1(save_file):

    df = load_data()

    plt.figure(figsize=(3.5, 3))
    ax = plt.gca()

    sns.barplot(x="sample", y="value", data=df, ci="sd", ax=ax, order=['complete', 'common DAPs', 'subsampled to H1'])
    ax.set_ylim([0, 100])
    ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax.set_xticklabels(['all\nDAPs', 'common\nDAPs', 'subsampled\nto H1'])
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.4)
    plt.tight_layout()

    plt.savefig(save_file)
    plt.close()




