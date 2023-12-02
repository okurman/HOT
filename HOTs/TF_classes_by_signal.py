import gzip
import warnings

import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, kruskal, spearmanr, sem

warnings.filterwarnings('ignore')

import os
import sys
sys.path.append(os.environ["HOT_CODE"])
from os.path import join
import numpy as np
from pybedtools import BedTool
import pandas
import matplotlib.pyplot as plt
from data_prep.basic import load_metadata

import seaborn as sns
from pathlib import Path

import os
DATA_PATH = Path(os.environ["HOT_DATA"])

get_loci_files = lambda x: [join(DATA_PATH, "log_bins/%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
TF_CLASS_FILE = DATA_PATH/"src_files/human_transcription_factors.txt"


def extract_TF_signals_binwise(from_file=True):

    save_file = join(DATA_PATH, "src_files/TF_binwise_signal_values.txt")

    if from_file:
        M = np.loadtxt(save_file)
        return M

    a = load_metadata()["HepG2"]

    tfs = sorted(a.keys())
    tf_ids = [a[_] for _ in tfs]
    tf_id2ind = {_: ind for ind, _ in enumerate(tf_ids)}

    M = np.zeros((len(tfs), 14))
    M_count = np.zeros((len(tfs), 14))

    loci = get_loci_files("HepG2")

    for bin_no, bin_loci in enumerate(loci):
        print(bin_no)
        locus = BedTool(bin_loci)
        for r in locus:
            r_tf_ids = r.fields[8].split(",")
            r_tf_signal_values = np.asarray(r.fields[9].split(","), dtype=float)
            for _tf_id, _tf_signal in zip(r_tf_ids, r_tf_signal_values):
                tf_ind = tf_id2ind[_tf_id]
                M[tf_ind, bin_no] += _tf_signal
                M_count[tf_ind, bin_no] += 1

    M /= M_count

    print(save_file)
    np.savetxt(save_file, M)


def extract_TF_signals_binwise_df(normalize=True):

    M = extract_TF_signals_binwise()

    a = load_metadata()
    tfs = sorted(a["HepG2"].keys())

    if normalize:
        M /= M.max(axis=1, keepdims=True)

    data = []
    for tf_ind, tf in enumerate(tfs):
        for bin in range(14):
            data.append([tf, bin, M[tf_ind, bin]])
    df = pd.DataFrame(data=data, columns=["tf", "bin", "value"])

    return df


def get_tf_signal_HOT_ratio(M=None):

    if M is None:
        M = extract_TF_signals_binwise()

    a = load_metadata()
    tfs = np.asarray(sorted(a["HepG2"].keys()))

    hot_M = M[:, -4:].mean(axis=1)
    non_hot_M = M[:, :-4].mean(axis=1)
    ratio = hot_M / non_hot_M

    ix = np.where(~np.isnan(ratio))

    ratio = np.log2(ratio[ix])
    tfs = tfs[ix]

    ix = np.argsort(ratio)[::-1]

    tf2ratio = {tfs[i]: ratio[i] for i in ix}

    return tf2ratio


def extract_TF_signals_per_tf(from_file=True):

    if not from_file:

        a = load_metadata()

        id2tf = {v: k for (k, v) in a["HepG2"].items()}

        hot_files = [join(DATA_PATH, "HOTs/HepG2_HOTs.proms.bed.gz"),
                     join(DATA_PATH, "HOTs/HepG2_HOTs.noproms.bed.gz")]

        data = []
        for hot_file, locus_type in zip(hot_files, ["HOT promoter", "HOT enhancer"]):
            for l in open(hot_file):
                parts = l.rstrip().split("\t")
                _ids = parts[8].split(",")
                _signal_values = parts[9].split(",")

                for _id, _v in zip(_ids, _signal_values):
                    data.append([locus_type, id2tf[_id], _v])

        df = pandas.DataFrame(data=data,
                              columns=["locus_type", "tf", "value"])


    else:

        save_file = join(DATA_PATH, "src_files/chipseq_signal_values_by_tfs.csv.gz")
        df = pandas.read_csv(save_file)

    return df


def plot_tf_signals_HOTs_4panels(save_file, df=None, all=False):

    if df is None:

        df = extract_TF_signals_per_tf(from_file=True)
        df = df[["tf", "value"]].\
            groupby(by=["tf"], as_index=False).\
            agg(mean=("value", np.mean), se=("value", sem)).\
            reset_index(drop=True).sort_values("mean", ascending=False).reset_index(drop=True)

    tf_stats = pandas.read_csv(DATA_PATH/"src_files/tf_stats_summary_table.txt", sep="\t")

    tf2hot_ratio = get_tf_signal_HOT_ratio()

    ratio_tfs = list(tf2hot_ratio.keys())

    df = df[df.tf.isin(ratio_tfs)].reset_index(drop=True)

    ratios = [tf2hot_ratio[tf] if tf in tf2hot_ratio else None for tf in df.tf]

    df["signal_ratio"] = ratios
    df = df.sort_values("signal_ratio", ascending=False).reset_index(drop=True)

    if not all:
        top = df.head(20)
        top["status"] = "HOT-specific"
        bottom = df.tail(20)
        bottom["status"] = "non-HOT-specific"
        df = pd.concat((top, bottom)).reset_index(drop=True)

    tf2hots_perc = {}
    for row in tf_stats.iterrows():
        tf = row[1]["tf"]
        hot_frac = row[1]["all_fracs"]
        enh_frac = row[1]["enh_fracs"]
        prom_frac = row[1]["prom_fracs"]
        tf2hots_perc[tf] = [hot_frac, enh_frac, prom_frac]

    perc_array = [tf2hots_perc[tf][0] for tf in df["tf"].values]
    df["perc"] = perc_array

    if all:
        fig, axs = plt.subplots(4, 1, figsize=(20, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1, 1, 2]})
    else:
        fig, axs = plt.subplots(4, 1, figsize=(5, 6), sharex=True, gridspec_kw={'height_ratios': [1, 1, 1, 1.5]})

    g = sns.barplot(x="tf", y="mean", data=df, ax=axs[0], hue="status", dodge=False, palette="Set1", alpha=0.5)
    g.legend(ncol=2, fontsize=7, frameon=False)

    axs[0].errorbar(x=np.arange(df.shape[0]), y=df["mean"], yerr=df["se"], ls="none")
    axs[0].set_xlabel("")
    axs[0].set_ylabel("ChIP-seq \n signal value")
    plt.setp(axs[0].get_xticklabels(), visible=False)

    axs[1].set_xlabel("")
    g = sns.barplot(x="tf", y="perc", data=df, ax=axs[1], hue="status", dodge=False, palette="Set1", alpha=0.5)
    g.legend_.remove()
    axs[1].set_ylim([0, 100])
    axs[1].grid(axis="y")
    axs[1].set_ylabel("HOT loci (%)")

    axs[2].set_xlabel("")
    g = sns.barplot(x="tf", y="signal_ratio", data=df, ax=axs[2], hue="status", dodge=False, palette="Set1", alpha=0.5)
    g.legend_.remove()
    axs[2].grid(axis="y")
    axs[2].set_ylabel("FC\n(HOT/non-HOT)")

    selected_tfs = df.tf.to_list()

    df_sig_binwise = extract_TF_signals_binwise_df(normalize=True)
    pivot = df_sig_binwise.pivot_table(index="bin", columns="tf", values="value")
    pivot = pivot.reindex(selected_tfs, axis=1)
    pivot = pivot.reindex(np.arange(14)[::-1], axis=0)

    im = axs[3].imshow(pivot.to_numpy(), cmap="Greys", aspect="auto")
    axs[3].set_xticklabels([_[:10] for _ in df["tf"].values], rotation=90, size=3 if all else 8)
    axs[3].set_yticks(np.arange(len(HEPG2_XTICK_LABELS))-0.5, HEPG2_XTICK_LABELS[::-1], size=7)
    axs[3].set_ylabel("# of ChIP-seq peaks")

    cbar_ax = fig.add_axes([0.99, 0.15, 0.01, 0.2])
    cbar = fig.colorbar(im, cax=cbar_ax)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(5)
    cbar.set_label("Norm. ChIP-seq signal", fontsize=5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.05)

    fig.savefig(save_file)
    plt.close()


def extract_TF_signals_particular():

    a = load_metadata()
    data = []

    ctcf_id = a["HepG2"]["CTCF"]
    med_id = a["HepG2"]["MED1"]
    p300_id = a["HepG2"]["EP300"]
    pol2_id = a["HepG2"]["POLR2A"]

    hot_file = join(DATA_PATH, "HOTs/HepG2_HOTs.proms.bed.gz")
    ctcf_values, med_values, p300_values, pol2_values = [], [], [], []
    for l in gzip.open(hot_file, "rt"):
        parts = l.rstrip().split("\t")
        _ids = parts[8].split(",")
        _signal_values = parts[9].split(",")

        for _id, _v in zip(_ids, _signal_values):
            if _id == ctcf_id:
                ctcf_values.append(float(_v))
            elif _id == med_id:
                med_values.append(float(_v))
            elif _id == p300_id:
                p300_values.append(float(_v))
            elif _id == pol2_id:
                pol2_values.append(float(_v))

    data.append([ctcf_values, med_values, p300_values, pol2_values])

    hot_file = join(DATA_PATH, "HOTs/HepG2_HOTs.noproms.bed.gz")
    ctcf_values, med_values, p300_values, pol2_values = [], [], [], []
    for l in gzip.open(hot_file, "rt"):
        parts = l.rstrip().split("\t")
        _ids = parts[8].split(",")
        _signal_values = parts[9].split(",")

        for _id, _v in zip(_ids, _signal_values):
            if _id == ctcf_id:
                ctcf_values.append(float(_v))
            elif _id == med_id:
                med_values.append(float(_v))
            elif _id == p300_id:
                p300_values.append(float(_v))
            elif _id == pol2_id:
                pol2_values.append(float(_v))

    data.append([ctcf_values, med_values, p300_values, pol2_values])

    return data


def extract_TF_signals():

    a = load_metadata()

    data = []

    ss_tfs = []
    nss_tfs = []

    for l in open(TF_CLASS_FILE):
        if not l.strip():
            continue

        parts = l.strip().split("\t")

        if parts[1] == "Yes":
            ss_tfs.append(parts[0])
        else:
            nss_tfs.append(parts[0])

    ss_tf_ids = []
    nss_tf_ids = []

    for tf, tf_id in a["HepG2"].items():
        if tf in ss_tfs:
            ss_tf_ids.append(tf_id)
        elif tf in nss_tfs:
            nss_tf_ids.append(tf_id)

    hot_file = DATA_PATH/"HOTs/HepG2_HOTs.proms.bed.gz"
    ss_values = []
    nss_values = []
    for l in gzip.open(hot_file, "rt"):
        parts = l.rstrip().split("\t")
        _ids = parts[8].split(",")
        _signal_values = parts[9].split(",")

        for _id, _v in zip(_ids, _signal_values):
            if _id in ss_tf_ids:
                ss_values.append(float(_v))
            elif _id in nss_tf_ids:
                nss_values.append(float(_v))

    ss_values, nss_values = np.asarray(ss_values), np.asarray(nss_values)
    data.append([ss_values, nss_values])

    hot_file = DATA_PATH / "HOTs/HepG2_HOTs.noproms.bed.gz"
    ss_values = []
    nss_values = []
    for l in gzip.open(hot_file, "rt"):
        parts = l.rstrip().split("\t")
        _ids = parts[8].split(",")
        _signal_values = parts[9].split(",")

        for _id, _v in zip(_ids, _signal_values):
            if _id in ss_tf_ids:
                ss_values.append(float(_v))
            elif _id in nss_tf_ids:
                nss_values.append(float(_v))

    ss_values, nss_values = np.asarray(ss_values), np.asarray(nss_values)
    data.append([ss_values, nss_values])

    return data


def merge_values():

    overall_data = extract_TF_signals()
    tfs_data = extract_TF_signals_particular()

    data = []
    for _ in overall_data[0][0]:
        data.append(["HOT promoters", "ssTFs", _])
    for _ in overall_data[0][1]:
        data.append(["HOT promoters", "nssTFs", _])
    for _ in overall_data[1][0]:
        data.append(["HOT enhancers", "ssTFs", _])
    for _ in overall_data[1][1]:
        data.append(["HOT enhancers", "nssTFs", _])

    for _ in tfs_data[0][0]:
        data.append(["HOT promoters", "CTCF", _])
    for _ in tfs_data[0][1]:
        data.append(["HOT promoters", "Mediator", _])
    for _ in tfs_data[0][2]:
        data.append(["HOT promoters", "P300", _])
    for _ in tfs_data[0][3]:
        data.append(["HOT promoters", "POLR2", _])

    for _ in tfs_data[1][0]:
        data.append(["HOT enhancers", "CTCF", _])
    for _ in tfs_data[1][1]:
        data.append(["HOT enhancers", "Mediator", _])
    for _ in tfs_data[1][2]:
        data.append(["HOT enhancers", "P300", _])
    for _ in tfs_data[1][3]:
        data.append(["HOT enhancers", "POLR2", _])

    df = pandas.DataFrame(data=data,
                          columns=["locus_type", "tf_category", "value"])

    return df


def plot_signals(save_file):

    df_small = merge_values()

    plt.figure(figsize=(5, 3))
    ax = plt.gca()

    hue_order = ["ssTFs", "nssTFs", "POLR2", "CTCF", "P300"]

    g = sns.boxplot(x="locus_type", y="value", hue="tf_category", hue_order=hue_order, showfliers=False, data=df_small,
                    ax=ax)
    g.legend(loc="upper left", ncol=2, prop={'size': 9})
    ax.set_xlabel("")
    ax.set_ylabel("Signal value")
    ax.grid(which="major", axis="y", alpha=0.2)
    plt.tight_layout()

    plt.savefig(save_file)

