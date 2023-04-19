import warnings

import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, kruskal, spearmanr, sem

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

from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

TF_CLASS_FILE = "ROOT_DIR/chipseq_files/human_transcription_factors.txt"

BINS_DIR = "ROOT_DIR/definitions/peak_8bp_v2/log_bins/"
get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed" % (x, i)) for i in range(14)]
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']


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

    hot_file = join(LOCI_DIR, "HepG2_HOTs.proms.bed")
    ss_values = []
    nss_values = []
    for l in open(hot_file):
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

    hot_file = join(LOCI_DIR, "HepG2_HOTs.noproms.bed")
    ss_values = []
    nss_values = []
    for l in open(hot_file):
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


def extract_TF_signals_binwise(from_file=True):

    save_file = join(LOCI_DIR, "TF_binwise_signal_values.txt")

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

        # for i in range(M.shape[0]):
        #     _arr = M[i, :]/np.max()
        #     _arr = (_arr - np.mean(_arr))/np.var(_arr)
        #     M[i, :] = _arr

    data = []
    for tf_ind, tf in enumerate(tfs):
        for bin in range(14):
            data.append([tf, bin, M[tf_ind, bin]])
    df = pd.DataFrame(data=data, columns=["tf", "bin", "value"])

    return df


def extract_TF_signals_pairwise(from_file=True, normalize=True):

    save_file = join(LOCI_DIR, "TF_pairwise_signal_values.txt")

    if from_file:
        M = np.loadtxt(save_file)
        a = load_metadata()
        tfs = sorted(a["HepG2"].keys())
        if normalize:
            M /= np.nanmax(M, axis=1, keepdims=True)
        M[np.isnan(M)] = 0
        M[np.diag_indices_from(M)] = 0
        data = []
        for tf_ind_1, tf_1 in enumerate(tfs):
            for tf_ind_2, tf_2 in enumerate(tfs):
                data.append([tf_1, tf_2, M[tf_ind_1, tf_ind_2]])
        df = pd.DataFrame(data=data, columns=["tf_1", "tf_2", "value"])

        return M, df

    a = load_metadata()["HepG2"]

    tfs = sorted(a.keys())
    tf_ids = [a[_] for _ in tfs]
    tf_id2ind = {_: ind for ind, _ in enumerate(tf_ids)}

    M = np.zeros((len(tfs), len(tfs)))
    M_count = np.zeros(M.shape)

    loci = get_loci_files("HepG2")

    for bin_no, bin_loci in enumerate(loci):
        print(bin_no)

        locus = BedTool(bin_loci)
        for r in locus:
            r_tf_ids = r.fields[8].split(",")
            r_tf_signal_values = np.asarray(r.fields[9].split(","), dtype=float)

            for i in range(len(r_tf_ids) - 1):
                tf_ind_i = tf_id2ind[r_tf_ids[i]]
                signal_i = r_tf_signal_values[i]
                for j in range(i+1, len(r_tf_ids)):
                    tf_ind_j = tf_id2ind[r_tf_ids[j]]

                    M[tf_ind_i, tf_ind_j] += signal_i
                    M_count[tf_ind_i, tf_ind_j] += 1

                    M[tf_ind_j, tf_ind_i] += signal_i
                    M_count[tf_ind_j, tf_ind_i] += 1

    M = np.divide(M, M_count)

    print(save_file)
    np.savetxt(save_file, M)


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


def extract_H3K27ac_signals():

    a = load_metadata()

    k27_file = "ROOT_DIR/definitions/regular_enhancers/encode_files/HepG2_H3K27ac_ENCFF749VEQ.bed.gz"
    k27ac = BedTool(k27_file)
    hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.PE.bed"))

    proms = BedTool([r for r in hots if r.fields[-1] == "P"])
    enhs = BedTool([r for r in hots if r.fields[-1] == "E"])

    prom_signals = np.asarray([_
                               for r in proms.intersect(k27ac, wo=True).groupby(c=11, o="collapse")
                               for _ in r.fields[-1].split(",")], dtype=float)

    enh_signals = np.asarray([_
                               for r in enhs.intersect(k27ac, wo=True).groupby(c=11, o="collapse")
                               for _ in r.fields[-1].split(",")], dtype=float)

    return prom_signals, enh_signals


def extract_TF_signals_per_tf(from_file=None):

    if not from_file:

        a = load_metadata()

        id2tf = {v: k for (k, v) in a["HepG2"].items()}

        hot_files = [join(LOCI_DIR, "HepG2_HOTs.proms.bed"),
                     join(LOCI_DIR, "HepG2_HOTs.noproms.bed")]

        data = []
        for hot_file, locus_type in zip(hot_files, ["HOT promoter", "HOT enhancer"]):
            print(hot_file)
            for l in open(hot_file):
                parts = l.rstrip().split("\t")
                _ids = parts[8].split(",")
                _signal_values = parts[9].split(",")

                for _id, _v in zip(_ids, _signal_values):
                    data.append([locus_type, id2tf[_id], _v])

        df = pandas.DataFrame(data=data,
                              columns=["locus_type", "tf", "value"])

        save_file = join(LOCI_DIR, "chipseq_signal_values_by_tfs.csv.gz")
        print(save_file)
        df.to_csv(save_file)

    else:

        save_file = join(LOCI_DIR, "chipseq_signal_values_by_tfs.csv.gz")
        df = pandas.read_csv(save_file)

    return df


def plot_TF_signals_binwise_heatmap():

    fig = plt.figure()

    df = extract_TF_signals_binwise_df(standardize=False)
    df = df.pivot("bin", "tf", "value")

    ax = sns.heatmap(df)

    save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_values_binwise.pdf")
    print(save_file)
    fig.savefig(save_file)
    plt.close()


def plot_tf_signals_PE(df):

    # df = extract_TF_signals_per_tf(from_file=True)
    # return df

    df_tf_hot_summary = pandas.read_table(
        "ROOT_DIR/definitions/peak_8bp_v2/HOTs/enrichment_files/tf_stats_summary_table.txt")
    df_tf_hot_summary = df_tf_hot_summary[df_tf_hot_summary["cell_line"] == "HepG2"]
    data_array = df_tf_hot_summary[["tf", "all_fracs", "enh_fracs", "prom_fracs"]].to_numpy()

    tf2hots_perc = {}
    tf2proms_perc = {}
    tf2enhs_perc = {}
    for _ in data_array:
        [tf, hots_perc, enh_fracs, prom_fracs] = _

        if hots_perc < 10:
            continue

        tf2hots_perc[tf] = hots_perc
        tf2proms_perc[tf] = prom_fracs
        tf2enhs_perc[tf] = enh_fracs

    df.value = df.value.astype(float)

    for _locus, frac_dict in zip(["HOT promoter", "HOT enhancer"], [tf2proms_perc, tf2enhs_perc]):
        _df = df[df.locus_type == _locus]
        sorted_tfs = _df.groupby(by=["tf"]).mean().sort_values(by=["value"], ascending=False)
        tfs_to_depict = [tf for tf in sorted_tfs.index.values if tf in tf2hots_perc][:100]

        fig, axs = plt.subplots(2, 1, figsize=(18, 10), sharex=True, gridspec_kw={'height_ratios': [4, 1]})

        perc_array = [frac_dict[tf] for tf in tfs_to_depict]

        sns.boxplot(x="tf", y="value", order=tfs_to_depict, data=_df, showfliers=False, ax=axs[0])
        axs[0].set_title(_locus+"s")
        axs[0].set_xlabel("")
        axs[0].set_ylabel("ChIP-seq signal value")

        axs[1].set_xlabel("")

        axs[1].bar(np.arange(len(tfs_to_depict)), perc_array)
        axs[1].set_ylim([0, 100])
        axs[1].grid(axis="y")
        axs[1].set_ylabel("% of loci")
        axs[1].set_xticklabels(tfs_to_depict, rotation=90, size=8)

        plt.setp(axs[0].get_xticklabels(), visible=False)

        plt.subplots_adjust(hspace=.0)

        save_file = join(PLOTS_DIR, "HepG2_" + _locus.replace(" ", "_") + "_tf_signal_values.pdf")
        print(save_file)
        fig.savefig(save_file)


def plot_tf_signals_HOTs_4panels(df=None, all=False):

    if df is None:

        df = extract_TF_signals_per_tf(from_file=True)
        df = df[["tf", "value"]].\
            groupby(by=["tf"], as_index=False).\
            agg(mean=("value", np.mean), se=("value", sem)).\
            reset_index(drop=True).sort_values("mean", ascending=False).reset_index(drop=True)

        return df

    tf_stats = pandas.read_table(
        "ROOT_DIR/definitions/peak_8bp_v2/HOTs/enrichment_files/tf_stats_summary_table.txt")

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
    if not all:
        save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_values_ratio_sorted_4panels_top_bottom_20.pdf")
    else:
        save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_values_ratio_sorted_4panels_all.pdf")

    a = BedTool("enhancer_file.hg38.bed").liftover(hg38_to_hg19)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.05)

    print(save_file)
    fig.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_tf_signals_binwise_all():

    # if df is None:
    #     df = extract_TF_signals_per_tf(from_file=True)
    #     df = df[["tf", "value"]]. \
    #         groupby(by=["tf"], as_index=False). \
    #         agg(mean=("value", np.mean), se=("value", sem)). \
    #         reset_index(drop=True).sort_values("mean", ascending=False).reset_index(drop=True)
    #
    #     return df

    df_sig_binwise = extract_TF_signals_binwise_df(standardize=True)

    pivot = df_sig_binwise.pivot_table(index="bin", columns="tf", values="value")
    pivot = pivot.reindex(np.arange(14)[::-1], axis=0)

    g = sns.clustermap(pivot, cmap="vlag", row_cluster=False, metric="seuclidean",
                       xticklabels=True, dendrogram_ratio=0.15)

    tf_labels = [_.get_text() for _ in g.ax_heatmap.axes.get_xticklabels()]
    g.ax_heatmap.axes.set_xticklabels(tf_labels, size=1)

    hm = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height * 0.15])
    g.ax_cbar.set_position((0, 0.7, .03, 0.1))

    g.ax_heatmap.set_yticks(np.arange(len(HEPG2_XTICK_LABELS)) - 0.5, HEPG2_XTICK_LABELS[::-1], size=7)

    save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_binwise_all.pdf")
    print(save_file)
    plt.savefig(save_file)
    plt.close()

    save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_binwise_all_labels.txt")
    with open(save_file, "w") as of:
        [of.write("%s\n" % tf) for tf in tf_labels]


def plot_tf_signals_pairwise_all():

    M, df_pairwise = extract_TF_signals_pairwise(normalize=True)

    pivot = df_pairwise.pivot_table(index="tf_1", columns="tf_2", values="value")

    g = sns.clustermap(pivot, cmap="vlag", metric="sqeuclidean", xticklabels=True,
                       yticklabels=True, row_cluster=False, dendrogram_ratio=0.14)

    g.ax_row_dendrogram.set_visible(False)

    tf_labels = [_.get_text() for _ in g.ax_heatmap.axes.get_yticklabels()]
    g.ax_heatmap.axes.set_yticklabels(tf_labels, size=1)
    tf_labels = [_.get_text() for _ in g.ax_heatmap.axes.get_xticklabels()]
    g.ax_heatmap.axes.set_xticklabels(tf_labels, size=1)

    g.ax_heatmap.axes.set_xlabel("DAP (n=544)")
    g.ax_heatmap.axes.set_ylabel("DAP (n=544)")

    # hm = g.ax_heatmap.get_position()
    # g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height * 0.15])
    # g.ax_cbar.set_position((0, 0.7, .03, 0.1))

    save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_pairwise_all_normalized.pdf")
    print(save_file)
    plt.savefig(save_file)
    plt.close()

    save_file = join(PLOTS_DIR, "HepG2_HOTs_tf_signal_pairwise_all_labels_normalized.txt")
    with open(save_file, "w") as of:
        [of.write("%s\n" % tf) for tf in tf_labels]


def plot_tf_signal_perc_correlation(df):

    # df = extract_TF_signals_per_tf(from_file=True)
    # return df

    df_tf_hot_summary = pandas.read_table(
        "ROOT_DIR/definitions/peak_8bp_v2/HOTs/enrichment_files/tf_stats_summary_table.txt")
    df_tf_hot_summary = df_tf_hot_summary[df_tf_hot_summary["cell_line"] == "HepG2"]
    data_array = df_tf_hot_summary[["tf", "all_fracs", "enh_fracs", "prom_fracs"]].to_numpy()

    tf2hots_perc = {}
    tf2proms_perc = {}
    tf2enhs_perc = {}
    for _ in data_array:
        [tf, hots_perc, enh_fracs, prom_fracs] = _

        if hots_perc < 5:
            continue

        tf2hots_perc[tf] = hots_perc
        tf2proms_perc[tf] = prom_fracs
        tf2enhs_perc[tf] = enh_fracs

    tf2avg_signal = df[["tf", "value"]].groupby(by="tf").mean().to_dict()["value"]

    data = [[hots_perc, tf2avg_signal[tf], tf] for tf, hots_perc in tf2hots_perc.items()]
    data = np.asarray(data)

    plt.figure(figsize=(5, 3.2))
    x = np.asarray(data[:, 0], dtype=np.float)
    y = np.asarray(data[:, 1], dtype=np.float)
    g = sns.regplot(x, y, line_kws={"color": "red"})
    g.set_ylim([0, 200])
    g.set_xlabel("% of HOT loci")
    g.set_ylabel("Signal value")

    sp = spearmanr(x, y)

    rho = sp.correlation
    pval = sp.pvalue

    # # title = "spearman's rho= %.1f -log10(P)= %d" % (rho, -1*np.log10(pval))
    # # plt.title(title)
    plt.tight_layout()

    save_file = join(PLOTS_DIR, "HepG2_HOT_perc_signal_scatterplot.pdf")
    print(save_file)
    plt.savefig(save_file)


def extract_TF_signals_particular():

    a = load_metadata()
    data = []

    ctcf_id = a["HepG2"]["CTCF"]
    med_id = a["HepG2"]["MED1"]
    p300_id = a["HepG2"]["EP300"]
    pol2_id = a["HepG2"]["POLR2A"]

    hot_file = join(LOCI_DIR, "HepG2_HOTs.proms.bed")
    ctcf_values, med_values, p300_values, pol2_values = [], [], [], []
    for l in open(hot_file):
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

    hot_file = join(LOCI_DIR, "HepG2_HOTs.noproms.bed")
    ctcf_values, med_values, p300_values, pol2_values = [], [], [], []
    for l in open(hot_file):
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

    prom_signals, enh_signals = extract_H3K27ac_signals()

    for _ in prom_signals:
        data.append(["HOT promoters", "H3K27ac", _])
    for _ in enh_signals:
        data.append(["HOT enhancers", "H3K27ac", _])

    df = pandas.DataFrame(data=data,
                          columns=["locus_type", "tf_category", "value"])

    return df


def plot_signals(df_small):

    plt.figure(figsize=(5, 3))
    ax = plt.gca()

    hue_order = ["ssTFs", "nssTFs", "POLR2", "CTCF", "P300", "H3K27ac"]

    g = sns.boxplot(x="locus_type", y="value", hue="tf_category", hue_order=hue_order, showfliers=False, data=df_small,
                    ax=ax)
    g.legend(loc="upper left", ncol=2, prop={'size': 9})
    ax.set_xlabel("")
    ax.set_ylabel("Signal value")
    ax.grid(which="major", axis="y", alpha=0.2)
    plt.tight_layout()

    save_file = join(PLOTS_DIR, "HepG2_ssTF_nssTF_extra_signals.pdf")
    plt.savefig(save_file)
    print(save_file)
