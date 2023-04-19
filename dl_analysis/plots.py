
from os.path import join

import numpy as np
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

WORK_DIR = ROOT_DIR
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']


def load_data(fname="auPRC.txt"):

    lines = []

    aucs = np.loadtxt(join(WORK_DIR, f"dl_runs/HepG2_400_14/{fname}"))
    lines += [["HepG2", 400, i, auc] for i, auc in enumerate(aucs)]
    aucs = np.loadtxt(join(WORK_DIR, f"dl_runs/K562_400_14/{fname}"))
    lines += [["K562", 400, i, auc] for i, auc in enumerate(aucs)]

    aucs = np.loadtxt(join(WORK_DIR, f"dl_runs/HepG2_1000_14/{fname}"))
    lines += [["HepG2", 1000, i, auc] for i, auc in enumerate(aucs)]
    aucs = np.loadtxt(join(WORK_DIR, f"dl_runs/K562_1000_14/{fname}"))
    lines += [["K562", 1000, i, auc] for i, auc in enumerate(aucs)]

    df = pandas.DataFrame(columns=["cl", "length", "category", "auc"],
                          data=lines)

    return df


def get_correlation():

    df = load_data(33)

    aucs = df[(df["cl"] == "HepG2") & (df["length"] == 400)]["auc"].values
    tfs = np.asarray(HEPG2_XTICK_LABELS, dtype=int)
    tfs = tfs[:-1] + (tfs[1:] - tfs[:1])/2

    print(pearsonr(aucs, tfs))
    print(spearmanr(aucs, tfs))


def plot_aucs():

    # df = load_data()
    df = load_data(14)

    # plt.figure()
    fig, axs = plt.subplots(1, 2, figsize=(7, 4))

    for ax, cl in zip(axs, ["HepG2", "K562"]):

        categories = df["category"].max()+1

        _df = df[df["cl"] == cl]
        g = sns.barplot(x="category", y="auc", hue="length", hue_order=[400, 1000], data=_df, palette="Set1", ax=ax)

        ax.grid(which="major", axis="y", alpha=0.5)
        g.get_legend().remove()
        ax.set_ylim([0, 1.05])

        ax.set_xticklabels(["COLD", "MILD", "HOT"])
        ax.set_xlabel("")
        ax.set_title(cl)

        for i in range(categories):
            _auc = _df[(_df["cl"] == cl) & (_df["category"] == i) & (_df["length"] == 400)]["auc"].values[0]
            _auc_str = str(round(_auc, 2))[1:]
            ax.text(i - 0.19, _auc + 0.03, _auc_str, color='black', ha="center", size=7)
            _auc = _df[(_df["cl"] == cl) & (_df["category"] == i) & (_df["length"] == 1000)]["auc"].values[0]
            _auc_str = str(round(_auc, 2))[1:]
            ax.text(i + 0.19, _auc + 0.03, _auc_str, color='black', ha="center", size=7)

    axs[0].set_ylabel("auROC")
    axs[1].set_ylabel("")

    handles, labels = axs[0].get_legend_handles_labels()
    plt.subplots_adjust(right=0.84)
    fig.legend(handles, labels, ncol=1, bbox_to_anchor=(1, 0.5), loc="center right", title="Input length")

    # save_file = join(WORK_DIR, "plots/aucROCs_3class.pdf")
    # plt.savefig(save_file)
    # print(save_file)
    # plt.tight_layout()
    plt.show()


def plot_aucs_14():

    df = load_data("auc.txt")
    df = df[df["category"] < 13]

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))

    for ax, cl in zip(axs, ["HepG2", "K562"]):

        _df = df[df["cl"] == cl]
        sns.barplot(x="category", y="auc", hue="length", hue_order=[400, 1000], data=_df, palette="Set1", ax=ax)

        ax.grid(which="both", axis="y", alpha=0.5)
        ax.set_ylim([0, 1.05])

        x_range = np.arange(len(PERC_XTICK_LABELS)-1)
        x_ticks = x_range - 0.5

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(HEPG2_XTICK_LABELS[:-1] if cl == "HepG2" else K562_XTICK_LABELS[:-1])
        ax.set_xlabel("DAP ChIP-seq peaks")

        ticks = ax.get_xticklabels()
        for t in ticks[5:]:
            t.set_rotation(30)
            t.set_fontsize(9)

        # ax.set_xlabel("")
        ax.set_title(cl)

    axs[0].legend(ncol=2)
    axs[1].legend(ncol=2)

    axs[0].set_ylabel("auROC")
    axs[1].set_ylabel("auROC")

    plt.tight_layout()
    save_file = join(WORK_DIR, "plots/aucROCs_14class.pdf")
    print(save_file)
    plt.savefig(save_file)

    plt.show()


def plot_auPRCs_14():

    df = load_data()
    df = df[df["category"]<13]

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))

    for ax, cl in zip(axs, ["HepG2", "K562"]):

        _df = df[df["cl"] == cl]
        sns.barplot(x="category", y="auc", hue="length", hue_order=[400, 1000], data=_df, palette="Set1", ax=ax)

        ax.grid(which="both", axis="y", alpha=0.5)
        ax.set_ylim([0, 1.05])

        x_range = np.arange(len(PERC_XTICK_LABELS)-1)
        x_ticks = x_range - 0.5

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(HEPG2_XTICK_LABELS[:-1] if cl == "HepG2" else K562_XTICK_LABELS[:-1])
        ax.set_xlabel("DAP ChIP-seq peaks")

        ticks = ax.get_xticklabels()
        for t in ticks[5:]:
            t.set_rotation(30)
            t.set_fontsize(9)

        # ax.set_xlabel("")
        ax.set_title(cl)

    axs[0].legend(ncol=2)
    axs[1].legend(ncol=2)

    axs[0].set_ylabel("auPRC")
    axs[1].set_ylabel("auPRC")

    plt.tight_layout()
    save_file = join(WORK_DIR, "plots/auPRCs_14class.pdf")
    print(save_file)
    plt.savefig(save_file)

    plt.show()


def plot_aurocs_14():

    df = load_data(14)

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))

    for ax, cl in zip(axs, ["HepG2", "K562"]):

        _df = df[df["cl"] == cl]
        sns.barplot(x="category", y="auc", hue="length", hue_order=[400, 1000], data=_df, palette="Set1", ax=ax)

        ax.grid(which="both", axis="y", alpha=0.5)
        ax.set_ylim([0, 1.05])

        x_range = np.arange(len(PERC_XTICK_LABELS))
        x_ticks = x_range - 0.5

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(PERC_XTICK_LABELS)
        ax.set_xlabel("TF ChIP-seq peaks")

        ticks = ax.get_xticklabels()
        for t in ticks[5:]:
            t.set_rotation(30)
            t.set_fontsize(9)

        ax.set_xlabel("")
        ax.set_title(cl)

    axs[0].legend(ncol=2)
    axs[1].legend(ncol=2)

    axs[0].set_ylabel("auROC")
    axs[1].set_ylabel("auROC")

    plt.tight_layout()
    save_file = join(WORK_DIR, "plots/aucROCs_14class.pdf")
    print(save_file)
    plt.savefig(save_file)

    plt.show()


def plot_aucs_14_cl():

    df = load_data(14)

    # fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))

    plt.close()

    for cl in ["HepG2", "K562"]:

        plt.figure(figsize=(5, 3.5))
        ax = plt.gca()

        _df = df[df["cl"] == cl]
        sns.barplot(x="category", y="auc", hue="length", hue_order=[400, 1000], data=_df, palette="Set1", ax=ax)

        ax.grid(which="both", axis="y", alpha=0.5)
        ax.set_ylim([0, 1.05])

        x_range = np.arange(len(PERC_XTICK_LABELS))
        x_ticks = x_range - 0.5

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(HEPG2_XTICK_LABELS if cl == "HepG2" else K562_XTICK_LABELS)
        ax.set_xlabel("TF ChIP-seq peaks")

        ticks = ax.get_xticklabels()
        for t in ticks[5:]:
            t.set_rotation(30)
            t.set_fontsize(9)

        ax.set_xlabel("DAPs")
        # ax.set_title(cl)

        ax.legend(ncol=2, title="length (bp)")
        ax.set_ylabel("auROC")

        plt.tight_layout()
        save_file = join(WORK_DIR, "plots/aucROCs_14class_%s.pdf" % cl)
        print(save_file)
        plt.savefig(save_file, bbox_inch="tight")
        plt.close()


def plot_aucs_14_400():

    df = load_data(classes=14)

    fig = plt.figure(figsize=(6, 4))
    ax = plt.gca()
    sns.barplot(x="category", y="auc", hue="cl", data=df, palette="Set1", ax=ax)

    ax.grid(which="major", axis="y", alpha=0.5)
    ax.set_ylim([0, 1.05])

    ax.set_ylabel("auROC")
    ax.set_xlabel("")
    ax.legend(title="")
    ax.set_title("Input length=400bps")

    x_range = np.arange(len(PERC_XTICK_LABELS))
    x_ticks = x_range - 0.5

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(PERC_XTICK_LABELS)
    ax.set_xlabel("TF ChIP-seq peaks")

    ticks = ax.get_xticklabels()
    for t in ticks[5:]:
        t.set_rotation(30)
        t.set_fontsize(9)

    save_file = join(WORK_DIR, "plots/aucROCs_14class.pdf")
    print(save_file)
    plt.savefig(save_file)

    plt.show()

