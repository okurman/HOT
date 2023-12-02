
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import seaborn as sns
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

from pathlib import Path
import os
DATA_PATH = Path(os.environ["HOT_DATA"])


def plot_common_variant_densities(save_file):

    data_file = DATA_PATH/"src_files/variant_density_values.txt"

    df_density = pd.read_csv(data_file, sep="\t")

    fig, axs = plt.subplots(3, 2, figsize=(6, 7))

    x_labels = ["HOT prom", "reg. prom", "HOT enh", "reg. enh"]
    order = ["HOT_prom", "RP", "HOT_enh", "RE"]

    _df = df_density[df_density["variant"] == "INDELs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 0])
    axs[0, 0].set_xlabel("")
    axs[0, 0].set_ylabel("var/10kb", fontsize=9)
    axs[0, 0].set_title("INDELs")
    axs[0, 0].set_xticks([])
    axs[0, 0].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[0, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "SNPs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 1])
    axs[0, 1].set_xlabel("")
    axs[0, 1].set_ylabel("var/10kb", fontsize=9)
    axs[0, 1].set_title("SNPs")
    axs[0, 1].set_xticks([])
    axs[0, 1].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[0, 1].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "eQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[1, 0])
    axs[1, 0].set_xlabel("")
    axs[1, 0].set_ylabel("var/10kb", fontsize=9)
    axs[1, 0].set_title("eQTLs")
    axs[1, 0].set_xticks([])
    axs[1, 0].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[1, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "caQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[1, 1])
    axs[1, 1].set_xlabel("")
    axs[1, 1].set_ylabel("var/10kb")
    axs[1, 1].set_title("caQTLs")
    axs[1, 1].set_xticks([])
    axs[1, 1].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[1, 1].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "raQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[2, 0])
    axs[2, 0].set_xlabel("")
    axs[2, 0].set_ylabel("var/10kb")
    axs[2, 0].set_title("raQTLs")
    axs[2, 0].grid(axis="y", alpha=0.3)
    axs[2, 0].set_xticklabels(x_labels, rotation=45, ha="right")
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[2, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "GWAS"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[2, 1])
    axs[2, 1].set_xlabel("")
    axs[2, 1].set_ylabel("var/10kb")
    axs[2, 1].set_title("GWAS")
    axs[2, 1].grid(axis="y", alpha=0.3)
    axs[2, 1].set_xticklabels(x_labels, rotation=45, ha="right")
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[2, 1].set_ylim([0, y_upper])

    plt.tight_layout()
    plt.savefig(save_file)
    plt.close()


def variant_stats_for_text(df_density, df_counts):

    _vars = ["INDELs", 'SNPs', 'eQTLs', 'caQTLs', "raQTLs", "GWAS", "ClinVar"]

    for _var in _vars:
        print(_var)

        _df = df_density[df_density["variant"] == _var]
        regions = ['HOT_enh', 'RE', 'HOT_prom', 'RP']
        counts = [_df[_df["region"] == _]["value"].values[0] for _ in regions]
        print("\t".join(regions))
        print("\t".join(["%.1f" % (100000 * _) for _ in counts]))
        folds = [counts[0] / counts[1], counts[2] / counts[3]]
        print("fold enhs \t fold proms")
        print("\t".join([str(_) for _ in folds]))

        _df = df_counts[df_counts["variant"] == _var]
        _df = df_counts[df_counts["variant"] == _var]
        v1 = _df[_df["region"] == "HOT_prom"].value
        v2 = _df[_df["region"] == "RP"].value
        p_prom = mannwhitneyu(v1, v2).pvalue
        v1 = _df[_df["region"] == "HOT_enh"].value
        v2 = _df[_df["region"] == "RE"].value
        p_enh = mannwhitneyu(v1, v2).pvalue
        print(p_enh, p_prom)

        print()
