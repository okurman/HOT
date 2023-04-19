import subprocess
import warnings
from collections import defaultdict

import pandas as pd

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

# matplotlib.style.use('ggplot')
import seaborn as sns

from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"


def clusterplot_by_percentage_full():

    data_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(data_file)
    df = df[df["cell_line"] == "HepG2"].reset_index(drop=True)

    return df

    df = df[["tf", "enh_fracs", "prom_fracs"]]
    # df = df[["tf", "all_fracs", "enh_fracs", "prom_fracs"]]
    # df = df[["tf", "pk_fracs_enh", "pk_fracs_prom"]]

    tfs = df.pop("tf").values

    g = sns.clustermap(df, cmap="vlag", col_cluster=False, metric="euclidean", yticklabels=True, vmin=0, vmax=100,
                       dendrogram_ratio=0.15)
    hm = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width * 0.15, hm.height])
    g.ax_cbar.set_position((0, 0.7, .03, 0.1))

    tf_labels = [tfs[int(_.get_text())] for _ in g.ax_heatmap.axes.get_yticklabels()]
    g.ax_heatmap.axes.set_yticklabels(tf_labels, size=1)
    g.ax_heatmap.axes.set_xticklabels(["HOT enh (%)", "HOT prom (%)"], rotation=45, ha="right")
    reor_ind = g.dendrogram_row.reordered_ind
    for i in reor_ind:
        print(i, tfs[i], df.loc[i].values)

    save_file = join(PLOTS_DIR, "clustermap_tfs_HOT_enh_proms_big.pdf")
    plt.savefig(save_file, bbox_inches='tight')


def clusterplot_by_percentage_essentials():

    data_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(data_file)
    df = df[df["cell_line"] == "HepG2"].reset_index(drop=True)

    df = df[["tf", "enh_fracs", "prom_fracs"]]

    essential_tfs = ["ARID4B", "MAX", "SAP130", "KDM1A", "FOXP1", "RCOR2", "ZNF217", "TEAD1", "SOX5", "NR2F6", "NFIL3", "HNF4A",
         "FOXA2", "CEBPA", "FOXO1", "PPARG", "MIXL1", "FOXA1", "CEBPG", "ZGPAT", "MED1", "TFAP4", "ZFX", "EGR1",
         "LIN54", "ZNF574", "HDAC1", "TBX2", "HOXA3", "KAT8", "THAP11", "KLF16", "PATZ1", "ERF", "ZNF331", "LCORL",
         "IRF2", "SKI", "ISL2", "ZBTB7B", "POGZ", "IKZF5", "HNF1B", "FOSL2", "TCF7L2", "LCOR", "FOXP4", "BCL6", "RXRB",
         "RARA", "ELF3", "GATAD2A", "SMAD4", "HDAC2", "FOXA3", "SOX6", "PAXIP1", "ZNF687", "KMT2A", "E2F4", "NONO",
         "HNRNPLL", "RBM39", "POLR2A", "POLR2AphosphoS5", "RBFOX2", "ARID4A", "TAF1", "PHF8", "POLR2G", "HMGXB4", "ZFY",
         "GABPB1", "ASH2L", "KDM2A", "MNX1", "UBTF", "GATAD1", "ZNF501", "DMAP1", "NR2C2", "DRAP1", "KMT2B", "YEATS4",
         "SPEN", "MAZ", "TFDP2"]

    df = df[df.tf.isin(essential_tfs)].reset_index(drop=True)

    tfs = df.pop("tf").values

    g = sns.clustermap(df, cmap="vlag", col_cluster=False, metric="seuclidean", yticklabels=True, vmin=0, vmax=100,
                       dendrogram_ratio=0.07, cbar_pos=(0, 0.06, 0.01, 0.05))
    hm = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width * 0.1, hm.height])
    tf_labels = [tfs[int(_.get_text())][:8] for _ in g.ax_heatmap.axes.get_yticklabels()]
    g.ax_heatmap.axes.set_yticklabels(tf_labels, size=8.5)

    reor_ind = g.dendrogram_row.reordered_ind
    for i in reor_ind:
        print(i, tfs[i], df.loc[i].values)

    g.ax_heatmap.axes.set_xticklabels([])

    save_file = join(PLOTS_DIR, "clustermap_tfs_HOT_enh_proms_smaller.pdf")
    print(save_file)
    plt.savefig(save_file)


def tf_cluster_stats():

    data_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(data_file)
    df = df[df["cell_line"] == "HepG2"].reset_index(drop=True)

    essential_tfs = ["ARID4B", "MAX", "SAP130", "KDM1A", "FOXP1", "RCOR2", "ZNF217", "TEAD1", "SOX5", "NR2F6", "NFIL3", "HNF4A",
         "FOXA2", "CEBPA", "FOXO1", "PPARG", "MIXL1", "FOXA1", "CEBPG", "ZGPAT", "MED1", "TFAP4", "ZFX", "EGR1",
         "LIN54", "ZNF574", "HDAC1", "TBX2", "HOXA3", "KAT8", "THAP11", "KLF16", "PATZ1", "ERF", "ZNF331", "LCORL",
         "IRF2", "SKI", "ISL2", "ZBTB7B", "POGZ", "IKZF5", "HNF1B", "FOSL2", "TCF7L2", "LCOR", "FOXP4", "BCL6", "RXRB",
         "RARA", "ELF3", "GATAD2A", "SMAD4", "HDAC2", "FOXA3", "SOX6", "PAXIP1", "ZNF687", "KMT2A", "E2F4", "NONO",
         "HNRNPLL", "RBM39", "POLR2A", "POLR2AphosphoS5", "RBFOX2", "ARID4A", "TAF1", "PHF8", "POLR2G", "HMGXB4", "ZFY",
         "GABPB1", "ASH2L", "KDM2A", "MNX1", "UBTF", "GATAD1", "ZNF501", "DMAP1", "NR2C2", "DRAP1", "KMT2B", "YEATS4",
         "SPEN", "MAZ", "TFDP2"]

    aux_df = df[~df.tf.isin(essential_tfs)]
    ess_df = df[df.tf.isin(essential_tfs)]

    print("Auxiliary DAPs. HOT: %f, HOT-enh: %f, HOT-prom: %f" % (aux_df["all_fracs"].median(), aux_df["enh_fracs"].median(), aux_df["prom_fracs"].median()))
    print("Essential DAPs. HOT: %f, HOT-enh: %f, HOT-prom: %f" % (ess_df["all_fracs"].median(), ess_df["enh_fracs"].median(), ess_df["prom_fracs"].median()))

    cl_list = [
        ["ZNF687", "ARID4B", "MAX", "SAP130"],
        ["KMT2A", "E2F4", "NONO", "HNRNPLL", "RBM39", "POLR2A", "POLR2AphosphoS5", "RBFOX2", "ARID4A", "TAF1", "ZFY",
         "GABPB1", "PHF8", "POLR2G", "NR2C2", "DRAP1", "YEATS4", "HMGXB4", "KMT2B", "TFDP2", "MAZ", "SPEN", "ASH2L",
         "KDM2A", "MNX1", "UBTF", "GATAD1", "ZNF501", "DMAP1"],
        ["NR2F6", "TEAD1", "SOX5", "NFIL3", "HNF4A", "FOXA2", "PPARG", "MIXL1", "FOXA1", "CEBPG", "FOXO1", "ZNF217", "CEBPA", "KDM1A", "FOXP1", "RCOR2"],
        ["ZGPAT", "MED1", "TFAP4", "ZFX", "EGR1", "LIN54", "ZNF574", "HDAC1", "TBX2", "THAP11", "KAT8", "HOXA3"],
        ["KLF16", "PATZ1", "ERF", "ZNF331", "LCORL", "IRF2", "SKI", "ISL2", "ZBTB7B", "POGZ", "IKZF5", "HNF1B", "FOSL2", "TCF7L2", "RXRB", "RARA", "LCOR", "FOXP4", "BCL6", "GATAD2A","SOX6", "PAXIP1", "ELF3", "SMAD4", "HDAC2", "FOXA3"]
    ]

    for i, cl_tfs in enumerate(cl_list):
        cl_df = df[df.tf.isin(cl_tfs)]
        print("CL %d: num=%d,  HOT=%f, HOT-enh:%f, HOT-prom:%f" % (i+1, cl_df.shape[0], cl_df["all_fracs"].mean(), cl_df["enh_fracs"].mean(), cl_df["prom_fracs"].mean()))

    # return df, aux_df, ess_df


def extract_TF_peaks_HOT_percentages():

    peaks_dir = "ROOT_DIR/chipseq_files/peaks"
    a = load_metadata()

    hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.PE.bed"))
    proms = BedTool([r for r in hots if r.fields[-1] == "P"])
    enhs = BedTool([r for r in hots if r.fields[-1] == "E"])

    # tfs = ["EP300", "CTCF", "POLR2A", "POLR2G", "HNF4A", "TCF7L2", "MAX", "ARID4A", "FOXA2", "FOXA1", "GATAD2A", "MED1",
    #        "FOXA3", "SOX6", "ZNF687", "SAP130", "ARID4B", "SOX6", "TFAP4", "KAT8", "HOXA3"]

    tf2peaks = {}
    tf2percs = {}

    cnt = 1
    for tf, tf_id in a["HepG2"].items():
        print(cnt, len(a["HepG2"]), tf)
        # tf_id = a["HepG2"][tf]
        tf_bed = BedTool(join(peaks_dir, "%s_hg19_peak_8bp.bed" % tf_id))
        tf_peaks = tf_bed.count()

        try:
            prom_peaks = tf_bed.intersect(proms, wo=True).groupby(c=4, o="collapse").count()
            enh_peaks = tf_bed.intersect(enhs, wo=True).groupby(c=4, o="collapse").count()
        except:
            continue
        prom_perc = 100 * prom_peaks/tf_peaks
        enh_perc = 100 * enh_peaks/tf_peaks

        tf2peaks[tf] = tf_peaks
        tf2percs[tf] = [prom_perc, enh_perc]

        cnt += 1

    return tf2peaks, tf2percs


def pie_chart_TFs_HOTs(tf):

    data_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(data_file)
    df = df[(df["cell_line"] == "HepG2") & (df["tf"] == tf)].reset_index(drop=True)

    prom_frac = df["pk_fracs_prom"].values[0]
    enh_frac = df["pk_fracs_enh"].values[0]
    num_peaks = df["num_peaks"].values[0]
    prom_frac = np.round(prom_frac)
    enh_frac = np.round(enh_frac)

    colors = ['#1f77b4', '#ff7f0e', "#808080"]
    percs = [prom_frac, enh_frac,  100 - prom_frac - enh_frac]
    labels = ["HOT prom", "HOT enh", "rest"]
    # plt.figure(figsize=(2.4, 2.4))
    plt.figure(figsize=(2.4, 2.4))
    ax = plt.gca()

    n = ax.pie(percs, autopct='%1.0f%%', counterclock=False, colors=colors, startangle=90)
    n[0][2].set_alpha(0.5)
    # for i in range(len(n[0])):
    #     n[0][i].set_alpha(0.5)
    ax.axis('equal')

    num_str = str(num_peaks)
    num_str = ",".join([num_str[i:i + 3] for i in range(0, len(num_str), 3)])
    num_str = num_str[::-1]

    if tf == "EP300":
        tf = "P300"
    if tf == "RAD21":
        tf = "Cohesin"

    # if tf == "RAD21":
    #     ax.set_title("%s" % tf, y=0.85)
    # else:
    ax.set_title("%s" % tf)

    ax.text(0, -1.35, "peaks=%s" % num_str, va='bottom', ha='center')
    plt.tight_layout()

    save_file = join(PLOTS_DIR, "TF_pie_charts/%s.pdf" % tf)
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    # plt.close()


def pie_chart_TFs_HOTs_iterate(tf2peaks, tf2percs):

    # tfs = ["EP300", "CTCF", "POLR2A", "POLR2G", "HNF4A", "TCF7L2"]

    for tf in tf2peaks.keys():
        pie_chart_TFs_HOTs(tf, tf2peaks, tf2percs)
        # return
