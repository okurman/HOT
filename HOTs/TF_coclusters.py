

import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import warnings
warnings.filterwarnings('ignore')
import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
DATA_PATH = Path(os.environ["HOT_DATA"])
LOCI_DIR = DATA_PATH/"HOTs/"


def clusterplot_by_percentage_full(save_file):

    data_file = DATA_PATH/"src_files/tf_stats_summary_table.txt"
    df = pandas.read_table(data_file)
    df = df[df["cell_line"] == "HepG2"].reset_index(drop=True)

    df = df[["tf", "enh_fracs", "prom_fracs"]]

    tfs = df.pop("tf").values

    g = sns.clustermap(df, cmap="vlag", col_cluster=False, metric="euclidean", yticklabels=True, vmin=0, vmax=100,
                       dendrogram_ratio=0.15)
    hm = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width * 0.15, hm.height])
    g.ax_cbar.set_position((0, 0.7, .03, 0.1))

    tf_labels = [tfs[int(_.get_text())] for _ in g.ax_heatmap.axes.get_yticklabels()]
    g.ax_heatmap.axes.set_yticklabels(tf_labels, size=1)
    g.ax_heatmap.axes.set_xticklabels(["HOT enh (%)", "HOT prom (%)"], rotation=45, ha="right")

    plt.savefig(save_file, bbox_inches='tight')
    plt.close("all")


def clusterplot_by_percentage_essentials(save_file):

    data_file = DATA_PATH/"src_files/tf_stats_summary_table.txt"
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

    g.ax_heatmap.axes.set_xticklabels([])

    plt.savefig(save_file)
    plt.close("all")

