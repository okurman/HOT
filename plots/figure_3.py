#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import sys
sys.path.append("../")
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import numpy as np
from os.path import join
from HOTs import HiC_analysis
from collections import defaultdict
import pandas as pd

DATA_PATH = Path("../data/data")
BINS_DIR = DATA_PATH / "log_bins"
PLOTS_DIR = DATA_PATH / "plots/figure_3"
PLOTS_DIR.mkdir(exist_ok=True)

HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]


def plot_binwise_contacts_encode(M, save_file):

	# M, _ = HiC_analysis.extract_contact_lists_encode()

	# return M

	cl = "HepG2"

	# contacts_file = join(PROJECT_DIR, "looping/encode/%s_log_binwise_contact_counts.txt" % cl)
	# M = np.loadtxt(contacts_file)

	files = get_loci_files(cl)

	locus_counts = np.asarray([BedTool(f).count() for f in files])

	M_freq = np.zeros_like(M)
	for i in range(len(files)):
		for j in range(i, len(files)):
			M_freq[i, j] = M[i, j] / (locus_counts[i] * locus_counts[j])

	M_norm = M_freq / np.max(M_freq)
	M_norm = M_norm + M_norm.T - np.diag(np.diag(M_norm))

	cmap = sns.cm.rocket_r
	plt.close("all")
	g = sns.heatmap(M_norm, cmap=cmap)

	tick_labels = HEPG2_XTICK_LABELS
	ticks = np.arange(len(tick_labels))

	g.yaxis.set_ticks(ticks)
	g.yaxis.set_ticklabels(tick_labels)
	g.tick_params(axis="y", rotation=0)

	g.xaxis.set_ticks(ticks)
	g.xaxis.set_ticklabels(tick_labels)
	g.tick_params(axis="x", rotation=45)

	plt.xlabel("DAPs")
	plt.ylabel("DAPs")
	plt.tight_layout()

	print(save_file)
	plt.savefig(save_file)


def plot_loci_contact_correlations(save_file):

	G = HiC_analysis.load_graph()
	hots = BedTool(str(DATA_PATH/"HOTs/HepG2_HOTs.bed.gz"))

	nodes_bed = BedTool("\n".join([n.replace("-", "\t") for n in G.nodes]), from_string=True)

	nodes_bed_loci = nodes_bed.intersect(hots, wo=True).groupby(c=[9], o="collapse")

	degree2tfs = defaultdict(list)
	degree2hots = defaultdict(list)
	hot_nodes = []
	for r in nodes_bed_loci:
		node = "%s-%d-%d" % (r.chrom, r.start, r.stop)
		tfs = [int(_) for _ in r.fields[3].split(",")]
		degree = len(G[node].values())
		degree2tfs[degree].append(max(tfs))
		degree2hots[degree].append(len(tfs))
		hot_nodes.append(node)

	data = []
	for d, tfs in degree2tfs.items():
		for tf in tfs:
			data.append([d, tf])
	df_tfs = pd.DataFrame(data=data, columns=["d", "tfs"])

	data = []
	for d, hots in degree2hots.items():
		for hot in hots:
			data.append([d, hot])
	df_hots = pd.DataFrame(data=data, columns=["d", "tfs"])

	fig, axs = plt.subplots(4, 1, gridspec_kw={'height_ratios': [.5, 1, 3, 3]}, sharex=True, figsize=(5, 7))

	d2count = defaultdict(int)
	degrees = sorted((d for n, d in G.degree()), reverse=True)
	degrees = list(np.unique(degrees, return_counts=True))
	for k, v in zip(degrees[0], degrees[1]):
		d2count[k] = v

	hot_d2count = defaultdict(int)
	hot_degrees = sorted((len(G[n].values()) for n in hot_nodes), reverse=True)
	hot_degrees = list(np.unique(hot_degrees, return_counts=True))
	for k, v in zip(hot_degrees[0], hot_degrees[1]):
		hot_d2count[k] = v

	hot_degrees_num = len(hot_degrees[0])
	degrees_num = len(degrees[0])

	if degrees_num - hot_degrees_num > 0:
		for d in range(len(hot_degrees[0]) + 1, len(degrees[0]) + 1):
			hot_degrees[0] = np.append(hot_degrees[0], d)
			hot_degrees[1] = np.append(hot_degrees[1], 0)

	for x, contacts in zip(degrees[0], degrees[1]):
		contacts = str(contacts)
		axs[0].text(x - 1, 0, contacts, horizontalalignment='center')
	axs[0].axis("off")

	x_range = np.asarray(degrees[0]) - 1
	hot_perc = [int(100 * hot_degrees[1][i] / degrees[1][i]) for i, d in enumerate(degrees[0])]
	axs[1].bar(x_range, hot_perc)
	hot_perc = [str(_) + "%" for _ in hot_perc]
	for x, perc in zip(x_range, hot_perc):
		axs[1].text(x, 8, perc, color="white", horizontalalignment='center')
	axs[1].tick_params('x', labelbottom=False)
	axs[1].grid(axis="y", which="both")
	axs[1].set_ylabel("% with HOTs")

	sns.boxplot(x="d", y="tfs", data=df_hots, color="green", boxprops=dict(alpha=.5), ax=axs[2], showfliers=False)
	axs[2].set_ylabel("# of HOTs in anchors")
	axs[2].set_xlabel("")
	axs[2].tick_params('x', labelbottom=False)

	sns.boxplot(x="d", y="tfs", data=df_tfs, color="red", boxprops=dict(alpha=.5), ax=axs[3], showfliers=False)
	axs[3].set_ylabel("DAPs bound to HOTs")

	axs[3].set_xlabel("# of chromatin contacts")

	fig.align_labels()
	plt.tight_layout()
	plt.savefig(save_file)

	plt.close("all")


if __name__ == "__main__":

	save_file = PLOTS_DIR/"Figure3_a.pdf"
	plot_binwise_contacts_encode(save_file)

