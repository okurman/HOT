#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
sys.path.append(os.environ["HOT_CODE"])

import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path

DATA_PATH = Path(os.environ["HOT_DATA"])
LOCI_PATH = DATA_PATH / "log_bins"
PLOTS_DIR = DATA_PATH / "plots/figure_8"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)


from HOTs import variant_analysis


def plot_enriched_traits(save_file):

	fname = DATA_PATH / "src_files/GWAS_enrichments.txt.gz"
	df = pd.read_csv(fname, sep="\t")

	ind = (df["logneg_pval_hot_enh"] > 3) & (df["logneg_pval_hot_prom"] > 3)
	df = df[ind].sort_values(by="enr_hot", ascending=False).reset_index(drop=True)

	traits = df["trait"].values

	data = []
	names = ["hot_prom", "rp", "hot_enh", "re", "dhs"]
	label_names = ["HOT prom", "reg. prom", "HOT enh", "reg. enh", "DHS"]

	for ind, row in df.iterrows():

		for y, name in enumerate(names):

			pvalue = row["logneg_pval_%s" % name]
			enr = row["enr_%s" % name]

			if pvalue < 3: continue

			_res = [ind, y, enr, min(pvalue, 30)]
			data.append(_res)

	df_plot = pd.DataFrame(data=data, columns=["trait", "y", "fold-enrichment", "-log10(p-value)"])

	plt.figure(figsize=(6, 6))
	ax = plt.gca()

	g = sns.scatterplot(x="trait", y="y", size="fold-enrichment", hue="-log10(p-value)", data=df_plot, ax=ax,
						sizes=(20, 200))

	handles, labels = g.get_legend_handles_labels()

	pval_hl = handles[1:7], labels[1:7]
	enr_hl = handles[8:], labels[8:]

	leg1 = ax.legend(*pval_hl, bbox_to_anchor=(0.05, 1.01), ncol=6, title="-log10(p-value)", loc="lower left",
					 handletextpad=0, columnspacing=0, title_fontsize=8, fontsize=7)
	leg2 = ax.legend(*enr_hl, bbox_to_anchor=(0.51, 1.01), ncol=6, title="fold-enrichment", loc="lower left",
					 handletextpad=0.2, columnspacing=0.5, title_fontsize=8, fontsize=7)
	ax.add_artist(leg1)

	ax.grid(axis="both", alpha=0.2, which="both")

	ax.yaxis.tick_right()
	ax.set_yticks([0, 1, 2, 3, 4])
	ax.set_ylim([-1, 5])
	ax.set_yticklabels(label_names)

	ax.set_xticks(np.arange(len(traits)))
	ax.set_xticklabels(traits, rotation=90)

	ax.set_xlabel("")
	ax.set_ylabel("")

	plt.tight_layout()
	plt.savefig(save_file)
	plt.close()


def plot_enriched_traits_adjusted(save_file):

	fname = DATA_PATH / "src_files/GWAS_enrichments.p_adj.txt.gz"
	df = pd.read_csv(fname, sep="\t")

	# at least 3 decimal points
	ind = (df["logneg_pval_hot_enh"] > 3) & (df["logneg_pval_hot_prom"] > 3)
	df = df[ind].sort_values(by="enr_hot", ascending=False).reset_index(drop=True)

	# shortened texts for long trait names
	traits = [
		"HDL cholesterol",
		"Sex hormone-binding globulin levels BMI",
		"Serum metabolite levels",
		"Sex hormone-binding globulin levels",
		"Low density lipoprotein chol. levels",
		"Triglyceride levels",
		"Blood protein levels"]

	data = []
	names = ["hot_prom", "rp", "hot_enh", "re", "dhs"]
	label_names = ["HOT prom", "reg. prom", "HOT enh", "reg. enh", "DHS"]

	for ind, row in df.iterrows():
		for y, name in enumerate(names):
			pvalue = row["logneg_pval_%s" % name]
			enr = row["enr_%s" % name]

			if pvalue < 3:
				continue

			_res = [ind, y, enr, min(pvalue, 30)]
			data.append(_res)

	df_plot = pd.DataFrame(data=data, columns=["trait", "y", "fold-enrichment", "-log10(p-value)"])

	plt.figure(figsize=(3.3, 5))
	ax = plt.gca()

	g = sns.scatterplot(x="trait", y="y", size="fold-enrichment", hue="-log10(p-value)", data=df_plot, ax=ax,
						sizes=(20, 200))

	handles, labels = g.get_legend_handles_labels()

	pval_hl = handles[1:7], labels[1:7]
	enr_hl = handles[8:], labels[8:]

	leg1 = ax.legend(*pval_hl, bbox_to_anchor=(1, 1.01), ncol=1, title="-log10(adj. P)", loc="upper left",
					 handletextpad=0, columnspacing=0, title_fontsize=7, fontsize=5)
	leg2 = ax.legend(*enr_hl, bbox_to_anchor=(1, 0.55), ncol=1, title="fold-enrichment", loc="upper left",
					 handletextpad=1, columnspacing=0, labelspacing=1.5, title_fontsize=7, fontsize=5)
	ax.add_artist(leg1)

	ax.grid(axis="both", alpha=0.2, which="both")

	ax.set_yticks([0, 1, 2, 3, 4])
	ax.set_ylim([-0.5, 4.5])
	ax.set_xlim([-0.5, len(traits) - 0.5])
	ax.set_yticklabels(label_names)

	ax.set_xticks(np.arange(len(traits)))
	ax.set_xticklabels(traits, rotation=90)

	ax.set_xlabel("")
	ax.set_ylabel("")

	plt.tight_layout()
	plt.savefig(save_file)
	plt.close()


if __name__ == "__main__":

	save_file = PLOTS_DIR / "Figure8_a.pdf"
	print(save_file)
	variant_analysis.plot_common_variant_densities(save_file)

	save_file = PLOTS_DIR / "Figure8_b.pdf"
	print(save_file)
	plot_enriched_traits_adjusted(save_file)

	save_file = PLOTS_DIR / "Figure8_b.unadjusted_pvalues.pdf"
	print(save_file)
	plot_enriched_traits(save_file)
