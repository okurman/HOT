#!/usr/bin/env python
# -*- coding: utf-8 -*-

# os.environ["HOT_CODE"] = "/Users/okurman/Projects/HOT/code/github/"
# os.environ["HOT_DATA"] = "/Users/okurman/Projects/HOT/data_reproducible/data"
# SYS_CODE_PATH = "/Users/okurman/Projects/HOT/code/github/plots"
# SYS_DATA_PATH = "/Users/okurman/Projects/HOT/data_reproducible/data"


import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")

from pathlib import Path
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import plots_data_factory
from pybedtools import BedTool
import pandas as pd
import seaborn as sns

import os
DATA_PATH = Path(os.environ["HOT_DATA"])
PLOTS_DIR = DATA_PATH / "plots/figure_5"
PHASTCONS_DIR = DATA_PATH/"phastCons"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)


def plot_phastcons_distribution_log_bins(df, species="vertebrate"):
	
	save_file = PLOTS_DIR / f"Figure5_a.pdf"
	print(save_file)
	
	# line_parser = lambda x: x.split("\t")[3]
	# file_amender = lambda x: x.replace(".gz", f".{species}.phastcons.gz")
	# df = plots_data_factory.load_bins_to_df(line_parser=line_parser, file_amender=file_amender, src_dir=PHASTCONS_DIR)

	# return df
	
	plots_data_factory.plot_categorical_worker(df, y_label="phastCons", plot_type="lineplot", save_file=save_file,
							leg_loc="upper left", palette="colorblind", figsize=(6, 3.4))


def load_categorized_phastcons_data(species="vertebrate"):

	hots = BedTool(str(PHASTCONS_DIR/f"HepG2_HOTs.proms.bed.{species}.phastcons.gz"))
	regular_enhancers = BedTool(str(PHASTCONS_DIR/f"HepG2_enhancers_DHS_H3K27ac.bed.{species}.phastcons.gz"))
	exons = BedTool(str(PHASTCONS_DIR/f"knownGene.exons.merged.bed.{species}.phastcons.gz"))

	data = []
	for name, scores_bed in zip(["exons", "reg.enh", "HOT"],
								[exons, regular_enhancers, hots]):
		data += [[name, 0 if r.fields[3] == "-" else float(r.fields[3])] for r in scores_bed]

	df = pd.DataFrame(data=data, columns=["region", "value"])

	return df


def plot_3panels_plots(df, species="vertebrate"):

	save_file = PLOTS_DIR / f"Figure5_b.pdf"
	print(save_file)

	# df = load_categorized_phastcons_data(species)

	plt.figure(figsize=(1.8, 3))
	ax = plt.gca()

	hot_scores = df[(df["region"] == "HOT")]["value"]
	enh_scores = df[(df["region"] == "regular enhancers")]["value"]
	exon_scores = df[(df["region"] == "exons")]["value"]

	print(mannwhitneyu(hot_scores, enh_scores))
	print(mannwhitneyu(exon_scores, enh_scores))

	hot_ratio = hot_scores.mean() / enh_scores.mean()
	exon_ratio = exon_scores.mean() / enh_scores.mean()

	data = [["reg enh", 1],
			["HOT", hot_ratio],
			["exons", exon_ratio]]

	df = pd.DataFrame(columns=["region", "value"], data=data)

	b = sns.barplot(data=df, x="region", y="value", ax=ax, hue="region")

	ax.text(1, hot_ratio, f"***", ha='center', va='bottom', color="k")
	ax.text(2, exon_ratio, f"***", ha='center', va='bottom', color="k")
	
	ax.grid(axis="y", alpha=0.3)
	ax.set_ylabel(f"phastCons")
	ax.set_xlabel("")
	ax.set_ylim([0, 5])

	ax.set_xticklabels(["reg.enh", "HOT", "exons"], rotation=45, ha="right")

	plt.tight_layout()
	plt.savefig(save_file)
	plt.close()




def pack_classification_results():

	data_file = DATA_PATH / "src_files/classification_results.txt"
	df = pd.read_csv(data_file, sep="\t")

	df_cnn = df[(df["method"] == "seq - CNN") & (df["len"] == 1000)]
	df_lsgkm = df[(df["method"] == "seq - SVM") & (df["kernel"] == "gapped-kmer") & (df["len"] == 1000)]

	df_feat = df[(df["method"].isin(["feat - LogReg", "feat - SVM"])) & (df["kernel"].isin(["linear", "-"]))]
	df_feat["method"] = df_feat["method"] + " - " + df_feat["feature"]

	df = pd.concat((df_cnn, df_lsgkm, df_feat)).reset_index(drop=True)

	return df


def pack_classification_results_brief():

	data_file = DATA_PATH / "src_files/classification_results.txt"
	df = pd.read_csv(data_file, sep="\t")

	df_cnn = df[(df["method"] == "seq - CNN") & (df["len"] == 1000)]
	df_lsgkm = df[(df["method"] == "seq - SVM") & (df["kernel"] == "gapped-kmer") & (df["len"] == 1000)]

	df_feat = df[(df["method"]=="feat - SVM") & (df["kernel"]=="linear") & (df["feature"]=="all")]

	df = pd.concat((df_cnn, df_lsgkm, df_feat)).reset_index(drop=True)

	return df


def plot_combined_aucs(save_file):

	df = pack_classification_results()

	method_order = ['seq - CNN', 'seq - SVM', 'feat - LogReg - all', 'feat - LogReg - GC', 'feat - LogReg - CpG',
					'feat - SVM - all', 'feat - SVM - GC', 'feat - SVM - CpG']
	x_tick_labels = ['seq - CNN', 'seq - SVM', 'feat - LogReg', 'GC - LogReg', 'CpG - LogReg', 'feat - SVM', 'GC - SVM',
					 'CpG - SVM']

	plt.close("all")
	fig = plt.figure(constrained_layout=True, figsize=(8, 6))

	subfigs = fig.subfigures(nrows=2, ncols=1)
	for subfig, cl in zip(subfigs, ["HepG2", "K562"]):
		subfig.suptitle(cl)

		axs = subfig.subplots(nrows=1, ncols=3)
		for ax, control in zip(axs, ["dhs", "proms", "re"]):
			_df = df[(df["cl"] == cl) & (df["control"] == control)]

			ax.set_title(f'controls={control}')
			g = sns.barplot(data=_df, x="method", y="value", hue="metric", ax=ax, order=method_order, palette="Set1")

			ax.set_ylim([0, 1])
			ax.set_ylabel("AUC" if control == "dhs" else "")
			ax.set_xlabel("")

			ax.grid(axis="y", which="both", alpha=0.6)
			ax.set_xticks(range(len(method_order)))
			ax.set_xticklabels(x_tick_labels, size=9, rotation=45, ha="right")

			g.legend(ncols=2,
					 frameon=False,
					 markerscale=0.2,
					 fontsize=9,
					 columnspacing=0.2,
					 labelspacing=0.1,
					 handletextpad=0.1,
					 borderaxespad=0.1,
					 loc="upper right")

	plt.savefig(save_file, bbox_inches='tight')
	plt.show()


def plot_combined_aucs_brief():
	
	save_file = PLOTS_DIR / f"Figure5_c_brief.pdf"
	print(save_file)

	df = pack_classification_results_brief()
	
	method_order = ['seq - CNN', 'seq - SVM', 'feat - SVM']
	x_tick_labels = method_order

	plt.close("all")
	fig, axs = plt.subplots(1, 2, figsize=(3.5, 3))

	df.loc[df["control"]=="dhs", "control"] = "DHS"
	df.loc[df["control"]=="proms", "control"] = "reg.prom"
	df.loc[df["control"]=="re", "control"] = "reg.enh"

	for ax, cl in zip(axs, ["HepG2", "K562"]):

		_df = df[(df["cl"]==cl) * (df["metric"]=="auROC")]

		g = sns.stripplot(data=_df, x="method", y="value", hue="control", order=method_order, palette="Set2",
			jitter=False, s=10, ax=ax)

		# g = sns.stripplot(data=_df, x="method", y="value", hue="control", order=method_order, palette="Set1",
		# 	jitter=False, s=10, ax=ax, legend=False)

		ax.grid(axis="y", which="both", alpha=0.6)
		ax.set_xticks(range(len(method_order)))
		ax.set_xticklabels(x_tick_labels, size=9, rotation=45, ha="right")
		ax.set_xlabel("")
		ax.set_ylabel("auROC")
		ax.set_ylim([0.4, 1])
		ax.set_title(cl)
		ax.legend(title="HOT vs.",
			# markerscale=0.2,
			fontsize=9,
			columnspacing=0.2,
			labelspacing=0.1,
			handletextpad=0.1,
			borderaxespad=0.1)

	# axs[0].get_legend().remove()
	# sns.move_legend(g, "upper left", bbox_to_anchor=(1.015, 1.015), title="")

	plt.tight_layout()

	plt.savefig(save_file, bbox_inches='tight')
	plt.close("all")

	return

	subfigs = fig.subfigures(nrows=2, ncols=1)
	for subfig, cl in zip(subfigs, ["HepG2", "K562"]):
		subfig.suptitle(cl)

		axs = subfig.subplots(nrows=1, ncols=3)
		for ax, control in zip(axs, ["dhs", "proms", "re"]):
			_df = df[(df["cl"] == cl) & (df["control"] == control)]

			ax.set_title(f'controls={control}')
			g = sns.barplot(data=_df, x="method", y="value", hue="metric", ax=ax, order=method_order, palette="Set1")

			ax.set_ylim([0, 1])
			ax.set_ylabel("AUC" if control == "dhs" else "")
			ax.set_xlabel("")

			ax.grid(axis="y", which="both", alpha=0.6)
			ax.set_xticks(range(len(method_order)))
			ax.set_xticklabels(x_tick_labels, size=9, rotation=45, ha="right")

			g.legend(ncols=2,
					 frameon=False,
					 # markerscale=0.2,
					 # fontsize=9,
					 # columnspacing=0.2,
					 # labelspacing=0.1,
					 # handletextpad=0.1,
					 borderaxespad=0.1)

	plt.savefig(save_file, bbox_inches='tight')
	plt.show()



if __name__ == "__main__":

	species = sys.argv[1] if len(sys.argv) > 1 else "vertebrate"

	save_file = PLOTS_DIR / f"Figure5_a.pdf"
	print(save_file)
	plot_phastcons_distribution_log_bins(save_file, species=species)

	save_file = PLOTS_DIR / f"Figure5_b.pdf"
	print(save_file)
	plot_3panels_plots(save_file, species=species)

	save_file = PLOTS_DIR / f"Figure5_c.pdf"
	print(save_file)
	plot_combined_aucs(save_file)

