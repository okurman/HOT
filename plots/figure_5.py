#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import sys

from scipy.stats import mannwhitneyu

sys.path.append("../")
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plots_data_factory
from pybedtools import BedTool
import pandas as pd
import seaborn as sns

DATA_PATH = Path("../data/data")
PLOTS_DIR = DATA_PATH / "plots/figure_5"
PLOTS_DIR.mkdir(exist_ok=True)


def plot_phastcons_distribution_log_bins(save_file, species="vertebrate"):

	line_parser = lambda x: x.split("\t")[3]
	file_amender = lambda x: x.replace(".gz", f".{species}.phastcons.gz")
	df = plots_data_factory.load_bins_to_df(line_parser=line_parser, file_amender=file_amender)

	plots_data_factory.plot_categorical_worker(df, y_label="phastCons", plot_type="lineplot", save_file=save_file,
							leg_loc="upper left", palette="colorblind", figsize=(5, 3.8))


def load_categorized_phastcons_data(species="vertebrate"):

	hots = BedTool(str(DATA_PATH/f"HOTs/HepG2_HOTs.proms.bed.{species}.phastcons.gz"))
	regular_enhancers = BedTool(str(DATA_PATH/f"src_files/HepG2_enhancers_DHS_H3K27ac.bed.{species}.phastcons.gz"))
	exons = BedTool(str(DATA_PATH/f"src_files/hg19_files/knownGene.exons.merged.bed.{species}.phastcons.gz"))

	data = []
	for name, scores_bed in zip(["exons", "regular enhancers", "HOT"],
								[exons, regular_enhancers, hots]):
		data += [[name, 0 if r.fields[3] == "-" else float(r.fields[3])] for r in scores_bed]

	df = pd.DataFrame(data=data, columns=["region", "value"])

	return df


def plot_3panels_plots(save_file, species="vertebrate"):

	df = load_categorized_phastcons_data(species)

	plt.figure(figsize=(3, 3))
	ax = plt.gca()

	hot_scores = df[(df["region"] == "HOT")]["value"]
	enh_scores = df[(df["region"] == "regular enhancers")]["value"]
	exon_scores = df[(df["region"] == "exons")]["value"]

	print(mannwhitneyu(hot_scores, enh_scores))
	print(mannwhitneyu(exon_scores, enh_scores))

	data = [["reg enh", 1],
			["HOT", hot_scores.mean() / enh_scores.mean()],
			["exons", exon_scores.mean() / enh_scores.mean()]]

	df = pd.DataFrame(columns=["region", "value"], data=data)

	sns.barplot(data=df, x="region", y="value", ax=ax)
	ax.grid(axis="y", alpha=0.3)
	ax.set_ylabel(f"phastCons {species}")
	ax.set_xlabel("")
	ax.set_ylim([0, 5])

	plt.tight_layout()
	ax.set_xticklabels(["regular\nenhancers", "HOT", "exons"])

	plt.savefig(save_file)
	plt.close()


if __name__ == "__main__":

	species = sys.argv[1] if len(sys.argv) > 1 else "vertebrate"

	save_file = PLOTS_DIR / f"Figure3_a.{species}.pdf"
	print(save_file)
	plot_phastcons_distribution_log_bins(save_file, species=species)

	save_file = PLOTS_DIR / f"Figure5_b.{species}.pdf"
	print(save_file)
	plot_3panels_plots(save_file, species=species)