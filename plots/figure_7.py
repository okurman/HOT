#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")

from pybedtools import BedTool
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import upsetplot
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path

import os
DATA_PATH = Path(os.environ["HOT_DATA"])
LOCI_PATH = DATA_PATH / "log_bins"
PLOTS_DIR = DATA_PATH / "plots/figure_7"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)

from HOTs import H1_subsamples


def overbinders_upset_plots(save_file, common_tfs=False):
	if common_tfs:

		hepg2 = BedTool(str(DATA_PATH / "src_files/HOTs_with_common_tfs/HepG2_HOTs.bed.gz"))
		k562 = BedTool(str(DATA_PATH / "src_files/HOTs_with_common_tfs/K562_HOTs.bed.gz"))
		h1 = BedTool(str(DATA_PATH / "src_files/HOTs_with_common_tfs/H1_HOTs.bed.gz"))

	else:

		hepg2 = BedTool(str(DATA_PATH / "HOTs/HepG2_HOTs.bed.gz"))
		k562 = BedTool(str(DATA_PATH / "HOTs/K562_HOTs.bed.gz"))
		h1 = BedTool(str(DATA_PATH / "HOTs/H1_HOTs.bed.gz"))

	multi_ix = pd.MultiIndex.from_tuples([(True, False, False),
										  (False, True, False),
										  (False, False, True),
										  (True, True, False),
										  (True, False, True),
										  (False, True, True),
										  (True, True, True)], names=["HepG2", "K562", "H1"])

	data = []
	data.append(hepg2.intersect(k562.cat(h1, postmerge=False), v=True).count())
	data.append(k562.intersect(hepg2.cat(h1, postmerge=False), v=True).count())
	data.append(h1.intersect(hepg2.cat(k562, postmerge=False), v=True).count())

	data.append(hepg2.intersect(k562, wa=True).intersect(h1, v=True).count())
	data.append(hepg2.intersect(h1, wa=True).intersect(k562, v=True).count())
	data.append(k562.intersect(h1, wa=True).intersect(hepg2, v=True).count())
	data.append(hepg2.intersect(k562, wa=True).intersect(h1, wa=True).count())

	df = pd.Series(data, index=multi_ix)

	fig = plt.figure(figsize=(4, 3))
	upsetplot.plot(df, fig=fig)
	plt.savefig(save_file, bbox_inches='tight')


def conservation_scores_H1(save_file):

	data = []
	for cl in ["HepG2", "K562", "H1"]:
		for i in range(10, 14):
			f = str(DATA_PATH / "phastCons" / f"{cl}_400_loci.{i}.bed.noprom.vertebrate.phastcons.gz")
			scores = np.loadtxt(f, usecols=(3))
			[data.append([cl, "enh", _s]) for _s in scores]

			f = str(DATA_PATH / "phastCons" / f"{cl}_400_loci.{i}.bed.prom.vertebrate.phastcons.gz")
			scores = np.loadtxt(f, usecols=(3))
			[data.append([cl, "prom", _s]) for _s in scores]
	df = pd.DataFrame(data=data, columns=["cell_line", "category", "score"])
	# return df

	for cl in ["HepG2", "K562", "H1"]:
		prom_scores = df[(df["cell_line"]==cl) & (df["category"]=="prom")]["score"].values
		enh_scores = df[(df["cell_line"]==cl) & (df["category"]=="enh")]["score"].values
		significance = mannwhitneyu(prom_scores, enh_scores)
		ratio = np.mean(prom_scores)/np.mean(enh_scores)
		print(f"cell line: {cl}, prom/enh= {ratio}, significance={significance}")

	plt.close()
	plt.figure(figsize=(4, 3))

	ax = plt.gca()
	g = sns.boxplot(x="cell_line", y="score", hue="category", data=df, showfliers=False, ax=ax)

	l = g.legend(loc="upper center", ncol=2)
	l.get_texts()[0].set_text('enhancers')
	l.get_texts()[1].set_text('promoters')
	ax.set_ylim([0, 1])
	ax.grid(axis="y", alpha=0.2)
	ax.set_xlabel("")
	ax.set_ylabel("phastCons score")

	plt.tight_layout()
	plt.savefig(save_file)


if __name__ == "__main__":

	save_file = PLOTS_DIR/"Figure7_a.pdf"
	print(save_file)
	overbinders_upset_plots(save_file, common_tfs=False)

	save_file = PLOTS_DIR / "Figure7_b.pdf"
	print(save_file)
	overbinders_upset_plots(save_file, common_tfs=True)

	save_file = PLOTS_DIR / "Figure7_c.pdf"
	print(save_file)
	H1_subsamples.plot_subsample_h1(save_file)

	save_file = PLOTS_DIR / "Figure7_d.pdf"
	print(save_file)
	conservation_scores_H1(save_file)