#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
sys.path.append(os.environ["HOT_CODE"])

import tempfile
from os.path import join
import numpy as np
import pandas
from pybedtools import BedTool

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from upsetplot import from_memberships
from upsetplot import plot as upset_plot

import warnings
warnings.filterwarnings('ignore')

DATA_PATH = Path(os.environ["HOT_DATA"])
BINS_DIR = DATA_PATH / "log_bins"

PLOTS_DIR = DATA_PATH / "plots/figure_1"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']


def plot_bins_loci_stats_2(save_file):

	tmp_file = tempfile.NamedTemporaryFile()
	columns = ["cell_line", "bin", "category", "value"]
	with open(tmp_file.name, "w") as outf:
		outf.write(",".join(columns) + "\n")

		for cell_line in ["HepG2", "K562"]:
			print(cell_line)
			prom_files = get_prom_files(cell_line)
			enh_files = get_enh_files(cell_line)
			for cat, files in zip(["promoters", "enhancers"], [prom_files, enh_files]):
				for bin, file in enumerate(files):
					bin_bed = BedTool(file)
					bin_count = bin_bed.count()
					out_line = "%s,%d,%s,%f\n" % (cell_line, bin, cat, bin_count)
					outf.write(out_line)
	df = pandas.read_csv(tmp_file.name)

	f = plt.figure(figsize=(11, 4))
	ax1 = f.add_subplot(121)
	ax2 = f.add_subplot(122)
	for ax, cell_line in zip([ax1, ax2], ["HepG2", "K562"]):
		x_tick_labels = HEPG2_XTICK_LABELS if cell_line == "HepG2" else K562_XTICK_LABELS
		x_range = np.arange(len(x_tick_labels))
		x_ticks = x_range - 0.5

		sns.barplot(x="bin", y="value", hue="category", data=df[df["cell_line"] == cell_line], ax=ax,
					hue_order=["enhancers", "promoters"], palette="colorblind")

		ax.set_xticks(x_ticks)
		ax.set_xticklabels(x_tick_labels, rotation=45)
		ax.set_xlabel("DAPs")
		ax.legend(loc="best")
		ax.set_title(cell_line)
		ax.set_yscale("log")
		ax.set_ylabel("# of loci")

	ax1.axvspan(9.5, 14, alpha=0.1, color="black")
	ax1.set_xlim([-1, 14])
	ax2.axvspan(9.5, 14, alpha=0.1, color="black")
	ax2.set_xlim([-1, 14])

	plt.tight_layout()
	plt.savefig(save_file)


def plot_enrichments_jointplot_size(save_file, cl="HepG2"):

	table_file = DATA_PATH / "src_files/tf_stats_summary_table.txt"
	df = pandas.read_table(table_file)
	df = df[df["cell_line"] == cl]

	cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
	g = sns.jointplot(
		data=df,
		x="all_fracs",
		y="pk_fracs_all",
		palette=cmap,
		height=6
	)

	g.ax_joint.cla()
	sc = sns.scatterplot(data=df, x="all_fracs", y="pk_fracs_all", hue="num_peaks", size="num_peaks", sizes=(10, 150),
						 ax=g.ax_joint)

	h, l = sc.get_legend_handles_labels()
	l = [_[:2] + "K" for _ in l]
	sc.legend(h, l, title="Total peaks")

	g.set_axis_labels("% of overlapping HOT loci", "% of ChIP-seq peaks in HOT loci")
	g.ax_joint.grid(alpha=0.4)
	g.ax_joint.set_xlim([-5, 100])
	g.ax_joint.set_ylim([0, 100])

	plt.savefig(save_file)


def hots_prom_enh_barplot(save_file, cl="HepG2"):

	introns_file = str(DATA_PATH / "src_files/hg19_files/knownGene.introns.merged.bed.gz")

	data = []
	hot_proms = BedTool(str(DATA_PATH / f"HOTs/{cl}_HOTs.proms.bed.gz"))
	hot_enhs = BedTool(str(DATA_PATH / f"HOTs/{cl}_HOTs.noproms.bed.gz"))
	cnt_proms = hot_proms.count()
	cnt_intronic = hot_enhs.intersect(introns_file, wa=True).groupby(c=4, o="collapse").count()
	cnt_intergenic = hot_enhs.intersect(introns_file, v=True, wa=True).groupby(c=4, o="collapse").count()
	data.append([cnt_proms, cnt_intronic, cnt_intergenic])

	data = np.array(data[0])
	data = 100 * (data / np.sum(data))

	labels = ["promoters", "intronic", "intergenic"]
	plt.figure(figsize=(3, 4))
	ax = plt.gca()

	bottoms = np.cumsum([0, data[0], data[1]])
	gs = [ax.bar([0], data[i], label=labels[i], bottom=bottoms[i]) for i in range(3)]
	[ax.bar_label(g, label_type="center", fmt='%d %%', fontweight="bold", color="white") for g in gs]

	ax.set_ylim([0, 100])
	ax.legend(bbox_to_anchor=(1.05, 0.5), ncol=1)
	ax.set_xticks([0])
	ax.set_xticklabels(["HOT loci"])
	ax.spines['right'].set_visible(False)

	plt.tight_layout()
	plt.savefig(save_file)


def load_ATAC_data(cl="HepG2"):
	# there's no ATAC-seq data for H1 cells in ENCODE
	_map = {"HepG2": "ENCFF439EIO", "K562": "ENCFF333TAT"}
	encode_id = _map[cl]

	hot_enhs = BedTool(join(DATA_PATH, f"HOTs/{cl}_HOTs.noproms.bed.gz"))
	hot_proms = BedTool(join(DATA_PATH, f"HOTs/{cl}_HOTs.proms.bed.gz"))

	data = []

	input_file = DATA_PATH / f"src_files/peaks/{encode_id}_hg19.bed.gz"
	ca_regions = BedTool(str(input_file))
	for reg, reg_name in zip([hot_proms, hot_enhs], ["HOT promoters", "HOT enhancers"]):
		cov = sum([r.length for r in reg.intersect(ca_regions).sort().merge()])
		total = sum([r.length for r in reg])
		fr_ca = cov / total

		data.append([reg_name, fr_ca])

	df = pandas.DataFrame(data=data, columns=["region", "value"])

	df["value"] = np.round(100 * df["value"])

	return df


def plot_ATAC(save_file, cl="HepG2"):

	df = load_ATAC_data(cl)

	plt.figure(figsize=(3, 3.1))
	ax = plt.gca()

	color = sns.color_palette("tab10")[0]

	sns.barplot(x="region", y="value", data=df, ax=ax, order=["HOT promoters", "HOT enhancers"], color=color)

	ax.set_xticks([0, 1])
	ax.set_xticklabels(["HOT\npromoters", "HOT\nenhancers"])
	ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
	ax.grid(axis="y", alpha=0.2)
	# ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False)

	ax.set_xlabel("")
	ax.set_ylabel("ATAC-seq")

	plt.tight_layout()
	plt.savefig(save_file, bbox_inches='tight')
	plt.close()


def draw_UpSet_SE_HOT_RE_HM(save_file):

	re_file = DATA_PATH/ "src_files/HepG2_enhancers_DHS_H3K27ac.bed.gz"
	se_file = DATA_PATH/ "src_files/HepG2_superenhancers.bed.gz"
	hot_enhs = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.noproms.bed.gz"))
	hot_proms = BedTool(join(DATA_PATH, "HOTs/HepG2_HOTs.proms.bed.gz"))

	k27ac = BedTool(join(DATA_PATH, "src_files/peaks/ENCFF392KDI_hg19.bed.gz"))
	k4me1 = BedTool(join(DATA_PATH, "src_files/peaks/ENCFF413EGR_hg19.bed.gz"))

	re = BedTool(re_file)
	se = BedTool(se_file).sort()

	data = from_memberships(
		[['HOT prom'],
		 ['HOT enh '],
		 ['Super-enh'],
		 ['Reg enh '],
		 ['HOT prom', 'H3K27ac '],
		 ['HOT prom', 'H3K4me1'],
		 ['HOT enh ', 'H3K27ac '],
		 ['HOT enh ', 'H3K4me1'],
		 ['HOT prom', 'Super-enh'],
		 ['HOT enh ', 'Super-enh'],
		 ['HOT enh ', 'Super-enh', 'Reg enh '],
		 ['Reg enh ', 'Super-enh'],
		 ],
		data=[
			hot_proms.count(),
			hot_enhs.count(),
			se.count(),
			re.count(),
			hot_proms.intersect(k27ac, wa=True).merge(d=-1).count(),
			hot_proms.intersect(k4me1, wa=True).merge(d=-1).count(),
			hot_enhs.intersect(k27ac, wa=True).merge(d=-1).count(),
			hot_enhs.intersect(k4me1, wa=True).merge(d=-1).count(),
			hot_proms.intersect(se, wa=True).merge(d=-1).count(),
			hot_enhs.intersect(se, wa=True).merge(d=-1).count(),
			hot_enhs.intersect(se, wa=True).merge(d=-1).intersect(re).sort().merge(d=-1).count(),
			re.intersect(se, wa=True).merge(d=-1).count(),
		]
	)

	fig = plt.figure()

	upset_plot(data, fig=fig)

	fig.savefig(save_file, bbox_inches='tight')
	return data


if __name__ == "__main__":

	save_file = join(PLOTS_DIR, "Figure1_a.pdf")
	print(save_file)
	plot_bins_loci_stats_2(save_file)

	for cl in ["HepG2", "K562"]:
		save_file = PLOTS_DIR / f"Figure1_b.{cl}.pdf"
		print(save_file)
		plot_enrichments_jointplot_size(save_file, cl)

	for cl in ["H1", "HepG2", "K562"]:
		save_file = join(PLOTS_DIR, f"Figure1_c.{cl}.pdf")
		print(save_file)
		hots_prom_enh_barplot(save_file, cl)

	for cl in ["HepG2", "K562"]:
		save_file = join(PLOTS_DIR, f"Figure1_d.{cl}.pdf")
		print(save_file)
		plot_ATAC(save_file, cl)

	save_file = join(PLOTS_DIR, "Figure1_e.pdf")
	print(save_file)
	draw_UpSet_SE_HOT_RE_HM(save_file)
