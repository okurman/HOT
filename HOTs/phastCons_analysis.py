#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import warnings
from collections import Counter, defaultdict

import pandas
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, wilcoxon

warnings.filterwarnings('ignore')

import sys
sys.path.append("../")
from os.path import join

import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt
from pathlib import Path

DATA_PATH = Path("../data/data")
PLOTS_DIR = DATA_PATH/"plots"
LOCI_DIR = DATA_PATH/"HOTs"
BINS_DIR = DATA_PATH/"log_bins"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

from data_prep.basic import load_metadata
from data_prep import phastCons


def extract_conservation_HOTs(score="phastCons", species="primates"):

	from lib import phastCons

	input_files = [join(LOCI_DIR, "HepG2_HOTs.proms.bed"), join(LOCI_DIR, "HepG2_HOTs.noproms.bed"),
				   join(LOCI_DIR, "K562_HOTs.proms.bed"), join(LOCI_DIR, "K562_HOTs.noproms.bed")]

	input_beds = [BedTool(_) for _ in input_files]
	save_files = [_.replace(".bed", ".%s.%s.bed.gz" % (score, species)) for _ in input_files]
	save_files = [gzip.open(_, "wt") for _ in save_files]

	chroms = ["chr%d" % _ for _ in range(1, 23)]
	chroms += ["chrX", "chrY"]

	for chrom in chroms[::-1]:
		print("chrom: " + chrom)
		pos2phylop = phastCons.load_scores(chrom=chrom, score=score, species=species)[chrom]

		for input_bed, of in zip(input_beds, save_files):
			chrom_bed = BedTool([r for r in input_bed if r.chrom == chrom]).sort()
			for r in chrom_bed:
				scores = [pos2phylop[i] for i in range(r.start, r.stop)]
				scores = np.asarray([_ for _ in scores if _ is not None])
				scores_str = ",".join([str(_) for _ in scores])
				out_line = "%s\t%d\t%d\t%f\t%s\n" % (r.chrom, r.start, r.stop, np.sum(scores) / r.length, scores_str)
				of.write(out_line)
				of.flush()


def extract_phastCons_RE(score="phastCons", species="primates"):

	print("Extracting %s for regular enhancers" % score, species)

	from lib import phastCons

	enh_dir = "ROOT_DIR/definitions/regular_enhancers/"
	save_dir = "ROOT_DIR/chipseq_files/"

	enh_files = [join(enh_dir, "HepG2_enhancers_DHS_H3K27ac.bed"), join(enh_dir, "K562_enhancers_DHS_H3K27ac.bed")]
	input_beds = [BedTool(_) for _ in enh_files]

	save_files = [join(save_dir, "HepG2_enhancers_DHS_H3K27ac.%s.%s.bed.gz" % (score, species)),
				join(save_dir, "K562_enhancers_DHS_H3K27ac.%s.%s.bed.gz" % (score, species))]
	save_files = [gzip.open(_, "wt") for _ in save_files]

	chroms = ["chr%d" % _ for _ in range(1, 23)]
	chroms += ["chrX", "chrY"]

	for chrom in chroms[::-1]:
		print("chrom: " + chrom)
		pos2phylop = phastCons.load_scores(chrom=chrom, score=score, species=species)[chrom]

		for input_bed, of in zip(input_beds, save_files):
			chrom_bed = BedTool([r for r in input_bed if r.chrom == chrom]).sort()
			for r in chrom_bed:
				scores = [pos2phylop[i] for i in range(r.start, r.stop)]
				scores = np.asarray([_ for _ in scores if _ is not None])
				scores_str = ",".join([str(_) for _ in scores])
				out_line = "%s\t%d\t%d\t%f\t%s\n" % (r.chrom, r.start, r.stop, np.sum(scores) / r.length, scores_str)
				of.write(out_line)
				of.flush()


def load_score_data_genes(species="primates"):

	hg19_dir = "/net/intdev/devdcode/common/genomes/hg19/annotations/genes/conservation_scores/"

	hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed")).merge(d=-1)
	rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True).merge(d=-1)

	if species == "vertebrates":

		exon_scores = BedTool(join(hg19_dir, "knownGene.exons.merged.phastCons.bed.gz"))
		intron_scores = BedTool(join(hg19_dir, "knownGene.introns.merged.phastCons.bed.gz"))
		prom_scores = BedTool(join(hg19_dir, "promoters.merged.phastCons.bed.gz")).intersect(rp, v=True)

	else:

		exon_scores = BedTool(join(hg19_dir, f"knownGene.exons.merged.phastCons.{species}.bed.gz"))
		intron_scores = BedTool(join(hg19_dir, f"knownGene.introns.merged.phastCons.{species}.bed.gz"))
		prom_scores = BedTool(join(hg19_dir, f"promoters.merged.phastCons.{species}.bed.gz")).intersect(rp, v=True)

	print("Extracting to data pool")
	data = []
	for name, scores_bed in zip(["exons", "introns", "regular promoters"],
								[exon_scores, intron_scores, prom_scores]):
		print(name)
		data += [[name, 0 if r.fields[3] == "-" else float(r.fields[3])] for r in scores_bed]

	df = pd.DataFrame(data=data, columns=["region", "value"])

	return df


def load_score_data_cell_types(species="primates"):

	hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed")).merge(d=-1)
	reg_enh_dir = "ROOT_DIR/chipseq_files/"
	if species == "vertebrates":
		re_scores = BedTool(join(reg_enh_dir, "HepG2_enhancers_DHS_H3K27ac.phastCons.bed.gz")).intersect(hot_enhs, v=True)
	else:
		re_scores = BedTool(join(reg_enh_dir, f"HepG2_enhancers_DHS_H3K27ac.phastCons.{species}.bed.gz")).intersect(hot_enhs, v=True)

	f1 = join(LOCI_DIR, "HepG2_HOTs.noproms.phastCons.bed.gz")
	f2 = join(LOCI_DIR, "HepG2_HOTs.proms.phastCons.bed.gz")
	hot_scores = BedTool(f1).cat(f2, postmerge=False)

	hot_enh_scores = hot_scores.intersect(hot_enhs, wa=True)
	hot_prom_scores = hot_scores.intersect(hot_enhs, wa=True, v=True)

	regions = [re_scores, hot_scores, hot_enh_scores, hot_prom_scores]
	names = ["regular enhancers", "HOT", "HOT enh", "HOT prom"]

	data = []
	for name, regions in zip(names, regions):
		data += [[name, float(r.fields[3])] for r in regions]

	df = pd.DataFrame(data=data, columns=["region", "value"])

	return df


def load_phastcons_data():

	from overbinders.log_bins import plots_data_factory

	line_parser = lambda x: x.split("\t")[3]
	file_amender = lambda x: x + ".phastcons"
	df = plots_data_factory.load_bins_to_df(line_parser=line_parser, file_amender=file_amender)
	hot_scores = df[(df["cell_line"] == "HepG2") & (df["bin"] > 10)].reset_index(drop=True)["value"].values

	f = "ROOT_DIR/chipseq_files/HepG2_regular_enhancers.phastCons.bed.gz"
	reg_enh_scores = np.loadtxt(f, usecols=(3,))
	ix = np.isnan(reg_enh_scores)
	reg_enh_scores = reg_enh_scores[~ix]

	f = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/pre/knownGene_exons.liftOver_hg19.phastCons.bed.gz"
	exon_scores = np.loadtxt(f, usecols=(3,))
	ix = np.isnan(exon_scores)
	exon_scores = exon_scores[~ix]

	data = []
	for cl, _data in zip(["HOT", "regular enhancers", "exons"], [hot_scores, reg_enh_scores, exon_scores]):
		for _ in _data:
			row = [cl, _]
			data.append(row)

	df = pd.DataFrame(data=data, columns=["cat", "value"])

	return df


def plot_mixed_plots(df_genes, df_cl):

	# df_genes, df_cl = load_score_data_genes(), load_score_data_cell_types()
	# return df_genes, df_cl

	fig, axs = plt.subplots(2, 1, figsize=(4, 4))

	cl_order = ['HOT enh', 'HOT prom', 'exons', 'introns', 'promoters', 'reg enh']

	_df_cl = df_cl[(df_cl["score"] == "phastCons") & (df_cl["cl"] == "HepG2")][["region", "score", "value"]]
	_df_genes = df_genes[df_genes["score"] == "phastCons"]
	_df = pd.concat([_df_genes, _df_cl])
	# sns.boxplot(data=_df, x="region", order=cl_order, y="value", ax=axs[0], showfliers=False)
	sns.barplot(data=_df, x="region", order=cl_order, y="value", ax=axs[0])

	for ax in [axs[0]]:
		ax.set_xlabel("")
		ax.set_xticklabels([])
		ax.set_ylabel("phastCons")
		ax.grid(axis="y", which="both", alpha=0.5)

	_df_cl = df_cl[(df_cl["score"] == "phyloP") & (df_cl["cl"] == "HepG2")][["region", "score", "value"]]
	_df_genes = df_genes[df_genes["score"] == "phyloP"]
	_df = pd.concat([_df_genes, _df_cl])
	# sns.boxplot(data=_df, x="region", order=cl_order, y="value", ax=axs[1], showfliers=False)
	sns.barplot(data=_df, x="region", order=cl_order, y="value", ax=axs[1])
	for ax in [axs[1]]:
		ax.set_xlabel("")
		ax.set_ylabel("phyloP")
		ax.grid(axis="y", which="both", alpha=0.5)

	axs[1].set_xticklabels(cl_order, rotation=45, ha="right")

	plt.tight_layout()

	# save_file = join(PLOTS_DIR, "conservation_scores_merged_boxplot.pdf")
	save_file = join(PLOTS_DIR, "conservation_scores_merged_barplot.pdf")
	print(save_file)
	plt.savefig(save_file)
	plt.close()


def plot_3panels_plots(species="primates"):
	
	df_genes, df_cl = load_score_data_genes(species=species), load_score_data_cell_types(species=species)
	# return df_genes, df_cl

	plt.figure(figsize=(3, 3))
	ax = plt.gca()

	hot_scores = df_cl[(df_cl["region"] == "HOT prom")]["value"]
	enh_scores = df_cl[(df_cl["region"] == "regular enhancers")]["value"]
	exon_scores = df_genes[(df_genes["region"] == "exons")]["value"]

	print(mannwhitneyu(hot_scores, enh_scores))
	print(mannwhitneyu(exon_scores, enh_scores))

	data = [["reg enh", 1],
			["HOT", hot_scores.mean()/enh_scores.mean()],
			["exons", exon_scores.mean()/enh_scores.mean()]]

	print(data)

	df = pd.DataFrame(columns=["region", "value"], data=data)

	sns.barplot(data=df, x="region", y="value", ax=ax)
	ax.grid(axis="y", alpha=0.3)
	ax.set_ylabel(f"phastCons {species}")
	ax.set_xlabel("")
	ax.set_ylim([0, 5])

	plt.tight_layout()
	ax.set_xticklabels(["regular\nenhancers", "HOT", "exons"])
	save_file = join(PLOTS_DIR, f"phastCons_3panels_HOTprom.{species}.pdf")
	print(save_file)
	plt.savefig(save_file)
	plt.close()


def plot_conservation_scores(df):

	# df = load_phastcons_data()

	reg_enh = df[df["cat"] == "regular enhancers"]["value"].values
	exon = df[df["cat"] == "exons"]["value"].values
	hot = df[df["cat"] == "HOT"]["value"].values

	data = [["HOT", np.mean(hot)/np.mean(reg_enh)],
			["exons", np.mean(exon)/np.mean(reg_enh)],
			["regular enhancers", 1]]
	plot_df = pd.DataFrame(data=data, columns=["cat", "value"])

	print(mannwhitneyu(reg_enh, hot))
	print(mannwhitneyu(reg_enh, exon))

	plt.figure(figsize=(3.5, 3))
	ax = plt.gca()

	sns.barplot(data=plot_df, order=["regular enhancers", "exons", "HOT"], x="cat", y="value", ax=ax)
	ax.set_ylim([0.5, 10**3])
	ax.set_yscale("log")
	ax.set_xlabel("")
	ax.set_ylabel("phastCons score ratio")
	ax.set_xticklabels(["regular\nenhancers", "exons", "HOT"])

	ax.bar_label(ax.containers[0], fmt='x%d')
	ax.grid(axis="y", alpha=0.5)

	plt.tight_layout()
	save_file = join(PLOTS_DIR, "phastCons.barplot.pdf")
	plt.savefig(save_file, bbox_inches='tight')
	plt.close()


if __name__ == "__main__":

	if sys.argv[1] == "1":
		print("HOTs", sys.argv[2])
		extract_conservation_HOTs(species=sys.argv[2])
	if sys.argv[1] == "2":
		extract_phastCons_RE(species=sys.argv[2])
	# extract_phastCons_RE_RP()
