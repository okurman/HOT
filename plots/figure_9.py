#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys

CODE_PATH = os.environ["HOT_CODE"]
if CODE_PATH not in sys.path:
	sys.path.append(CODE_PATH)

import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
from os.path import join
from scipy.stats import hypergeom, fisher_exact, binomtest, ttest_ind

DATA_PATH = Path(os.environ["HOT_DATA"])
PLOTS_DIR = DATA_PATH / "plots/figure_9"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)

from data_prep.basic import load_metadata
a = load_metadata()


RBP_tfs = ["AGO1", "AGO2", "CCAR2", "CHD2", "FIP1L1", "FUBP1", "GTF3A", "HNRNPH1", "HNRNPL", "HNRNPLL", "HNRNPUL1", 
			"JRK", "MATR3", "MBNL3", "NCOA5", "NONO", "PCBP1", "PCBP2", "PHF5A", "PRPF4", "PTBP1", "RBM22", "RBM39", 
			"REPIN1", "SAFB2", "SMYD3", "SNRNP70", "SRSF1", "SRSF4", "SRSF9", "SSRP1", "SYNCRIP", "TAF15", "TARDBP", 
			"TOE1", "U2AF1", "U2AF2", "ZC3H13", "ZC3H4", "ZCCHC11", "ZFP36", "ZFP36L1", "ZFP36L2", "ZMAT3", "ZNF207", 
			"ZNF326", "ZNF598", "ZNF768"]


def load_signal_values(score_type="raw"):

	file = join(CODE_PATH, "plots/data_files/HepG2_ChIP_signal_values.txt.gz")
	df = pd.read_csv(file, sep="\t")

	df = df[df["score_type"]==score_type]
	df = df[["cl", "tf", "hot", "mean", "median"]].set_index(["cl", "tf", "hot"]).unstack("hot")
	df.columns = [f"{t[1]}_{t[0]}" for t in df.columns.to_flat_index()]

	return df


def parse_CDCODE():

	tfs = set(a["HepG2"].keys())

	src_file = DATA_PATH/"src_files/LLPS_proteins/CD_CODE.csv"

	df = pd.read_csv(src_file)
	df = df[df["Species"] == "Homo sapiens"].reset_index(drop=True)

	df = df[(df["Biomolecular condensate count"] > 0) | (df["Synthetic condensate count"] > 0)]
	
	return set(df["Gene Name"].values)


def load_enr_dataset(cl="HepG2"):

	file = join(CODE_PATH, "plots/data_files/DAP_enrichment_values.txt.gz")
	df = pd.read_csv(file, sep="\t")
	df = df[df["cl"]==cl]

	df["pval_logneg"] = df["pval_fisher"].apply(lambda x: -1*np.log10(x) if x>0 else 340)

	df = df.set_index(["cl", "tf", "bkg"])
	df = df.unstack("bkg")
	df.columns = [f"{t[1]}_{t[0]}" for t in df.columns.to_flat_index()]

	table_file = DATA_PATH/"src_files/tf_stats_summary_table.txt"
	df_pre = pd.read_csv(table_file, sep="\t")
	df_pre = df_pre.rename(columns={"cell_line": "cl"}).set_index(["cl", "tf"])[["all_fracs", "pk_fracs_all", "num_peaks"]]
	df = df.merge(df_pre, left_index=True, right_index=True).reset_index()
	
	llps_tfs = parse_CDCODE()
	df["LLPS"] = df["tf"].apply(lambda x: True if x in llps_tfs else False)
	df["RBP"] = df["tf"].apply(lambda x: True if x in RBP_tfs else False)

	df_signals = load_signal_values(score_type="std")

	df = df.set_index(["cl", "tf"])
	df = df.merge(df_signals, left_index=True, right_index=True).reset_index()

	df["s_frac_mean"] = np.log2(df["1_mean"]/df["0_mean"])
	df["s_frac_median"] = np.log2(df["1_median"]/df["0_median"])

	return df


def plot_piechart_all_LLPS():

	df = load_enr_dataset()

	plt.close("all")
	save_file = join(PLOTS_DIR, "Figure9_a.pdf")
	print(save_file)

	plt.figure(figsize=(2, 2))
	ax = plt.gca()

	nums = [df.shape[0] - df["LLPS"].sum(), df["LLPS"].sum()]
	labels = ["non-LLPS", "LLPS"]
	ax.pie(nums, labels=labels, autopct='%1.0f%%')

	plt.tight_layout()
	plt.savefig(save_file, bbox_inches="tight")


def get_perc_LLPS():

	df = load_enr_dataset()
	
	bin_col = "all_fracs"
	df = df[["LLPS", bin_col]]

	bins = np.arange(0, 100, 15)
	labels = bins[1:]

	df["bins"] = pd.cut(df[bin_col], bins=bins, labels=labels)
	df["bins"] = df["bins"]

	df = df.groupby("bins").agg({"LLPS": ["sum", "count"]}).droplevel(axis=1, level=0).reset_index()
	df["perc"] = 100 * df["sum"]/df["count"]

	return df


def plot_perc_binned_LLPS_barplot():

	plt.close("all")
	save_file = join(PLOTS_DIR, "Figure9_b.pdf")
	print(save_file)
	
	df = get_perc_LLPS()
	df = df[df["perc"]>0]

	fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3.9, 2.4))

	x_range = np.arange(df.shape[0])
	p = axs[0].bar(x_range, df["count"], width=0.9, color=sns.color_palette()[0])
	axs[0].set_ylim([1, 10**3/3])
	axs[0].set_yscale("log")
	axs[0].bar_label(p, label_type='center')
	axs[0].set_ylabel("DAPs")
	axs[0].grid(alpha=0.4, axis="y")

	perc_vals = df["perc"].values

	p = axs[1].bar(x_range, perc_vals, label="LLPS", color=sns.color_palette()[1], width=0.9)
	perc_fmt = lambda x: f"{int(np.round(x))}%"
	axs[1].bar_label(p, label_type='center', fmt=perc_fmt)
	axs[1].set_ylabel(f"LLPS (%)")
	axs[1].grid(alpha=0.4, axis="y")
	
	x_ticks = np.asarray(list(x_range) + [x_range[-1] + 1]) - 0.5
	_fmt = lambda x: f"{x}%"
	x_tick_labels = map(_fmt, [0] + list(df["bins"].values))
	axs[1].set_xticks(x_ticks)
	axs[1].set_xticklabels(x_tick_labels)
	axs[1].set_xlabel(f"% of HOT loci a DAP is present")

	plt.tight_layout()
	plt.savefig(save_file)


def plot_signal_distr(df=None):
	
	plt.close("all")

	if df is None:
		df = load_enr_dataset()
	
	save_file = join(PLOTS_DIR, "Figure9_c.pdf")
	print(save_file)
	
	plt.figure(figsize=(1.8, 2))
	ax = plt.gca()

	sns.boxenplot(data=df, x="LLPS", order=[True, False], y="1_mean", hue="LLPS", legend=False, ax=ax, showfliers=False)
	
	_df = df.groupby("LLPS")["1_mean"].mean().reset_index()
	sns.stripplot(data=_df, x="LLPS", order=[True, False], y="1_mean", zorder=10, ax=ax,
		color='red', linewidth=0.5, jitter=False, edgecolor='lightgray')

	ax.set_ylim([0, 4.5])
	ax.set_xticklabels(["LLPS", "non-LLPS"])
	ax.set_ylabel("ChIP-seq signal")
	ax.set_xlabel("")

	ax.grid(alpha=0.3, axis="y")

	y = 4
	h = 0.1

	llps_sig_vals = df[df["LLPS"]==True]["1_mean"].values
	nonllps_sig_vals = df[df["LLPS"]==False]["1_mean"].values
	pval = ttest_ind(llps_sig_vals, nonllps_sig_vals).pvalue
	
	print(_df)
	print(pval)

	ax.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, c="k")
	# ax.text(0.5, y+h+0.05, f"pvalue={np.round(pval, 3)}", ha='center', va='bottom', color="k")
	ax.text(0.5, y+h-0.1, f"***", ha='center', va='bottom', color="k")
	
	plt.tight_layout()
	
	plt.savefig(save_file, bbox_inches="tight")
	plt.close()


def load_mobidb():

	file = join(DATA_PATH, "src_files/HUMAN_9606_idmapping.GeneCards.dat")
	tf2uniprot = {l.split("\t")[2].strip():l.split("\t")[0] for l in open(file)}

	# MobiDB-lite
	file = join(DATA_PATH, "src_files/LLPS_proteins/mobidb_mobidbLite.tsv")
	df = pd.read_csv(file, sep='\t')
	uniprot2mobidb = {}
	for acc, perc in df[["acc", "content_fraction"]].values:
		uniprot2mobidb[acc] = perc

	tfs = a["HepG2"].keys()

	tf2alias = {
		"KIAA2018": "USF3",
		"ZCCHC11": "TUT4",
		"ZNF788": "ZNF788P",
		"ZNF720": "KRBOX5",
		"POLR2AphosphoS2": "POLR2A",
		"ARNTL": "BMAL1",
		"RNF219": "OBI1",
		"POLR2AphosphoS5": "POLR2A"}

	tf2mobidb = {}
	for tf in tfs:
		uniprot_id = tf2uniprot[tf] if tf in tf2uniprot else tf2uniprot[tf2alias[tf]]
		
		tf2mobidb[tf] = uniprot2mobidb[uniprot_id] if uniprot_id in uniprot2mobidb else np.NaN
		
	df = load_enr_dataset()
	df = df[["tf", "LLPS"]]

	df["mobidb"] = np.asarray([tf2mobidb[_] for _ in df["tf"].values])

	return df


def plot_signal_IDR_distr(df=None):
	
	plt.close("all")
	
	if df is None:
		df_enr = load_enr_dataset()
		df_enr = df_enr.set_index("tf")
		df_dis = load_mobidb().drop("LLPS", axis=1)
		df_dis = df_dis.set_index("tf")

		df = df_enr.merge(df_dis, left_index=True, right_index=True)

		df["mobidb"] = 100 * df["mobidb"]
	

	save_file = join(PLOTS_DIR, "Figure9_d.pdf")
	print(save_file)
	
	plt.figure(figsize=(2.1, 2))
	ax = plt.gca()

	# idr_col = "disprot"
	idr_col = "mobidb"

	df = df[~df[idr_col].isnull()]

	# sns.boxenplot(data=df, x="LLPS", order=[True, False], y=idr_col, hue="LLPS", legend=False, ax=ax, showfliers=False)
	sns.boxenplot(data=df, x="LLPS", order=[True, False], y=idr_col, hue="LLPS", legend=False, ax=ax, showfliers=False)

	_df = df.groupby("LLPS")[idr_col].mean().reset_index()
	sns.stripplot(data=_df, x="LLPS", order=[True, False], y=idr_col, zorder=10, ax=ax,
		color='red', linewidth=0.5, jitter=False, edgecolor='lightgray')
	ax.grid(alpha=0.3, axis="y")

	ax.set_ylim([0, 77])
	ax.set_xticks([0, 1])
	ax.set_xticklabels(["LLPS", "non-LLPS"])
	ax.set_ylabel(f"% of IDRs in DAPs")
	ax.set_xlabel("")

	y = 70
	h = 0.2

	llps_idr_vals = df[df["LLPS"]==True][idr_col].values
	nonllps_idr_vals = df[df["LLPS"]==False][idr_col].values
	pval = ttest_ind(llps_idr_vals, nonllps_idr_vals).pvalue
	print(pval)
	ax.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, c="k")
	# ax.text(0.5, y+h+0.05, f"pvalue={np.round(pval, 2)}", ha='center', va='bottom', color="k")
	ax.text(0.5, y-2.5, f"**", ha='center', va='bottom', color="k")

	plt.tight_layout()
	
	plt.savefig(save_file, bbox_inches="tight")
	plt.close()


def plot_RBP_enr():

	plt.close("all")
	df = load_enr_dataset()

	save_file = join(PLOTS_DIR, "Figure9_e.pdf")
	print(save_file)
	
	plt.figure(figsize=(1.8, 2))
	ax = plt.gca()

	sns.boxenplot(data=df, x="RBP", order=[True, False], y="DHS_cnt_enr", ax=ax, hue="RBP",
		legend=False, hue_order=[True, False], showfliers=False, palette="Set3")
	
	_df = df.groupby("RBP")["DHS_cnt_enr"].mean().reset_index()
	sns.stripplot(data=_df, x="RBP", order=[True, False], y="DHS_cnt_enr", zorder=10, ax=ax,
		color='red', linewidth=0.5, jitter=False, edgecolor='lightgray')

	ax.set_ylim([0, 2.6])
	ax.set_xticklabels(["RBPs", "DAPs"])
	ax.set_ylabel("log2(FC)")
	ax.set_xlabel("")

	ax.grid(alpha=0.3, axis="y")

	y = 2.3
	h = 0.05

	v1 = df[(df["RBP"]==True) & (df["LLPS"]==True)]["DHS_cnt_enr"].values
	v2 = df[(df["RBP"]==False) & (df["LLPS"]==True)]["DHS_cnt_enr"].values
	pval = ttest_ind(v1, v2).pvalue
	
	print(_df)
	print(pval)

	ax.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, c="k")
	ax.text(0.5, y+h-0.07, "**", ha='center', va='bottom', color="k")

	plt.tight_layout()
	
	plt.savefig(save_file, bbox_inches="tight")
	plt.close()


if __name__ == "__main__":

	plot_piechart_all_LLPS()
	
	plot_perc_binned_LLPS_barplot()

	plot_signal_distr()

	plot_signal_IDR_distr()

	plot_RBP_enr()
