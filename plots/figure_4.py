#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append("../")

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

import warnings
warnings.filterwarnings('ignore')

import plots_data_factory
from HOTs import TF_classes_by_signal


DATA_PATH = Path("../data/data")
BINS_DIR = DATA_PATH / "log_bins"
PLOTS_DIR = DATA_PATH / "plots/figure_4"
PLOTS_DIR.mkdir(exist_ok=True)

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']


def plot_chipseq_signals(save_file):

	plt.close("all")
	cl = "HepG2"

	line_parser = lambda x: np.mean([float(_) for _ in x.rstrip().split("\t")[-1].split(",")])
	df = plots_data_factory.load_bins_to_df(line_parser=line_parser)

	plots_data_factory.plot_categorical_single_cell_line(df,
									cell_line=cl,
									y_label="Signal value",
									plot_type="boxplot",
									save_file=save_file,
									leg_loc="upper left",
									figsize=(4.3, 2.8))


if __name__ == "__main__":

	save_file = PLOTS_DIR / "Figure4_a.pdf"
	print(save_file)
	plot_chipseq_signals(save_file)

	save_file = PLOTS_DIR / "Figure4_b_c.pdf"
	print(save_file)
	TF_classes_by_signal.plot_tf_signals_HOTs_4panels(save_file)

	save_file = PLOTS_DIR / "Figure4_d.pdf"
	print(save_file)
	TF_classes_by_signal.plot_signals(save_file)
