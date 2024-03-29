#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
from HOTs import PCA_analysis
from HOTs import TF_coclusters
import warnings
warnings.filterwarnings('ignore')
import os

DATA_PATH = Path(os.environ["HOT_DATA"])
PLOTS_DIR = DATA_PATH / "plots/figure_2"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)


def plot_PCA_figures():

	pca_matrix, info_mat, explained_variances = PCA_analysis.load_PCA()

	data = [pca_matrix, info_mat]

	save_file = PLOTS_DIR/"Figure2_a.pdf"
	PCA_analysis.plot_pca_PE(data=data, save_file=save_file)

	save_file = PLOTS_DIR / "Figure2_b.pdf"
	PCA_analysis.plot_pca_param(data=data, tf="p300", save_file=save_file)

	save_file = PLOTS_DIR / "Figure2_c.pdf"
	PCA_analysis.plot_pca_param(data=data, tf="ctcf", y_ind=3, save_file=save_file)


def plot_TF_clusters():

	save_file = PLOTS_DIR / "Figure2_d_full.pdf"
	TF_coclusters.clusterplot_by_percentage_full(save_file)

	save_file = PLOTS_DIR / "Figure2_d.pdf"
	TF_coclusters.clusterplot_by_percentage_essentials(save_file)


if __name__ == "__main__":

	plot_PCA_figures()

	plot_TF_clusters()
