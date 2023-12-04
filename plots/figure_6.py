#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use("Agg")

from pathlib import Path
import os
DATA_PATH = Path(os.environ["HOT_DATA"])
PLOTS_DIR = DATA_PATH / "plots/figure_6"
PLOTS_DIR.mkdir(exist_ok=True, parents=True)

from HOTs import HK_genes_analysis


if __name__ == "__main__":

	save_file = PLOTS_DIR/"Figure6_a.pdf"
	print(save_file)
	HK_genes_analysis.plot_fractions_merged(save_file)

	save_file = PLOTS_DIR/"Figure6_b.pdf"
	print(save_file)
	HK_genes_analysis.plot_tau_single_cl(save_file)
