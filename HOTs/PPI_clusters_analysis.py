#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gzip
import time
import warnings
from collections import Counter, defaultdict
from random import shuffle

import pandas
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, wilcoxon, ttest_1samp

warnings.filterwarnings('ignore')

import sys
import os
from os.path import join

import numpy as np
import configparser

###############################################################

config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
import matplotlib.pyplot as plt

PROJECT_DIR = "ROOT_DIR/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

from overbinders.data_prep.basic import load_metadata
import requests
import json

cl_list = [
	["ZNF687", "ARID4B", "MAX", "SAP130"],
	["KMT2A", "E2F4", "NONO", "HNRNPLL", "RBM39", "POLR2A", "POLR2AphosphoS5", "RBFOX2", "ARID4A", "TAF1", "ZFY",
	 "GABPB1", "PHF8", "POLR2G", "NR2C2", "DRAP1", "YEATS4", "HMGXB4", "KMT2B", "TFDP2", "MAZ", "SPEN", "ASH2L",
	 "KDM2A", "MNX1", "UBTF", "GATAD1", "ZNF501", "DMAP1"],
	["NR2F6", "TEAD1", "SOX5", "NFIL3", "HNF4A", "FOXA2", "PPARG", "MIXL1", "FOXA1", "CEBPG", "FOXO1", "ZNF217",
	 "CEBPA", "KDM1A", "FOXP1", "RCOR2"],
	["ZGPAT", "MED1", "TFAP4", "ZFX", "EGR1", "LIN54", "ZNF574", "HDAC1", "TBX2", "THAP11", "KAT8", "HOXA3"],
	["KLF16", "PATZ1", "ERF", "ZNF331", "LCORL", "IRF2", "SKI", "ISL2", "ZBTB7B", "POGZ", "IKZF5", "HNF1B", "FOSL2",
	 "TCF7L2", "RXRB", "RARA", "LCOR", "FOXP4", "BCL6", "GATAD2A", "SOX6", "PAXIP1", "ELF3", "SMAD4", "HDAC2", "FOXA3"]
]

string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "ppi_enrichment"

all_tfs = list(load_metadata()["HepG2"].keys())

from pathlib import Path


def get_enrichment(cl_no=1):

	save_dir = Path(os.path.join(PLOTS_DIR, f"PPI_enrichments/{cl_no}"))
	save_dir.mkdir(exist_ok=True)

	params = {
		"identifiers": "%0d".join([_ for _ in cl_list[cl_no - 1]]),
		"species": 9606,  # species NCBI identifier
		"caller_identity": "hudaiber@nih.gov"
	}

	request_url = "/".join([string_api_url, output_format, method])

	response = requests.post(request_url, data=params)

	save_file = save_dir/"main_PPI.txt"
	print(save_file)
	with open(save_file, "w") as of:
		of.write(response.text)

	for i in range(100):

		shuffle(all_tfs)

		rand_tfs = all_tfs[:len(cl_list[cl_no-1])]

		params["identifiers"] = "%0d".join(rand_tfs)
		request_url = "/".join([string_api_url, output_format, method])
		response = requests.post(request_url, data=params)

		save_file = save_dir/f"rand_{i+1}_PPI.txt"
		print(save_file)
		with open(save_file, "w") as of:
			of.write(response.text)

		time.sleep(2)


def get_enrichment_stats():

	src_dir = Path(os.path.join(PLOTS_DIR, "PPI_enrichments/"))

	for i in range(1, 5):

		main_pvalue = float(open(src_dir/f"{i}/main_PPI.txt").readlines()[1].split()[-1])
		rand_pvalues = [float(open(f).readlines()[-1].split()[-1]) for f in (src_dir/f"{i}").glob("rand_*_PPI.txt")]

		print(i, main_pvalue, ttest_1samp(rand_pvalues, main_pvalue))

		# print(main_pvalue, rand_pvalues)


