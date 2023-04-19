import warnings
warnings.filterwarnings('ignore')
import sys
import os
from os.path import join
import configparser

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
import pandas
import matplotlib.pyplot as plt
from overbinders.data_prep.basic import load_metadata
from matplotlib_venn import venn3
import numpy as np

PROJECT_DIR = ROOT_DIR

LOCI_DIR = PROJECT_DIR+ "/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = PROJECT_DIR+ "/plots/HOTs/"


def extract_remarker_HOTs():

	loci_file = ROOT_DIR + "/chipseq_files/loci_v2/H1_400_loci.bed"
	save_file = ROOT_DIR +  "/alternate_definitions/remaker/H1.bed"
	total_tfs = 47

	with open(loci_file, "r") as inf, open(save_file, "w") as of:

		for l in inf:

			parts = l.split("\t")

			try:

				if int(parts[5]) > (total_tfs/5):
					of.write(l)

			except:
				continue


def exract_overlaps():

	# hot_file = os.path.join(LOCI_DIR, "HepG2_HOTs.bed")
	# remaker_file = os.path.join(PROJECT_DIR, "alternate_definitions/remaker/HepG2.bed")
	# wr_file = os.path.join(PROJECT_DIR, "alternate_definitions/wreckzycka/maphot_hs_selection_reg_hg_simP05_hot.shuf.bed")
	# save_file = os.path.join(PROJECT_DIR, "alternate_definitions/Venn_HepG2.pdf")

	hot_file = os.path.join(LOCI_DIR, "K562_HOTs.bed")
	remaker_file = os.path.join(PROJECT_DIR, "alternate_definitions/remaker/K562.bed")
	wr_file = os.path.join(PROJECT_DIR,
						   "alternate_definitions/wreckzycka/maphot_hs_selection_reg_k5_simP05_hot.shuf.bed")
	save_file = os.path.join(PROJECT_DIR, "alternate_definitions/Venn_K562.pdf")

	files = [remaker_file, wr_file]

	beds = [BedTool("\n".join(["\t".join(l.split("\t")[:3]) for l in open(hot_file)]), from_string=True).sort()] + \
		   [BedTool(_).sort() for _ in files]

	subsets = []

	n = sum([r.length for r in (beds[0].subtract((beds[1].cat(beds[2])))).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[1].subtract((beds[2].cat(beds[0])))).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[0].intersect(beds[1]).subtract(beds[2])).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[2].subtract(beds[1].cat(beds[0]))).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[0].intersect(beds[2]).subtract(beds[1])).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[1].intersect(beds[2]).subtract(beds[0])).sort().merge()])
	subsets.append(n)

	n = sum([r.length for r in (beds[1].intersect(beds[2]).intersect(beds[0])).sort().merge()])
	subsets.append(n)

	plt.close("all")
	plt.figure(figsize=(4, 4))

	v = venn3(subsets=subsets, set_labels=("HOT", "Remaker et al.", "Wreckzycka et al."))
	for _id in ['100', '110', '010', '101', '111', '011', '001']:
		v.get_label_by_id(_id).set_text("")

	total_hots = sum([r.length for r in beds[0].sort().merge()])

	percs = np.asarray([subsets[0], subsets[2], subsets[6], subsets[4]])/total_hots
	percs = np.round(100*percs)
	print(sum(percs))
	percs = ["%d%%" % _ for _ in percs]

	# v.get_label_by_id('100').set_text(perc_1)

	for _id, perc in zip(['100', '110', '101', '111'], percs):
		v.get_label_by_id(_id).set_text(perc)

	print(save_file)
	plt.savefig(save_file)




