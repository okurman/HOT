#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings

warnings.filterwarnings('ignore')

import pandas
import seaborn
import sys

sys.path.append("../")

import os
from os.path import join
import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt
from pathlib import Path
import networkx
import gzip

DATA_PATH = Path("../data/data/")
PLOTS_DIR = DATA_PATH / "plots"
LOCI_DIR = DATA_PATH / "HOTs"
BINS_DIR = DATA_PATH / "log_bins"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]


def extract_contact_lists_encode():

	cl = "HepG2"

	# FitHiChIP file was generated from
	contact_file = DATA_PATH/"src_files/FitHiChIP.interactions_FitHiC_Q0.0001.hg19.bed.gz"

	contacts = BedTool(contact_file)
	reverse_contacts = []
	for r in contacts:
		parts = r.fields
		new_parts = [parts[0]] + parts[4:6] + [parts[0]] + parts[1:3] + parts[6:]
		reverse_contacts.append(new_parts)
	reverse_contacts = BedTool(reverse_contacts)
	contacts = contacts.cat(reverse_contacts, postmerge=False)

	files = get_loci_files(cl)

	counts = np.zeros((len(files), len(files)))
	freqs = np.zeros((len(files), len(files)))

	locus_counts = [BedTool(f).count() for f in files]

	for i in range(len(files)):
		print(f"Bin no: {i}")

		l1_bed = BedTool(files[i])
		l1_contacts = set(["\t".join(r.fields[10:16]) for r in l1_bed.intersect(contacts, wo=True)])

		target_bed = BedTool("\n".join(["\t".join(r.split("\t")[3:]) for r in l1_contacts]),
							 from_string=True)

		for j in range(i, len(files)):
			l2_bed = BedTool(files[j])
			count = target_bed.intersect(l2_bed, wo=True).sort().merge().count()
			counts[i, j] = count
			freqs[i, j] = np.log10(count / (locus_counts[i] * locus_counts[j]))

	return counts, freqs


def load_graph():

	loops_file = DATA_PATH/"src_files/ENCFF050EKS.bedpe.hg19.gz"

	G = networkx.Graph()

	for line in gzip.open(loops_file, "rt"):

		if line.startswith("#"):
			continue

		parts = line.rstrip().split("\t")

		n1 = "%s-%s-%s" % (parts[0], parts[1], parts[2])
		n2 = "%s-%s-%s" % (parts[3], parts[4], parts[5])

		G.add_edge(n1, n2)

	return G



if __name__ == "__main__":
	pass
