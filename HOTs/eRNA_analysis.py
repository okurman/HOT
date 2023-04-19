
import warnings
from collections import defaultdict, Counter
import h5py
warnings.filterwarnings('ignore')

import sys
import os
from os.path import join

import numpy as np
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
from lib import phastCons, phyloP
import seaborn as sns
from upsetplot import from_memberships
from upsetplot import plot as upset_plot


from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

eRNA_DIR = "ROOT_DIR/eRNAs/"


def extract_enrichments():

    # re_file = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers.bed.gz"
    # dhs = BedTool(join(LOCI_DIR, "variant_analysis/HepG2_DHS.bed"))
    # hots_file = join(LOCI_DIR, "HepG2_HOTs.PE.bed")

    re_file = "ROOT_DIR/definitions/regular_enhancers/K562_enhancers.bed.gz"
    dhs = BedTool("ROOT_DIR/chipseq_files/DHS_K562.bed.gz")
    hots_file = join(LOCI_DIR, "K562_HOTs.PE.bed")

    re = BedTool(re_file)
    re_total = np.sum([r.length for r in re])

    dhs_total = np.sum([r.length for r in dhs])

    hots = BedTool(hots_file)
    hot_proms = BedTool([r for r in hots if r.fields[-1] == "P"])
    prom_total = np.sum([r.length for r in hot_proms])
    hot_enhs = BedTool([r for r in hots if r.fields[-1] == "E"])
    enh_total = np.sum([r.length for r in hot_enhs])

    # eRNAs = BedTool("ROOT_DIR/eRNAs/CAGE/HepG2_ENCFF933JJT.bed.gz")
    # eRNAs = BedTool("ROOT_DIR/eRNAs/CAGE/ENCFF318EPW.bed.gz")
    # eRNAs = BedTool("ROOT_DIR/eRNAs/FANTOM/HepG2_eRNAs.bed")

    eRNAs = BedTool("ROOT_DIR/eRNAs/PRO-cap/ENCFF902VME.hg19.bed")

    dhs_eRNAs = dhs.intersect(eRNAs, wa=True)
    # dhs_eRNAs_count = sum([r.length for r in dhs_eRNAs])
    dhs_eRNAs_count = dhs_eRNAs.count()

    re_eRNAs = re.intersect(eRNAs, wa=True)
    # re_eRNAs_count = sum([r.length for r in re_eRNAs])
    re_eRNAs_count = re_eRNAs.count()

    enh_eRNAs = hot_enhs.intersect(eRNAs, wa=True)
    # enh_eRNAs_count = sum([r.length for r in enh_eRNAs])
    enh_eRNAs_count = enh_eRNAs.count()

    prom_eRNAs = hot_enhs.intersect(eRNAs, wa=True)
    # prom_eRNAs_count = sum([r.length for r in prom_eRNAs])
    prom_eRNAs_count = prom_eRNAs.count()

    dhs_density = dhs_eRNAs_count / dhs.count()
    re_density = re_eRNAs_count / re.count()
    enh_density = enh_eRNAs_count / hot_enhs.count()
    prom_density = prom_eRNAs_count / hot_proms.count()

    print(re_density/dhs_density, enh_density/dhs_density, prom_density/dhs_density)











if __name__ == "__main__":

    pass






