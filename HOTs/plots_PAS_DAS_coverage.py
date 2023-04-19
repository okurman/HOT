
import warnings
from collections import Counter, defaultdict

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

INTRONS_FILE ="ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownGene_introns.liftOver_hg19.bed"
PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed"

PROJECT_DIR = "ROOT_DIR/"
BINS_DIR = "ROOT_DIR/definitions/peak_8bp_v2/log_bins/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom" % (x, i)) for i in range(14)]

K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"


def extract_PAS_coverage():

    d = "ROOT_DIR/definitions/peak_8bp_v2/PAS_DAS/"

    for cell_line in ["HepG2", "K562"]:
        print(cell_line)

        hots_loci = BedTool(join(LOCI_DIR, "%s_HOTs.PE.bed" % cell_line))

        pas_loci = BedTool(join(d, "%s_PAS_DAS.bed" % cell_line))
        motifs_loci = BedTool(join(PROJECT_DIR, "chipseq_files/%s_motif_concat_fimo.bed.gz" % cell_line)).intersect(
            hots_loci).sort().merge()

        total_hot_area = np.sum([r.length for r in hots_loci])
        empty_hot_area = np.sum([r.length for r in hots_loci.subtract(pas_loci).subtract(motifs_loci)])

        hot_pas_area = np.sum([r.length for r in hots_loci.intersect(pas_loci).sort().merge()])
        hot_motif_area = np.sum([r.length for r in motifs_loci])

        pas_ovp_motifs_area = sum([r.length for r in pas_loci.intersect(motifs_loci).sort().merge()])

        print([total_hot_area, empty_hot_area, hot_pas_area, hot_motif_area, pas_ovp_motifs_area])


def plot_venn():

    from matplotlib_venn import venn3

    # (100, 010, 110, 001, 101, 011, 111)

    plt.figure()

    values = (19,  0, 50,  0, 6.1,  0, 24)
    v = venn3(subsets=(values), set_labels=["HOTs", "PAS/DAS", "motifs"], set_colors=["r", "g", "b"])
    # for _id in ['100', '110', '010', '101', '111', '011', '001']:
    #     v.get_label_by_id(_id).set_text("")
    for idx, subset in enumerate(v.subset_labels):
        try:
            v.subset_labels[idx].set_visible(False)
        except:
            continue

    save_file = join(PLOTS_DIR, "Venn_HepG2_PAS_DAS_motifs.pdf")
    save_file = join(PLOTS_DIR, "Venn_HepG2_PAS_DAS_motifs_nolabels.pdf")
    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)


