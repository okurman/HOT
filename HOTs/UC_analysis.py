import sys
import os
from os.path import join

import warnings
warnings.filterwarnings('ignore')
import configparser
###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
from overbinders.data_prep.basic import load_metadata
import numpy as np

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

# PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_promoters.hg19.bed.gz"
# TSS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.hg19.bed.gz"

TSS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/TSS.bed.gz"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"

SE_DIR = "/panfs/pan1/devdcode/sanjar/superenhancers/"
UC_FILE = "/panfs/pan1/devdcode/common/UCNEbase/hg19_UCNE_coord.bed"

from lib.TSS import get_tx_id2gene, get_enst2gene, get_enst2ensg


def load_TSS():

    enst2ensg = get_enst2ensg()
    tss_bed = BedTool(TSS_FILE)
    bed_lines = []
    for r in tss_bed:
        p = r.fields
        out_line = "\t".join(p[:4] + ["-" if p[4] == "-" else enst2ensg[p[4]]] + [p[5]])
        bed_lines.append(out_line)

    return BedTool("\n".join(bed_lines), from_string=True)


def extract_loci(cl="HepG2"):

    PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_promoters.hg19.bed.gz"
    TSS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.hg19.bed.gz"

    save_file = join(LOCI_DIR, "UC_%s.bed" % cl)
    print(save_file)

    hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl))

    ucs = BedTool(UC_FILE).intersect(hots, wo=True)
    ucs = ucs.groupby(c=[7, 8, 10], o=["min", "max", "collapse"])

    print("Total UCs: %d" % ucs.count())

    uc_tss = ucs.closest(TSS_FILE, d=True).groupby(g=[1, 2, 3, 4, 5, 6], c=[12, 13], o=["distinct"])

    columns = "#chrom\tUC start\tUC stop\tHOT start\tHOT stop\tDAPs\tGene\tdistance"
    uc_tss.saveas(save_file, trackline=columns)


def extract_loci_v2(cl="HepG2"):

    hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl))

    ucs = BedTool(UC_FILE).intersect(hots, wo=True)
    ucs = ucs.groupby(c=[7, 8, 10], o=["min", "max", "collapse"])

    print("Total UCs: %d" % ucs.count())

    columns = "#chrom\tUC start\tUC stop\tHOT start\tHOT stop\tDAPs\tGene ID\tGene name\tdistance"

    tss = load_TSS().groupby(c=[5, 6], o="distinct")
    window_tss = tss.window(ucs, w=50000).groupby(c=[4, 5], o="distinct")
    tss_closest_ucs = window_tss.closest(ucs, d=True)

    res_lines = []
    sel_ix = [5, 6, 7, 8, 9, 10, 3, 4, 11]
    for r in tss_closest_ucs:
        p = r.fields
        out_line = "\t".join([p[i] for i in sel_ix])
        res_lines.append(out_line)

    ucs_tss = BedTool("\n".join(res_lines), from_string=True).groupby(g=[1, 2, 3, 4, 5, 6, 7, 8], c=9, o="min")
    save_file = join(LOCI_DIR, "UC_%s.bed" % cl)
    print(save_file)
    ucs_tss.saveas(save_file, trackline=columns)


if __name__ == "__main__":

    pass






