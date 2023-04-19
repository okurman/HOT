
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

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

# PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed"
PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_promoters.hg19.bed.gz"

# from lib.TSS import get_ENST2gene_symbols
# enst2gene = get_ENST2gene_symbols()


def extract_loci(cl="HepG2"):

    ###################################################################################
    ###################################################################################
    ####    extract blacklist-free loci
    blacklist_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/blacklisted_regions/hg19-blacklist.v2.bed"
    orig_dir = "ROOT_DIR/definitions/peak_8bp_v2"
    orig_loci = ["%s_400_loci.last_4.bed", "%s_400_loci.last_4.bed.prom", "%s_400_loci.last_4.bed.noprom"]
    orig_loci = [join(orig_dir, _ % cl) for _ in orig_loci]
    filtered_loci = ["%s_HOTs.bed", "%s_HOTs.proms.bed", "%s_HOTs.noproms.bed"]
    filtered_loci = [join(LOCI_DIR, _ % cl) for _ in filtered_loci]

    orig_loci = [BedTool(_) for _ in orig_loci]
    for _orig_loci, _filtered_loci in zip(orig_loci, filtered_loci):
        _orig_loci.intersect(blacklist_file, v=True, wa=True).sort().saveas(_filtered_loci)
    ####################################################################################

    ####################################################################################
    ##### extract nearest TSS
    input_files = [join(LOCI_DIR, "%s_HOTs.proms.bed" % cl), join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)]

    for input_file in input_files:
        loci = BedTool(input_file).merge(d=-1)
        save_file = input_file.replace(".bed", ".nearest_promoter.bed")

        grp = loci.closest(PROMOTERS_FILE, d=True).sort().groupby(c=[8, 9], o="distinct")

        out_list = []
        for r in grp:
            parts = r.fields
            genes = ",".join([_ for _ in set([enst2gene[_] if _ in enst2gene else _ for _ in parts[3].split(",")])])
            parts[3] = genes
            out_line = "\t".join(parts)
            out_list.append(out_line)

        BedTool("\n".join(out_list), from_string=True).saveas(save_file)


def extract_loci_new_proms(cl="HepG2"):
    ###################################################################################
    ###################################################################################
    ####    extract blacklist-free loci
    blacklist_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/blacklisted_regions/hg19-blacklist.v2.bed"
    orig_dir = "ROOT_DIR/definitions/peak_8bp_v2"
    orig_loci = ["%s_400_loci.last_4.bed", "%s_400_loci.last_4.bed.prom", "%s_400_loci.last_4.bed.noprom"]
    orig_loci = [join(orig_dir, _ % cl) for _ in orig_loci]
    filtered_loci = ["%s_HOTs.bed", "%s_HOTs.proms.bed", "%s_HOTs.noproms.bed"]
    filtered_loci = [join(LOCI_DIR, _ % cl) for _ in filtered_loci]

    orig_loci = [BedTool(_) for _ in orig_loci]
    for _orig_loci, _filtered_loci in zip(orig_loci, filtered_loci):
        _orig_loci.intersect(blacklist_file, v=True, wa=True).sort().saveas(_filtered_loci)
    ####################################################################################

    ####################################################################################
    ##### extract nearest TSS
    input_files = [join(LOCI_DIR, "%s_HOTs.proms.bed" % cl), join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)]

    for input_file in input_files:
        loci = BedTool(input_file).merge(d=-1)
        save_file = input_file.replace(".bed", ".nearest_promoter.bed")

        grp = loci.closest(PROMOTERS_FILE, d=True).sort().groupby(c=[7,8], o="distinct")
        grp.saveas(save_file)

        # out_list = []
        # for r in grp:
        #     parts = r.fields
        #     # genes = ",".join([_ for _ in set([enst2gene[_] if _ in enst2gene else _ for _ in parts[3].split(",")])])
        #     # genes = ",".join([_ for _ in set([enst2gene[_] if _ in enst2gene else _ for _ in parts[3].split(",")])])
        #     # parts[3] = genes
        #     out_line = "\t".join(parts)
        #     out_list.append(out_line)

        # BedTool("\n".join(out_list), from_string=True).saveas(save_file)








if __name__ == "__main__":

    pass






