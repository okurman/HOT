
from collections import Counter, defaultdict
import seaborn as sns
import pandas
from scipy.stats import mannwhitneyu, kruskal, kstest
import warnings
warnings.filterwarnings('ignore')

import sys
import os
from os.path import join
import numpy as np


from pybedtools import BedTool

import matplotlib.pyplot as plt

INTRONS_FILE ="ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownGene_introns.liftOver_hg19.bed"
PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.bed"
PRIM_PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.primary_TSS.liftOver_hg19.bed"
HG19_FASTA = "/panfs/pan1/devdcode/common/genomes/hg19/hg19.fa"

PROJECT_DIR = "ROOT_DIR/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"

hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed")).merge(d=-1)
hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed")).merge(d=-1)


def get_HOT_overlaps(test_regions_file=None):

    test_regions_file = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/cpgIslandExt.bed"

    prom_ovp = hot_proms.intersect(test_regions_file, wa=True).count()/hot_proms.count()
    enh_ovp = hot_enhs.intersect(test_regions_file, wa=True).count()/hot_enhs.count()

    print("HOT proms fraction: %f" % prom_ovp)
    print("HOT enhs fraction: %f" % enh_ovp)


def get_GC():

    print("Extracting nuc")
    prom_ovp = hot_proms.nucleotide_content(fi=HG19_FASTA)
    enh_ovp = hot_enhs.nucleotide_content(fi=HG19_FASTA)

    print("Loading to data")
    data = []
    data += [["prom", float(r.fields[4])] for r in prom_ovp]
    data += [["enh", float(r.fields[4])] for r in enh_ovp]

    df = pandas.DataFrame(data=data, columns=["region", "gc"])
    return df

