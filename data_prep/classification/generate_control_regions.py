#!/usr/bin/env python

import warnings
warnings.filterwarnings('ignore')
import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import random
import tempfile
from os.path import join
from pybedtools import BedTool
from pybedtools import featurefuncs
from pathlib import Path

DATA_PATH = Path(os.environ["HOT_DATA"])
BINS_DIR = DATA_PATH/"log_bins"
HOTS_DIR = DATA_PATH/"HOTs"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

ALL_DHS_FILE = DATA_PATH / "src_files/all_dhs_merged.bed.gz"
HG19_FILE = DATA_PATH / "src_files/hg19_files/hg19.genome"
PROMOTERS_FILE = DATA_PATH / "src_files/hg19_files/promoters.merged.bed.gz"

SAVE_DIR = DATA_PATH / "classification/datasets/control_regions"
SAVE_DIR.mkdir(exist_ok=True, parents=True)

random.seed(1024)


def generate_controls(cl, name):

    # ["hots", "dhs_controls", "proms_controls", "re_controls"]

    hots_file = DATA_PATH / f"HOTs/{cl}_HOTs.bed.gz"
    hots = BedTool(hots_file)
    hots_count = hots.count()

    seq_lens = [400, 1000]

    if name == "hots":

        for seq_len in seq_lens:
            rs = hots.each(featurefuncs.midpoint).slop(g=HG19_FILE, l=seq_len / 2, r=(seq_len / 2) - 1).sort()
            save_file = join(SAVE_DIR, f"{cl}_{name}_{seq_len}.bed.gz")
            rs.saveas(save_file)
        return

    if name == "dhs_controls":

        dhs = list(BedTool(ALL_DHS_FILE).intersect(hots, v=True, wa=True))
        dhs = random.sample(dhs, hots_count * 10)
        dhs = BedTool(dhs)
        controls = dhs

    elif name == "proms_controls":
        proms = BedTool(str(PROMOTERS_FILE))
        proms = proms.intersect(hots, v=True, wa=True).each(featurefuncs.midpoint).merge(d=-1)
        controls = proms

    elif name == "re_controls":
        fname = DATA_PATH / f"src_files/{cl}_enhancers_DHS_H3K27ac.bed.gz"
        regular_enhancers = BedTool(fname)
        regular_enhancers = regular_enhancers.intersect(hots, v=True, wa=True).each(featurefuncs.midpoint).merge(d=-1)
        controls = regular_enhancers

    for seq_len in seq_lens:

        save_file = join(SAVE_DIR, f"{cl}_{name}_{seq_len}.bed.gz")
        controls = controls.each(featurefuncs.midpoint).slop(g=HG19_FILE, l=seq_len / 2, r=(seq_len / 2) - 1).sort()
        controls.saveas(save_file)


def generate_controls_all(cl="HepG2"):

    hots_file = DATA_PATH/f"HOTs/{cl}_HOTs.bed.gz"
    hots = BedTool(hots_file)
    hots_count = hots.count()

    # DHS-based controls
    dhs = list(BedTool(ALL_DHS_FILE).intersect(hots, v=True, wa=True))
    dhs = random.sample(dhs, hots_count * 10)
    dhs = BedTool(dhs)

    # controls from promoters and regular enhancers
    proms = BedTool(str(PROMOTERS_FILE))
    proms = proms.intersect(hots, v=True, wa=True).each(featurefuncs.midpoint).merge(d=-1)

    fname = DATA_PATH/f"src_files/{cl}_enhancers_DHS_H3K27ac.bed.gz"
    regular_enhancers = BedTool(fname)
    regular_enhancers = regular_enhancers.intersect(hots, v=True, wa=True).each(featurefuncs.midpoint).merge(d=-1)

    bed_pool = [hots, dhs, proms, regular_enhancers]
    names_pool = ["hots", "dhs_controls", "proms_controls", "re_controls"]

    for bed, name in zip(bed_pool, names_pool):
        first_count = bed.count()
        for seq_len in [400, 1000]:
            final_regions = bed.each(featurefuncs.midpoint).slop(g=HG19_FILE, l=seq_len / 2, r=(seq_len / 2) - 1).sort()
            save_file = join(SAVE_DIR, f"{cl}_{name}_{seq_len}.bed.gz")
            print(first_count, final_regions.count(), save_file)
            final_regions.saveas(save_file)


if __name__ == "__main__":

    cl = sys.argv[1]
    generate_controls_all(cl)

    # for cl in ["HepG2", "K562"]:
    #     generate_controls_test(cl)
