import warnings
warnings.filterwarnings('ignore')

from collections import defaultdict
import sys
import os
import tempfile
from os.path import join
import numpy as np
import configparser

###############################################################
import pandas

config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool

import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.style.use('ggplot')
import seaborn as sns
from overbinders.data_prep.LLPS import get_llps_tfs
from overbinders.data_prep.basic import load_metadata
from overbinders.log_bins import extract_data

INTRONS_FILE ="ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownGene_introns.liftOver_hg19.bed"
PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.alternative_TSS.liftOver_hg19.gene_symbols.bed"

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
HIC_DIR = "/panfs/pan1/devdcode/sanjar/HiC/encode/"


def describe_TADs():

    data = []
    labels = []
    cell_lines = ["HepG2", "K562"]

    for cl in cell_lines:

        tads_file = join(HIC_DIR, "%s_TADs.bedpe.gz" % cl)
        tads = BedTool(tads_file)

        tad_lens = np.asarray([r.length for r in tads])

        _label = "%s median: %d kbps" % (cl, int(np.median(tad_lens)/1000))
        labels.append(_label)

        data.append(tad_lens)

    plt.close()
    plt.boxplot(data, labels=labels)
    plt.title("TAD sizes")
    plt.legend()
    ax = plt.gca()
    ax.set_xticklabels(cell_lines)
    plt.show()


def load_TADs_with_HOTs():

    tmp_file = tempfile.NamedTemporaryFile()
    columns = ["cell_line", "category", "value"]
    outf = open(tmp_file.name, "w")
    outf.write(",".join(columns) + "\n")
    out_fmt = "%s,%s,%d\n"
    cell_lines = ["HepG2", "K562"]
    for cl in cell_lines:

        tads_file = join(HIC_DIR, "%s_TADs.bedpe.gz" % cl)
        tads = BedTool(tads_file)
        tad_lens = np.asarray([r.length for r in tads])
        [outf.write(out_fmt % (cl, "all_TADs", _)) for _ in tad_lens]

        hots_file = join(LOCI_DIR, "%s_HOTs.bed" % cl)
        hots_prom_file = join(LOCI_DIR, "%s_HOTs.proms.bed" % cl)
        hots_enhs_file = join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)

        hot_tads = tads.intersect(hots_file, wa=True).groupby(c=[4, 5], o="collapse")
        hot_tad_lens = np.asarray([r.length for r in hot_tads])
        [outf.write(out_fmt % (cl, "HOT TADs", _)) for _ in hot_tad_lens]

        hot_prom_enh_tads = tads.intersect(hots_prom_file, wa=True).\
                                 groupby(c=3, o="collapse").\
                                 intersect(hots_enhs_file, wa=True).\
                                 groupby(c=3, o="collapse")
        hot_prom_enh_tad_lens = np.asarray([r.length for r in hot_prom_enh_tads])
        [outf.write(out_fmt % (cl, "HOT_PE TADs", _)) for _ in hot_prom_enh_tad_lens]

        hot_enhs_tads = tads.intersect(hots_enhs_file, wa=True).\
                             groupby(c=3, o="collapse").\
                             intersect(hots_prom_file, v=True).\
                             groupby(c=3, o="collapse")
        [outf.write(out_fmt % (cl, "HOT_E TADs", r.length)) for r in hot_enhs_tads]
    df = pandas.read_csv(tmp_file.name)

    _order = ["all_TADs", "HOT TADs", "HOT_PE TADs", "HOT_E TADs"]
    plt.figure()
    sns.boxplot(x="cell_line", y="value", hue="category", data=df, showfliers=False, hue_order=_order)
    plt.legend()
    plt.xlabel("")
    plt.title("HOTs in TADs")
    plt.ylabel("Size (bps)")
    plt.tight_layout()
    save_file = join(PLOTS_DIR, "TAD_sizes.pdf")
    plt.savefig(save_file)
    plt.show()

    return

    plt.figure()
    sns.barplot(x="cell_line", y="value", hue="category",
                data=df.groupby(["cell_line", "category"]).size().reset_index(name='value'), hue_order=_order)
    # plt.legend([], [], frameon=False)
    plt.legend()
    plt.title("HOTs in TADs")
    plt.ylabel("# of TADs")
    plt.xlabel("")
    plt.tight_layout()
    save_file = join(PLOTS_DIR, "TAD_numbers.pdf")
    plt.savefig(save_file)
    plt.show()


def TADs_gene_stats(df):

    # tmp_file = tempfile.NamedTemporaryFile()
    # columns = ["cell_line", "category", "value"]
    # outf = open(tmp_file.name, "w")
    # outf.write(",".join(columns) + "\n")
    # out_fmt = "%s,%s,%d\n"
    # cell_lines = ["HepG2", "K562"]
    # for cl in cell_lines:
    #     tads_file = join(HIC_DIR, "%s_TADs.bed.gz" % cl)
    #     hots_file = join(LOCI_DIR, "%s_HOTs.bed" % cl)
    #     all_tads = BedTool(tads_file)
    #     tad_w_proms = all_tads.intersect(PROMOTERS_FILE, wo=True).groupby(c=[9, 9], o=["count_distinct", "distinct"])
    #     all_tads_gene_nos = [int(r.fields[3]) for r in tad_w_proms] + [0] * (all_tads.count() - tad_w_proms.count())
    #     [outf.write(out_fmt % (cl, "All TADs", _)) for _ in all_tads_gene_nos]
    #     gene_tads_gene_nos = [int(r.fields[3]) for r in tad_w_proms]
    #     [outf.write(out_fmt % (cl, "Gene TADs", _)) for _ in gene_tads_gene_nos]
    #     hot_tads = all_tads.intersect(hots_file, wa=True).intersect(PROMOTERS_FILE, wo=True).groupby(c=[9, 9], o=["count_distinct", "distinct"])
    #     hot_tads_gene_nos = [int(r.fields[3]) for r in hot_tads]
    #     [outf.write(out_fmt % (cl, "HOT TADs", _)) for _ in hot_tads_gene_nos]
    # df = pandas.read_csv(tmp_file.name)

    _order = ["All TADs", "Gene TADs", "HOT TADs"]
    plt.figure()
    sns.boxplot(x="cell_line", y="value", hue="category", data=df, hue_order=_order, showfliers=False)
    plt.legend()
    plt.xlabel("")
    plt.title("Genes in TADs")
    plt.ylabel("# of genes")
    plt.tight_layout()
    save_file = join(PLOTS_DIR, "TAD_genes.pdf")
    plt.savefig(save_file)
    plt.show()


def extract_contacts(cl="HepG2"):

    # cl = "K562"

    contacts_file = join(HIC_DIR, "%s_chromatin_interactions.bed.gz" % cl)

    hots_file = join(LOCI_DIR, "%s_HOTs.PE.gene_symbol.bed" % cl)

    tads = BedTool(join(HIC_DIR, "%s_TADs.bed.gz" % cl))

    src_id2trg_id = {}
    src_id2hot_ids = defaultdict(list)
    trg_id2hot_ids = defaultdict(list)

    src_hots_overlaps = BedTool(contacts_file).intersect(hots_file, wo=True)

    trg_loci = []

    for r in src_hots_overlaps:
        src_id = "%s-%d-%d" % (r.chrom, r.start, r.stop)
        trg_id = "%s-%s-%s" % (r.fields[3], r.fields[4], r.fields[5])
        hot_id = "%s-%s-%s-%s-%s" % (r.fields[6], r.fields[7], r.fields[8], r.fields[9], r.fields[10])

        src_id2trg_id[src_id] = trg_id
        src_id2hot_ids[src_id].append(hot_id)

        trg_loci.append(trg_id.replace("-", "\t"))

    trg_loci = BedTool("\n".join(set(trg_loci)), from_string=True)

    trg_hots_overlaps = trg_loci.intersect(hots_file, wo=True)

    for r in trg_hots_overlaps:
        trg_id = "%s-%s-%s" % (r.fields[0], r.fields[1], r.fields[2])
        hot_id = "%s-%s-%s-%s-%s" % (r.fields[3], r.fields[4], r.fields[5], r.fields[6], r.fields[7])
        trg_id2hot_ids[trg_id].append(hot_id)

    lines = []

    for src_id, trg_id in src_id2trg_id.items():
        if trg_id not in trg_id2hot_ids:
            continue

        [_, src_from, src_to] = src_id.split("-")
        src_from, src_to = int(src_from), int(src_to)
        [_, trg_from, trg_to] = trg_id.split("-")
        trg_from, trg_to = int(trg_from), int(trg_to)

        src_hot_ids = set(src_id2hot_ids[src_id])
        trg_hot_ids = set(trg_id2hot_ids[trg_id])

        src_hot_symbols = set([symbol for _id in src_hot_ids for symbol in _id.split("-")[-1].split(",")])
        src_hot_symbols = ",".join(src_hot_symbols)
        trg_hot_symbols = set([symbol for _id in trg_hot_ids for symbol in _id.split("-")[-1].split(",")])
        trg_hot_symbols = ",".join(trg_hot_symbols)

        _str = src_id.replace("-", "\t") + "\n" + trg_id.replace("-", "\t")
        _bed = BedTool(_str, from_string=True)

        tad_intersects = tads.intersect(_bed, wo=True)
        total_tads = tad_intersects.count()
        if total_tads > 0:
            tads_groupby = tad_intersects.groupby(c=[5, 5, 6], o=["count", "collapse", "collapse"])
            tad_mode = ",".join([r.fields[3] for r in tads_groupby])
        else:
            tad_mode = "0"

        out_line = src_id.replace("-", "\t") + "\t" + \
                   trg_id.replace("-", "\t") + "\t" + \
                   str(trg_from-src_to) + "\t" + \
                   str(len(src_hot_ids)) + "\t" + \
                   str(len(trg_hot_ids)) + "\t" + \
                   tad_mode + "\t" + \
                   ";".join(src_hot_ids) + "\t" + \
                   ";".join(trg_hot_ids) + "\t" + \
                   src_hot_symbols + "\t" + \
                   trg_hot_symbols + "\n"

        lines.append(out_line)

    out_file = join(LOCI_DIR, "%s_contact_hots.txt" % cl)
    of = open(out_file, "w")
    header = ["chr_1", "start_1", "stop_1", "chr_2", "start_2", "stop_2",
              "distance",
              "n_HOTs_1",
              "n_HOTs_2",
              "TAD_mode",
              "HOTs_1",
              "HOTs_2",
              "gene_symbols_1",
              "gene_symbols_2"]
    header = "#" + "\t".join(header) + "\n"
    of.write(header)
    [of.write(l) for l in lines]


def extract_contacts_v2():

    contacts_file = join(HIC_DIR, "ENCFF128OJZ.bedpe.gz")

    hots_file = join(LOCI_DIR, "HepG2_HOTs.bed")

    src_id2trg_id = {}
    src_id2hot_ids = defaultdict(list)
    trg_id2hot_ids = defaultdict(list)

    src_hots_overlaps = BedTool(contacts_file).intersect(hots_file, wo=True)

    return src_hots_overlaps

    trg_loci = []

    for r in src_hots_overlaps:
        src_id = "%s-%d-%d" % (r.chrom, r.start, r.stop)
        trg_id = "%s-%s-%s" % (r.fields[3], r.fields[4], r.fields[5])
        hot_id = "%s-%s-%s-%s-%s" % (r.fields[6], r.fields[7], r.fields[8], r.fields[9], r.fields[10])

        src_id2trg_id[src_id] = trg_id
        src_id2hot_ids[src_id].append(hot_id)

        trg_loci.append(trg_id.replace("-", "\t"))

    trg_loci = BedTool("\n".join(set(trg_loci)), from_string=True)

    trg_hots_overlaps = trg_loci.intersect(hots_file, wo=True)

    for r in trg_hots_overlaps:
        trg_id = "%s-%s-%s" % (r.fields[0], r.fields[1], r.fields[2])
        hot_id = "%s-%s-%s-%s-%s" % (r.fields[3], r.fields[4], r.fields[5], r.fields[6], r.fields[7])
        trg_id2hot_ids[trg_id].append(hot_id)

    lines = []

    for src_id, trg_id in src_id2trg_id.items():
        if trg_id not in trg_id2hot_ids:
            continue

        [_, src_from, src_to] = src_id.split("-")
        src_from, src_to = int(src_from), int(src_to)
        [_, trg_from, trg_to] = trg_id.split("-")
        trg_from, trg_to = int(trg_from), int(trg_to)

        src_hot_ids = set(src_id2hot_ids[src_id])
        trg_hot_ids = set(trg_id2hot_ids[trg_id])

        src_hot_symbols = set([symbol for _id in src_hot_ids for symbol in _id.split("-")[-1].split(",")])
        src_hot_symbols = ",".join(src_hot_symbols)
        trg_hot_symbols = set([symbol for _id in trg_hot_ids for symbol in _id.split("-")[-1].split(",")])
        trg_hot_symbols = ",".join(trg_hot_symbols)

        _str = src_id.replace("-", "\t") + "\n" + trg_id.replace("-", "\t")
        _bed = BedTool(_str, from_string=True)

        tad_intersects = tads.intersect(_bed, wo=True)
        total_tads = tad_intersects.count()
        if total_tads > 0:
            tads_groupby = tad_intersects.groupby(c=[5, 5, 6], o=["count", "collapse", "collapse"])
            tad_mode = ",".join([r.fields[3] for r in tads_groupby])
        else:
            tad_mode = "0"

        out_line = src_id.replace("-", "\t") + "\t" + \
                   trg_id.replace("-", "\t") + "\t" + \
                   str(trg_from - src_to) + "\t" + \
                   str(len(src_hot_ids)) + "\t" + \
                   str(len(trg_hot_ids)) + "\t" + \
                   tad_mode + "\t" + \
                   ";".join(src_hot_ids) + "\t" + \
                   ";".join(trg_hot_ids) + "\t" + \
                   src_hot_symbols + "\t" + \
                   trg_hot_symbols + "\n"

        lines.append(out_line)

    out_file = join(LOCI_DIR, "%s_contact_hots.txt" % cl)
    of = open(out_file, "w")
    header = ["chr_1", "start_1", "stop_1", "chr_2", "start_2", "stop_2",
              "distance",
              "n_HOTs_1",
              "n_HOTs_2",
              "TAD_mode",
              "HOTs_1",
              "HOTs_2",
              "gene_symbols_1",
              "gene_symbols_2"]
    header = "#" + "\t".join(header) + "\n"
    of.write(header)
    [of.write(l) for l in lines]

