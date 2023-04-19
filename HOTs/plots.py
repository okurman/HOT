
import warnings
from collections import Counter, defaultdict
import seaborn as sns
import pandas
from scipy.stats import mannwhitneyu, kruskal, kstest

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

INTRONS_FILE ="/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/knownGene.introns.merged.bed.gz"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"
# PRIM_PROMOTERS_FILE = "ROOT_DIR/data/genomes/hg19/annotations/genes/gencode/knownCanonical.primary_TSS.liftOver_hg19.bed"

PROJECT_DIR = "ROOT_DIR/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"
LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"


def hots_prom_enh_barplot(data=None):

    if not data:
        data = []
        for cl in ["HepG2"]:
        # for cl in ["HepG2", "K562"]:
            hot_proms = BedTool(join(LOCI_DIR, "%s_HOTs.proms.bed" % cl))
            hot_enhs = BedTool(join(LOCI_DIR, "%s_HOTs.enhs.bed" % cl))
            cnt_proms = hot_proms.count()
            cnt_intronic = hot_enhs.intersect(INTRONS_FILE, wa=True).groupby(c=4, o="collapse").count()
            cnt_intergenic = hot_enhs.intersect(INTRONS_FILE, v=True, wa=True).groupby(c=4, o="collapse").count()
            data.append([cnt_proms, cnt_intronic, cnt_intergenic])

        return data

    data = np.array(data[0])
    data = 100 * (data/np.sum(data))

    labels = ["promoters", "intronic", "intergenic"]
    plt.figure(figsize=(3, 4))
    ax = plt.gca()

    bottoms = np.cumsum([0, data[0], data[1]])
    gs = [ax.bar([0], data[i], label=labels[i], bottom=bottoms[i]) for i in range(3)]
    [ax.bar_label(g, label_type="center", fmt='%d %%', fontweight="bold", color="white") for g in gs]

    ax.set_ylim([0, 100])
    ax.legend(bbox_to_anchor=(1.05, 0.5), ncol=1)
    ax.set_xticks([0])
    ax.set_xticklabels(["HOT loci"])
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "HepG2_introns_intergenic_barchart_tall_v2.pdf")
    print(save_file)
    plt.savefig(save_file)


def hots_prom_enh_barplot_in_legend(data=None):

    if not data:
        data = []
        for cl in ["HepG2"]:
            # for cl in ["HepG2", "K562"]:
            hot_proms = BedTool(join(LOCI_DIR, "%s_HOTs.proms.bed" % cl))
            hot_enhs = BedTool(join(LOCI_DIR, "%s_HOTs.enhs.bed" % cl))
            cnt_proms = hot_proms.count()
            cnt_intronic = hot_enhs.intersect(INTRONS_FILE, wa=True).groupby(c=4, o="collapse").count()
            cnt_intergenic = hot_enhs.intersect(INTRONS_FILE, v=True, wa=True).groupby(c=4, o="collapse").count()
            data.append([cnt_proms, cnt_intronic, cnt_intergenic])

        # return data

    data = np.array(data[0])
    data = 100 * (data / np.sum(data))
    data = np.round(data)

    labels = ["promoters", "intronic", "intergenic"]
    plt.figure(figsize=(1.8, 3))
    ax = plt.gca()

    bottoms = np.cumsum([0, data[0], data[1]])
    gs = [ax.bar([0], data[i], label=labels[i], bottom=bottoms[i]) for i in range(3)]
    [ax.bar_label(g, label_type="center", fmt='%d %%\n{0}'.format(l), fontweight="bold", color="white") for g, l in zip(gs, labels)]

    # ax.set_ylim([0, 100])
    ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
    # ax.legend(bbox_to_anchor=(1.05, 0.5), ncol=1)
    ax.set_xticks([0])
    ax.set_xticklabels(["HOT loci"])
    # ax.spines['right'].set_visible(False)

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "HepG2_introns_intergenic_barchart_tall_in_legend_v2.pdf")
    print(save_file)
    plt.savefig(save_file)


def hot_gene_distance():

    cl = "HepG2"
    hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl)).sort()
    enhs = BedTool(join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)).sort()

    closest = enhs.closest(BedTool(PROMOTERS_FILE).sort(), d=True)

    distal = [r for r in closest if int(r.fields[-1]) > 50000]

    print(len(distal), hots.count(), closest.count())


def hots_piechart_DHS():

    data = []
    for cl in ["HepG2", "K562"]:
        hot_proms = BedTool(join(LOCI_DIR, "%s_HOTs.proms.bed" % cl))
        hot_enhs = BedTool(join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl))
        dhs_file = "ROOT_DIR/chipseq_files/DHS_%s.bed.gz" % cl

        dhs = BedTool(dhs_file)
        proms = BedTool(hot_proms)
        enhs = BedTool(hot_enhs)

        dhs_cov = sum([r.length for r in dhs])
        dhs_prom_cov = sum([r.length for r in dhs.intersect(proms)])
        dhs_enh_cov = sum([r.length for r in dhs.intersect(enhs)])

        data.append([dhs_cov-dhs_prom_cov-dhs_enh_cov, dhs_prom_cov, dhs_enh_cov])

    labels = ["rest", "HOT prom", "HOT enh"]
    colors = ["grey", "orange", "blue"]
    fig1, [ax1, ax2] = plt.subplots(1, 2)
    ax1.pie(data[0], labels=labels, autopct='%2.1f%%', startangle=90, counterclock=False)
    ax1.axis('equal')
    ax1.set_title("HepG2")
    ax2.pie(data[1], labels=labels, autopct='%2.1f%%', startangle=90, counterclock=False)
    ax2.axis('equal')
    ax2.set_title("K562")

    plt.tight_layout()
    save_file = join(PLOTS_DIR, "piechart_DHS_intersections.pdf")
    plt.savefig(save_file)
    plt.show()


def cluster_HOTs():

    BINS = 6
    # CLUSTER_SIZE = 20000
    # save_file = join(PLOTS_DIR, "clusters_20kbps.pdf")
    # CLUSTER_SIZE = 10000
    # save_file = join(PLOTS_DIR, "clusters_10kbps.pdf")
    CLUSTER_SIZE = 5000
    save_file = join(PLOTS_DIR, "clusters_5kbps.pdf")

    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(8, 4))

    for ax, cl in zip([ax1, ax2], ["HepG2", "K562"]):

        hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl)).sort()
        merged_slops = hots.slop(b=(CLUSTER_SIZE-400)/2, genome="hg19").merge()
        clusters = merged_slops.intersect(hots, wo=True).groupby(c=4, o="count")

        bins_proms = []
        bins_enhs = []

        for i in range(BINS):

            if i < BINS - 1:
                bin_rs = [r for r in clusters if int(r.fields[3]) == i+1]
            else:
                bin_rs = [r for r in clusters if int(r.fields[3]) >= i]

            bin_bed = BedTool(bin_rs)
            cnt = bin_bed.intersect(PROMOTERS_FILE, wo=True).groupby(c=4, o="collapse").count()
            bins_proms.append(cnt)
            cnt = bin_bed.intersect(PROMOTERS_FILE, v=True).groupby(c=4, o="collapse").count()
            bins_enhs.append(cnt)

        bins_enhs = np.array(bins_enhs)
        bins_proms = np.array(bins_proms)
        totals = bins_enhs + bins_proms
        prom_perc = 100 * bins_proms / totals
        prom_perc = np.asarray(prom_perc, dtype=int)

        x_range = np.arange(len(bins_enhs))

        ax.bar(x_range, bins_enhs, label="enhancers")
        ax.bar(x_range, bins_proms, bottom=bins_enhs, label="promoters")

        for i in range(len(bins_enhs)):
            ax.text(i, bins_enhs[i] + bins_proms[i]/2, "%d%%" % prom_perc[i], ha="center")

        x_tick_labels = ["%d" % (_ + 1) for _ in range(BINS - 1)] + [">%d" % BINS]

        ax.set_xticks(x_range)
        ax.set_xticklabels(x_tick_labels)
        ax.set_title(cl)
        ax.legend()
        ax.set_xlabel("# of HOTs in clusters")
        ax.set_ylabel("# clusters")

    plt.tight_layout()

    plt.savefig(save_file)

    plt.show()


def load_chromHMM_data(cl):

    if cl == "HepG2":
        chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmHepg2HMM.bed"
    else:
        chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmK562HMM.bed"

    hots_file = join(LOCI_DIR, "%s_HOTs.PE.bed" % cl)

    hots = BedTool(hots_file)
    proms = [r for r in hots if r.fields[3] == "P"]
    enhs = [r for r in hots if r.fields[3] == "E"]

    enhs = BedTool(enhs).intersect(chromhmm_file, wo=True).groupby(c=8, o="collapse")
    proms = BedTool(proms).intersect(chromhmm_file, wo=True).groupby(c=8, o="collapse")

    e_states = Counter([r.fields[13] for r in enhs])
    p_states = Counter([r.fields[13] for r in proms])

    tmp = defaultdict(int)

    for k, v in p_states.items():
        _k = " ".join(k.split("_")[1:])
        tmp[_k] += v
    p_states = tmp

    tmp = defaultdict(int)
    for k, v in d_states.items():
        _k = " ".join(k.split("_")[1:])
        tmp[_k] += v
    d_states = tmp

    main_keys = ['Active Promoter', 'Weak Promoter', 'Strong Enhancer', 'Weak Enhancer']

    tmp = defaultdict(int)
    for k, v in p_states.items():
        if k in main_keys:
            tmp[k] = v
        else:
            tmp["Rest"] += v
    p_states = tmp

    tmp = defaultdict(int)
    for k, v in d_states.items():
        if k in main_keys:
            tmp[k] = v
        else:
            tmp["Rest"] += v
    d_states = tmp

    return p_states, d_states


def load_DHS_ATAC_data():

    eids = ["ENCFF439EIO", "ENCFF897NME"]
    names = ["ATAC-seq", "DHS"]

    hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed"))
    hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed"))

    data = []

    for eid, name in zip(eids, names):
        input_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/%s_hg19.bed.gz" % eid
        ca_regions = BedTool(input_file)

        for reg, reg_name in zip([hot_proms, hot_enhs], ["HOT promoters", "HOT enhancers"]):

            cov = sum([r.length for r in reg.intersect(ca_regions).sort().merge()])
            total = sum([r.length for r in reg])
            fr_ca = cov / total

            # regions_count = reg.count()
            # fr_ca = reg.intersect(ca_regions, f=0.5, wa=True).merge(d=-1).count()/regions_count

            data.append([name, reg_name, fr_ca])

    df = pandas.DataFrame(data=data, columns=["ca_data", "region", "value"])

    df["value"] = np.round(100 * df["value"])

    return df


def plot_DHS_ATAC(df):

    # df = load_DHS_ATAC_data()
    save_file = join(PLOTS_DIR, "HepG2_HOT_DHS_ATAC_cov.pdf")
    plt.figure(figsize=(3, 3.1))
    ax = plt.gca()

    sns.barplot(x="region", y="value", hue="ca_data", data=df, ax=ax, order=["HOT promoters", "HOT enhancers"])

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["HOT\npromoters", "HOT\nenhancers"])
    ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
    ax.grid(axis="y", alpha=0.2)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False)

    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_ATAC(df):

    # df = load_DHS_ATAC_data()
    df = df[df["ca_data"] == "ATAC-seq"]
    save_file = join(PLOTS_DIR, "HepG2_HOT_ATAC_cov.pdf")
    plt.figure(figsize=(3, 3.1))
    ax = plt.gca()
    
    color = sns.color_palette("tab10")[0]

    sns.barplot(x="region", y="value", data=df, ax=ax, order=["HOT promoters", "HOT enhancers"], color=color)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["HOT\npromoters", "HOT\nenhancers"])
    ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
    ax.grid(axis="y", alpha=0.2)
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False)

    ax.set_xlabel("")
    ax.set_ylabel("ATAC-seq")

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_DHS_fractions():

    # dhs_file = "ROOT_DIR/chipseq_files/DHS_HepG2/ENCFF628UDH_hg19.bed"

    data_type = "DHS"
    """Two DHS experiments on ENCODE. Default files are: ENCFF897NME ENCFF748QCZ """
    # eid = "ENCFF897NME"
    eid = "ENCFF748QCZ"
    dhs_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/%s_hg19.bed.gz" % eid
    save_file = join(PLOTS_DIR, "HepG2_HOT_DHS_fractions_tall_%s_cov.pdf" % eid)

    # data_type = "ATAC-seq"
    # f=0.5
    # """Two ATAC experiments on ENCODE. Default files are: ENCFF439EIO ENCFF913MQB """
    # # eid = "ENCFF439EIO"
    # eid = "ENCFF913MQB"
    # dhs_file = "/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/%s_hg19.bed.gz" % eid
    # save_file = join(PLOTS_DIR, "HepG2_HOT_ATAC_fractions_tall_%s_f%s.pdf" % (eid, str(f)))

    hot_enh_file = join(LOCI_DIR, "HepG2_HOTs.enhs.bed")
    hot_pr_file = join(LOCI_DIR, "HepG2_HOTs.proms.bed")

    hot_en = BedTool(hot_enh_file)
    hot_pr = BedTool(hot_pr_file)
    dhs = BedTool(dhs_file)

    hot_en_count = hot_en.count()
    hot_en_non_dhs_count = hot_en.intersect(dhs, f=f, v=True).count()

    hot_pr_count = hot_pr.count()
    hot_pr_non_dhs_count = hot_pr.intersect(dhs, f=f, v=True).count()

    enh_dhs_frac = (hot_en_count - hot_en_non_dhs_count) / hot_en_count
    pr_dhs_frac = (hot_pr_count - hot_pr_non_dhs_count) / hot_pr_count

    dhs_percs = np.ceil([enh_dhs_frac*100, pr_dhs_frac*100])
    non_dhs_percs = np.floor([(1 - enh_dhs_frac) * 100, (1 - pr_dhs_frac) * 100])

    plt.figure(figsize=(3, 4))
    ax = plt.gca()

    g1 = plt.bar([0, 1], dhs_percs, label=data_type)
    g2 = plt.bar([0, 1], non_dhs_percs, bottom=dhs_percs, label="non %s" % data_type)

    plt.bar_label(g1, label_type="center", fmt='%d %%', fontweight="bold", color="white")
    plt.bar_label(g2, label_type="center", fmt='%d %%', fontweight="bold", color="white")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["HOT\nenhancers", "HOT\npromoters"])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, frameon=False)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)
    plt.close()



def conservation_scores_H1(df):

    # _loci_dir = "ROOT_DIR/definitions/peak_8bp_v2/log_bins/"
    # data = []
    # for cl in ["HepG2", "K562", "H1"]:
    #     for i in range(10, 14):
    #         _l = join(_loci_dir, "%s_400_loci.%d.bed.noprom.phastcons" % (cl, i))
    #         scores = np.loadtxt(_l, usecols=(3))
    #         [data.append([cl, "enh", _s]) for _s in scores]
    #
    #         _l = join(_loci_dir, "%s_400_loci.%d.bed.prom.phastcons" % (cl, i))
    #         scores = np.loadtxt(_l, usecols=(3))
    #         [data.append([cl, "prom", _s]) for _s in scores]
    # df = pandas.DataFrame(data=data, columns=["cell_line", "category", "score"])
    # return df

    # for cl in ["HepG2", "K562", "H1"]:
    #     prom_scores = df[(df["cell_line"]==cl) & (df["category"]=="prom")]["score"].values
    #     enh_scores = df[(df["cell_line"]==cl) & (df["category"]=="enh")]["score"].values
    #
    #     print(cl)
    #     print(mannwhitneyu(prom_scores, enh_scores))
    #     print(kruskal(prom_scores, enh_scores))
    #     print(np.mean(prom_scores)/np.mean(enh_scores))
    #     print("\n")
    # return

    plt.close()
    plt.figure(figsize=(4, 3))

    ax = plt.gca()
    g = sns.boxplot(x="cell_line", y="score", hue="category", data=df, showfliers=False, ax=ax)

    # l = g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2)
    l = g.legend(loc="upper center", ncol=2)
    l.get_texts()[0].set_text('enhancers')
    l.get_texts()[1].set_text('promoters')
    ax.set_ylim([0, 1])
    # ax.set_ymargin(0.5)
    ax.grid(axis="y", alpha=0.2)
    ax.set_xlabel("")
    ax.set_ylabel("phastCons score")

    plt.tight_layout()
    save_file = join(PLOTS_DIR, "phastcons_3CL_HOTs.pdf")
    print(save_file)
    plt.savefig(save_file)



