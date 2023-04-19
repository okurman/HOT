import warnings
from collections import Counter, defaultdict
from copy import deepcopy
import pandas
import seaborn as sns
from scipy.stats import mannwhitneyu, binom_test, hypergeom

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

PROJECT_DIR = "ROOT_DIR/"
BINS_DIR = "ROOT_DIR/definitions/peak_8bp_v2/log_bins/"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
GWAS_FILE = "/panfs/pan1/devdcode/common/GWAS/gwas_catalog_v1.0.2-associations_e93_r2019-01-31.tsv"
LD_FILE = "/panfs/pan1/devdcode/common/GWAS/gwasSNP08.proxy"

from overbinders.data_prep.basic import load_metadata


def load_SNP2proxy_map():

    snp2ld_snps = defaultdict(list)

    for line in open(LD_FILE):
        parts = line.split("\t")
        if parts[3] == parts[4]:
            continue
        snp2ld_snps[parts[4]].append("\t".join(parts[:3]))

    return snp2ld_snps


def extract_traitswise_snps(snp2ld_snps=None):

    if snp2ld_snps:
        snp2ld_snps = load_SNP2proxy_map()

    trait2snps = defaultdict(set)
    snp_id2snp_bed = {}
    for l in open(GWAS_FILE):
        if l.startswith("DATE ADDED TO CATALOG"):
            continue
        parts = l.split("\t")
        trait = parts[7]

        if "x" in parts[21]:
            continue

        for snp_id, chrom, pos in zip(parts[21].split(";"), parts[11].split(";"), parts[12].split(";")):
            snp_id = snp_id.strip()
            chrom = "chr"+chrom.strip()
            pos = pos.strip()

            if not pos:
                continue

            snp_bed = "%s\t%s\t%s" % (chrom, pos, pos)
            snp_id2snp_bed[snp_id] = snp_bed
            trait2snps[trait].update([snp_id])

    trait2bed = dict()

    for trait, gwas_snps in trait2snps.items():
        bed_list = []
        for gwas_snp in gwas_snps:

            bed_list.append(snp_id2snp_bed[gwas_snp])
            if gwas_snp in snp2ld_snps:
                for ld_snp in snp2ld_snps[gwas_snp]:
                    bed_list.append(ld_snp)

        bed_str = "\n".join(set(bed_list))

        trait2bed[trait] = BedTool(bed_str, from_string=True)

    return trait2bed


def write_trait_beds_to_files(trait2bed):

    _dir = join(LOCI_DIR, "GWAS_analysis/traitwise_files/")

    trait2id = {}
    id2trait = {}
    cnt = 1
    for _trait, bed in trait2bed.items():
        trait2id[_trait] = cnt
        id2trait[cnt] = _trait

        save_file = join(_dir, "%d.bed.gz" % cnt)
        bed.saveas(save_file)
        cnt += 1

    map_file = join(LOCI_DIR, "GWAS_analysis/trait_no_map.txt")
    with open(map_file, "w") as of:
        for i in range(1, cnt):
            out_line = "%d\t%s\n" % (i, id2trait[i])
            of.write(out_line)


def extract_trait_overlap_SNPs(trait2bed):

    map_file = join(LOCI_DIR, "GWAS_analysis/trait_no_map.txt")
    id2trait = {}
    trait2id = {}
    for l in open(map_file):
        parts = l.strip().split("\t")
        id = int(parts[0])
        trait = parts[2]
        id2trait[id] = trait
        trait2id[trait] = id

    save_dir = join(LOCI_DIR, "GWAS_analysis/traitwise_files/")

    hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.PE.bed"))
    hot_enhs = BedTool([r for r in hots if r.fields[3] == "E"])
    hot_proms = BedTool([r for r in hots if r.fields[3] == "P"])
    dhs = BedTool(join(LOCI_DIR, "variant_analysis/HepG2_DHS.bed"))
    all_dhs = BedTool("ROOT_DIR/Projects/enhancers/data/roadmap/all_dhs_merged.bed")

    total = len(trait2bed)
    cnt = 1
    for trait, bed in trait2bed.items():

        print(cnt, total)
        cnt += 1

        trait_id = trait2id[trait]

        enh_bed = bed.intersect(hot_enhs, wa=True)
        if enh_bed.count():
            enh_bed = enh_bed.sort().merge()
            save_file = join(save_dir, "%d_enh.bed.gz" % trait_id)
            enh_bed.saveas(save_file)

        prom_bed = bed.intersect(hot_proms, wa=True)
        if prom_bed.count():
            prom_bed = prom_bed.sort().merge()
            save_file = join(save_dir, "%d_prom.bed.gz" % trait_id)
            prom_bed.saveas(save_file)

        dhs_bed = bed.intersect(dhs, wa=True)
        if dhs_bed.count():
            dhs_bed = dhs_bed.sort().merge()
            save_file = join(save_dir, "%d_dhs.bed.gz" % trait_id)
            dhs_bed.saveas(save_file)

        all_dhs_bed = bed.intersect(all_dhs, wa=True)
        if all_dhs_bed.count():
            dhs_bed = all_dhs_bed.sort().merge()
            save_file = join(save_dir, "%d_all_dhs.bed.gz" % trait_id)
            dhs_bed.saveas(save_file)


def extract_traitwise_counts():

    dhs = BedTool(join(LOCI_DIR, "variant_analysis/HepG2_DHS.bed"))
    dhs_total = np.sum([r.length for r in dhs])
    hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.PE.bed"))
    hot_enhs = BedTool([r for r in hots if r.fields[3] == "E"])
    enh_total = np.sum([r.length for r in hot_enhs])
    hot_proms = BedTool([r for r in hots if r.fields[3] == "P"])
    prom_total = np.sum([r.length for r in hot_proms])

    map_file = join(LOCI_DIR, "GWAS_analysis/trait_no_map.txt")
    id2trait = {}
    trait2id = {}
    for l in open(map_file):
        parts = l.strip().split("\t")
        id = parts[0]
        trait = parts[1]
        id2trait[id] = trait
        trait2id[trait] = id

    data = []

    _dir = join(LOCI_DIR, "GWAS_analysis/")
    for i in range(3123):

        _id = str(i + 1)

        trait_file = join(_dir, "traitwise_files/%s.bed.gz" % _id)
        prom_file = join(_dir, "traitwise_overlaps/%s_prom.bed.gz" % _id)
        enh_file = join(_dir, "traitwise_overlaps/%s_enh.bed.gz" % _id)
        dhs_file = join(_dir, "traitwise_overlaps/%s_dhs.bed.gz" % _id)
        all_dhs_file = join(_dir, "traitwise_overlaps/%s_all_dhs.bed.gz" % _id)

        if not (os.path.exists(prom_file) and os.path.exists(enh_file)):
            continue

        trait_count = BedTool(trait_file).count()
        prom_count = BedTool(prom_file).count()
        enh_count = BedTool(enh_file).count()
        dhs_count = BedTool(dhs_file).count() if os.path.exists(dhs_file) else 0
        all_dhs_count = BedTool(all_dhs_file).count() if os.path.exists(all_dhs_file) else 0

        data.append([_id, id2trait[_id], trait_count, all_dhs_count, dhs_count, prom_count, enh_count])

    df = pandas.DataFrame(data=data,
                columns=["trait_id", "trait", "trait_count", "all_dhs_count", "dhs_count", "prom_count", "enh_count"])

    return df


def plot_enriched_traits(df):

    # df = extract_fold_enrichments()
    # return df

    print(df.count())

    df = df[(df["fold_prom"]>1.2) & (df["fold_enh"]>1.2)].sort_values(by="fold_prom", ascending=False)[["trait", "fold_prom", "fold_enh"]]
    print(df.count())
    return

    traits = df["trait"].values
    traits = [_ if len(_) < 50 else _[:47]+"..." for _ in traits]
    fold_prom = df["fold_prom"].values
    fold_enh = df["fold_enh"].values

    # plt.figure(figsize=(8, 5))
    plt.figure(figsize=(5, 8))
    ax = plt.gca()

    x_range = np.arange(len(traits))
    colors = sns.color_palette("tab10")

    x_range = x_range[::-1]
    # fold_enh = fold_enh[::-1]
    # fold_prom = fold_prom[::-1]
    # traits = traits[::-1]

    label = False
    for y, x in zip(x_range, fold_enh):
        if not label:
            plt.scatter(x, y, c=colors[1], label="HOT enhancers")
            label=True
        else:
            plt.scatter(x, y, c=colors[1])
    # plt.scatter(x_range, fold_enh, c=colors[1], label="HOT enhancers")
    plt.plot(fold_enh, x_range, 'go--', c=colors[1], alpha=0.3)

    label = False
    for y, x in zip(x_range, fold_prom):
        if not label:
            plt.scatter(x, y, c=colors[2], label="HOT promoters")
            label = True
        else:
            plt.scatter(x, y, c=colors[2])
    plt.plot(fold_prom, x_range,  c=colors[2], alpha=0.3)
    ax.xaxis.tick_top()

    ax.yaxis.tick_right()

    # plt.scatter(x_range, fold_prom, c=colors[2], label="HOT promoters")
    # plt.plot(x_range, fold_prom, c=colors[2], alpha=0.3)

    # ax.set_xlabel("fold enrichment of \n GWAS SNPs", fontsize=8)

    ax.set_xticks(np.arange(1, 4))
    ax.tick_params(axis='x', which='major', labelsize=8)
    ax.set_yticks(x_range)
    ax.margins(y=0.01)
    ax.set_yticklabels(traits, fontsize=11)
    # plt.legend(fontsize=8, loc="lower right")
    ax.legend().set_visible(False)
    plt.subplots_adjust(right=0.9)

    plt.grid(axis="both", alpha=0.2)

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "GWAS_fold_enrichment_15SNPs_vertical.pdf")
    print(save_file)
    plt.savefig(save_file)


def get_GWAS_counts(counts_map=None):

    if not counts_map:

        hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.PE.bed"))
        hot_enhs = BedTool([r for r in hots if r.fields[3] == "E"])
        hot_proms = BedTool([r for r in hots if r.fields[3] == "P"])
        dhs = BedTool(join(LOCI_DIR, "variant_analysis/HepG2_DHS.bed"))
        all_dhs = BedTool("ROOT_DIR/Projects/enhancers/data/roadmap/all_dhs_merged.bed")

        counts_map = {"hots_count": hots.intersect(LD_FILE).count(),
                      "hots_len": sum([r.length for r in hots]),
                      "prom_count": hot_proms.intersect(LD_FILE).count(),
                      "prom_len": sum([r.length for r in hot_proms]),
                      "enh_count": hot_enhs.intersect(LD_FILE).count(),
                      "enh_len": sum([r.length for r in hot_enhs]),
                      "dhs_count": dhs.intersect(LD_FILE).count(),
                      "dhs_len": sum([r.length for r in dhs]),
                      "all_dhs_count": all_dhs.intersect(LD_FILE).count(),
                      "all_dhs_len": sum([r.length for r in all_dhs])}

        return counts_map

    dhs_density = counts_map["dhs_count"]/counts_map["dhs_len"]
    all_dhs_density = counts_map["all_dhs_count"]/counts_map["all_dhs_len"]
    hots_density = counts_map["hots_count"]/counts_map["hots_len"]
    prom_density = counts_map["prom_count"]/counts_map["prom_len"]
    enh_density = counts_map["enh_count"]/counts_map["enh_len"]

    dhs_p = binom_test(counts_map["dhs_count"], counts_map["all_dhs_count"], counts_map["dhs_len"]/counts_map["all_dhs_len"], alternative="greater")
    hots_p = binom_test(counts_map["hots_count"], counts_map["all_dhs_count"], counts_map["hots_len"]/counts_map["all_dhs_len"], alternative="greater")
    enh_p = binom_test(counts_map["enh_count"], counts_map["all_dhs_count"], counts_map["enh_len"]/counts_map["all_dhs_len"], alternative="greater")
    prom_p = binom_test(counts_map["prom_count"], counts_map["all_dhs_count"], counts_map["prom_len"]/counts_map["all_dhs_len"], alternative="greater")

    print("all_dhs_density", np.round(all_dhs_density*1000, 3))
    print("dhs_density", np.round(dhs_density*1000, 3), dhs_p)
    print("hots_density", np.round(hots_density*1000, 3), hots_p)
    print("prom_density", np.round(prom_density*1000, 3), prom_p)
    print("enh_density", np.round(enh_density*1000, 3), enh_p)

    return counts_map


def add_pvalues(df_orig, cm):

    df = deepcopy(df_orig)

    # df = extract_traitwise_stats()
    # cm = get_GWAS_counts()

    # df = df_orig[(df_orig["fold_prom"]>1) & (df_orig["fold_enh"]>1)].reset_index()

    dhs_prob_uni = cm["dhs_len"]/cm["all_dhs_len"]
    prom_prob_uni = cm["prom_len"]/cm["all_dhs_len"]
    enh_prob_uni = cm["enh_len"]/cm["all_dhs_len"]
    hot_prob_uni = (cm["enh_len"] + cm["prom_len"]) / cm["all_dhs_len"]

    prom_prob = cm["prom_len"] / cm["dhs_len"]
    enh_prob = cm["enh_len"] / cm["dhs_len"]
    hot_prob = (cm["enh_len"] + cm["prom_len"]) / cm["dhs_len"]

    def _get_pvalues(df, column, p):
        column_vals = df[column].to_numpy()
        trait_counts = df["trait_count"].to_numpy()
        pvalues = np.zeros(df.shape[0])
        pvalues[:] = np.nan
        ix = np.argwhere((column_vals > 0))[:, 0]
        pvalues[ix] = np.asarray([binom_test(column_vals[_], trait_counts[_], p) for _ in ix])
        return pvalues

    df["dhs_p"] = _get_pvalues(df, "dhs_count", dhs_prob_uni)
    df["enh_p"] = _get_pvalues(df, "enh_count", enh_prob_uni)
    df["prom_p"] = _get_pvalues(df, "prom_count", prom_prob_uni)

    df = df[(df["dhs_p"] < 30**-4) & (df["enh_p"] < 30**-4) & (df["prom_p"] < 30**-4)].reset_index(drop=True)

    df["prom_fold"] = (df["prom_count"]/cm["prom_len"]) / (df["dhs_count"]/cm["dhs_len"])
    df["enh_fold"] = (df["enh_count"]/cm["enh_len"]) / (df["dhs_count"]/cm["dhs_len"])

    return df

    # return a, b, dhs_prob_uni
    # binom_test(df["dhs_count"], df["trait_count"], dhs_prob_uni)
    hot_p = binom_test((df["enh_count"]+df["prom_count"]), df["trait_count"], hot_prob_uni)
    enh_p = binom_test(df["enh_count"], df["trait_count"], enh_prob_uni)
    prom_p = binom_test(df["prom_count"], df["trait_count"], prom_prob_uni)

    return dhs_p

    # dhs_pval_list = []
    prom_pval_list = []
    enh_pval_list = []



    for _, row in df.iterrows():

        print(_)

        hot_pval = binom_test((row.enh_count + row.prom_count), (row.enh_count + row.prom_count + row.dhs_count),
                              hot_prob, alternative="less")

        try:
            hot_pval = binom_test((row.enh_count + row.prom_count), (row.enh_count + row.prom_count + row.dhs_count),
                                  hot_prob, alternative="less")
        except:

            print("Error:", (row.enh_count + row.prom_count), row.dhs_count)
            print(row)
            return

        tr_hot_prob = (row.enh_count + row.prom_count) / (row.dhs_count + row.enh_count + row.prom_count)
        print("hot", hot_pval, tr_hot_prob, hot_prob)

        # # dhs_pval = binom_test(row.dhs_count, all_dhs_count, dhs_prob)
        # enh_pval = binom_test(row.enh_count, row.dhs_count, enh_prob)
        # tr_enh_prob = row.enh_count/row.dhs_count
        # print("enh", enh_pval, tr_enh_prob, enh_prob)
        # tr_prom_prom = row.prom_count/row.dhs_count
        # prom_pval = binom_test(row.prom_count, row.dhs_count, prom_prob)
        # print("prom", prom_pval, tr_prom_prom, prom_prob)

        print()
        # dhs_pval_list.append(dhs_pval)
        # prom_pval_list.append(prom_pval)
        # enh_pval_list.append(enh_pval)

        if _ == 100:
            break

    # df["dhs_pval"] = np.asarray(dhs_pval_list)
    # df["prom_pval"] = np.asarray(prom_pval_list)
    # df["enh_pval"] = np.asarray(enh_pval_list)

    return df








