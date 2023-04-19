import gzip
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
import pandas
import seaborn as sns
import pandas as pd
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
GWAS_FILE = "/panfs/pan1/devdcode/common/GWAS_v1.0.2/gwas_catalog_v1.0.2-associations_e0_r2022-11-29.tsv"
LD_FILE = "/panfs/pan1/devdcode/common/GWAS_v1.0.2/LD/plink_output/gwas.ld.gz"
LIFTOVER_CHAIN_FILE = "/panfs/pan1/devdcode/common/progs/liftover/hg38ToHg19.over.chain"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"

HEPG2_DHS_FILE = "/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF897NME_hg19.bed.gz"
HEPG2_RE_FILE = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers.bed"
HEPG2_HOT_FILE = join(LOCI_DIR, "HepG2_HOTs.bed")
HEPG2_HOT_PROM_FILE = join(LOCI_DIR, "HepG2_HOTs.proms.bed")
HEPG2_HOT_ENH_FILE = join(LOCI_DIR, "HepG2_HOTs.enhs.bed")


def load_SNP2proxy_map():

    snp2ld_snps = defaultdict(list)

    for line in gzip.open(LD_FILE, "rt"):

        if line.startswith("CHR_A"):
            continue

        parts = line.split()
        snp = parts[2]
        ld_snp = parts[5]
        if snp == ld_snp:
            continue
        snp2ld_snps[snp].append(ld_snp)

    return snp2ld_snps


def load_gwas():

    gwas = pd.read_table(GWAS_FILE)
    gwas[gwas["CONTEXT"].notnull()].reset_index(drop=True)

    return gwas


def load_snp2bed():

    snp2bed = {}
    snp2traits = defaultdict(list)
    # LD file coordinates are in hg19
    for line in gzip.open(LD_FILE, "rt"):
        if line.strip().startswith("CHR_A"):
            continue
        parts = line.split()
        chr = parts[0]
        pos1 = int(parts[1])
        pos2 = int(parts[4])
        snp1 = parts[2]
        snp2 = parts[5]
        if snp1 not in snp2bed:
            snp2bed[snp1] = "chr%s\t%d\t%d\t%s" % (chr, pos1-1, pos1, snp1)
        if snp2 not in snp2bed:
            snp2bed[snp2] = "chr%s\t%d\t%d\t%s" % (chr, pos2-1, pos2, snp2)

    gwas = load_gwas()
    bed_list = []
    for parts in gwas[["CHR_ID", "CHR_POS", "SNPS", "DISEASE/TRAIT"]].values:
        [chr, pos, snps, trait] = parts
        pos = str(pos)
        chr = str(chr)

        if " x " in snps:
            snps = snps.split(" x ")
            pos = pos.split(" x ") if " x " in pos else [pos, pos]
            chr = chr.split(" x ") if " x " in chr else [chr, chr]

            for [_chr, _pos, _snp] in zip(chr, pos, snps):
                if _snp not in snp2bed:
                    try:
                        _bed = "chr%s\t%s\t%s\t%s" % (_chr, str(int(_pos) - 1), _pos, _snp)
                    except:
                        continue
                    bed_list.append(_bed)
                    snp2traits[_snp].append(trait)
            continue

        if ";" in snps or "," in snps or "/" in snps:
            continue

        try:
            _bed = "chr%s\t%s\t%s\t%s" % (chr, str(int(pos) - 1), pos, snps)
            if snps not in snp2bed:
                bed_list.append(_bed)
            snp2traits[snps].append(trait)
        except:
            continue

    hg19_bed = BedTool("\n".join(bed_list), from_string=True).liftover(LIFTOVER_CHAIN_FILE, liftover_args="-bedPlus=4")
    for r in hg19_bed:
        snp = r.fields[-1]
        snp2bed[snp] = str(r)

    return snp2bed, snp2traits


def extract_traitswise_snps(snp2bed=None, snp2traits=None, snp2ld_snps=None):

    print("loading data")

    if not snp2bed:
        snp2bed, snp2traits = load_snp2bed()

    if not snp2ld_snps:
        snp2ld_snps = load_SNP2proxy_map()

    trait2snps = defaultdict(set)
    for snp, traits in snp2traits.items():
        for trait in traits:
            trait2snps[trait].update([snp])

    for trait, snps in trait2snps.items():
        snps.update([_ for snp in snps for _ in snp2ld_snps[snp]])
        trait2snps[trait] = snps

    print("Starting")

    _dir = join(LOCI_DIR, "GWAS_analysis/traitwise_files/")

    id2trait = {}
    cnt = 1
    for trait, snps in sorted(trait2snps.items(), key=lambda x: len(x[1]), reverse=True):

        print(cnt, len(trait2snps))

        id2trait[cnt] = trait

        trait_bed = [snp2bed[snp] for snp in snps if snp in snp2bed]
        trait_bed = BedTool("\n".join(trait_bed), from_string=True).groupby(c=4, o="distinct")
        save_file = join(_dir, "%d.bed.gz" % cnt)
        trait_bed.saveas(save_file)

        cnt += 1

    map_file = join(LOCI_DIR, "GWAS_analysis/trait_no_map.txt")
    with open(map_file, "w") as of:
        for i in range(1, cnt):
            out_line = "%d\t%d\t%s\n" % (i, len(trait2snps[id2trait[i]]), id2trait[i])
            of.write(out_line)


def extract_trait_overlap_SNPs():

    map_file = join(LOCI_DIR, "GWAS_analysis/trait_no_map.txt")
    id2trait = {}
    trait2id = {}
    for l in open(map_file):
        parts = l.strip().split("\t")
        id = int(parts[0])
        trait = parts[2]
        id2trait[id] = trait
        trait2id[trait] = id

    traits_dir = join(LOCI_DIR, "GWAS_analysis/traitwise_files/")

    hot_enhs = BedTool(HEPG2_HOT_ENH_FILE)
    hot_proms = BedTool(HEPG2_HOT_PROM_FILE)
    dhs = BedTool(HEPG2_DHS_FILE)
    all_dhs = BedTool("ROOT_DIR/Projects/enhancers/data/roadmap/all_dhs_merged.bed")
    re = BedTool(HEPG2_RE_FILE).intersect(hot_enhs, v=True)
    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True)

    header = ["trait_id", "trait", "trait_snps", "all_dhs", "dhs", "hot enh", "hot prom", "re", "rp"]
    save_file = join(LOCI_DIR, "GWAS_analysis/trait_overlaps.txt")

    with open(save_file, "w") as of:

        of.write("\t".join(header) + "\n")

        total = len(id2trait)
        for trait, id in trait2id.items():

            print(id, total)

            trait_bed = BedTool(join(traits_dir, "%d.bed.gz" % id))

            trait_count = trait_bed.count()
            hot_enhs_count = trait_bed.intersect(hot_enhs, wo=True).count()
            hot_prom_count = trait_bed.intersect(hot_proms, wo=True).count()
            dhs_count = trait_bed.intersect(dhs, wo=True).count()
            all_dhs_count = trait_bed.intersect(all_dhs, wo=True).count()
            re_count = trait_bed.intersect(re, wo=True).count()
            rp_count = trait_bed.intersect(rp, wo=True).count()

            parts = [id, trait, trait_count, all_dhs_count, dhs_count, hot_enhs_count, hot_prom_count, re_count, rp_count]

            out_line = "\t".join([str(_) for _ in parts]) + "\n"
            of.write(out_line)

            of.flush()


def get_gwas_snps():

    snp2bed, snp2traits = load_snp2bed()
    snp2ld_snps = load_SNP2proxy_map()

    gwas_snps = list(snp2traits.keys())
    gwas_snps += [_ for snp in gwas_snps for _ in snp2ld_snps[snp]]
    gwas_snps_bed = [snp2bed[_] for _ in gwas_snps if _ in snp2bed]
    gwas_snps_bed = BedTool("\n".join(gwas_snps_bed), from_string=True)

    return gwas_snps_bed


def extract_global_overlaps(snp2bed=None, snp2traits=None, snp2ld_snps=None):

    if not snp2bed:
        snp2bed, snp2traits = load_snp2bed()
        snp2ld_snps = load_SNP2proxy_map()

        return snp2bed, snp2traits, snp2ld_snps

    def extract_variant_overlap(input_regions, variant_regions):
        ovp = input_regions.intersect(variant_regions, wo=True).groupby(c=4, o="count")
        var_counts = [int(r.fields[3]) for r in ovp]
        var_counts += [0 for _ in range(input_regions.count() - ovp.count())]
        var_counts = np.asarray(var_counts)
        density = np.sum(var_counts) / np.sum([r.length for r in input_regions])
        return var_counts, density

    print("Loading beds")
    hot_enhs = BedTool(HEPG2_HOT_ENH_FILE)
    hot_proms = BedTool(HEPG2_HOT_PROM_FILE)
    dhs = BedTool(HEPG2_DHS_FILE)
    all_dhs = BedTool("ROOT_DIR/Projects/enhancers/data/roadmap/all_dhs_merged.bed")
    re = BedTool(HEPG2_RE_FILE).intersect(hot_enhs, v=True)
    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True)

    gwas_snps = list(snp2traits.keys())
    gwas_snps += [_ for snp in gwas_snps for _ in snp2ld_snps[snp]]
    gwas_snps_bed = [snp2bed[_] for _ in gwas_snps if _ in snp2bed]
    gwas_snps_bed = BedTool("\n".join(gwas_snps_bed), from_string=True)

    regions = [all_dhs, dhs, re, rp, hot_proms, hot_enhs]
    region_names = ["all_DHS", "DHS", "RE", "RP", "HOT_prom", "HOT_enh"]

    density_data = []
    counts_data = []
    print("overlaps")
    for regions, name in zip(regions, region_names):
        print(name)
        counts, _density = extract_variant_overlap(regions, gwas_snps_bed)
        density_data.append([name, _density])
        counts_data += [[name, _] for _ in counts]

    df_density = pd.DataFrame(data=density_data, columns=["region", "value"])
    df_counts = pd.DataFrame(data=counts_data, columns=["region", "value"])

    return df_density, df_counts


def plot_common_global_densities(df_density):

    save_file = join(PLOTS_DIR, "HepG2_GWAS_global.pdf")

    plt.figure(figsize=(3, 3))
    ax = plt.gca()

    x_labels = ["HOT prom", "HOT enh", "reg. prom", "reg. enh"]
    order = ["HOT_prom", "HOT_enh", "RP", "RE"]

    df_density["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=df_density, ax=ax)

    ax.grid(axis="y", alpha=0.3)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def get_traitwise_enrichments(adjust=False):

    print("loading counts")

    hot_enhs = BedTool(HEPG2_HOT_ENH_FILE)
    hot_enhs_cnt = hot_enhs.count() * 400
    hot_proms = BedTool(HEPG2_HOT_PROM_FILE)
    hot_proms_cnt = hot_proms.count() * 400
    dhs = BedTool(HEPG2_DHS_FILE)
    dhs_cnt = sum([r.length for r in dhs])
    all_dhs = BedTool("ROOT_DIR/Projects/enhancers/data/roadmap/all_dhs_merged.bed")
    all_dhs_cnt = sum([r.length for r in all_dhs])

    re = BedTool(HEPG2_RE_FILE).intersect(hot_enhs, v=True)
    re_cnt = sum([r.length for r in re])
    rp = BedTool(PROMOTERS_FILE).intersect(hot_proms, v=True)
    rp_cnt = sum([r.length for r in rp])

    print("Calculating enrichments")

    trait_stats_file = "ROOT_DIR/definitions/peak_8bp_v2/HOTs/GWAS_analysis/trait_overlaps.txt"

    df = pd.read_table(trait_stats_file)
    trait_num = df.shape[0]
    data = []

    for i, row in df.iterrows():

        target_snps = [row["dhs"],
                         row["hot enh"] + row["hot prom"],
                         row["hot enh"],
                         row["hot prom"],
                         row["re"],
                         row["rp"]]
        target_lens = [dhs_cnt,
                       hot_enhs_cnt+hot_proms_cnt,
                       hot_enhs_cnt,
                       hot_proms_cnt,
                       re_cnt,
                       rp_cnt]

        pval_list = []
        enr_list = []

        for target_snps, target_length in zip(target_snps, target_lens):

            p = target_length/all_dhs_cnt
            try:
                pval = binom_test(target_snps, row["all_dhs"], p)
                if adjust:
                    pval *= trait_num
                pval = -1 * np.log10(pval)
            except:
                pval = np.NAN

            try:
                enr = (target_snps / target_length) / (row["all_dhs"] / all_dhs_cnt)
            except:
                enr = np.NAN

            pval_list.append(pval)
            enr_list.append(enr)

        _data = [row["trait_id"], row["trait"]] + pval_list + enr_list

        data.append(_data)

    names = ["dhs", "hot", "hot_enh", "hot_prom", "re", "rp"]
    columns = ["trait_id", "trait"] + ["pval_"+_ for _ in names] + ["enr_"+_ for _ in names]

    df_pval = pd.DataFrame(data=data, columns=columns)

    return df_pval


def plot_enriched_traits(df):

    fname = "GWAS_fold_enrichment.pdf"

    # df = get_traitwise_enrichments()
    # return df

    # ind = (df["pval_hot_enh"] > 3) & (df["pval_hot_prom"] > 3) & (df["pval_dhs"] > 3)
    ind = (df["pval_hot_enh"] > 3) & (df["pval_hot_prom"] > 3)
    df = df[ind].sort_values(by="enr_hot", ascending=False).reset_index(drop=True)

    # traits = df["trait"].values
    # traits = [_ if len(_) < 40 else _[:37]+"..." for _ in traits]

    traits = ["Esophageal cancer (squamous cell)",
            "LDL cholesterol levels",
            "HDL cholesterol",
            "Sex hormone-binding globulin levels BMI",
            "Serum metabolite levels",
            "Sex hormone-binding globulin levels",
            "Red cell distribution width",
            "Non-HDL cholesterol levels",
            "Alanine aminotransferase levels",
            "LDL cholesterol",
            "Low density lipoprotein chol. levels",
            "Triglyceride levels",
            "Total cholesterol levels",
            "Total testosterone levels",
            "Hemoglobin concentration",
            "Hemoglobin",
            "Height",
            "Hematocrit",
            "Type 2 diabetes",
            "Blood protein levels",
            "Systolic blood pressure",
            "Monocyte count",
            "Diastolic blood pressure",
            "Red blood cell count"]

    data = []
    names = ["hot_prom", "rp", "hot_enh", "re", "dhs"]
    label_names = ["HOT prom", "reg. prom", "HOT enh", "reg. enh", "DHS"]

    for ind, row in df.iterrows():

        for y, name in enumerate(names):

            pvalue = row["pval_%s" % name]
            enr = row["enr_%s" % name]

            if pvalue < 3:
                continue

            _res = [ind, y, enr, min(pvalue, 30)]
            # print(name, _res)
            data.append(_res)

    df_plot = pd.DataFrame(data=data, columns=["trait", "y", "fold-enrichment", "-log10(p-value)"])

    plt.figure(figsize=(6, 5))
    ax = plt.gca()

    g = sns.scatterplot(x="trait", y="y", size="fold-enrichment", hue="-log10(p-value)", data=df_plot, ax=ax, sizes=(20, 200))

    handles, labels = g.get_legend_handles_labels()

    # print(handles)
    # print(labels)

    pval_hl = handles[1:7], labels[1:7]
    enr_hl = handles[8:], labels[8:]

    leg1 = ax.legend(*pval_hl, bbox_to_anchor=(0.05, 1.01), ncol=6, title="-log10(p-value)", loc="lower left", handletextpad=0, columnspacing=0, title_fontsize=8, fontsize=7)
    leg2 = ax.legend(*enr_hl, bbox_to_anchor=(0.51, 1.01), ncol=6, title="fold-enrichment", loc="lower left", handletextpad=0.2, columnspacing=0.5, title_fontsize=8, fontsize=7)
    ax.add_artist(leg1)

    ax.grid(axis="both", alpha=0.2, which="both")

    ax.yaxis.tick_right()
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylim([-1, 5])
    ax.set_yticklabels(label_names)

    ax.set_xticks(np.arange(len(traits)))
    ax.set_xticklabels(traits, rotation=90)

    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, fname)
    print(save_file)
    plt.savefig(save_file)
    plt.close()


def plot_enriched_traits_adjusted(df):

    fname = "GWAS_fold_enrichment_adj.pdf"

    # df = get_traitwise_enrichments(adjust=True)
    # return df

    ind = (df["pval_hot_enh"] > 3) & (df["pval_hot_prom"] > 3)
    df = df[ind].sort_values(by="enr_hot", ascending=False).reset_index(drop=True)

    # traits = df["trait"].values
    # traits = [_ if len(_) < 40 else _[:37] + "..." for _ in traits]

    traits = [
            "HDL cholesterol",
            "Sex hormone-binding globulin levels BMI",
            "Serum metabolite levels",
            "Sex hormone-binding globulin levels",
            "Low density lipoprotein chol. levels",
            "Triglyceride levels",
            "Blood protein levels"]

    data = []
    names = ["hot_prom", "rp", "hot_enh", "re", "dhs"]
    label_names = ["HOT prom", "reg. prom", "HOT enh", "reg. enh", "DHS"]

    for ind, row in df.iterrows():
        for y, name in enumerate(names):
            pvalue = row["pval_%s" % name]
            enr = row["enr_%s" % name]

            if pvalue < 3:
                continue

            _res = [ind, y, enr, min(pvalue, 30)]
            data.append(_res)

    df_plot = pd.DataFrame(data=data, columns=["trait", "y", "fold-enrichment", "-log10(p-value)"])

    plt.figure(figsize=(3.3, 5))
    ax = plt.gca()

    g = sns.scatterplot(x="trait", y="y", size="fold-enrichment", hue="-log10(p-value)", data=df_plot, ax=ax,
                        sizes=(20, 200))

    handles, labels = g.get_legend_handles_labels()

    pval_hl = handles[1:7], labels[1:7]
    enr_hl = handles[8:], labels[8:]

    leg1 = ax.legend(*pval_hl, bbox_to_anchor=(1, 1.01), ncol=1, title="-log10(adj. P)", loc="upper left",
                     handletextpad=0, columnspacing=0, title_fontsize=7, fontsize=5)
    leg2 = ax.legend(*enr_hl, bbox_to_anchor=(1, 0.55), ncol=1, title="fold-enrichment", loc="upper left",
                     handletextpad=1, columnspacing=0, labelspacing=1.5, title_fontsize=7, fontsize=5)
    ax.add_artist(leg1)

    ax.grid(axis="both", alpha=0.2, which="both")

    # ax.yaxis.tick_right()
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylim([-0.5, 4.5])
    ax.set_xlim([-0.5, len(traits) - 0.5])
    ax.set_yticklabels(label_names)

    ax.set_xticks(np.arange(len(traits)))
    ax.set_xticklabels(traits, rotation=90)

    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, fname)
    print(save_file)
    plt.savefig(save_file, bbox_inch="tight")
    plt.close()


if __name__ == "__main__":

    extract_trait_overlap_SNPs()

    # extract_traitswise_snps()





