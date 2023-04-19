import gzip
import warnings
from collections import Counter, defaultdict
from copy import deepcopy
import pandas
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

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

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom" % (x, i)) for i in range(14)]

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PROMOTERS_FILE = "/panfs/pan1/devdcode/common/genomes/hg19/annotations/genes/promoters.merged.bed.gz"
SNPS_FILE="/panfs/pan1/devdcode/common/gnomAD/v2_hg19/common_SNPs/gnomad.genomes.r2.1.1.sites.MAF_5.PASS.SNVs.bed.gz"
INDELS_FILE="/panfs/pan1/devdcode/common/gnomAD/v2_hg19/common_SNPs/gnomad.genomes.r2.1.1.sites.MAF_5.PASS.INDELs.bed.gz"
CLINVAR_FILE="/panfs/pan1/devdcode/common/ClinVar/clinvar_20210404.6f.bed"

HEPG2_DHS_FILE = "/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF897NME_hg19.bed.gz"
HEPG2_RE_FILE = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers.bed"
HEPG2_HOT_FILE = join(LOCI_DIR, "HepG2_HOTs.bed")
HEPG2_HOT_PROM_FILE = join(LOCI_DIR, "HepG2_HOTs.proms.bed")
HEPG2_HOT_ENH_FILE = join(LOCI_DIR, "HepG2_HOTs.enhs.bed")
CAQTL_FILE = "ROOT_DIR/caQTLs/HepG2_caQTLs.bed"
EQTL_FILE = "ROOT_DIR/eQTLs/Liver.eQTLs.bed.gz"
RAQTL_FILE = "ROOT_DIR/raQLTs/hepg2.sign.id.LP190708.bed"
LIFTOVER_CHAIN_FILE = "/panfs/pan1/devdcode/common/progs/liftover/hg38ToHg19.over.chain"

from overbinders.data_prep.basic import load_metadata


def extract_common_SNPs():

    snps_file = "/panfs/pan1/devdcode/common/gnomAD/gnomad.genomes.r2.1.1.sites.MAF_5.PASS.SNVs.bed.gz"
    indels_file = "/panfs/pan1/devdcode/common/gnomAD/gnomad.genomes.r2.1.1.sites.MAF_5.PASS.INDELs.bed.gz"

    cl = "HepG2"

    hots_file = join(LOCI_DIR, "%s_HOTs.PE.bed" % cl)
    save_file_snps = join(LOCI_DIR, "variant_analysis/%s_HOTs.PE.SNPs.bed" % cl)
    save_file_indels = join(LOCI_DIR, "variant_analysis/%s_HOTs.PE.INDELs.bed" % cl)

    print("INDELs")
    BedTool(hots_file). \
        intersect(indels_file, wo=True). \
        groupby(g=[1, 2, 3, 4], c=[9, 9], o=["count", "collapse"]). \
        sort(). \
        saveas(save_file_indels)

    print("SNPs")
    BedTool(hots_file).\
        intersect(snps_file, wo=True).\
        groupby(g=[1, 2, 3, 4], c=[9, 9], o=["count", "collapse"]).\
        sort().\
        saveas(save_file_snps)


def extract_common_eQTLs():

    liver_eqtls = "ROOT_DIR/eQTLs/Liver.eQTLs.bed.gz"
    blood_eqtls = "ROOT_DIR/eQTLs/Whole_Blood.eQTLs.bed.gz"

    cl = "HepG2"
    cl = "K562"

    eqtl_file = liver_eqtls if cl == "HepG2" else blood_eqtls

    hots_file = join(LOCI_DIR, "%s_HOTs.PE.bed" % cl)
    save_file = join(LOCI_DIR, "%s_HOTs.PE.eQTLs.bed" % cl)

    print("eQTLs:", cl)
    BedTool(hots_file). \
        intersect(eqtl_file, wo=True). \
        groupby(g=[1, 2, 3, 4], c=[8, 8], o=["count", "collapse"]). \
        sort(). \
        saveas(save_file)


def extract_variant_overlap(input_regions, variant_regions):

    ovp = input_regions.intersect(variant_regions, wo=True).groupby(c=4, o="count")

    var_counts = [int(r.fields[3]) for r in ovp]
    var_counts += [0 for _ in range(input_regions.count() - ovp.count())]
    var_counts = np.asarray(var_counts)
    density = np.sum(var_counts)/np.sum([r.length for r in input_regions])

    return var_counts, density


def load_common_variants_data_density_otf():

    snps = BedTool(SNPS_FILE)
    indels = BedTool(INDELS_FILE)
    caqtls = BedTool(CAQTL_FILE)
    eqtls = BedTool(EQTL_FILE)
    # eqtls = load_eQTL_bed()
    print("eQTL v7")
    raqtls = BedTool(RAQTL_FILE)
    clinvars = BedTool(CLINVAR_FILE)
    from GWAS_analysis_v2 import get_gwas_snps
    gwas = get_gwas_snps()

    dhs = BedTool(HEPG2_DHS_FILE)
    hots = BedTool(HEPG2_HOT_FILE)
    hot_proms = BedTool(HEPG2_HOT_PROM_FILE)
    hot_enhs = BedTool(HEPG2_HOT_ENH_FILE)

    re = BedTool(HEPG2_RE_FILE).intersect(hots, v=True)
    rp = BedTool(PROMOTERS_FILE).intersect(hots, v=True)

    regions = [dhs, re, rp, hots, hot_proms, hot_enhs]
    region_names = ["DHS", "RE", "RP", "HOT", "HOT_prom", "HOT_enh"]

    variants = [snps, indels, caqtls, eqtls, raqtls, clinvars, gwas]
    variant_names = ["SNPs", "INDELs", "caQTLs", "eQTLs", "raQTLs", "ClinVar", "GWAS"]

    density_data = []
    counts_data = []
    for regions, name in zip(regions, region_names):
        print(name)
        for vars, var_name in zip(variants, variant_names):

            counts, _density = extract_variant_overlap(regions, vars)
            print(name, var_name, _density)
            density_data.append([name, var_name, _density])
            counts_data += [[name, var_name, _] for _ in counts]

    df_density = pd.DataFrame(data=density_data, columns=["region", "variant", "value"])
    df_counts = pd.DataFrame(data=counts_data, columns=["region", "variant", "value"])

    return df_density, df_counts


def load_eQTLs_data_slopes():

    eqtls = BedTool(EQTL_FILE)

    dhs = BedTool(HEPG2_DHS_FILE)
    hots = BedTool(HEPG2_HOT_FILE)
    hot_proms = BedTool(HEPG2_HOT_PROM_FILE)
    hot_enhs = BedTool(HEPG2_HOT_ENH_FILE)

    re = BedTool(HEPG2_RE_FILE).intersect(hots, v=True)
    rp = BedTool(PROMOTERS_FILE).intersect(hots, v=True)

    regions = [dhs, re, rp, hots, hot_proms, hot_enhs]
    region_names = ["DHS", "RE", "RP", "HOT", "HOT_prom", "HOT_enh"]

    data = []

    for regions, name in zip(regions, region_names):
        print(name)
        data += [[name, float(r.fields[3]), float(r.fields[4])] for r in eqtls.intersect(regions, wa=True)]

    df = pd.DataFrame(data=data, columns=["region", "slope", "p-value"]).reset_index(drop=True)
    df["p-value"] = -1*np.log10(df["p-value"])
    df["abs_slope"] = abs(df["slope"])

    return df


def plot_eQTL_slopes(df):

    # df = load_eQTLs_data_slopes()

    save_file = join(PLOTS_DIR, "HepG2_eQTL_slopes_barplot.pdf")

    fig, axs = plt.subplots(1, 3, figsize=(7, 3))

    order = ["HOT_prom", "RP", "HOT_enh", "RE"]
    x_labels = ["HOT prom", "reg. prom", "HOT enh", "reg. enh"]

    sns.barplot(x="region", y="p-value", data=df, ax=axs[1], order=order)
    axs[1].set_xlabel("")
    axs[1].set_ylabel("slope")
    axs[1].set_title("eQTL effect sizes")
    axs[1].set_xticks(np.arange(len(x_labels)), x_labels, rotation=45, ha="right")
    axs[1].set_ylim([0, 0.1])
    axs[1].grid(axis="y", alpha=0.3)

    # _df = df[df["metric"] == "slopes"]
    # _df["value"] = np.abs(_df["value"])
    #
    # sns.barplot(x="category", y="value", data=_df, ax=axs[2])
    # axs[2].set_xlabel("")
    # axs[2].set_ylabel("abs(slope)")
    # axs[2].set_title("eQTL effect sizes")
    # axs[2].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    # axs[2].set_ylim([0, 0.7])
    # axs[2].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)

    # print("Counts")
    # df_t = df[df["metric"]=="counts"]
    # values = [df_t[(df_t["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")
    #
    # print("Slopes")
    # df_t = df[df["metric"] == "slopes"]
    # values = [df_t[(df_t["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")
    #
    # print("Slopes abs")
    #
    # values = [_df[(_df["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")

    """
    Counts
    RE - HOT enh
    MannwhitneyuResult(statistic=215065104.0, pvalue=0.00994432902385867)
    RE - HOT prom
    MannwhitneyuResult(statistic=213586125.0, pvalue=2.2219806721520244e-135)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=139185698.5, pvalue=5.684495189884156e-86)
    
    Slopes
    RE - HOT enh
    MannwhitneyuResult(statistic=6327900.0, pvalue=0.3276121301487832)
    RE - HOT prom
    MannwhitneyuResult(statistic=12446159.5, pvalue=0.013241044066134673)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=8915177.0, pvalue=0.07618516217789674)
    
    Slopes abs
    RE - HOT enh
    MannwhitneyuResult(statistic=6071384.5, pvalue=0.0003722676931697666)
    RE - HOT prom
    MannwhitneyuResult(statistic=12523192.5, pvalue=0.04474305085429737)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=8487941.5, pvalue=1.867644240287438e-07)
    """


def plot_eQTL_slopes_pos_neg(df):

    # df = load_eQTLs_data_slopes()

    save_file = join(PLOTS_DIR, "HepG2_eQTL_slopes_barplot_pos_neg.pdf")

    fig, axs = plt.subplots(1, 4, figsize=(8, 3))

    order = ["RE", "HOT enh", "HOT prom"]
    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="value", data=df[df["metric"] == "counts"], ax=axs[0], order=order)
    axs[0].set_xlabel("")
    axs[0].set_ylabel("variants/locus")
    axs[0].set_title("eQTLs")
    axs[0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0].grid(axis="y", alpha=0.3)
    axs[0].set_ylim([0, 0.45])

    _df = df[df["metric"] == "slopes"]
    _df["value"] = np.abs(_df["value"])

    g = sns.barplot(x="category", y="value", data=_df, ax=axs[1], order=order)
    g.bar_label(g.containers[0], fmt="%.3f", padding=5)
    axs[1].set_xlabel("")
    axs[1].set_ylabel("abs(slope)")
    axs[1].set_title("all")
    axs[1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1].set_ylim([0, 0.7])
    axs[1].grid(axis="y", alpha=0.3)

    _df = df[(df["metric"] == "slopes") & (df["value"] > 0)][["category", "value"]]

    g = sns.barplot(x="category", y="value", data=_df, ax=axs[2], order=order)
    g.bar_label(g.containers[0], fmt="%.3f", padding=5)
    axs[2].set_xlabel("")
    axs[2].set_ylabel("slope")
    axs[2].set_title("positives")
    axs[2].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[2].set_ylim([0, 0.75])
    axs[2].grid(axis="y", alpha=0.3)

    _df = df[(df["metric"] == "slopes") & (df["value"] < 0)][["category", "value"]]
    _df["value"] = np.abs(_df["value"])
    g = sns.barplot(x="category", y="value", data=_df, ax=axs[3], order=order)
    g.bar_label(g.containers[0], fmt="%.3f", padding=5)
    axs[3].set_xlabel("")
    axs[3].set_ylabel("abs(slope)")
    axs[3].set_title("negatives")
    axs[3].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[3].set_ylim([0, 0.6])
    axs[3].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)

    # print("Counts")
    # df_t = df[df["metric"]=="counts"]
    # values = [df_t[(df_t["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")
    #
    # print("Slopes")
    # df_t = df[df["metric"] == "slopes"]
    # values = [df_t[(df_t["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")
    #
    # print("Slopes abs")
    #
    # values = [_df[(_df["category"] == _type)]["value"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    # print("RE - HOT enh")
    # print(mannwhitneyu(values[0], values[1]))
    # print("RE - HOT prom")
    # print(mannwhitneyu(values[0], values[2]))
    # print("HOT enh - HOT prom")
    # print(mannwhitneyu(values[1], values[2]))
    # print("")

    """
    Counts
    RE - HOT enh
    MannwhitneyuResult(statistic=215065104.0, pvalue=0.00994432902385867)
    RE - HOT prom
    MannwhitneyuResult(statistic=213586125.0, pvalue=2.2219806721520244e-135)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=139185698.5, pvalue=5.684495189884156e-86)

    Slopes
    RE - HOT enh
    MannwhitneyuResult(statistic=6327900.0, pvalue=0.3276121301487832)
    RE - HOT prom
    MannwhitneyuResult(statistic=12446159.5, pvalue=0.013241044066134673)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=8915177.0, pvalue=0.07618516217789674)
    
    Slopes abs
    RE - HOT enh
    MannwhitneyuResult(statistic=6071384.5, pvalue=0.0003722676931697666)
    RE - HOT prom
    MannwhitneyuResult(statistic=12523192.5, pvalue=0.04474305085429737)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=8487941.5, pvalue=1.867644240287438e-07)
    """


def plot_common_SNPs_INDELs(df):

    # df = load_common_variants_data()
    # return df

    # for _var in ["INDELs", "SNPs"]:
    #     print(_var)
    #     values = [df[(df["variant"] == _var) & (df["category"] == _type)]["density"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    #     print("RE - HOT enh")
    #     print(mannwhitneyu(values[0], values[1]))
    #     print("RE - HOT prom")
    #     print(mannwhitneyu(values[0], values[2]))
    #     print("HOT enh - HOT prom")
    #     print(mannwhitneyu(values[1], values[2]))
    #     print("")
    """
    INDELs
    RE - HOT enh
    MannwhitneyuResult(statistic=214966530.5, pvalue=0.022695190947905713)
    RE - HOT prom
    MannwhitneyuResult(statistic=226892747.0, pvalue=2.887459071364543e-07)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=147859075.0, pvalue=0.003456518681646981)

    SNPs
    RE - HOT enh
    MannwhitneyuResult(statistic=207584159.5, pvalue=1.1576407467242702e-14)
    RE - HOT prom
    MannwhitneyuResult(statistic=284556722.0, pvalue=1.2012872896843295e-273)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=175818718.0, pvalue=0.0)

    """

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_barplot.pdf")

    fig, axs = plt.subplots(1, 2, figsize=(5, 3))

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="density", data=df[df["variant"] == "INDELs"], ax=axs[0])
    axs[0].set_xlabel("")
    axs[0].set_ylabel("variants/locus")
    axs[0].set_title("INDELs")
    axs[0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0].grid(axis="y", alpha=0.3)
    axs[0].set_ylim([0, 0.27])

    sns.barplot(x="category", y="density", data=df[df["variant"] == "SNPs"], ax=axs[1])
    axs[1].set_xlabel("")
    axs[1].set_ylabel("variants/locus")
    axs[1].set_title("SNPs")
    axs[1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1].set_ylim([0, 1.25])
    axs[1].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)

    print("Average values:\n")
    print(df.groupby(by=["variant", "category"]).mean())


def plot_common_SNPs_INDELs_eQTLs(df):

    # df = load_common_variants_data()
    # return df

    # for _var in ["INDELs", "SNPs", "eQTLs"]:
    #     print(_var)
    #     values = [df[(df["variant"] == _var) & (df["category"] == _type)]["density"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    #     print("RE - HOT enh")
    #     print(mannwhitneyu(values[0], values[1]))
    #     print("RE - HOT prom")
    #     print(mannwhitneyu(values[0], values[2]))
    #     print("HOT enh - HOT prom")
    #     print(mannwhitneyu(values[1], values[2]))
    #     print("")
    # return
    """
    INDELs
    RE - HOT enh
    MannwhitneyuResult(statistic=214966530.5, pvalue=0.022695190947905713)
    RE - HOT prom
    MannwhitneyuResult(statistic=226892747.0, pvalue=2.887459071364543e-07)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=147859075.0, pvalue=0.003456518681646981)
    
    SNPs
    RE - HOT enh
    MannwhitneyuResult(statistic=207584159.5, pvalue=1.1576407467242702e-14)
    RE - HOT prom
    MannwhitneyuResult(statistic=228291954.0, pvalue=0.028638531531496033)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=144656819.5, pvalue=5.288398817698053e-08)
    
    eQTLs
    RE - HOT enh
    MannwhitneyuResult(statistic=215065104.0, pvalue=0.00994432902385867)
    RE - HOT prom
    MannwhitneyuResult(statistic=213586125.0, pvalue=2.2219806721520244e-135)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=139185698.5, pvalue=5.684495189884156e-86)     

    """

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_barplot.pdf")

    fig, axs = plt.subplots(1, 3, figsize=(6, 3))

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="density", data=df[df["variant"] == "INDELs"], ax=axs[0])
    axs[0].set_xlabel("")
    axs[0].set_ylabel("variants/locus")
    axs[0].set_title("INDELs")
    axs[0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0].grid(axis="y", alpha=0.3)
    axs[0].set_ylim([0, 0.27])

    sns.barplot(x="category", y="density", data=df[df["variant"] == "SNPs"], ax=axs[1])
    axs[1].set_xlabel("")
    axs[1].set_ylabel("variants/locus")
    axs[1].set_title("SNPs")
    axs[1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1].set_ylim([0, 1.25])
    axs[1].grid(axis="y", alpha=0.3)

    sns.barplot(x="category", y="density", data=df[df["variant"] == "eQTLs"], ax=axs[2])
    axs[2].set_xlabel("")
    axs[2].set_ylabel("variants/locus")
    axs[2].set_title("eQTLs")
    axs[2].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[2].set_ylim([0, 0.45])
    axs[2].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)

    print("Average values:\n")
    print(df.groupby(by=["variant", "category"]).mean())


def plot_clinvar(df):

    save_file = join(PLOTS_DIR, "HepG2_clinvar_barplot.pdf")

    plt.figure(figsize=(2, 2))
    ax = plt.gca()

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="value", data=df, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("variants/locus")
    # ax.set_title("INDELs")
    ax.set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.3)
    # ax.set_ylim([0, 0.27])

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)


def plot_common_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoFs(df_variants, df_eqtls, df_ca, df_ra, df_lofs):

    # df_variants = load_common_variants_data()
    # df_eqtls = load_eQTLs_data_slopes()
    # df_ca = load_caQTLs_data()
    # df_ra = load_raQTLs_data()
    # df_lofs, _ = load_gnomAD_LoF_data()
    # return df_variants, df_eqtls, df_ca, df_ra, df_lofs

    # for _var in ["INDELs", "SNPs", "eQTLs"]:
    #     print(_var)
    #     values = [df[(df["variant"] == _var) & (df["category"] == _type)]["density"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    #     print("RE - HOT enh")
    #     print(mannwhitneyu(values[0], values[1]))
    #     print("RE - HOT prom")
    #     print(mannwhitneyu(values[0], values[2]))
    #     print("HOT enh - HOT prom")
    #     print(mannwhitneyu(values[1], values[2]))
    #     print("")
    # return
    """

    """

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoF_barplot_2x3.pdf")

    # fig, axs = plt.subplots(1, 6, figsize=(11, 3))
    fig, axs = plt.subplots(3, 2, figsize=(4, 7))

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="density", data=df_variants[df_variants["variant"] == "INDELs"], ax=axs[0, 0])
    axs[0, 0].set_xlabel("")
    axs[0, 0].set_ylabel("variants/locus")
    axs[0, 0].set_title("INDELs")
    # axs[0, 0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0, 0].set_xticks([])
    axs[0, 0].grid(axis="y", alpha=0.3)
    axs[0, 0].set_ylim([0, 0.27])

    sns.barplot(x="category", y="density", data=df_variants[df_variants["variant"] == "SNPs"], ax=axs[0, 1])
    axs[0, 1].set_xlabel("")
    # axs[0, 1].set_ylabel("variants/locus")
    axs[0, 1].set_ylabel("")
    axs[0, 1].set_title("SNPs")
    # axs[0, 1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0, 1].set_xticks([])
    axs[0, 1].set_ylim([0, 1.25])
    axs[0, 1].grid(axis="y", alpha=0.3)

    _df = df_eqtls[df_eqtls["metric"] == "counts"]
    sns.barplot(x="category", y="value", data=_df, ax=axs[1, 0])
    axs[1, 0].set_xlabel("")
    axs[1, 0].set_ylabel("variants/locus")
    axs[1, 0].set_title("eQTLs")
    # axs[1, 0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1, 0].set_xticks([])
    axs[1, 0].set_ylim([0, 0.45])
    axs[1, 0].grid(axis="y", alpha=0.3)

    sns.barplot(x="category", y="value", data=df_lofs, ax=axs[1, 1])
    axs[1, 1].set_xlabel("")
    # axs[1, 1].set_ylabel("variants/locus")
    axs[1, 1].set_ylabel("")
    axs[1, 1].set_title("LoF")
    axs[1, 1].set_xticks([])
    axs[1, 1].set_ylim([0, 0.73])
    axs[1, 1].grid(axis="y", alpha=0.3)

    sns.barplot(x="category", y="value", data=df_ra, ax=axs[2, 0])
    axs[2, 0].set_xlabel("")
    axs[2, 0].set_ylabel("variants/locus")
    axs[2, 0].set_title("raQTLs")
    axs[2, 0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[2, 0].set_ylim([0, 0.03])
    axs[2, 0].grid(axis="y", alpha=0.3)

    sns.barplot(x="category", y="value", data=df_ca, ax=axs[2, 1])
    axs[2, 1].set_xlabel("")
    # axs[2, 1].set_ylabel("variants/locus")
    axs[2, 1].set_ylabel("")
    axs[2, 1].set_title("caQTLs")
    axs[2, 1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[2, 1].set_ylim([0, 0.01])
    axs[2, 1].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)


def plot_common_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoFs_density_v1(df_variants, df_eqtls, df_caqtls, df_raqtls, df_clinvar):

    df_variants, _ = load_common_variants_data_density()
    # df_eqtls = load_eQTLs_data_slopes()
    # df_clinvar = load_raQTLs_data(var="ClinVar", density=True)
    # df_raqtls = load_raQTLs_data(density=True)
    # df_caqtls = load_caQTLs_data(density=True)
    # return df_variants, df_eqtls, df_caqtls, df_raqtls, df_clinvar

    # for _var in ["INDELs", "SNPs", "eQTLs"]:
    #     print(_var)
    #     values = [df[(df["variant"] == _var) & (df["category"] == _type)]["density"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    #     print("RE - HOT enh")
    #     print(mannwhitneyu(values[0], values[1]))
    #     print("RE - HOT prom")
    #     print(mannwhitneyu(values[0], values[2]))
    #     print("HOT enh - HOT prom")
    #     print(mannwhitneyu(values[1], values[2]))
    #     print("")
    # return

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoF_barplot_2x3_density.pdf")

    fig, axs = plt.subplots(3, 2, figsize=(6, 7))

    x_labels = ["HOT prom", "HOT enh", "reg. prom", "reg. enh"]
    order = ["HOT prom", "HOT enh", "RP", "RE"]

    _df = df_variants[df_variants["variant"] == "INDELs"]
    _df["density"] *= 10000
    sns.barplot(x="category", order=order, y="density", data=_df, ax=axs[0, 0])
    axs[0, 0].set_xlabel("")
    axs[0, 0].set_ylabel("var/10kb", fontsize=9)
    axs[0, 0].set_title("INDELs")
    axs[0, 0].set_xticks([])
    axs[0, 0].grid(axis="y", alpha=0.3)
    # axs[0, 0].set_ylim([0, 0.27])

    _df = df_variants[df_variants["variant"] == "SNPs"]
    _df["density"] *= 1000
    sns.barplot(x="category", order=order, y="density", data=_df, ax=axs[0, 1])
    axs[0, 1].set_xlabel("")
    axs[0, 1].set_ylabel("var/1kb", fontsize=9)
    axs[0, 1].set_title("SNPs")
    axs[0, 1].set_xticks([])
    axs[0, 1].grid(axis="y", alpha=0.3)

    _df = df_eqtls[df_eqtls["metric"] == "density"]
    _df["value"] *= 10000
    sns.barplot(x="category", order=order, y="value", data=_df, ax=axs[1, 0])
    axs[1, 0].set_xlabel("")
    axs[1, 0].set_ylabel("var/10kb", fontsize=9)
    axs[1, 0].set_title("eQTLs")
    axs[1, 0].set_xticks([])
    axs[1, 0].grid(axis="y", alpha=0.3)

    _df = deepcopy(df_caqtls)
    _df["value"] *= 10000
    sns.barplot(x="category", y="value", order=order, data=_df, ax=axs[1, 1])
    axs[1, 1].set_xlabel("")
    axs[1, 1].set_ylabel("var/10kb")
    axs[1, 1].set_title("caQTLs")
    axs[1, 1].set_xticks([])
    axs[1, 1].grid(axis="y", alpha=0.3)

    _df = deepcopy(df_raqtls)
    _df["value"] *= 10000
    sns.barplot(x="category", y="value", order=order, data=_df, ax=axs[2, 0])
    axs[2, 0].set_xlabel("")
    axs[2, 0].set_ylabel("var/10kb")
    axs[2, 0].set_title("raQTLs")
    axs[2, 0].grid(axis="y", alpha=0.3)
    axs[2, 0].set_xticklabels(x_labels, rotation=45, ha="right")

    _df = deepcopy(df_clinvar)
    _df["value"] *= 10000
    sns.barplot(x="category", order=order, y="value", data=_df, ax=axs[2, 1])
    axs[2, 1].set_xlabel("")
    axs[2, 1].set_ylabel("var/10kb")
    axs[2, 1].set_title("ClinVar")
    # axs[2, 1].set_xticks([])
    axs[2, 1].grid(axis="y", alpha=0.3)
    axs[2, 1].set_xticklabels(x_labels, rotation=45, ha="right")

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_common_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoFs_density(df_density, df_counts):

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_caQTLs_raQTLs_LoF_barplot_2x3_density_v2.pdf")

    fig, axs = plt.subplots(3, 2, figsize=(6, 7))

    x_labels = ["HOT prom", "HOT enh", "reg. prom", "reg. enh"]
    order = ["HOT_prom", "HOT_enh", "RP", "RE"]

    _df = df_density[df_density["variant"] == "INDELs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 0])
    axs[0, 0].set_xlabel("")
    axs[0, 0].set_ylabel("var/10kb", fontsize=9)
    axs[0, 0].set_title("INDELs")
    axs[0, 0].set_xticks([])
    axs[0, 0].grid(axis="y", alpha=0.3)
    # axs[0, 0].set_ylim([0, 0.27])

    _df = df_density[df_density["variant"] == "SNPs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 1])
    axs[0, 1].set_xlabel("")
    axs[0, 1].set_ylabel("var/10kb", fontsize=9)
    axs[0, 1].set_title("SNPs")
    axs[0, 1].set_xticks([])
    axs[0, 1].grid(axis="y", alpha=0.3)

    _df = df_density[df_density["variant"] == "eQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[1, 0])
    axs[1, 0].set_xlabel("")
    axs[1, 0].set_ylabel("var/10kb", fontsize=9)
    axs[1, 0].set_title("eQTLs")
    axs[1, 0].set_xticks([])
    axs[1, 0].grid(axis="y", alpha=0.3)

    _df = df_density[df_density["variant"] == "caQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[1, 1])
    axs[1, 1].set_xlabel("")
    axs[1, 1].set_ylabel("var/10kb")
    axs[1, 1].set_title("caQTLs")
    axs[1, 1].set_xticks([])
    axs[1, 1].grid(axis="y", alpha=0.3)

    _df = df_density[df_density["variant"] == "raQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[2, 0])
    axs[2, 0].set_xlabel("")
    axs[2, 0].set_ylabel("var/10kb")
    axs[2, 0].set_title("raQTLs")
    axs[2, 0].grid(axis="y", alpha=0.3)
    axs[2, 0].set_xticklabels(x_labels, rotation=45, ha="right")

    _df = df_density[df_density["variant"] == "ClinVar"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[2, 1])
    axs[2, 1].set_xlabel("")
    axs[2, 1].set_ylabel("var/10kb")
    axs[2, 1].set_title("ClinVar")
    # axs[2, 1].set_xticks([])
    axs[2, 1].grid(axis="y", alpha=0.3)
    axs[2, 1].set_xticklabels(x_labels, rotation=45, ha="right")

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_common_SNPs_INDELs_eQTLs_caQTLs_raQTLs_GWAS_density(df_density):

    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_caQTLs_raQTLs_GWAS_barplot.pdf")

    fig, axs = plt.subplots(3, 2, figsize=(6, 7))

    x_labels = ["HOT prom", "reg. prom", "HOT enh", "reg. enh"]
    order = ["HOT_prom", "RP", "HOT_enh", "RE"]

    _df = df_density[df_density["variant"] == "INDELs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 0])
    axs[0, 0].set_xlabel("")
    axs[0, 0].set_ylabel("var/10kb", fontsize=9)
    axs[0, 0].set_title("INDELs")
    axs[0, 0].set_xticks([])
    axs[0, 0].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[0, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "SNPs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[0, 1])
    axs[0, 1].set_xlabel("")
    axs[0, 1].set_ylabel("var/10kb", fontsize=9)
    axs[0, 1].set_title("SNPs")
    axs[0, 1].set_xticks([])
    axs[0, 1].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[0, 1].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "eQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[1, 0])
    axs[1, 0].set_xlabel("")
    axs[1, 0].set_ylabel("var/10kb", fontsize=9)
    axs[1, 0].set_title("eQTLs")
    axs[1, 0].set_xticks([])
    axs[1, 0].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[1, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "caQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[1, 1])
    axs[1, 1].set_xlabel("")
    axs[1, 1].set_ylabel("var/10kb")
    axs[1, 1].set_title("caQTLs")
    axs[1, 1].set_xticks([])
    axs[1, 1].grid(axis="y", alpha=0.3)
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[1, 1].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "raQTLs"]
    _df["value"] *= 10000
    sns.barplot(x="region", y="value", order=order, data=_df, ax=axs[2, 0])
    axs[2, 0].set_xlabel("")
    axs[2, 0].set_ylabel("var/10kb")
    axs[2, 0].set_title("raQTLs")
    axs[2, 0].grid(axis="y", alpha=0.3)
    axs[2, 0].set_xticklabels(x_labels, rotation=45, ha="right")
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[2, 0].set_ylim([0, y_upper])

    _df = df_density[df_density["variant"] == "GWAS"]
    _df["value"] *= 10000
    sns.barplot(x="region", order=order, y="value", data=_df, ax=axs[2, 1])
    axs[2, 1].set_xlabel("")
    axs[2, 1].set_ylabel("var/10kb")
    axs[2, 1].set_title("GWAS")
    axs[2, 1].grid(axis="y", alpha=0.3)
    axs[2, 1].set_xticklabels(x_labels, rotation=45, ha="right")
    y_upper = _df["value"].values.max()
    y_upper += 0.4 * y_upper
    axs[2, 1].set_ylim([0, y_upper])

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close()


def plot_common_SNPs_INDELs_eQTLs_slopes(df_variants, df_eqtls):

    # df_variants = load_common_variants_data()
    # df_eqtls = load_eQTLs_data_slopes()
    # return df_variants, df_eqtls

    # for _var in ["INDELs", "SNPs", "eQTLs"]:
    #     print(_var)
    #     values = [df[(df["variant"] == _var) & (df["category"] == _type)]["density"].values for _type in ["RE", "HOT enh", "HOT prom"]]
    #     print("RE - HOT enh")
    #     print(mannwhitneyu(values[0], values[1]))
    #     print("RE - HOT prom")
    #     print(mannwhitneyu(values[0], values[2]))
    #     print("HOT enh - HOT prom")
    #     print(mannwhitneyu(values[1], values[2]))
    #     print("")
    # return
    """
    INDELs
    RE - HOT enh
    MannwhitneyuResult(statistic=214966530.5, pvalue=0.022695190947905713)
    RE - HOT prom
    MannwhitneyuResult(statistic=226892747.0, pvalue=2.887459071364543e-07)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=147859075.0, pvalue=0.003456518681646981)

    SNPs
    RE - HOT enh
    MannwhitneyuResult(statistic=207584159.5, pvalue=1.1576407467242702e-14)
    RE - HOT prom
    MannwhitneyuResult(statistic=228291954.0, pvalue=0.028638531531496033)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=144656819.5, pvalue=5.288398817698053e-08)

    eQTLs
    RE - HOT enh
    MannwhitneyuResult(statistic=215065104.0, pvalue=0.00994432902385867)
    RE - HOT prom
    MannwhitneyuResult(statistic=213586125.0, pvalue=2.2219806721520244e-135)
    HOT enh - HOT prom
    MannwhitneyuResult(statistic=139185698.5, pvalue=5.684495189884156e-86)

    """

    # save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_eQTLs_barplot.pdf")

    fig, axs = plt.subplots(1, 4, figsize=(7, 3))

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="density", data=df_variants[df_variants["variant"] == "INDELs"], ax=axs[0])
    axs[0].set_xlabel("")
    axs[0].set_ylabel("variants/locus")
    axs[0].set_title("INDELs")
    axs[0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0].grid(axis="y", alpha=0.3)
    axs[0].set_ylim([0, 0.27])

    sns.barplot(x="category", y="density", data=df_variants[df_variants["variant"] == "SNPs"], ax=axs[1])
    axs[1].set_xlabel("")
    axs[1].set_ylabel("variants/locus")
    axs[1].set_title("SNPs")
    axs[1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1].set_ylim([0, 1.25])
    axs[1].grid(axis="y", alpha=0.3)

    _df = df_eqtls[df_eqtls["metric"] == "counts"]
    sns.barplot(x="category", y="value", data=_df, ax=axs[2])
    axs[2].set_xlabel("")
    axs[2].set_ylabel("variants/locus")
    axs[2].set_title("eQTLs")
    axs[2].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[2].set_ylim([0, 0.45])
    axs[2].grid(axis="y", alpha=0.3)

    _df = df_eqtls[df_eqtls["metric"] == "slopes"]
    sns.barplot(x="category", y="value", data=_df, ax=axs[3])
    axs[3].set_xlabel("")
    axs[3].set_ylabel("abs(eQTL slopes)")
    axs[3].set_title("eQTLs")
    axs[3].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    # axs[3].set_ylim([0, 0.45])
    axs[3].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    print(save_file)
    plt.savefig(save_file)

    print("Average values:\n")
    print(df.groupby(by=["variant", "category"]).mean())


def plot_common_SNPs_INDELs_fold(df):

    # df = load_common_variants_data_density
    # return df

    dhs_snps_density = float(df[(df["category"] == "DHS") & (df["variant"] == "SNPs")]["density"])
    dhs_indels_density = float(df[(df["category"] == "DHS") & (df["variant"] == "INDELs")]["density"])

    data=[]
    for cat in ["RE", "HOT enh", "HOT prom"]:
        for _var, factor in zip(["SNPs", "INDELs"], [dhs_snps_density, dhs_indels_density]):
            density = float(df[(df["category"] == cat) & (df["variant"] == _var)]["density"])
            fold = density / factor
            data.append([cat, _var, fold])

    df = pandas.DataFrame(data=data, columns=["category", "variant", "fold"])

    fig, axs = plt.subplots(1, 2, figsize=(5, 3))

    x_labels = ["enhancers", "HOT enh", "HOT prom"]

    sns.barplot(x="category", y="fold", data=df[df["variant"] == "INDELs"], ax=axs[0])
    axs[0].set_xlabel("")
    axs[0].set_ylabel("fold vs. DHS")
    axs[0].set_title("INDELs")
    axs[0].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[0].grid(axis="y", alpha=0.3)
    axs[0].set_ylim([0, 1.2])

    sns.barplot(x="category", y="fold", data=df[df["variant"] == "SNPs"], ax=axs[1])
    axs[1].set_xlabel("")
    axs[1].set_ylabel("fold vs. DHS")
    axs[1].set_title("SNPs")
    axs[1].set_xticks(np.arange(3), x_labels, rotation=45, ha="right")
    axs[1].grid(axis="y", alpha=0.3)
    axs[1].set_ylim([0, 1.2])

    plt.tight_layout()
    save_file = join(PLOTS_DIR, "HepG2_SNPs_INDELs_fold_vs_DHS_barplot.pdf")
    print(save_file)
    plt.savefig(save_file)
    #
    # print("Average values:\n")
    # print(df.groupby(by=["variant", "category"]).mean())


def variant_stats_for_text(df_density, df_counts):

    _vars = ["INDELs", 'SNPs', 'eQTLs', 'caQTLs', "raQTLs", "GWAS", "ClinVar"]

    for _var in _vars:
        print(_var)

        _df = df_density[df_density["variant"] == _var]
        regions = ['HOT_enh', 'RE', 'HOT_prom', 'RP']
        counts = [_df[_df["region"] == _]["value"].values[0] for _ in regions]
        # counts = ["%.1f" % (100000 * _) for _ in counts]
        print("\t".join(regions))
        print("\t".join(["%.1f" % (100000 * _) for _ in counts]))
        folds = [counts[0] / counts[1], counts[2] / counts[3]]
        print("fold enhs \t fold proms")
        print("\t".join([str(_) for _ in folds]))

        _df = df_counts[df_counts["variant"] == _var]
        _df = df_counts[df_counts["variant"] == _var]
        v1 = _df[_df["region"] == "HOT_prom"].value
        v2 = _df[_df["region"] == "RP"].value
        p_prom = mannwhitneyu(v1, v2).pvalue
        v1 = _df[_df["region"] == "HOT_enh"].value
        v2 = _df[_df["region"] == "RE"].value
        p_enh = mannwhitneyu(v1, v2).pvalue
        print(p_enh, p_prom)

        print()
