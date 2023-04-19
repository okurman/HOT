
import warnings
from collections import defaultdict, Counter
import h5py
warnings.filterwarnings('ignore')

import sys
import os
from os.path import join

import numpy as np
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
from lib import phastCons, phyloP
import seaborn as sns
from upsetplot import from_memberships
from upsetplot import plot as upset_plot

from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"

SE_DIR = "/panfs/pan1/devdcode/sanjar/superenhancers/"


def extract_conservation_for_loci(cons_type="phastcons"):

    for ch_no in range(22, 0, -1):
        chrom = "chr%d" % ch_no
        print(chrom)
        if cons_type == "phastcons":
            pos2score_def = phastCons.load_scores(chrom=chrom)[chrom]
            fname_fmt = "%s.phastcons.bed"
        else:
            pos2score_def = phyloP.load_phylop_scores(chrom=chrom)[chrom]
            fname_fmt = "%s.phylop.bed"

        pos2score = {k: v for k, v in pos2score_def.items()}

        for cell_line in ["HepG2", "K562"]:
            print(cell_line)

            locus_file = join(SE_DIR, "%s.bed" % cell_line)
            scores_file = join(SE_DIR, fname_fmt % cell_line)

            chrom_bed = [r for r in BedTool(locus_file) if r.chrom == chrom]
            with open(scores_file, "a+") as outf:
                for r in chrom_bed:
                    scores = [pos2score[p] if p in pos2score else "None" for p in range(r.start, r.stop)]

                    out_line = "%s\t%d\t%d\t%s\n" % (r.chrom, r.start, r.stop, ",".join([str(_) for _ in scores]))
                    outf.write(out_line)
                    outf.flush()


def HOTs_in_SE_conservation():

    HepG2_SE_file = join(SE_DIR, "HepG2.phastcons.bed")
    K562_SE_file = join(SE_DIR, "K562.phastcons.bed")
    save_file = join(PROJECT_DIR, "chipseq_files/SE_HOT_scores.phastcons.h5")
    print(save_file)

    # HepG2_SE_file = join(SE_DIR, "HepG2.phylop.bed")
    # K562_SE_file = join(SE_DIR, "K562.phylop.bed")
    # save_file = join(PROJECT_DIR, "chipseq_files/SE_HOT_scores.phylop.h5")

    HepG2_SE = BedTool(HepG2_SE_file)
    K562_SE = BedTool(K562_SE_file)

    data = []

    for cl, SE, other_SE in zip(["HepG2", "K562"], [HepG2_SE, K562_SE], [K562_SE, HepG2_SE]):

        print(cl)
        hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl))

        SE_HOT_ovp = SE.intersect(hots, wo=True).groupby(c=[6, 7], o="collapse")
        enh_id2scores = {"%s-%d-%d" % (r.chrom, r.start, r.stop):
                             [float(_) if not _ == "None" else np.nan for _ in r.fields[3].rstrip().split(",")]
                              for r in SE}
        se_scores = np.asarray([_ for scores in enh_id2scores.values() for _ in scores])

        se_hot_scores = []
        for r in SE_HOT_ovp:

            enh_id = "%s-%d-%d" % (r.chrom, r.start, r.stop)
            r_scores = enh_id2scores[enh_id]

            hots_starts = r.fields[3].split(",")
            hots_stops = r.fields[4].split(",")

            for _start, _stop in zip(hots_starts, hots_stops):
                _start = int(_start) - r.start
                _stop = int(_stop) - r.start

                se_hot_scores += [r_scores[p] for p in range(_start, min(_stop, r.length))]

        se_hot_scores = np.asarray(se_hot_scores)
        data.append([se_scores, se_hot_scores])

    with h5py.File(save_file, "w") as of:
        of["HepG2/SE"] = data[0][0]
        of["HepG2/SE_HOT"] = data[0][1]
        of["K562/SE"] = data[1][0]
        of["K562/SE_HOT"] = data[1][1]


def plot_HOTs_in_SE_conservation():

    # scores_file = join(PROJECT_DIR, "chipseq_files/SE_HOT_scores.phylop.h5")
    # save_file = join(PLOTS_DIR, "SE_HOT_phylop_scores.pdf")

    scores_file = join(PROJECT_DIR, "chipseq_files/SE_HOT_scores.phastcons.h5")
    save_file = join(PLOTS_DIR, "SE_HOT_phastcons_scores.pdf")

    data = []

    with h5py.File(scores_file, "r") as of:
        data.append([of["HepG2/SE"].value, of["HepG2/SE_HOT"].value])
        data.append([of["K562/SE"].value, of["K562/SE_HOT"].value])

    list = [["HepG2", "SE", np.nanmean(data[0][0])],
            ["HepG2", "SE_HOT", np.nanmean(data[0][1])],
            ["K562", "SE", np.nanmean(data[1][0])],
            ["K562", "SE_HOT", np.nanmean(data[1][1])]]

    columns = ["cell_line", "category", "value"]

    df = pandas.DataFrame(list, columns=columns)

    plt.close()
    sns.barplot(x="cell_line", y="value", hue="category", data=df)
    plt.legend()
    plt.xlabel("")
    # plt.ylabel("PhyloP scores")
    plt.ylabel("phastCons scores")
    plt.title("Conservation scores")
    plt.tight_layout()
    plt.savefig(save_file)
    plt.show()


def HOTs_in_SE_fractions():

    # hot_enh_file = join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)
    # hot_prom_file = join(LOCI_DIR, "%s_HOTs.proms.bed" % cl)

    HepG2_DHS = BedTool("ROOT_DIR/chipseq_files/DHS_HepG2.old_version.bed.gz")
    K562_DHS = BedTool("ROOT_DIR/chipseq_files/DHS_K562.old_version.bed.gz")

    HepG2_enhancers = BedTool(join(PROJECT_DIR, "enhancer_definitions/E118_enhancers.bed"))
    K562_enhancers = BedTool(join(PROJECT_DIR, "enhancer_definitions/E123_enhancers.bed"))

    HepG2_SE_file = join(SE_DIR, "HepG2.bed")
    HepG2_SE = BedTool(HepG2_SE_file)
    K562_SE_file = join(SE_DIR, "K562.bed")
    K562_SE = BedTool(K562_SE_file)

    data = []

    for cl, SE, other_SE in zip(["HepG2", "K562"], [HepG2_SE, K562_SE], [K562_SE, HepG2_SE]):

        hots = BedTool(join(LOCI_DIR, "%s_HOTs.bed" % cl))

        total_SE_cov = sum([r.length for r in SE])
        print(cl, total_SE_cov)
        total_SE_HOTs_cov = sum([r.length for r in SE.intersect(hots)])
        total_SE_density = total_SE_HOTs_cov / total_SE_cov

        UBQ_SE = HepG2_SE.intersect(other_SE, wa=True)
        UBQ_cov = sum([r.length for r in UBQ_SE])
        UBQ_HOTs_cov = sum([r.length for r in UBQ_SE.intersect(hots)])
        UBQ_density = UBQ_HOTs_cov/UBQ_cov

        OVP_SE = SE.intersect(other_SE)
        OVP_cov = sum([r.length for r in OVP_SE])
        OVP_HOTs_cov = sum([r.length for r in OVP_SE.intersect(hots)])
        OVP_density = OVP_HOTs_cov / OVP_cov

        data.append([total_SE_density, UBQ_density, OVP_density])

    return data


def plot_HOTs_in_SE_fractions(data=None):

    if not data:
        data = HOTs_in_SE_fractions()

    list = [["HepG2", "SE", data[0][0]],
            ["HepG2", "SE OVP", data[0][2]],
            ["K562", "SE", data[1][0]],
            ["K562", "SE OVP", data[1][2]]]

    columns = ["cell_line", "category", "value"]

    df = pandas.DataFrame(list, columns=columns)

    plt.close()
    sns.barplot(x="cell_line", y="value", hue="category", data=df)
    plt.legend()
    plt.xlabel("")
    plt.ylabel("Density of HOT regions")
    plt.title("HOTs in Superenhancers")
    plt.tight_layout()
    plt.show()

    save_file = join(PLOTS_DIR, "HOTs_in_SE_densities.pdf")
    plt.savefig(save_file)


def vista_overlaps():

    vistas = BedTool("/panfs/pan1/devdcode/common/vista/database.bed")
    all_cls = np.asarray([len(r.fields[4].split(";")) for r in vistas])

    HepG2_enhancers = BedTool(join(PROJECT_DIR, "enhancer_definitions/E118_enhancers.bed"))
    K562_enhancers = BedTool(join(PROJECT_DIR, "enhancer_definitions/E123_enhancers.bed"))

    columns = ["cell_line", "category", "value"]
    list = []

    for cl, enhancers in zip(["HepG2", "K562"], [HepG2_enhancers, K562_enhancers]):

        hots = join(LOCI_DIR, "%s_HOTs.bed" % cl)

        enh_cls = np.asarray([len(r.fields[4].split(";")) for r in vistas.intersect(enhancers, wa=True)])
        [list.append([cl, "regular enhancers", _]) for _ in enh_cls]
        hot_cls = np.asarray([len(r.fields[4].split(";")) for r in vistas.intersect(hots, wa=True)])
        [list.append([cl, "HOT enhancers", _]) for _ in hot_cls]

    df = pandas.DataFrame(list, columns=columns)

    plt.figure()
    ax = plt.gca()

    g = sns.violinplot(x="cell_line", y="value", hue="category", data=df, ax=ax, split="True")
    g.set_ylim([-0.2, 10])

    plt.legend(loc="upper center")
    plt.xlabel("")
    plt.ylabel("# of active tissues")
    plt.title("Vista enhancers")

    plt.show()

    save_file = join(PLOTS_DIR, "VISTA_tissues_violinplot.pdf")
    plt.savefig(save_file)


def vista_overlaps_tissues(data=None):

    vistas = BedTool("/panfs/pan1/devdcode/common/vista/database.bed")

    if not data:
        data = []
        for cl in ["HepG2", "K562"]:
            hots = join(LOCI_DIR, "%s_HOTs.bed" % cl)
            tissues = [_.split("[")[0] for r in vistas.intersect(hots, wa=True) for _ in r.fields[-1].split(";")]
            tissue2count = Counter(tissues)
            data.append(tissue2count)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    for cl, ax, _map in zip(["HepG2", "K562"], axs, data):

        sorted_map = sorted(_map.items(), key=lambda x: x[1], reverse=True)
        tissues = [_[0].split("(")[0] for _ in sorted_map]
        counts = [_[1] for _ in sorted_map]

        x_range = np.arange(len(tissues))

        ax.set_title(cl)
        ax.bar(x_range, counts)
        ax.set_xticks(x_range)
        ax.set_xticklabels(tissues, rotation=45, ha="right")
        ax.set_ylabel("# of HOT enhancers")

    plt.suptitle("VISTA enhancer active tissues")
    plt.tight_layout()

    save_file = join(PLOTS_DIR, "VISTA_tissues_counts.pdf")
    plt.savefig(save_file)


def draw_UpSet_SE_HOT_RE_HM():

    re_file = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac.bed"
    se_file = "/panfs/pan1/devdcode/sanjar/superenhancers/HepG2.bed"
    hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed"))
    hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed"))
    k27ac = BedTool("/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF392KDI_hg19.bed.gz")
    k4me1 = BedTool("/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF413EGR_hg19.bed.gz")

    save_file = join(PLOTS_DIR, "HepG2_HOTs_SE_k27ac_UpSet_HM.pdf")

    re = BedTool(re_file)
    se = BedTool(se_file).sort()

    data = from_memberships(
        [['HOT prom'],
         ['HOT enh '],
         ['Super-enh'],
         ['Reg enh '],
         ['HOT prom', 'H3K27ac '],
         ['HOT prom', 'H3K4me1'],
         ['HOT enh ', 'H3K27ac '],
         ['HOT enh ', 'H3K4me1'],
         ['HOT prom', 'Super-enh'],
         ['HOT enh ', 'Super-enh'],
         ['HOT enh ', 'Super-enh', 'Reg enh '],
         ['Reg enh ', 'Super-enh'],
         ],
        data=[
            hot_proms.count(),
            hot_enhs.count(),
            se.count(),
            re.count(),
            hot_proms.intersect(k27ac, wa=True).merge(d=-1).count(),
            hot_proms.intersect(k4me1, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(k27ac, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(k4me1, wa=True).merge(d=-1).count(),
            hot_proms.intersect(se, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(se, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(se, wa=True).merge(d=-1).intersect(re).sort().merge(d=-1).count(),
            re.intersect(se, wa=True).merge(d=-1).count(),
        ]
    )

    fig = plt.figure()

    upset_plot(data, fig=fig)

    print(save_file)
    fig.savefig(save_file, bbox_inches='tight')
    return data


def draw_UpSet_SE_HOT_RE():

    re_file = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac.bed"
    # re_file = "ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac_OR_H3K4me1.bed"
    se_file = "/panfs/pan1/devdcode/sanjar/superenhancers/HepG2.bed"

    hot_enhs = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed"))
    hot_proms = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed"))
    save_file = join(PLOTS_DIR, "HepG2_HOTs_SE_k27ac_UpSet_v2.pdf")
    # save_file = join(PLOTS_DIR, "HepG2_HOTs_SE_k27ac_or_k4me1_UpSet.pdf")

    # re_file = "ROOT_DIR/definitions/regular_enhancers/K562_enhancers.bed.gz"
    # se_file = "/panfs/pan1/devdcode/sanjar/superenhancers/K562.bed"
    # hot_file = join(LOCI_DIR, "K562_HOTs.PE.bed")
    # save_file = join(PLOTS_DIR, "K562_HOTs_SE_UpSet.pdf")

    re = BedTool(re_file).intersect().sort()
    se = BedTool(se_file).sort()

    data = from_memberships(
        [['HOT promoters      '],
         ['HOT enhancers      '],
         ['Super-enhancers   '],
         ['Regular enhancers'],
         ['HOT promoters      ', 'Super-enhancers   '],
         ['HOT enhancers      ', 'Super-enhancers   '],
         ['HOT enhancers      ', 'Super-enhancers   ', 'Regular enhancers'],
         ['Regular enhancers', 'Super-enhancers   '],
         ],
        data=[
            hot_proms.count(),
            hot_enhs.count(),
            se.count(),
            re.count(),
            hot_proms.intersect(se, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(se, wa=True).merge(d=-1).count(),
            hot_enhs.intersect(se, wa=True).merge(d=-1).intersect(re).sort().merge(d=-1).count(),
            re.intersect(se, wa=True).merge(d=-1).count(),
        ]
    )

    fig = plt.figure()

    upset_plot(data, fig=fig)

    print(save_file)
    fig.savefig(save_file, bbox_inches='tight')
    return data


def get_percentages():

    re = BedTool("ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac.bed")
    # re_2 = BedTool("ROOT_DIR/definitions/regular_enhancers/HepG2_enhancers_DHS_H3K27ac_OR_H3K4me1.bed")
    se = BedTool("/panfs/pan1/devdcode/sanjar/superenhancers/HepG2.bed").sort().merge()
    hots = BedTool(join(LOCI_DIR, "HepG2_HOTs.bed"))
    hots_count = hots.count()

    hot_enh = BedTool(join(LOCI_DIR, "HepG2_HOTs.enhs.bed"))
    hot_enh_count = hot_enh.count()
    hot_prom = BedTool(join(LOCI_DIR, "HepG2_HOTs.proms.bed"))
    hot_prom_count = hot_prom.count()

    hot_enh_se_count = hot_enh.intersect(se, wa=True).sort().merge(d=-1).count()
    hot_prom_se_count = hot_prom.intersect(se, wa=True).sort().merge(d=-1).count()

    print("HOT enh SE fraction: %f" % (hot_enh_se_count/hot_enh_count))
    print("HOT prom SE fraction: %f" % (hot_prom_se_count/hot_prom_count))

    hot_se = hots.intersect(se).merge()
    hot_se_cov = sum([r.length for r in hot_se])/sum([r.length for r in se])
    print("HOT SE cov fraction: %f" % hot_se_cov)

    k27ac = BedTool("/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF392KDI_hg19.bed.gz")
    k27ac_cov = hots.intersect(k27ac, wa=True).merge(d=-1).count()/hots_count
    print("HOTs K27ac fraction: %f" % k27ac_cov)

    k4me1 = BedTool("/panfs/pan1/devdcode/common/ENCODE_phase4/files_hg19/ENCFF413EGR_hg19.bed.gz")
    k4me1_cov = hots.intersect(k4me1, wa=True).merge(d=-1).count()/hots_count
    print("HOTs K4me1 fraction: %f" % k4me1_cov)

    both_hm = k27ac.cat(k4me1)
    both_cov = hots.intersect(both_hm, wa=True).merge(d=-1).count() / hots_count
    print("HOTs K4me1 AND K27ac fraction: %f" % both_cov)
















if __name__ == "__main__":

    pass






