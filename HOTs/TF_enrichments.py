import subprocess
import warnings
from collections import defaultdict

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

# matplotlib.style.use('ggplot')
import seaborn as sns

from overbinders.data_prep.basic import load_metadata

PROJECT_DIR = "ROOT_DIR/"

LOCI_DIR = "ROOT_DIR/definitions/peak_8bp_v2/HOTs"
PLOTS_DIR = "ROOT_DIR/plots/HOTs/"


def get_tfs_dhs_density():

    peaks_dir = "ROOT_DIR/chipseq_files/peaks"

    tf_id2density = {}
    a = load_metadata()

    for cl in ["HepG2", "K562"]:
        print(cl)
        dhs_file = "ROOT_DIR/chipseq_files/DHS_%s.bed.gz" % cl
        dhs = BedTool(dhs_file)

        dhs_cov = sum([r.length for r in dhs])

        all_tfs = len(a[cl])
        cnt = 1
        for tf, _id in a[cl].items():
            print(cnt, all_tfs, tf, _id)
            cnt += 1
            peaks_file = join(peaks_dir, "%s_hg19_peak_8bp.bed" % _id)

            peaks = BedTool(peaks_file)
            dhs_ovp = dhs.intersect(peaks)

            if dhs_ovp.count() == 0:
                peak_cov = 0
            else:
                peak_cov = sum([r.length for r in dhs_ovp])

            density = peak_cov/dhs_cov

            tf_id2density[_id] = density

    save_file = "ROOT_DIR/chipseq_files/tf_dhs_densities.txt"
    with open(save_file, "w") as of:
        for k,v in tf_id2density.items():
            out_line = "%s\t%s\n" % (k, str(v))
            of.write(out_line)


def tf_HOT_density():
    
    id2hot_density = defaultdict(float)
    id2proms_density = defaultdict(float)
    id2enhs_density = defaultdict(float)

    a = load_metadata()

    for cl in ["HepG2", "K562"]:

        cl_tf_ids = set(a[cl].values())

        hot_enh_file = join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)
        hot_prom_file = join(LOCI_DIR, "%s_HOTs.proms.bed" % cl)
        hot_file = join(LOCI_DIR, "%s_HOTs.bed" % cl)

        hot_enhs = BedTool(hot_enh_file)
        hot_proms = BedTool(hot_prom_file)
        hots = BedTool(hot_file)

        for _map, rs in zip([id2enhs_density, id2proms_density, id2hot_density],
                            [hot_enhs, hot_proms, hots]):
            rs_total = sum([r.length for r in rs])

            for r in rs:
                for tf_id in r.fields[8].split(","):
                    _map[tf_id] += 8

            for tf_id in cl_tf_ids:
                _map[tf_id] /= rs_total

    for _map, fname in zip([id2enhs_density, id2proms_density, id2hot_density],
                            ["tf_enhs_densities.txt", "tf_proms_densities.txt", "tf_hots_densities.txt"]):

        save_file = join(LOCI_DIR, fname)
        with open(save_file, "w") as of:
            for k,v in _map.items():
                out_line = "%s\t%s\n" % (k, str(v))
                of.write(out_line)


def tf_HOT_percentage():

    a = load_metadata()

    id2hot = defaultdict(float)
    id2proms = defaultdict(float)
    id2enhs = defaultdict(float)

    for cl in ["HepG2", "K562"]:

        cl_tf_ids = set(a[cl].values())

        hot_enh_file = join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)
        hot_prom_file = join(LOCI_DIR, "%s_HOTs.proms.bed" % cl)
        hot_file = join(LOCI_DIR, "%s_HOTs.bed" % cl)

        hot_enhs = BedTool(hot_enh_file)
        hot_proms = BedTool(hot_prom_file)
        hots = BedTool(hot_file)

        for _map, rs in zip([id2enhs, id2proms, id2hot],
                            [hot_enhs, hot_proms, hots]):
            rs_total = len(rs)
            for r in rs:
                for tf_id in set(r.fields[8].split(",")):
                    _map[tf_id] += 1

            for tf_id in cl_tf_ids:
                _map[tf_id] /= rs_total

    for _map, fname in zip([id2enhs, id2proms, id2hot],
                           ["tf_enhs_frac.txt", "tf_proms_frac.txt", "tf_hots_frac.txt"]):

        save_file = join(LOCI_DIR, fname)
        with open(save_file, "w") as of:
            for k, v in _map.items():
                out_line = "%s\t%s\n" % (k, str(v))
                of.write(out_line)


def tf_HOT_percentage_tf_wise():

    peaks_dir = "ROOT_DIR/chipseq_files/peaks/"

    a = load_metadata()

    # data_pool = []
    #
    # for cl in ["HepG2", "K562"]:
    #     print(cl)
    #     cl_tf_ids = set(a[cl].values())
    #
    #     hot_enh_file = join(LOCI_DIR, "%s_HOTs.noproms.bed" % cl)
    #     hot_prom_file = join(LOCI_DIR, "%s_HOTs.proms.bed" % cl)
    #     hot_file = join(LOCI_DIR, "%s_HOTs.bed" % cl)
    #
    #     hot_enhs = BedTool(hot_enh_file)
    #     hot_proms = BedTool(hot_prom_file)
    #     hots = BedTool(hot_file)
    #
    #     for cnt, tf_id in enumerate(cl_tf_ids):
    #         print(cnt, len(cl_tf_ids))
    #
    #         peak_file = join(peaks_dir, "%s_hg19_peak_8bp.bed" % tf_id)
    #         peaks = BedTool(peak_file)
    #         all_peaks = peaks.count()
    #         hot_peaks = peaks.intersect(hots, wa=True).count()
    #         enh_peaks = peaks.intersect(hot_enhs, wa=True).count()
    #         prom_peaks = peaks.intersect(hot_proms, wa=True).count()
    #
    #         data_pool.append([tf_id,
    #                           all_peaks,
    #                           100 * hot_peaks/all_peaks,
    #                           100 * enh_peaks/all_peaks,
    #                           100 * prom_peaks/all_peaks])

    # save_file = join(LOCI_DIR, "tf_perc_in_hots.txt")
    # with open(save_file, "w") as of:
    #     of.write("#tf_id\tall_peaks\tperc_all\tperc_enhs\tperc_proms\n")
    #     for _d in data_pool:
    #         out_line = "%s\t%d\t%2.2f\t%2.2f\t%2.2f\n" % (_d[0], _d[1], _d[2], _d[3], _d[4])
    #         of.write(out_line)

    tf_id2tf = {}
    for _map in a.values():
        for k, v in _map.items():
            tf_id2tf[v] = k

    f1 = join(LOCI_DIR, "tf_perc_in_hots.txt")
    f2 = join(LOCI_DIR, "tf_perc_in_hots_2.txt")

    lines = open(f1).readlines()[1:]
    data = []
    for line in lines:
        _d = line.split("\t")
        _d.insert(1, tf_id2tf[_d[0]])
        data.append(_d)
    with open(f2, "w") as of:
        of.write("#tf_id\ttf\tall_peaks\tperc_all\tperc_enhs\tperc_proms\n")
        for _d in data:
            out_line = "\t".join(_d)
            of.write(out_line)


def plot_enrichments():

    dhs_density_file = "ROOT_DIR/chipseq_files/tf_dhs_densities.txt"
    id2dhs = {}
    for l in open(dhs_density_file):
        parts = l.rstrip().split()
        id2dhs[parts[0]] = float(parts[1])

    a = load_metadata()
    tf_id2tf = {}
    for m in a.values():
        for k, v in m.items():
            tf_id2tf[v] = k

    fig, axs = plt.subplots(3, 2, figsize=(15, 10))

    for cl, axrow in zip(["HepG2", "K562"], [axs[:, 0], axs[:, 1]]):

        for fracs_file, density_file, category, ax in zip(["tf_enhs_frac.txt", "tf_proms_frac.txt", "tf_hots_frac.txt"],
                                                          ["tf_enhs_densities.txt", "tf_proms_densities.txt", "tf_hots_densities.txt"],
                                                          ["enhancers", "promoters", "all"],
                                                          axrow):

            fracs_file = join(LOCI_DIR, fracs_file)
            density_file = join(LOCI_DIR, density_file)
            id2fracs, id2density = defaultdict(float), defaultdict(float)
            for _map, fname in zip([id2fracs, id2density],
                                   [fracs_file, density_file]):
                for l in open(fname):
                    parts = l.rstrip().split()
                    _map[parts[0]] = float(parts[1])

            fracs = []
            enrichments = []
            tfs = []
            for tf, tf_id in a[cl].items():

                dhs_density = id2dhs[tf_id]
                frac = id2fracs[tf_id]

                hot_density = id2density[tf_id]

                if hot_density == 0 or dhs_density == 0:
                    print("skipping", tf, tf_id, dhs_density, hot_density)
                    continue

                fracs.append(frac)
                enrichments.append(np.log2(hot_density/dhs_density))
                tfs.append(tf)

            fracs = 100 * np.asarray(fracs)
            enrichments = np.asarray(enrichments)

            ax.scatter(fracs, enrichments, facecolors='none', edgecolors='grey')
            ax.set_title("%s %s HOTs" % (cl, category))
            ax.grid()
            ax.set_xlim([0, 100])
            ax.set_ylabel("log2(HOT/DHS)")

    axs[2, 1].set_xlabel("% of loci")
    axs[2, 0].set_xlabel("% of loci")

    plt.tight_layout()

    save_file = join(PLOTS_DIR, "tf_perc_enrichments.pdf")
    plt.savefig(save_file)
    plt.show()


def extract_enriched_tfs():

    save_dir = join(LOCI_DIR, "enrichment_files/")

    dhs_density_file = "ROOT_DIR/chipseq_files/tf_dhs_densities.txt"
    id2dhs = {}
    for l in open(dhs_density_file):
        parts = l.rstrip().split()
        id2dhs[parts[0]] = float(parts[1])

    a = load_metadata()
    tf_id2tf = {}
    for m in a.values():
        for k, v in m.items():
            tf_id2tf[v] = k

    for cl in ["HepG2", "K562"]:

        for fracs_file, density_file, category in zip(["tf_enhs_frac.txt", "tf_proms_frac.txt", "tf_hots_frac.txt"],
                                                      ["tf_enhs_densities.txt", "tf_proms_densities.txt", "tf_hots_densities.txt"],
                                                      ["enhancers", "promoters", "all"]):

            fracs_file = join(LOCI_DIR, fracs_file)
            density_file = join(LOCI_DIR, density_file)
            id2fracs, id2density = defaultdict(float), defaultdict(float)
            for _map, fname in zip([id2fracs, id2density],
                                   [fracs_file, density_file]):
                for l in open(fname):
                    parts = l.rstrip().split()
                    _map[parts[0]] = float(parts[1])

            fracs = []
            enrichments = []
            tfs = []
            for tf, tf_id in a[cl].items():

                dhs_density = id2dhs[tf_id]
                frac = id2fracs[tf_id]

                hot_density = id2density[tf_id]

                if hot_density == 0 or dhs_density == 0:
                    print("skipping", tf, tf_id, dhs_density, hot_density)
                    continue

                fracs.append(frac)
                enrichments.append(np.log2(hot_density / dhs_density))
                tfs.append(tf)

            fracs = 100 * np.asarray(fracs)
            enrichments = np.asarray(enrichments)

            save_file = join(save_dir, "%s_%s.txt" % (cl, category))
            print(save_file)
            with open(save_file, "w") as of:
                for tf, frac, enr in zip(tfs, fracs, enrichments):
                    out_line = "%s\t%f\t%f\n" % (tf, frac, enr)
                    of.write(out_line)


def merge_to_big_table():

    a = load_metadata()

    fname = join(LOCI_DIR, "tf_enhs_frac.txt")
    id2enh_fracs = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}
    fname = join(LOCI_DIR, "tf_proms_frac.txt")
    id2prom_fracs = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}
    fname = join(LOCI_DIR, "tf_hots_frac.txt")
    id2all_fracs = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}

    fname = "ROOT_DIR/chipseq_files/tf_dhs_densities.txt"
    id2dhs_dens = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}

    fname = join(LOCI_DIR, "tf_enhs_densities.txt")
    id2enh_enr = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}
    fname = join(LOCI_DIR, "tf_proms_densities.txt")
    id2prom_enr = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}
    fname = join(LOCI_DIR, "tf_hots_densities.txt")
    id2all_enr = {l.split()[0]: float(l.split()[1].rstrip()) for l in open(fname)}
    for _map in [id2enh_enr, id2prom_enr, id2all_enr]:
        for k,v in _map.items():
            try:
                _map[k] = np.log2(v/id2dhs_dens[k])
            except:
                _map[k] = np.NAN

    tf_fracs_peakwise = join(LOCI_DIR, "tf_perc_in_hots.txt")
    id2fracs_pk_enh = {}
    id2fracs_pk_prom = {}
    id2fracs_pk_all = {}
    id2num_peaks = {}
    for l in open(tf_fracs_peakwise).readlines()[1:]:
        parts = l.rstrip().split("\t")
        tf_id = parts[0]
        num_peaks = parts[2]
        frac_all = parts[3]
        frac_enh = parts[4]
        frac_prom = parts[5]

        id2fracs_pk_all[tf_id] = frac_all
        id2fracs_pk_enh[tf_id] = frac_enh
        id2fracs_pk_prom[tf_id] = frac_prom
        id2num_peaks[tf_id] = num_peaks

    save_file = join(LOCI_DIR, "tf_stats_summary_table.txt")

    data_pool = []
    for cl in ["HepG2", "K562"]:
        for tf, tf_id in a[cl].items():
            _d = [tf,
                  tf_id,
                  cl,
                  id2num_peaks[tf_id],
                  np.round(100 * id2all_fracs[tf_id], 2),
                  np.round(100 * id2enh_fracs[tf_id], 2),
                  np.round(100 * id2prom_fracs[tf_id], 2),
                  id2all_enr[tf_id],
                  id2enh_enr[tf_id],
                  id2prom_enr[tf_id],
                  id2fracs_pk_all[tf_id],
                  id2fracs_pk_enh[tf_id],
                  id2fracs_pk_prom[tf_id]]

            data_pool.append(_d)

    header = "tf\ttf_id\tcell_line\tnum_peaks\tall_fracs\tenh_fracs\tprom_fracs\tall_enr\tenh_enr\tprom_enr\tpk_fracs_all\tpk_fracs_enh\tpk_fracs_prom\n"
    with open(save_file, "w") as of:
        of.write(header)
        for _d in data_pool:
            out_line = "\t".join([str(_) for _ in _d]) + "\n"
            of.write(out_line)


def plot_enrichments_scatterplot():

    table_file = join(LOCI_DIR, "tf_stats_summary_table.txt")
    df = pandas.read_table(table_file)

    cl = "HepG2"
    # cl = "K562"

    df = df[df["cell_line"] == cl]

    x_arrays = ["all_fracs", "prom_fracs", "enh_fracs"]
    y_arrays = ["pk_fracs_all", "pk_fracs_prom", "pk_fracs_enh"]
    enr_arrays = ["all_enr", "prom_enr", "enh_enr"]
    fnames = ["all_hots.pdf", "prom_hots.pdf", "enh_hots.pdf"]
    titles = ["All HOTs", "Promoter HOTs", "Enhancer HOTs"]

    for _x, _y, enr, fname, title in zip(x_arrays, y_arrays, enr_arrays, fnames, titles):

        _df = df[df[enr] != -np.inf]

        cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
        g = sns.relplot(
            data=_df,
            x=_x, y=_y,
            hue=enr,
            palette=cmap
        )

        g.set(xlim=[-5, 100])
        g.set(ylim=[0, 100])
        g.ax.xaxis.grid(True, alpha=0.5)
        g.ax.yaxis.grid(True, alpha=0.5)
        g._legend.set_title("Enr")
        plt.xlabel("% of HOT loci")
        plt.ylabel("% of TFBSs")
        plt.suptitle(title + " " + cl)

        save_file = join(PLOTS_DIR, "enrichment_subplots/%s_%s" % (cl, fname))
        plt.savefig(save_file)
        plt.show()


def plot_enrichments_jointplot():

    table_file = join(LOCI_DIR, "tf_stats_summary_table.txt")
    df = pandas.read_table(table_file)

    cl = "HepG2"
    # cl = "K562"

    df = df[df["cell_line"] == cl]

    x_arrays = ["all_fracs", "prom_fracs", "enh_fracs"]
    y_arrays = ["pk_fracs_all", "pk_fracs_prom", "pk_fracs_enh"]
    enr_arrays = ["all_enr", "prom_enr", "enh_enr"]
    fnames = ["all_hots.pdf", "prom_hots.pdf", "enh_hots.pdf"]
    titles = ["All HOTs", "Promoter HOTs", "Enhancer HOTs"]

    for _x, _y, enr, fname, title in zip(x_arrays, y_arrays, enr_arrays, fnames, titles):
        _df = df[df[enr] != -np.inf]

        cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
        g = sns.jointplot(
            data=_df,
            x=_x, y=_y,
            palette=cmap
        )

        g.set_axis_labels("% of HOT loci", "% of TFBSs")
        g.ax_joint.grid(alpha=0.4)
        g.ax_joint.set_xlim([-5, 100])
        g.ax_joint.set_ylim([0, 100])
        plt.suptitle(title + " " + cl)

        save_file = join(PLOTS_DIR, "enrichment_subplots/jointplot_%s_%s" % (cl, fname))
        plt.savefig(save_file)
        plt.show()


def plot_enrichments_jointplot_size():

    table_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(table_file)

    cl = "HepG2"
    # cl = "K562"

    df = df[df["cell_line"] == cl]

    cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    g = sns.jointplot(
        data=df,
        x="all_fracs",
        y="pk_fracs_all",
        palette=cmap,
        height=6
    )

    g.ax_joint.cla()
    sc = sns.scatterplot(data=df, x="all_fracs", y="pk_fracs_all", hue="num_peaks", size="num_peaks", sizes=(10, 150), ax=g.ax_joint)

    h, l = sc.get_legend_handles_labels()
    l = [_[:2]+"K" for _ in l]
    sc.legend(h, l, title="Total peaks")

    g.set_axis_labels("% of overlapping HOT loci", "% of ChIP-seq peaks in HOT loci")
    g.ax_joint.grid(alpha=0.4)
    g.ax_joint.set_xlim([-5, 100])
    g.ax_joint.set_ylim([0, 100])

    save_file = join(PLOTS_DIR, "enrichment_subplots/jointplot_HepG2_TFs_hue_size.pdf")
    print(save_file)
    plt.savefig(save_file)



def plot_enrichments_jointplot_size_labels():

    table_file = join(LOCI_DIR, "enrichment_files/tf_stats_summary_table.txt")
    df = pandas.read_table(table_file)

    cl = "HepG2"
    # cl = "K562"

    df = df[df["cell_line"] == cl]

    cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    g = sns.jointplot(
        data=df,
        x="all_fracs",
        y="pk_fracs_all",
        palette=cmap,
        height=15
    )

    g.ax_joint.cla()
    sc = sns.scatterplot(data=df, x="all_fracs", y="pk_fracs_all", hue="num_peaks", size="num_peaks", sizes=(10, 150), ax=g.ax_joint)

    h, l = sc.get_legend_handles_labels()
    l = [_[:2]+"K" for _ in l]
    sc.legend(h, l, title="Total peaks")

    g.set_axis_labels("% of overlapping HOT loci", "% of ChIP-seq peaks in HOT loci")
    g.ax_joint.grid(alpha=0.4)
    g.ax_joint.set_xlim([-5, 100])
    g.ax_joint.set_ylim([0, 100])

    def label_point(x, y, val, ax):
        a = pandas.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x'] + .02, point['y'], str(point['val']))

    label_point(df.all_fracs, df.pk_fracs_all, df.tf, sc)

    save_file = join(PLOTS_DIR, "enrichment_subplots/jointplot_HepG2_TFs_hue_size_large_labels.pdf")
    print(save_file)
    plt.savefig(save_file)

