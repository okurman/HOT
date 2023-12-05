import gzip
import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import tempfile
from collections import defaultdict
from os.path import join
import numpy as np
from pybedtools import BedTool
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DATA_PATH = Path(os.environ["HOT_DATA"])
BINS_DIR = DATA_PATH/"log_bins"
PLOTS_DIR = DATA_PATH/"plots"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]

K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']



def plot_categorical_single_cell_line(df,
                                     y_label="Coverage",
                                     plot_type="boxplot",
                                     save_file=None,
                                     leg_loc="best",
                                     suptitle=None,
                                     single_category=False,
                                     cell_line="K562",
                                     figsize=(7,4)):

    f = plt.figure(figsize=figsize)
    ax = plt.gca()

    x_tick_labels = HEPG2_XTICK_LABELS if cell_line == "HepG2" else K562_XTICK_LABELS
    x_range = np.arange(len(x_tick_labels))
    x_ticks = x_range - 0.5

    if plot_type == "boxplot":
        if single_category:
            sns.boxplot(x="bin", y="value", data=df[df["cell_line"] == cell_line], showfliers=False, ax=ax, color=(0.20392156862745098, 0.5411764705882353, 0.7411764705882353))
        else:
            sns.boxplot(x="bin", y="value", data=df[df["cell_line"] == cell_line], showfliers=False, ax=ax, hue="category", hue_order=["enhancers", "promoters"])
    elif plot_type == "barplot":
        df_2 = df.groupby(["bin", "category"])["value"].mean().reset_index()
        sns.barplot(x="bin", y="value", hue="category", data=df_2, ax=ax, hue_order=["enhancers", "promoters"])

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels, size=8)
    ax.set_xlabel("DAPs")
    ax.legend(loc=leg_loc)
    ax.set_ylabel(y_label)
    # ax.set_title(cell_line)

    ax.axvspan(9.5, 14, alpha=0.1, color="black")
    ax.set_xlim([-1, 13.5])

    if suptitle:
        f.suptitle(suptitle)

    plt.tight_layout()
    plt.savefig(save_file)
    plt.close()


def plot_categorical_worker(df,
                            y_label="Coverage",
                            plot_type="boxplot",
                            save_file=None,
                            leg_loc="best",
                            suptitle=None,
                            single_category=False,
                            y_upper_limit=None,
                            y_lower_limit=None,
                            same_scale=True,
                            log_scale=False,
                            palette="colorblind",
                            figsize=(6, 4),
                            hline=None,
                            shadow=True):

    if single_category:
        df = df[df["category"] == "enhancers"]

    if plot_type == "lineplot":

        plt.figure(figsize=figsize)
        ax = plt.gca()

        sns.lineplot(x="bin", y="value", hue="cell_line", style="category", markers=True, data=df, ax=ax,
                     estimator="mean", ci=80,
                     palette=palette, style_order=["enhancers", "promoters"], hue_order=["HepG2", "K562"])

        x_range = np.arange(len(HEPG2_XTICK_LABELS))
        x_ticks = x_range - 0.5

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(PERC_XTICK_LABELS)
        ax.set_xlabel("DAPs")

        ticks = ax.get_xticklabels()
        for t in ticks[5:]:
            t.set_rotation(30)
            t.set_fontsize(9)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:3] + handles[4:], labels=labels[1:3] + labels[4:])

        ax.set_ylabel(y_label)

        if hline:
            plt.axhline(y=hline, linestyle="--")

        if y_upper_limit or y_lower_limit:
            ylims = ax.get_ylim()
            _lower = y_lower_limit if y_lower_limit is not None else ylims[0]
            _upper = y_upper_limit if y_upper_limit else ylims[1]
            ax.set_ylim([_lower, _upper])

        if shadow:
            ax.axvspan(9.5, 14, alpha=0.1, color="black")
            ax.set_xlim([-1, 13.5])

    else:

        f = plt.figure(figsize=figsize)
        ax1 = f.add_subplot(121)
        ax2 = f.add_subplot(122)
        for ax, cell_line in zip([ax1, ax2], ["HepG2", "K562"]):

            x_tick_labels = HEPG2_XTICK_LABELS if cell_line == "HepG2" else K562_XTICK_LABELS
            x_range = np.arange(len(x_tick_labels))
            x_ticks = x_range - 0.5

            if plot_type == "boxplot":
                if single_category:
                    sns.boxplot(x="bin", y="value", data=df[df["cell_line"] == cell_line], showfliers=False, ax=ax,
                                color=(0.20392156862745098, 0.5411764705882353, 0.7411764705882353))
                else:
                    sns.boxplot(x="bin", y="value", data=df[df["cell_line"] == cell_line], showfliers=False, ax=ax,
                                hue="category", hue_order=["enhancers", "promoters"], palette=palette)

            elif plot_type == "barplot":

                sns.barplot(x="bin", y="value", hue="category", data=df[df["cell_line"] == cell_line], ax=ax,
                            hue_order=["enhancers"] if single_category else ["enhancers", "promoters"], palette=palette)


            elif plot_type == "violin":

                sns.violinplot(x="bin", y="value", hue="category", data=df[df["cell_line"] == cell_line], ax=ax,
                               hue_order=["enhancers"] if single_category else ["enhancers", "promoters"], inner=None,
                               split=True)

            elif plot_type == "boxenplot":
                df_2 = df
                sns.boxenplot(x="bin", y="value", hue="category", data=df_2[df_2["cell_line"] == cell_line], ax=ax,
                              hue_order=["enhancers"] if single_category else ["enhancers", "promoters"])

            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_tick_labels, size=8)
            ax.set_xlabel("DAPs")
            ax.legend(loc=leg_loc)
            if single_category:
                ax.legend([])
            ax.set_ylabel(y_label)
            ax.set_title(cell_line)
            if y_upper_limit:
                ax.set_ylim([0, y_upper_limit])

        if same_scale:
            ylims_1 = ax1.get_ylim()
            ylims_2 = ax2.get_ylim()
            ylim = [min(ylims_1[0], ylims_2[0]), max(ylims_1[1], ylims_2[1])]
            if y_upper_limit:
                ylim = [ylim[0], y_upper_limit]
            ax1.set_ylim(ylim)
            ax2.set_ylim(ylim)

        if log_scale:
            ax1.set_yscale("log")
            ax2.set_yscale("log")

        if shadow:
            for ax in [ax1, ax2]:
                ax.axvspan(9.5, 14, alpha=0.1, color="black")
                ax.set_xlim([-1, 13.5])

    if suptitle:
        plt.suptitle(suptitle)

    plt.tight_layout()
    plt.savefig(save_file)


def load_bins_to_df(line_parser, file_amender=lambda x: x, zero_comp=False, coefficient=1, src_dir=None):

    tmp_file = tempfile.NamedTemporaryFile()
    columns = ["cell_line", "bin", "category", "value"]
    with open(tmp_file.name, "w") as outf:
        outf.write(",".join(columns) + "\n")
        for cell_line in ["HepG2", "K562"]:
            prom_files_loci = get_prom_files(cell_line)
            prom_files = list(map(file_amender, prom_files_loci))
            enh_files_loci = get_enh_files(cell_line)
            enh_files = list(map(file_amender, enh_files_loci))

            for bin in range(len(prom_files)):
                prom_file = prom_files[bin]
                enh_file = enh_files[bin]
                for cat, f in zip(["promoters", "enhancers"], [prom_file, enh_file]):

                    if src_dir:
                        f = Path(src_dir)/Path(f).name
                    
                    try:
                        lines = gzip.open(f, "rt").readlines()
                    except:
                        lines = open(f, "rt").readlines()

                    for line in lines:
                        value = line_parser(line)
                        if value is None:
                            continue

                        if isinstance(value, list):
                            for _ in value:
                                out_line = "%s,%d,%s,%s\n" % (cell_line, bin, cat, _)
                                outf.write(out_line)
                        else:
                            out_line = "%s,%d,%s,%s\n" % (cell_line, bin, cat, value)
                            outf.write(out_line)

            if zero_comp:
                for bin in range(len(prom_files)):
                    prom_file = prom_files[bin]
                    prom_file_loci = prom_files_loci[bin]
                    enh_file = enh_files[bin]
                    enh_file_loci = enh_files_loci[bin]

                    for cat, f, f_loci in zip(["promoter", "enhancer"],
                                              [prom_file, enh_file],
                                              [prom_file_loci, enh_file_loci]):

                        comp_loci = BedTool(f_loci).intersect(f, v=True).count()
                        out_line = "%s,%d,%s,0\n" % (cell_line, bin, cat)
                        [outf.write(out_line) for _ in range(comp_loci)]

    df = pd.read_csv(tmp_file)

    if coefficient != 1:
        df['value'] = df['value'].apply(lambda x: x * coefficient)

    return df


def load_peaks_motifs_to_df(file_amender=lambda x: x):

    tmp_file = tempfile.NamedTemporaryFile()
    columns = ["cell_line", "bin", "category", "bio", "value"]
    with open(tmp_file.name, "w") as outf:
        outf.write(",".join(columns) + "\n")
        for cell_line in ["HepG2", "K562"]:
            prom_files_loci = get_prom_files(cell_line)
            prom_files = list(map(file_amender, prom_files_loci))
            enh_files_loci = get_enh_files(cell_line)
            enh_files = list(map(file_amender, enh_files_loci))

            for bin in range(len(prom_files)):
                prom_file = prom_files[bin]
                enh_file = enh_files[bin]
                for bio, f in zip(["promoters", "enhancers"], [prom_file, enh_file]):
                    for line in open(f):
                        parts = line.split()
                        peaks = int(parts[3])

                        peak_ids = set(parts[5].split(","))
                        motif_ids = set(parts[10].split(",")).intersection(peak_ids)
                        motifs = len(motif_ids)

                        out_line = "%s,%d,%s,%s,%s\n" % (cell_line, bin, "peaks", bio, peaks)
                        outf.write(out_line)
                        out_line = "%s,%d,%s,%s,%s\n" % (cell_line, bin, "motifs", bio, motifs)
                        outf.write(out_line)

            for bin in range(len(prom_files)):
                prom_file = prom_files[bin]
                prom_file_loci = prom_files_loci[bin]
                enh_file = enh_files[bin]
                enh_file_loci = enh_files_loci[bin]

                for bio, f, f_loci in zip(["promoter", "enhancer"],
                                          [prom_file, enh_file],
                                          [prom_file_loci, enh_file_loci]):

                    comp_loci = BedTool(f_loci).intersect(f, v=True)
                    for r in comp_loci:
                        out_line = "%s,%d,%s,%s,%s\n" % (cell_line, bin, "peaks", bio, int(r.fields[5]))
                        outf.write(out_line)
                        out_line = "%s,%d,%s,%s,%s\n" % (cell_line, bin, "motifs", bio, 0)
                        outf.write(out_line)

    df = pd.read_csv(tmp_file)

    return df


def load_peaks_motifs_fractions_to_df(file_amender=lambda x: x):

    tmp_file = tempfile.NamedTemporaryFile()
    columns = ["cell_line", "bin", "category", "value"]
    with open(tmp_file.name, "w") as outf:
        outf.write(",".join(columns) + "\n")
        for cell_line in ["HepG2", "K562"]:
            prom_files_loci = get_prom_files(cell_line)
            prom_files = list(map(file_amender, prom_files_loci))
            enh_files_loci = get_enh_files(cell_line)
            enh_files = list(map(file_amender, enh_files_loci))

            for bin in range(len(prom_files)):
                prom_file = prom_files[bin]
                enh_file = enh_files[bin]
                for bio, f in zip(["promoters", "enhancers"], [prom_file, enh_file]):
                    for line in open(f):
                        parts = line.split()
                        peaks = int(parts[3])

                        peak_ids = set(parts[5].split(","))
                        motif_ids = set(parts[10].split(",")).intersection(peak_ids)
                        motifs = len(motif_ids)

                        out_line = "%s,%d,%s,%f\n" % (cell_line, bin, bio, motifs/peaks)
                        outf.write(out_line)

            for bin in range(len(prom_files)):
                prom_file = prom_files[bin]
                prom_file_loci = prom_files_loci[bin]
                enh_file = enh_files[bin]
                enh_file_loci = enh_files_loci[bin]

                for bio, f, f_loci in zip(["promoters", "enhancers"],
                                          [prom_file, enh_file],
                                          [prom_file_loci, enh_file_loci]):

                    comp_loci = BedTool(f_loci).intersect(f, v=True)
                    for r in comp_loci:
                        out_line = "%s,%d,%s,%f\n" % (cell_line, bin, bio, 0)
                        outf.write(out_line)

    df = pd.read_csv(tmp_file)

    return df


def load_df_from_overlaps(overlap_file, normalize=True, coverage=False, separate_files=None, coefficient=1):

    tmp_file = tempfile.NamedTemporaryFile()
    columns = ["cell_line", "bin", "category", "value"]
    with open(tmp_file.name, "w") as outf:
        outf.write(",".join(columns) + "\n")

        for cell_line in ["HepG2", "K562"]:
            print(cell_line)
            prom_files = get_prom_files(cell_line)
            enh_files = get_enh_files(cell_line)

            if separate_files:
                overlap_file = separate_files[0] if cell_line == "HepG2" else separate_files[1]

            for cat, files in zip(["promoters", "enhancers"], [prom_files, enh_files]):
                for bin, file in enumerate(files):
                    bin_bed = BedTool(file)
                    bin_count = bin_bed.count()
                    overlaps = bin_bed.intersect(overlap_file, wo=True)
                    _cnt = overlaps.count()
                    if _cnt == 0:
                        out_line = "%s,%d,%s,%f\n" % (cell_line, bin, cat, 0)
                    else:
                        if coverage:
                            metric = np.mean([float(r.fields[-1]) / 400 for r in overlaps])
                        else:
                            num_fields = len(overlaps[0].fields)
                            metric = overlaps.sort().groupby(c=[num_fields - 1], o="collapse").count()

                        if normalize:
                            metric /= bin_count
                        out_line = "%s,%d,%s,%f\n" % (cell_line, bin, cat, metric*coefficient)
                    outf.write(out_line)

    df = pd.read_csv(tmp_file.name)

    return df


def load_histone_overlap_maps():

    file_amender = lambda x: x + ".all_histones.bed"
    map_pool = []
    for cell_line in ["HepG2", "K562"]:
        cl_pool = [[], []]
        prom_files = get_prom_files(cell_line)
        enh_files = get_enh_files(cell_line)
        for i, files in enumerate([prom_files, enh_files]):
            counts = [BedTool(f).count() for f in files]
            files = list(map(file_amender, files))
            for count, f in zip(counts, files):
                h2c = defaultdict(int)
                for line in open(f):
                    hms = set(line.rstrip().split()[-1].split(","))
                    for hm in hms:
                        h2c[hm] += 1
                h2c = defaultdict(float, {k:(v/count) for (k,v) in h2c.items()})
                cl_pool[i].append(h2c)
        map_pool.append(cl_pool)

    return map_pool

