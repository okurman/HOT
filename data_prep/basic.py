#!//usr/local/bin/python

import subprocess
import sys
import os
import tempfile
from os.path import join

import numpy as np

from pybedtools import BedTool
from collections import defaultdict, Counter

PROJECT_DIR = "/net/intdev/devdcode/sanjar/overbinders"

exons = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/exon_regions.bed"
promoters = "/panfs/pan1/patternquest/data/genomes/hg19/annotations/from_wei/output_ENSG_promoter.UCSC.alternativetss"

# hg19 = "/panfs/pan1/patternquest/data/genomes/hg19/hg19.genome"
hg19_fa = "/data/Dcode/common/hg19.fa"
hg19 = "/data/Dcode/common/hg19.genome"



def load_metadata_backup():

    metadata_file = os.path.join(PROJECT_DIR, "chipseq_files/metadata_filtered.tsv")
    cell_line2tfs = defaultdict(list)
    exp_id2lines = defaultdict(list)

    for line in open(metadata_file):
        parts = line.strip().split("\t")
        exp_id = parts[6]
        exp_id2lines[exp_id].append(line)

    filtered_lines = []
    for exp_id, lines in exp_id2lines.items():

        if len(lines) == 1:
            filtered_lines += lines
        elif len(lines) == 2:

            if lines[0].split("\t")[4] == "optimal IDR thresholded peaks":
                filtered_lines.append(lines[0])
            elif lines[1].split("\t")[4] == "optimal IDR thresholded peaks":
                filtered_lines.append(lines[1])
            else:
                raise Exception("Check output type!")

        else:
            raise Exception("shouldn't be here!")

    cell_line2tf2lines = defaultdict(dict)
    for line in filtered_lines:

        parts = line.strip().split("\t")
        tf = parts[21].replace("-human", "")
        cell_line = parts[9]
        cell_line2tfs[cell_line].append(tf)

        if tf not in cell_line2tf2lines[cell_line]:
            cell_line2tf2lines[cell_line][tf] = []

        cell_line2tf2lines[cell_line][tf].append(line)

    # for cell_line, tfs in cell_line2tfs.items():
    #     if len(tfs) != len(set(tfs)):
    #         print cell_line, Counter(tfs)

    cell_line2tf2file_id = defaultdict(dict)

    with open(os.path.join(PROJECT_DIR, "chipseq_files/cell_line_tf_map_v1.txt"), "w") as of:
        of.write("#Cell_line\tTF\tfile_id\n")
        for cell_line in cell_line2tf2lines:
            for tf, lines in cell_line2tf2lines[cell_line].items():

                if len(lines) == 1:
                    most_recent = lines[0]
                else:
                    most_recent = sorted(lines, key=lambda x: x.split("\t")[28], reverse=True)[0]

                file_id = most_recent.split("\t")[0]
                of.write("%s\t%s\t%s\n" % (cell_line, tf, file_id))
                cell_line2tf2file_id[cell_line][tf] = file_id

    return cell_line2tf2file_id


def load_metadata(metadata_file=None):

    if not metadata_file:
        metadata_file = os.path.join(PROJECT_DIR, "chipseq_files/metadata_HepG2_K569_H1.txt")
    cell_line2tfs = defaultdict(list)
    exp_id2lines = defaultdict(list)

    for line in open(metadata_file):
        parts = line.strip().split("\t")
        if not parts[53] == "released":
            continue
        #if there are audit error messages
        if len(parts) > 56:
            continue
        exp_id = parts[6]
        exp_id2lines[exp_id].append(line)

    filtered_lines = []
    for exp_id, lines in exp_id2lines.items():

        if len(lines) == 1:
            filtered_lines += lines
            continue

        lines = [l for l in lines if l.split("\t")[52].startswith("ENCODE4")]
        if not lines:
            continue

        lines = [l for l in lines if l.split("\t")[4] == "conservative IDR thresholded peaks"]
        if not lines:
            continue

        assert len(lines) == 1
        filtered_lines.append(lines[0])

    [filtered_lines.append(l) for l in open(metadata_file) if l.startswith("ENCFF932FVX") or l.startswith("ENCFF180XUM")]

    cell_line2tf2lines = defaultdict(dict)
    for line in filtered_lines:

        parts = line.strip().split("\t")
        tf = parts[21].replace("-human", "")
        cell_line = parts[9]
        cell_line2tfs[cell_line].append(tf)

        if tf not in cell_line2tf2lines[cell_line]:
            cell_line2tf2lines[cell_line][tf] = []

        cell_line2tf2lines[cell_line][tf].append(line)

    cell_line2tf2file_id = defaultdict(dict)

    # with open(os.path.join(PROJECT_DIR, "chipseq_files/cell_line_tf_map.txt"), "w") as of:
    #     of.write("#Cell_line\tTF\tfile_id\tgenome_assembly\tlink\n")
    #     for cell_line in cell_line2tf2lines:
    #         for tf, lines in cell_line2tf2lines[cell_line].items():
    #
    #             if len(lines) == 1:
    #                 most_recent = lines[0]
    #             else:
    #                 most_recent = sorted(lines, key=lambda x: x.split("\t")[28], reverse=True)[0]
    #
    #             file_id = most_recent.split("\t")[0]
    #             assembly = most_recent.split("\t")[5]
    #             of.write("%s\t%s\t%s\t%s\n" % (cell_line, tf, file_id, assembly))
    #             cell_line2tf2file_id[cell_line][tf] = file_id


    for cell_line in cell_line2tf2lines:
        for tf, lines in cell_line2tf2lines[cell_line].items():
            if len(lines) == 1:
                most_recent = lines[0]
            else:
                most_recent = sorted(lines, key=lambda x: x.split("\t")[28], reverse=True)[0]
            file_id = most_recent.split("\t")[0]
            cell_line2tf2file_id[cell_line][tf] = file_id

    return cell_line2tf2file_id



def load_metadata2(metadata_file=None):

    if not metadata_file:
        metadata_file = os.path.join(PROJECT_DIR, "chipseq_files/metadata_HepG2_K569_H1.txt")
    cell_line2tfs = defaultdict(list)
    exp_id2lines = defaultdict(list)

    for line in open(metadata_file):
        parts = line.strip().split("\t")
        if not parts[54] == "released":
            continue
        #if there are audit error messages
        if len(parts) > 56:
            continue
        exp_id = parts[6]
        exp_id2lines[exp_id].append(line)

    filtered_lines = []
    for exp_id, lines in exp_id2lines.items():

        if len(lines) == 1:
            filtered_lines += lines
            continue

        lines = [l for l in lines if l.split("\t")[53].startswith("ENCODE4")]
        if not lines:
            continue

        lines = [l for l in lines if l.split("\t")[4] == "conservative IDR thresholded peaks"]
        if not lines:
            continue

        assert len(lines) == 1

        filtered_lines.append(lines[0])

    # [filtered_lines.append(l) for l in open(metadata_file) if l.startswith("ENCFF932FVX") or l.startswith("ENCFF180XUM")]

    cell_line2tf2lines = defaultdict(dict)
    for line in filtered_lines:

        parts = line.strip().split("\t")
        tf = parts[22].replace("-human", "")
        cell_line = parts[10]
        cell_line2tfs[cell_line].append(tf)

        if tf not in cell_line2tf2lines[cell_line]:
            cell_line2tf2lines[cell_line][tf] = []

        cell_line2tf2lines[cell_line][tf].append(line)

    cell_line2tf2file_id = defaultdict(dict)

    with open(os.path.join(PROJECT_DIR, "chipseq_files/cell_line_tf_map.txt"), "w") as of:
        of.write("#Cell_line\tTF\tfile_id\tgenome_assembly\tlink\n")
        for cell_line in cell_line2tf2lines:
            for tf, lines in cell_line2tf2lines[cell_line].items():

                if len(lines) == 1:
                    most_recent = lines[0]
                else:
                    most_recent = sorted(lines, key=lambda x: x.split("\t")[28], reverse=True)[0]

                file_id = most_recent.split("\t")[0]
                assembly = most_recent.split("\t")[5]
                of.write("%s\t%s\t%s\t%s\n" % (cell_line, tf, file_id, assembly))
                cell_line2tf2file_id[cell_line][tf] = file_id

    return cell_line2tf2file_id



def download_files():

    map_f = os.path.join(PROJECT_DIR, "chipseq_files/cell_line_tf_map.txt")
    file_ids = [l.split()[2] for l in open(map_f) if not l.startswith("#")]

    cmd_fmt = "wget https://www.encodeproject.org/files/%s/@@download/%s.bed.gz"

    cmd = "\n".join([cmd_fmt % (_id, _id) for _id in file_ids]) + "\n"
    with open(join(PROJECT_DIR, "chipseq_files/peaks/download.sh"), "w") as of:
        of.write(cmd)


def liftover():

    d = join(PROJECT_DIR, "chipseq_files/peaks")
    cmd_fmt = "/panfs/pan1/devdcode/common/progs/liftover/liftOver -bedPlus=6 %s /panfs/pan1/devdcode/common/progs/liftover/hg38ToHg19.over.chain %s %s"

    for f in os.listdir(d):
        print(f)
        nf = f.replace(".bed", "_hg19.bed")
        uf = f.replace(".bed", "_unmapped.bed")

        cmd = cmd_fmt % (join(d, f), join(d, nf), join(d, uf))
        subprocess.call(cmd, shell=True)



def get_pioneer_factors():

    """List of pioneer TFs is obtained from:
        - https://www.cell.com/trends/molecular-medicine/fulltext/S1471-4914(19)30019-X
        - https://www.jbc.org/article/S0021-9258(20)30933-9/fulltext
    """

    # tf_names = ["ascl1", "cebp", "ebf1", "esrrb", "foxa", "foxa1", "foxa2", "foxd3", "foxe1", "foxm1", "foxo1", "gata",
    #             "gata1", "gata2", "gata3", "gata4", "gata5", "gata6", "grhl1", "grhl2", "klf4", "mash1", "neurod1",
    #             "nrf1", "oct4", "p53", "pax7", "pbx1", "pou5f1", "pu.1", "pu.i", "sox2", "spi1"]
    # a = load_metadata()
    # for _tf in tf_names:
    #     print("\n%s  :" % _tf)
    #     for cl in ["HepG2", "K562"]:
    #         for k, v in a[cl].items():
    #             if _tf in k.lower():
    #                 print("%s\t%s\t%s" % (cl, k, v))

    hepg2_tfs = {"CEBPA": "ENCFF775MPD", "FOXA1": "ENCFF675TLS", "FOXA2": "ENCFF683KBO", "GATA4": "ENCFF210DMP",
                 "NRF1": "ENCFF003GNY"}

    k562_tfs = {"FOXA1": "ENCFF809YFY", "FOXM1": "ENCFF699ZII", "GATA1": "ENCFF632NQI", "GATA2": "ENCFF732HOE",
                "NEUROD1": "ENCFF581YKY", "NRF1": "ENCFF144PPR", "POU5F1": "ENCFF990CFV", "SPI1": "ENCFF664XPS",
                "SREBF1": "ENCFF985TSH"}

    return hepg2_tfs, k562_tfs


def define_enhancers():

    # dhs_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E116-DNase.macs2.narrowPeak.gz"
    # k27_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E116-H3K27ac.narrowPeak.gz"
    # k4_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E116-H3K4me1.narrowPeak.gz"
    # save_file = os.path.join(data_path, "enhancer_definitions/GM12787.ehancers.bed.gz")

    # dhs_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E118-DNase.macs2.narrowPeak.gz"
    # k27_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E118-H3K27ac.narrowPeak.gz"
    # k4_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E118-H3K4me1.narrowPeak.gz"
    # save_file = os.path.join(data_path, "enhancer_definitions/HepG2.ehancers.bed.gz")

    dhs_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E123-DNase.macs2.narrowPeak.gz"
    k27_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E123-H3K27ac.narrowPeak.gz"
    k4_file = "/panfs/pan1/patternquest/Projects/enhancers/data/roadmap/peaks/E123-H3K4me1.narrowPeak.gz"
    save_file = os.path.join(PROJECT_DIR, "enhancer_definitions/K562.ehancers.bed.gz")

    marked_dhs = BedTool(dhs_file).sort().intersect(k27_file, wa=True).intersect(k4_file, wa=True)

    bed_pool = []
    for r in marked_dhs:
        center = r.start + (r.length/2)
        bed_pool.append("%s\t%d\t%d" % (r.chrom, center - 200, center + 200))

    BedTool("\n".join(bed_pool), from_string=True).merge().subtract(exons).subtract(promoters).sort().saveas(save_file)


    # cell_line = "K562"

    cell_line2tf2file_id = load_metadata()
    chipseqs_dir = os.path.join(PROJECT_DIR, "chipseq_files")
    peaks_dir = os.path.join(PROJECT_DIR, "chipseq_files/peaks")
    tf2file = cell_line2tf2file_id[cell_line]

    # generated by tf_peak_centers()
    # all_peaks_file = os.path.join(chipseqs_dir, "%s_all_tf_centers.bed" % cell_line)
    # all_1k_bindigs_file = os.path.join(chipseqs_dir, "%s_1k_tfs.bed" % cell_line)
    all_peaks_file = os.path.join(chipseqs_dir, "%s_all_tf_centers_8bp.bed" % cell_line)
    all_1k_bindigs_file = os.path.join(chipseqs_dir, "%s_1k_tfs_8bp.bed" % cell_line)

    groupby_handle = open(all_1k_bindigs_file, "w")

    file2tf = {v: k for k, v in tf2file.items()}

    fname_list = tf2file.values()

    tmp_file = os.path.join(PROJECT_DIR, "tmp.%s.8bp.bed" % cell_line)
    tmp_tf_centers = os.path.join(PROJECT_DIR, "tmp.%s.all_tfs.8bp.bed" % cell_line)

    def exclude_tf(tf):
        with open(all_peaks_file) as inf, open(tmp_tf_centers, "w") as outf:
            for line in inf:
                if not line.split("\t")[3] == tf:
                    outf.write(line)

    cnt = 1
    for fname in fname_list:

        print(cnt, len(fname_list), fname)

        tf = file2tf[fname]
        exclude_tf(tf)

        fname = os.path.join(peaks_dir, "%s_simplified.bed.gz" % fname)
        tf_bed = BedTool(fname)

        with open(tmp_file, "w") as of:
            for r in tf_bed:
                center = r.start + r.length / 2
                out_line = "%s\t%d\t%d\t%d\t%d\t%s\n" % (
                r.chrom, center - 500, center + 500, r.start, r.stop, r.fields[3])
                of.write(out_line)

        intersect = BedTool(tmp_file).intersect(tmp_tf_centers, wo=True)
        grp = intersect.groupby(g=[1, 2, 3, 4, 5], c=[10, 11, 12], o=["collapse", "collapse", "collapse"])

        total_grp = grp.count()
        cnt_grp = 0

        for r in grp:

            center = r.start + r.length / 2

            tf_starts = r.fields[6].split(",") + [center - 4]
            tf_stops = r.fields[7].split(",") + [center + 4]
            tfs = r.fields[5] + "," + tf
            num_tfs = len(tf_stops)

            r_array = set(range(r.length))

            tf_pos_pool = []
            for (_start, _stop) in zip(tf_starts, tf_stops):
                _start, _stop = int(_start), int(_stop)
                for _pos in range(_start - r.start, _stop - r.start):
                    tf_pos_pool.append(_pos)

            int_cnt = len(r_array.intersection(tf_pos_pool))
            coverage = float(int_cnt)/r.length

            out_line = "%s\t%d\t%d\t%d\t%f\t%s\n" % (r.chrom, r.start, r.stop, num_tfs, coverage, tfs)
            groupby_handle.write(out_line)

            if cnt_grp % 10000 == 0:
                print("\t", cnt_grp, total_grp)

            cnt_grp += 1

        cnt += 1


def generate_control_regions(tf_id):

    root_dir = "/data/Dcode/sanjar/Projects/overbinders/data/"

    getfasta_cmd = "/home/hudaiber/progs/bedtools/bin/bedtools getfasta -fi %s -bed %s -fo %s"

    peaks_dir = join(root_dir, "peaks/")
    dhs_file = join(root_dir, "DHS_125_merged.bed")

    controls_dir = join(peaks_dir, "%s_controls" % tf_id)
    if not os.path.exists(controls_dir):
        os.mkdir(controls_dir)
    tmp_file = tempfile.NamedTemporaryFile()
    print(controls_dir)

    peak_file = join(peaks_dir, "%s.bed" % tf_id)
    peaks_bed = BedTool(peak_file)

    control_file_1 = join(controls_dir, "controls_1.bed")
    control_file_1_fa = join(controls_dir, "controls_1.fa")
    try:
        peaks_bed.shuffle(g=hg19, incl=dhs_file, excl=peak_file, maxTries=10).saveas(control_file_1)
    except:
        print("Exclude failed. Genertaing control without exlcusion:")
        print(control_file_1)
        peaks_bed.shuffle(g=hg19, incl=dhs_file).saveas(control_file_1)

    cmd = getfasta_cmd % (hg19_fa, control_file_1, control_file_1_fa)
    subprocess.call(cmd, shell=True)
    cmd = "cat %s %s > %s" % (peak_file, control_file_1, tmp_file.name)
    subprocess.call(cmd, shell=True)

    control_file_2 = join(controls_dir, "controls_2.bed")
    control_file_2_fa = join(controls_dir, "controls_2.fa")
    try:
        peaks_bed.shuffle(g=hg19, incl=dhs_file, excl=tmp_file.name, maxTries=10).saveas(control_file_2)
    except:
        print("Exclude failed. Genertaing control without exlcusion:")
        print(control_file_2)
        peaks_bed.shuffle(g=hg19, incl=dhs_file).saveas(control_file_2)
    cmd = getfasta_cmd % (hg19_fa, control_file_2, control_file_2_fa)
    subprocess.call(cmd, shell=True)
    cmd = "cat %s %s %s > %s" % (peak_file, control_file_1, control_file_2, tmp_file.name)
    subprocess.call(cmd, shell=True)

    control_file_3 = join(controls_dir, "controls_3.bed")
    control_file_3_fa = join(controls_dir, "controls_3.fa")
    try:
        peaks_bed.shuffle(g=hg19, incl=dhs_file, excl=tmp_file.name, maxTries=10).saveas(control_file_3)
    except:
        print("Exclude failed. Genertaing control without exlcusion:")
        print(control_file_3)
        peaks_bed.shuffle(g=hg19, incl=dhs_file).saveas(control_file_3)
    cmd = getfasta_cmd % (hg19_fa, control_file_3, control_file_3_fa)
    subprocess.call(cmd, shell=True)


def generate_control_regions_iterate():

    a = load_metadata()

    for cell_line in ["HepG2", "K562"]:

        tf_ids = a[cell_line].values()
        for tf_id in tf_ids:
            tf_id += "_hg19"
            generate_control_regions(tf_id)


def plot_encode_datasets():

    print("here")

    f = '/panfs/pan1/devdcode/sanjar/overbinders/chipseq_files/metadata_all_chipseqs.txt'
    a = load_metadata2(f)

    sorted_kv = sorted(a.items(), key=lambda x:len(x[1]), reverse=True)
    tfs = [len(_[1]) for _ in sorted_kv]
    labels = [_[0] for _ in sorted_kv]

    import matplotlib.pyplot as plt

    x = np.arange(10)
    y = tfs[:10]
    labels = labels[:10]

    plt.figure(figsize=(5, 4))
    bars = plt.bar(x, y, color="grey")
    bars[0].set_color("red")
    bars[1].set_color("red")
    bars[6].set_color("red")
    ax = plt.gca()
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)

    ax.set_title("ENCODE4 ChIP-seq datasets")
    # ax.set_xlabel("Cell-lines")
    ax.set_ylabel("Assayed TFs", fontdict={"size":12})

    plt.tight_layout()

    save_file = "/panfs/pan1/devdcode/sanjar/overbinders/plots/ENCODE4_TFs_barplot.pdf"
    plt.savefig(save_file)


if __name__ == "__main__":

    generate_control_regions(sys.argv[1])
