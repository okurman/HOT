
import warnings
warnings.filterwarnings('ignore')

import gzip
import tempfile
from pathlib import Path
import numpy as np
from pybedtools import BedTool
import subprocess as sp
from basic import load_metadata
import os

# DATA_PATH = Path("../data")
DATA_PATH = Path(os.environ["HOT_DATA"])

# TSS-related annotation files extracted from knownGene table from UCSC Genome Browser database.
# refer to the methods for details.
EXONS_FILE = DATA_PATH/"src_files/hg19_files/knownGene.exons.merged.bed.gz"
PROMOTERS_FILE = DATA_PATH/"src_files/hg19_files/promoters.merged.bed.gz"
INTRONS_FILE = DATA_PATH/"src_files/hg19_files/knownGene.introns.merged.bed.gz"
BLACKLIST_FILE = DATA_PATH/"src_files/hg19_files/hg19-blacklist.v2.bed.gz"
MATADATA_FILE = DATA_PATH/"src_files/metadata_HepG2_K569_H1.txt"

PEAKS_DIR = DATA_PATH/"src_files/peaks/"

HOTS_DIR = DATA_PATH/"HOTs"
HOTS_DIR.mkdir(exist_ok=True)
LOG_BINS_DIR = DATA_PATH/"log_bins"
LOG_BINS_DIR.mkdir(exist_ok=True)
LOCI_DIR = DATA_PATH/"loci"
LOCI_DIR.mkdir(exist_ok=True)


def run_locus_extraction(cell_line="HepG2"):

    print("Concatenating all ChIP-seq peaks (8bp)")
    concat_all_rep_loci(cell_line)
    print("extract_bound_loci")
    extract_bound_loci(cell_line)
    print("add_stats_columns")
    add_stats_columns(cell_line)
    print("extract_to_bins")
    extract_to_log_bins(cell_line)
    print("overlap_promoters")
    overlap_promoters(cell_line)
    print("Extract HOT loci")
    extract_HOT_loci(cell_line)


def concat_all_rep_loci(cell_line="H1"):

    cell_line2tf2file_id = load_metadata(MATADATA_FILE)
    all_files = cell_line2tf2file_id[cell_line].values()

    tmp_file = tempfile.NamedTemporaryFile()

    for i, d in enumerate(all_files):
        rep_file = PEAKS_DIR/("%s_hg19_peak_8bp.bed.gz" % d)
        cmd = "zcat %s >> %s" % (rep_file, tmp_file.name)
        sp.call(cmd, shell=True)

    save_file = LOCI_DIR/("%s_all_tfbs.bed.gz" % cell_line)

    print("Saving to:", save_file)
    BedTool(tmp_file.name).sort().groupby(g=[1, 2, 3], c=[4, 5], o="collapse").saveas(save_file)


def extract_bound_loci(cell_line="H1"):

    all_rep_loci_file = LOCI_DIR/("%s_all_tfbs.bed.gz" % cell_line)

    windows_file = DATA_PATH/"src_files/windows_400_hg19.bed.gz"
    windows_grp_file = LOCI_DIR/("%s_windows_400_groupby.bed.gz" % cell_line)

    print(f"Saving to: {windows_grp_file}")

    BedTool(windows_file).\
        intersect(str(all_rep_loci_file), wo=True).\
        groupby(g=[1, 2, 3], c=[5, 6, 7, 8], o="collapse").\
        saveas(windows_grp_file)


def add_stats_columns(cell_line="H1"):

    windows_grp_file = LOCI_DIR/("%s_windows_400_groupby.bed.gz" % cell_line)
    save_file = LOCI_DIR/("%s_400_loci.bed.gz" % cell_line)

    print(f"Saving to: {save_file}")

    with gzip.open(save_file, "wt") as of:

        header_line = "#chr\tstart\tstop\tcov\tuq_tfbs\tuq_tfs\ttf_starts\ttf_stops\tencode_ids\tsignal_values\n"
        of.write(header_line)

        out_fmt = "%s\t%d\t%d\t%f\t%d\t%d\t%s\t%s\t%s\t%s\n"

        for r in BedTool(windows_grp_file):

            total_range = set([_ for _ in range(r.start, r.stop)])

            tf_ranges = []
            tf_starts = [int(_) for _ in r.fields[3].split(",")]
            tf_stops = [int(_) for _ in r.fields[4].split(",")]
            for (_start, _stop) in zip(tf_starts, tf_stops):
                tf_ranges += [_ for _ in range(_start, _stop)]

            overlap_cnt = len(total_range.intersection(tf_ranges))
            cov = np.round(overlap_cnt / len(total_range), 5)

            tfs = len(tf_starts)
            tfs_uq = len(set(r.fields[5].split(",")))

            out_line = out_fmt % (r.chrom,
                                  r.start,
                                  r.stop,
                                  cov,
                                  tfs,
                                  tfs_uq,
                                  r.fields[3],
                                  r.fields[4],
                                  r.fields[5],
                                  r.fields[6])

            of.write(out_line)


def extract_to_log_bins(cell_line="H1"):

    locus_file = LOCI_DIR/("%s_400_loci.bed.gz" % cell_line)
    lines = gzip.open(locus_file, "rt").readlines()
    tfs_column = [int(_.split("\t")[5]) for _ in lines[1:]]

    # adding the edges to the last bins
    if cell_line == "H1":
        max_tfs = max(tfs_column) + 1
    else:
        max_tfs = max(tfs_column) + 10

    bin_edges = [1, 2, 3, 4] + list(np.logspace(np.log10(5), np.log10(max_tfs), 11, dtype=int))
    bin_edges = np.asarray(bin_edges)

    fnames = ["%s_400_loci.%d.bed.gz" % (cell_line, i) for i in range(bin_edges.shape[0])]

    lines_list = [[] for _ in range(np.size(bin_edges, 0))]

    # skip the header
    for l in lines[1:]:
        parts = l.split()
        tfs = int(parts[5])
        ind = np.argmax((bin_edges - tfs) > 0) - 1
        lines_list[ind].append(l)

    for _lines, _fname in zip(lines_list, fnames):
        if not _lines:
            continue
        print(_fname)
        with gzip.open(LOG_BINS_DIR/_fname, "wt") as of:
            for line in _lines:
                of.write(line)


def overlap_promoters(cell_line="H1"):

    fnames = LOG_BINS_DIR.glob("%s_400_loci.*.bed.gz" % cell_line)
    tmp_file = tempfile.NamedTemporaryFile()

    for loci_file in fnames:
        loci_bed = BedTool(loci_file)
        noprom_file = str(loci_file).replace(".gz", ".noprom.gz")
        prom_file = str(loci_file).replace(".gz", ".prom.gz")
        loci_bed.sort().intersect(str(PROMOTERS_FILE), v=True).sort().saveas(tmp_file.name)
        cmd = "uniq %s > %s" % (tmp_file.name, noprom_file)
        sp.call(cmd, shell=True)
        loci_bed.sort().intersect(str(PROMOTERS_FILE), wa=True).sort().saveas(tmp_file.name)
        cmd = "uniq %s > %s" % (tmp_file.name, prom_file)
        sp.call(cmd, shell=True)
        print(f"Bin promoter file:{prom_file}")
        print(f"Bin enhancer file: {noprom_file}")


def extract_HOT_loci(cl="H1"):

    # the last four bins
    hot_bins = [10, 11, 12, 13]
    _files = [LOG_BINS_DIR/f"{cl}_400_loci.{i}.bed.gz" for i in hot_bins]
    hot_loci = BedTool(_files[0]).\
        cat(BedTool(_files[1]), postmerge=False).\
        cat(BedTool(_files[2]), postmerge=False).\
        cat(BedTool(_files[3]), postmerge=False)

    hot_loci = hot_loci.intersect(str(BLACKLIST_FILE), v=True, wa=True).sort()

    save_file = HOTS_DIR/f"{cl}_HOTs.bed.gz"
    hot_loci.saveas(save_file)
    
    hot_proms = hot_loci.intersect(str(PROMOTERS_FILE), wa=True, u=True)
    save_file = HOTS_DIR / f"{cl}_HOTs.proms.bed.gz"
    hot_proms.saveas(save_file)

    hot_enhs = hot_loci.intersect(str(PROMOTERS_FILE), wa=True, v=True)
    save_file = HOTS_DIR / f"{cl}_HOTs.noproms.bed.gz"
    hot_enhs.saveas(save_file)

    hot_proms_nearest_tss = hot_proms.closest(str(PROMOTERS_FILE), d=True).sort().groupby(c=[15,16], o="distinct")
    save_file = HOTS_DIR / f"{cl}_HOTs.proms.nearest_tss.bed.gz"
    hot_proms_nearest_tss.saveas(save_file)

    hot_enhs_nearest_tss = hot_enhs.closest(str(PROMOTERS_FILE), d=True).sort().groupby(c=[15, 16], o="distinct")
    save_file = HOTS_DIR / f"{cl}_HOTs.noproms.nearest_tss.bed.gz"
    hot_enhs_nearest_tss.saveas(save_file)


if __name__ == "__main__":

    for cl in ["H1", "HepG2", "K562"]:
        print(f"Extracting loci for: {cl}")
        run_locus_extraction(cl)
        print("\n\n")

    print("Done!")