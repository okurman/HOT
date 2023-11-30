
import os
import sys
from collections import defaultdict
import gzip
from pathlib import Path
from os.path import join
from pybedtools import BedTool
import numpy as np
import gzip
import urllib.request

DATA_PATH = Path("../data/data/")
PHASTCONS_DIR = DATA_PATH / "phastCons"
BINS_DIR = DATA_PATH / "log_bins"

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.gz" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom.gz" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom.gz" % (x, i)) for i in range(14)]


def load_single_chromosome(score_file, pos2score):

    if score_file.endswith(".gz"):
        inf = gzip.open(score_file, "rt")
    else:
        inf = open(score_file, "r")

    for line in inf:
        if line.startswith("fixedStep"):
            local_cnt = 0
            pivot = int(line.split()[2].split("=")[1])
            continue
        pos2score[pivot + local_cnt] = float(line)
        local_cnt += 1

    inf.close()


def load_scores(chrom="chr8", score="phastCons", species="vertebrate", all_chromosomes=False):

    chrom2pos2score = dict()

    if score == "phastCons":
        if species == "vertebrate":
            file_fmt = os.path.join(PHASTCONS_DIR, "vertebrate", "%s.phastCons46way.wigFix.gz")
        elif species == "primates":
            file_fmt = os.path.join(PHASTCONS_DIR, "primates", "%s.phastCons46way.primates.wigFix.gz")
        elif species == "placentalMammals":
            file_fmt = os.path.join(PHASTCONS_DIR, "placentalMammals", "%s.phastCons46way.placental.wigFix.gz")
    elif score == "phyloP":
        pass
        # file_fmt = os.path.join(phyloP_dir, "%s.phyloP46way.wigFix.gz")

    else:
        raise ValueError("score value should be either phastCons or phyloP")

    if all_chromosomes:
        chroms = ["chr%d" % i for i in range(1, 23)]
    elif chrom:
        chroms = [chrom]
    else:
        chroms = []

    for chrom in chroms:

        score_file = file_fmt % chrom
        
        chrom2pos2score[chrom] = defaultdict(int)
        load_single_chromosome(score_file, chrom2pos2score[chrom])

    return chrom2pos2score


def extract_phastcons_for_bins(species="vertebrate"):

    for ch_no in range(22, 0, -1):
        chrom = "chr%d" % ch_no
        print(chrom)

        pos2phylop = load_scores(chrom=chrom, score="phastCons", species=species)[chrom]

        for cell_line in ["H1", "HepG2", "K562"]:
            print(cell_line)

            for loci_files in [get_prom_files(cell_line), get_enh_files(cell_line)]:

                phylop_files = [str(f).replace(".gz", f".{species}.phastcons.gz") for f in loci_files]

                for locus_file, phylop_file in zip(loci_files, phylop_files):

                    if os.path.exists(phylop_file):
                        continue

                    chrom_bed = [r for r in BedTool(locus_file) if r.chrom == chrom]

                    with gzip.open(phylop_file, "at") as outf:
                        for r in chrom_bed:
                            locus_scores = [pos2phylop[p] for p in range(r.start, r.stop)]
                            locus_scores = np.average([s for s in locus_scores])

                            tf_starts = [int(_) for _ in r.fields[6].split(",")]
                            tf_stops = [int(_) for _ in r.fields[7].split(",")]

                            tfbs_positions = []
                            for _start, _stop in set([(a, b) for a, b in zip(tf_starts, tf_stops)]):
                                tfbs_positions += range(_start, _stop)

                            tfbs_sc = np.average([pos2phylop[p] for p in set(tfbs_positions)])
                            out_line = "%s\t%d\t%d\t%f\t%f\n" % (r.chrom, r.start, r.stop, locus_scores, tfbs_sc)
                            outf.write(out_line)


def extract_phastcons_for_HOT_and_buddies(species="vertebrate"):

    src_files = [DATA_PATH/"src_files/HepG2_enhancers_DHS_H3K27ac.bed.gz",
                 DATA_PATH/"src_files/hg19_files/knownGene.exons.merged.bed.gz",
                 DATA_PATH/"HOTs/HepG2_HOTs.bed.gz",
                 DATA_PATH/"HOTs/HepG2_HOTs.proms.bed.gz",
                 DATA_PATH/"HOTs/HepG2_HOTs.noproms.bed.gz"]

    save_files = [str(f).replace(".gz", f".{species}.phastcons.gz") for f in src_files]
    save_files = [gzip.open(_, "wt") for _ in save_files]

    chroms = ["chr%d" % _ for _ in range(1, 23)]
    chroms += ["chrX", "chrY"]

    for chrom in chroms[::-1]:
        print("chrom: " + chrom)
        pos2phylop = load_scores(chrom=chrom, species=species)[chrom]

        for src_file, of in zip(src_files, save_files):

            input_bed = BedTool(str(src_file))
            chrom_bed = BedTool([r for r in input_bed if r.chrom == chrom]).sort()
            for r in chrom_bed:
                scores = [pos2phylop[i] for i in range(r.start, r.stop)]
                scores = np.asarray([_ for _ in scores if _ is not None])
                scores_str = ",".join([str(_) for _ in scores])
                out_line = "%s\t%d\t%d\t%f\t%s\n" % (r.chrom, r.start, r.stop, np.sum(scores) / r.length, scores_str)
                of.write(out_line)
                of.flush()


def check_files_exist(species="vertebrate"):

    chroms = ["chr%d" % _ for _ in range(1, 23)] + ["chrX", "chrY"]
    
    if species == "vertebrate":
        files = [PHASTCONS_DIR / species /f"{chrom}.phastCons46way.wigFix.gz" for chrom in chroms]
        url_prefix = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/"
    else:
        files = [PHASTCONS_DIR / species / f"{chrom}.phastCons46way.{species}.wigFix.gz" for chrom in chroms]
        url_prefix = f"https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/{species}/"

    if not all([f.exists() for f in files]):
        print(f" {species} phastCons files (or some of them) are not present. Proceeding to download from:")
        print(url_prefix)
        print("\n")

        species_dir = PHASTCONS_DIR/species
        species_dir.mkdir(exist_ok=True, parents=True)

        for file in files:
            if not file.exists():
                file_url = url_prefix + file.name
                print(f"Downloading: {file.name}. Saving to: {file}")
                urllib.request.urlretrieve(file_url, file)




if __name__ == "__main__":

    if len(sys.argv) > 1:
        species = sys.argv[1]
        if species not in ["placentalMammals", "vertebrate", "primates"]:
            print("Species need to be one of: placentalMammals, vertebrate, primates")
        sys.exit()

    else:
        species = "vertebrate"

    check_files_exist(species)

    # print(f"Extracting phastCons{species} scores for binned loci.")
    # extract_phastcons_for_bins(species)

    print(f"Extracting phastCons{species} scores for hg19 exons, regular enhancers and HOT loci .")
    extract_phastcons_for_HOT_and_buddies(species)