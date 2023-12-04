import gzip
import urllib
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import sys
from pathlib import Path
import subprocess as sp
import tempfile
import os
DATA_PATH = Path(os.environ["HOT_DATA"])


def download_hg19_fasta(save_file):

    url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
    print(f"Downloading hg19 fasta file from: {url}")
    temp = tempfile.NamedTemporaryFile()
    urllib.request.urlretrieve(url, temp.name)

    print(f"Uncompressing to: {save_file}")
    cmd = f"zcat {temp.name} > {save_file}"
    sp.call(cmd, shell=True)

    return


def load_metadata(metadata_file=None):

    if not metadata_file:
        metadata_file = DATA_PATH/"src_files/metadata_HepG2_K569_H1.txt"

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

    for cell_line in cell_line2tf2lines:
        for tf, lines in cell_line2tf2lines[cell_line].items():
            if len(lines) == 1:
                most_recent = lines[0]
            else:
                most_recent = sorted(lines, key=lambda x: x.split("\t")[28], reverse=True)[0]
            file_id = most_recent.split("\t")[0]
            cell_line2tf2file_id[cell_line][tf] = file_id

    return cell_line2tf2file_id


def get_chrom2seq(fasta_file):

    if not os.path.exists(fasta_file):
        print(f"Couldn't find the hg10 fasta file: {fasta_file}")
        print(f"Quitting!")
        sys.exit()

    chrom2seq = {}

    inf = gzip.open(fasta_file, "rt") if str(fasta_file).endswith(".gz") else open(fasta_file, "rt")

    for seq in SeqIO.parse(inf, "fasta"):
        chrom2seq[seq.description.split()[0]] = str(seq.seq)

    return chrom2seq


def seq2one_hot(seq):
    d = np.array(['A', 'C', 'G', 'T'])
    m = np.zeros((len(seq), 4), dtype=bool)
    seq = seq.upper()
    for i in range(len(seq)):
        m[i, :] = (d == seq[i])

    return m