#!//usr/local/bin/python

import subprocess
import os
import tempfile
from os.path import join
from pybedtools import BedTool
from collections import defaultdict


def load_metadata(metadata_file=None):

    if not metadata_file:
        metadata_file = "../data/data/src_files/metadata_HepG2_K569_H1.txt"

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

