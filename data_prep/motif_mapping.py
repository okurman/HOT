import glob
import shutil
import subprocess
import sys
import os
from collections import defaultdict, OrderedDict, Counter
from os.path import join

import networkx
import numpy as np
import configparser

###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
from lib import phyloP
from lib import phastCons

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

PROJECT_DIR = "/panfs/pan1/devdcode/sanjar/overbinders/"

from overbinders.data_prep.basic import load_metadata


def extract_motifs():

    a = load_metadata()

    f = "/panfs/pan1/devdcode/sanjar/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    m2id = {l.split("\t")[1].upper().rstrip(): l.split("\t")[0].upper()[1:] for l in open(f) if l.startswith(">")}

    all_tfs = list(set(list(a["HepG2"].keys()) + list(a["K5622"].keys())))

    # found_cnt = 0
    # for i, tf in enumerate(all_tfs):
    #     # _motifs = ["%s\t%s" %(k,v) for (k,v) in m2id.items() if tf in k]
    #     _motifs = [k for (k,v) in m2id.items() if tf in k]
    #     if _motifs:
    #         print(i, tf, _motifs)
    #         found_cnt += 1
    # print(len(all_tfs), found_cnt)

    tf2jaspar_names = {"REST": ['REST'], "MAFK": ['MAFK'], "NR5A1": ['NR5A1'], "NFYC": ['NFYC'],
                      "THRB": ['THRB', 'THRB(VAR.2)', 'THRB(VAR.3)'], "ISX": ['ISX'], "NRF1": ['NRF1'], "TBP": ['TBP'],
                      "RARA": ['RARA', 'RARA(VAR.2)'], "MAFF": ['MAFF'], "ELF1": ['ELF1'], "STAT6": ['STAT6'],
                      "NFAT5": ['NFAT5'], "SMAD3": ['SMAD3'], "JUND": ['JUND', 'JUND(VAR.2)'], "FOXK2": ['FOXK2'],
                      "TBX3": ['TBX3'], "ZNF384": ['ZNF384'], "ZNF24": ['ZNF24'], "TEAD1": ['TEAD1'],
                      "PITX1": ['PITX1'], "SP1": ['SP1'], "ATF2": ['ATF2'], "CREM": ['CREM'], "DLX6": ['DLX6'],
                      "ZBTB7A": ['ZBTB7A'], "TGIF2": ['TGIF2'], "USF2": ['USF2'], "MIXL1": ['MIXL1'],
                      "KLF12": ['KLF12'], "HNF4A": ['HNF4A', 'HNF4A(VAR.2)'], "CEBPA": ['CEBPA', ],
                      "ZKSCAN5": ['ZKSCAN5'], "CEBPB": ['CEBPB'], "RREB1": ['RREB1'], "SNAI1": ['SNAI1'],
                      "AHR": ['AHR::ARNT'], "FOXO1": ['FOXO1'], "ONECUT1": ['ONECUT1'], "ZBTB26": ['ZBTB26'],
                      "ZNF740": ['ZNF740'], "CREB3": ['CREB3'], "RFX5": ['RFX5'], "SIX1": ['SIX1'], "MEF2A": ['MEF2A'],
                      "MAX": ['MAX'], "ZSCAN29": ['ZSCAN29'], "PRDM15": ['PRDM15'], "MEF2D": ['MEF2D'],
                      "MAFG": ['MAFG'], "SOX5": ['SOX5'], "PPARG": ['PPARG'], "LBX2": ['LBX2'], "ARNT": ['ARNT'],
                      "ZBTB33": ['ZBTB33'], "RFX3": ['RFX3'], "BCL6": ['BCL6'], "IRF3": ['IRF3'], "SP2": ['SP2'],
                      "CTCF": ['CTCF', 'CTCFL'], "NFIA": ['NFIA'], "FOXP1": ['FOXP1'], "GABPA": ['GABPA'],
                      "NFE2L1": ['NFE2L1'], "KLF11": ['KLF11'], "TP53": ['TP53'], "ELF4": ['ELF4'], "SOX13": ['SOX13'],
                      "NFIL3": ['NFIL3'], "TFDP1": ['TFDP1'], "THAP11": ['THAP11'], "MNT": ['MNT'],
                      "ARID3A": ['ARID3A'], "HNF1A": ['HNF1A'], "ZFX": ['ZFX'], "BHLHE40": ['BHLHE40'],
                      "ZBTB14": ['ZBTB14'], "TCF7": ['TCF7'], "RXRB": ['RXRB', 'RXRB(VAR.2)'],
                      "DDIT3": ['DDIT3::CEBPA'], "ELF3": ['ELF3'], "STAT5B": ['STAT5B'], "ETV5": ['ETV5'],
                      "TEF": ['TEF'], "SRF": ['SRF'], "HNF1B": ['HNF1B'], "ZNF274": ['ZNF274'],
                      "CEBPG": ['CEBPG', 'CEBPG(VAR.2)'], "CEBPD": ['CEBPD'], "FOXA1": ['FOXA1'], "E2F4": ['E2F4'],
                      "ATF7": ['ATF7'], "ATF3": ['ATF3'], "RFX1": ['RFX1'], "GFI1": ['GFI1'], "NFKB2": ['NFKB2'],
                      "ZNF460": ['ZNF460'], "MYBL2": ['MYBL2'], "USF1": ['USF1'], "NFYA": ['NFYA'], "MTF1": ['MTF1'],
                      "NFATC3": ['NFATC3'], "JUN": ['JUN'], "HNF4G": ['HNF4G'], "ZNF382": ['ZNF382'],
                      "NFIC": ['NFIC', 'NFIC(VAR.2)'], "ONECUT2": ['ONECUT2'], "CREB1": ['CREB1'], "IRF2": ['IRF2'],
                      "ZNF652": ['ZNF652'], "TCF7L2": ['TCF7L2'], "FOXA3": ['FOXA3'], "ZBTB7B": ['ZBTB7B'],
                      "ZEB1": ['ZEB1'], "FOXK1": ['FOXK1'], "ATF6": ['ATF6'], "HLF": ['HLF'],
                      "RORA": ['RORA', 'RORA(VAR.2)'], "FOXA2": ['FOXA2'],
                      "RARG": ['RARG', 'RARG(VAR.2)', 'RARG(VAR.3)'], "ESRRA": ['ESRRA'], "HINFP": ['HINFP'],
                      "NRL": ['NRL'], "NFE2L2": ['NFE2L2'], "ZNF282": ['ZNF282'], "MNX1": ['MNX1'], "PBX2": ['PBX2'],
                      "SOX18": ['SOX18'], "KLF9": ['KLF9'], "ZNF317": ['ZNF317'], "ISL2": ['ISL2'],
                      "ZNF263": ['ZNF263'], "ATF1": ['ATF1'], "NR3C1": ['NR3C1'], "TCF12": ['TCF12', 'TCF12(VAR.2)'],
                      "NKX3-1": ['NKX3-1'], "TEAD4": ['TEAD4'], "TFE3": ['TFE3'], "YY1": ['YY1'], "E2F1": ['E2F1'],
                      "RELA": ['RELA'], "TFAP4": ['TFAP4', 'TFAP4(VAR.2)'], "ETV6": ['ETV6'], "IRF5": ['IRF5'],
                      "EGR1": ['EGR1'], "TBX2": ['TBX2'], "SMAD4": ['SMAD4'], "ZNF281": ['ZNF281'], "IRF1": ['IRF1'],
                      "GATA4": ['GATA4'], "GMEB1": ['GMEB1'], "ETS1": ['ETS1'], "HOXA5": ['HOXA5'], "MXI1": ['MXI1'],
                      "NFYB": ['NFYB'], "NR2C2": ['NR2C2', 'NR2C2(VAR.2)'], "ETV4": ['ETV4'], "E2F2": ['E2F2'],
                      "NR2F2": ['NR2F2'], "KLF13": ['KLF13'], "HSF2": ['HSF2'], "LIN54": ['LIN54'],
                      "ZKSCAN1": ['ZKSCAN1'], "GMEB2": ['GMEB2'], "MLX": ['MLX'], "FOSL2": ['FOSL2'],
                      "FOXJ3": ['FOXJ3'], "TCF3": ['TCF3'], "ARNTL": ['ARNTL'], "KLF16": ['KLF16'],
                      "NR2F1": ['NR2F1', 'NR2F1(VAR.2)', 'NR2F1(VAR.3)'], "GATA2": ['GATA2'], "NFE2": ['NFE2'],
                      "RBPJ": ['RBPJ'], "XBP1": ['XBP1'], "SOX6": ['SOX6'], "MAZ": ['MAZ'], "ZNF143": ['ZNF143'],
                      "ELK1": ['ELK1'], "RXRA": ['RXRA'], "NR2F6": ['NR2F6', 'NR2F6(VAR.2)', 'NR2F6(VAR.3)'],
                      "ERF": ['ERF']}

    # id2m = {v:k for (k,v) in m2id.items()}
    tf2jaspar_ids = {tf: [m2id[n] for n in names] for (tf, names) in tf2jaspar_names.items()}

    # for (k,v) in tf2jaspar_ids.items():
    #     print(k, v)

    return tf2jaspar_ids


def split_fimo_homer():

    a = load_metadata()

    tf2jaspar_ids = extract_motifs()

    hepg2_tfs = list(a["HepG2"].keys())
    hepg2_fimo_tfs = [tf for tf in hepg2_tfs if tf in tf2jaspar_ids]
    hepg2_homer_tfs = set(hepg2_tfs).difference(hepg2_fimo_tfs)
    k562_tfs = list(a["K562"].keys())
    k562_fimo_tfs = [tf for tf in k562_tfs if tf in tf2jaspar_ids]
    k562_homer_tfs = set(k562_tfs).difference(k562_fimo_tfs)

    hepg2_fimo_ids = [a["HepG2"][tf] for tf in hepg2_fimo_tfs]
    hepg2_homer_ids = [a["HepG2"][tf] for tf in hepg2_homer_tfs]
    k562_fimo_ids = [a["K562"][tf] for tf in k562_fimo_tfs]
    k562_homer_ids = [a["K562"][tf] for tf in k562_homer_tfs]

    split_map = {"HepG2": [hepg2_fimo_ids, hepg2_homer_ids],
                 "K562": [k562_fimo_ids, k562_homer_ids]}

    assert len(hepg2_tfs) == len(hepg2_fimo_tfs) + len(hepg2_homer_tfs)
    assert len(k562_tfs) == len(k562_fimo_tfs) + len(k562_homer_tfs)

    return split_map


def master_fimo_runs():

    split_map = split_fimo_homer()
    tf2jaspar_ids = extract_motifs()

    a = load_metadata()
    id2tf = {}
    for cl in ["HepG2", "K562"]:
        for k, v in a[cl].items():
            id2tf[v] = k

    peaks_dir = join(PROJECT_DIR, "chipseq_files/peaks")
    fimo_dir = join(PROJECT_DIR, "chipseq_files/JASPAR_fimo_runs/")
    motifs_dir = "/panfs/pan1/devdcode/sanjar/motif_databases/JASPAR/PFMs/"

    fimo_tf_ids = split_map["HepG2"][0] + split_map["K562"][0]

    print(len(fimo_tf_ids))

    def get_fasta(tf_id):
        hg19 = "/panfs/pan1/patternquest/data/genomes/hg19/hg19.fa"
        peak_file = join(peaks_dir, "%s_hg19.bed" % tf_id)
        tf_dir = join(fimo_dir, tf_id)
        save_file = join(tf_dir, "%s.fa" % tf_id)
        cmd = "bedtools getfasta -fi %s -bed %s -fo %s" % (hg19, peak_file, save_file)
        subprocess.call(cmd, shell=True)

    for i, tf_id in enumerate(fimo_tf_ids):

        tf_dir = join(fimo_dir, tf_id)
        if not os.path.exists(tf_dir):
            os.mkdir(tf_dir)

        print(i, len(fimo_tf_ids), tf_id)

        # get_fasta(tf_id)
        # continue
        # break

        for motif_id in tf2jaspar_ids[id2tf[tf_id]]:
            print(motif_id)
            motif_file = join(motifs_dir, "%s.meme" % motif_id)

            tf_motif_dir = join(tf_dir, motif_id)
            if os.path.exists(tf_motif_dir):
                shutil.rmtree(tf_motif_dir)

            cmd = "fimo -o %s %s %s >%s 2>%s" % (
                tf_motif_dir,
                motif_file,
                join(tf_dir, "%s.fa" % tf_id),
                join(tf_dir, "log_%s.out" % motif_id),
                join(tf_dir, "log_%s.err" % motif_id)
            )
            subprocess.call(cmd, shell=True)


def annotate_meme(prefix, meme_dirs):

    src_dir = join(PROJECT_DIR, "chipseq_files/meme/results/")
    peaks_dir = join(PROJECT_DIR, "chipseq_files/peaks/")

    report_file = join(PROJECT_DIR, "chipseq_files/meme/%s_annotated_peaks_stats_q5.txt" % prefix)
    rf = open(report_file, "w")

    # meme_dirs = os.listdir(src_dir)

    for cnt, d in enumerate(meme_dirs):

        concat_fimo_file = join(src_dir, d, "all_motifs_q5.bed")

        with open(concat_fimo_file, "w") as of:
            for root, dirs, files in os.walk(join(src_dir, d)):
                for f in files:
                    if not f == "fimo.tsv":
                        continue
                    fimo_file = join(root, files[0])
                    for l in open(fimo_file):
                        if l.startswith("motif_id") or l.startswith("#") or not l.strip():
                            continue
                        parts = l.split("\t")
                        qvalue = float(parts[8])
                        if qvalue > 0.05:
                            continue
                        out_line = "%s\t%s\t%s\t%s\t%s\n" % (parts[2], parts[3], parts[4], parts[7], parts[8])
                        of.write(out_line)

        peaks_file = join(peaks_dir, "%s.bed" % d)
        save_file = join(src_dir, d, "annotated_peaks_q5.bed")

        peaks_bed = BedTool(peaks_file).groupby(c=[3], o="collapse").sort()
        fimo_bed = BedTool(concat_fimo_file).sort()

        if fimo_bed.count() == 0:
            print(cnt, len(dirs), d, "None", "None", 0)
            continue

        grpb = peaks_bed.intersect(fimo_bed, F=1.0, wo=True).groupby(c=[6, 7, 8, 9], o="collapse")

        out_pool = []
        for r in grpb:
            [starts, stops, pvalues, qvalues] = [r.fields[i].split(",") for i in range(3, 7)]
            qvalues = list(map(lambda x: float(x), qvalues))
            ind = np.argsort(qvalues)[0]
            out_line = [r.chrom, r.start, r.stop, starts[ind], stops[ind], pvalues[ind], qvalues[ind]]
            out_line = "\t".join(map(lambda x: str(x), out_line))
            out_pool.append(out_line)

        annotated_peaks_bed = BedTool("\n".join(out_pool), from_string=True)
        total_peaks = peaks_bed.count()
        annotated_peaks = annotated_peaks_bed.count()
        print(cnt, len(meme_dirs), d, total_peaks, annotated_peaks, annotated_peaks/total_peaks)

        stats_terms = [d, total_peaks, annotated_peaks, annotated_peaks/total_peaks]
        stats_terms = "\t".join([str(_) for _ in stats_terms])+"\n"
        rf.write(stats_terms)
        rf.flush()

        annotated_peaks_bed.saveas(save_file)


def concat_annotated_peaks():

    a = load_metadata()
    for cell_line in ["HepG2", "K562"]:
        print(cell_line)
        _dirs = a[cell_line].values()
        save_file = join(PROJECT_DIR, "chipseq_files/meme", "%s_concat.bed" % cell_line)
        with open(save_file, "w") as outf:
            for cnt, _d in enumerate(_dirs):
                cur_file = join(PROJECT_DIR, "chipseq_files/meme/results/%s_hg19/annotated_peaks_q5.bed" % _d)
                if not os.path.exists(cur_file):
                    continue
                [outf.write(l.replace("\n", "\t%s\n" % _d)) for l in open(cur_file)]
                print(cnt, _d)


if __name__ == "__main__":

    inf = sys.argv[1]
    meme_dirs = [l.strip() for l in open(inf).readlines()]

    annotate_meme(os.path.basename(inf), meme_dirs)