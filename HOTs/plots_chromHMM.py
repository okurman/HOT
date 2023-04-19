import math
import warnings
from collections import Counter, defaultdict

from scipy.stats import mannwhitneyu, kstest, kruskal

warnings.filterwarnings('ignore')

import sys
import os
import tempfile
from os.path import join

import networkx
import numpy as np
import configparser

###############################################################
import pandas

config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('enhancers', 'code_path')
sys.path.append(code_path)
###############################################################

from pybedtools import BedTool
import matplotlib.pyplot as plt
import seaborn as sns


PROJECT_DIR = "ROOT_DIR/"
BINS_DIR = "ROOT_DIR/definitions/peak_8bp_v2/log_bins/"
PLOTS_DIR = "ROOT_DIR/plots/log_bins/binwise_pool/"

import plots_data_factory

get_loci_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed" % (x, i)) for i in range(14)]
get_prom_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.prom" % (x, i)) for i in range(14)]
get_enh_files = lambda x: [join(BINS_DIR, "%s_400_loci.%d.bed.noprom" % (x, i)) for i in range(14)]

K562_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '11', '16', '24', '37', '55', '82', '123', '184', '275']
HEPG2_XTICK_LABELS = ['1', '2', '3', '4', '5', '7', '12', '19', '31', '48', '77', '122', '192', '304', '480']
PERC_XTICK_LABELS = ['1', '2', '3', '4', '5', '2%', '3%', '5%', '8%', '12%', '18%', '28%', '42%', '65%', '100%']


def load_chromHMM_data(cl, chromhmm_file):

    # cl = "HepG2"
    # chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmHepg2HMM.bed"
    cl = "K562"
    chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmK562HMM.bed"

    distal_files = get_enh_files(cl)
    distals = BedTool(distal_files[-1]).cat(BedTool(distal_files[-2]), postmerge=False).\
        cat(BedTool(distal_files[-3]), postmerge=False).\
        cat(BedTool(distal_files[-4]), postmerge=False)

    prom_files = get_prom_files(cl)
    proms = BedTool(prom_files[-1]).cat(BedTool(prom_files[-2]), postmerge=False).\
        cat(BedTool(prom_files[-3]), postmerge=False).\
        cat(BedTool(prom_files[-4]), postmerge=False)

    distals = distals.intersect(chromhmm_file, wo=True)
    proms = proms.intersect(chromhmm_file, wo=True)

    d_states = Counter([r.fields[13] for r in distals])
    p_states = Counter([r.fields[13] for r in proms])

    tmp = defaultdict(int)
    for k, v in p_states.items():
        _k = " ".join(k.split("_")[1:])
        tmp[_k] += v
    p_states = tmp

    tmp = defaultdict(int)
    for k, v in d_states.items():
        _k = " ".join(k.split("_")[1:])
        tmp[_k] += v
    d_states = tmp

    main_keys = ['Active Promoter', 'Weak Promoter', 'Strong Enhancer', 'Weak Enhancer']

    tmp = defaultdict(int)
    for k, v in p_states.items():
        if k in main_keys:
            tmp[k] = v
        else:
            tmp["Rest"] += v
    p_states = tmp

    tmp = defaultdict(int)
    for k, v in d_states.items():
        if k in main_keys:
            tmp[k] = v
        else:
            tmp["Rest"] += v
    d_states = tmp

    return p_states, d_states


def plot_piecharts():

    # p_states, d_states = load_chromHMM_data()
    p_states, d_states = get_average_values()

    def _autopct(pct):
        return ('%.1f%%' % pct) if pct > 7 else ''

    fig1, (ax1, ax2) = plt.subplots(1, 2)

    p_sorted = sorted(p_states.items(), key=lambda x: x[1], reverse=True)
    d_sorted = sorted(d_states.items(), key=lambda x: x[1], reverse=True)

    color_map = {
        'Active Promoter': "royalblue",
        'Weak Promoter': "dodgerblue",
        'Strong Enhancer': "brown",
        'Weak Enhancer': "lightsalmon",
        'Rest': "grey"}

    sizes = [_[1] for _ in p_sorted]
    labels = [_[0] for _ in p_sorted]
    patches, texts, _ = ax1.pie(sizes, colors=[color_map[l] for l in labels], autopct=_autopct, counterclock=False,
                             startangle=90)
    ax1.legend(patches, labels, loc='upper center', bbox_to_anchor=(1.0, -0.12), ncol=3, fontsize=8)
    ax1.set_title("Promoters")

    sizes = [_[1] for _ in d_sorted]
    labels = [_[0] for _ in d_sorted]
    ax2.pie(sizes, colors=[color_map[l] for l in labels], autopct=_autopct, counterclock=False, startangle=90)
    ax2.set_title("Distal")

    # ax3.legend(patches, labels, loc="best")

    # save_file = join(PLOTS_DIR, "HepG2_chromHMM_pie.pdf")
    # save_file = join(PLOTS_DIR, "K562_chromHMM_pie.pdf")
    save_file = join(PLOTS_DIR, "Average_chromHMM_pie.pdf")
    print(save_file)
    plt.tight_layout()
    plt.savefig(save_file)
    # plt.show()


def get_average_values():

    cl = "K562"
    chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmK562HMM.bed"
    k562_p_states, k562_d_states = load_chromHMM_data(cl, chromhmm_file)

    cl = "HepG2"
    chromhmm_file = "/panfs/pan1/devdcode/common/chromHMM/wgEncodeBroadHmmHepg2HMM.bed"
    hepg2_p_states, hepg2_d_states = load_chromHMM_data(cl, chromhmm_file)

    p_states = {k:(v+k562_p_states[k])/2 for k,v in hepg2_p_states.items()}
    d_states = {k:(v+k562_d_states[k])/2 for k,v in hepg2_d_states.items()}

    return p_states, d_states