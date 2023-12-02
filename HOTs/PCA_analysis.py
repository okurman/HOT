
import warnings
warnings.filterwarnings('ignore')
import os
import sys
sys.path.append(os.environ["HOT_CODE"])
import pandas
import seaborn
from os.path import join
import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt
from pathlib import Path

import os
DATA_PATH = Path(os.environ["HOT_DATA"])
PLOTS_DIR = DATA_PATH/"plots/"
LOCI_DIR = DATA_PATH/"HOTs/"

from data_prep.basic import load_metadata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


columns = ["x", "y", "PE", "tfs", "ctcf", "rad21", "p300", "polr2"]


def get_loci_matrices():

    cl = "HepG2"
    loci_file = join(LOCI_DIR, f"{cl}_HOTs.bed.gz")
    locus_id2PE = dict()

    for _loci_type, _file in zip(["P", "E"],
                                 [join(LOCI_DIR, f"{cl}_HOTs.proms.bed.gz"),
                                  join(LOCI_DIR, f"{cl}_HOTs.noproms.bed.gz")]):

        locus_id2PE.update({"%s-%d-%d"%(r.chrom, r.start, r.stop): _loci_type for r in BedTool(_file)})

    loci = BedTool(loci_file)

    metadata_file = DATA_PATH/"src_files/metadata_HepG2_K569_H1.txt"
    a = load_metadata(metadata_file)

    tf_ids = list(a[cl].values())
    tf_id2ind = {tf_id: i for i, tf_id in enumerate(tf_ids)}

    M = np.zeros((loci.count(), len(tf_ids)))
    PE = np.zeros((loci.count(), 1))
    num_tfs = np.zeros((loci.count(), 1))

    has_ctcf = np.zeros((loci.count(), 1))
    ctcf_id = a[cl]["CTCF"]

    has_rad21 = np.zeros((loci.count(), 1))
    rad21_id = a[cl]["RAD21"]

    has_p300 = np.zeros((loci.count(), 1))
    p300_id = a[cl]["EP300"]

    has_polr2 = np.zeros((loci.count(), 1))
    polr2_ids = set([a[cl][_] for _ in a[cl] if "POLR2" in _])

    for cnt, l in enumerate(loci):
        _tf_ids = set(l.fields[8].split(","))

        for _tf in _tf_ids:
            _ind = tf_id2ind[_tf]
            M[cnt, _ind] = 1

        locus_id = "%s-%d-%d" % (l.chrom, l.start, l.stop)
        if locus_id2PE[locus_id] == "E":
            PE[cnt] = 1

        num_tfs[cnt] = int(l.fields[5])

        if ctcf_id in _tf_ids:
            has_ctcf[cnt] = 1

        if rad21_id in _tf_ids:
            has_rad21[cnt] = 1

        if p300_id in _tf_ids:
            has_p300[cnt] = 1

        if len(polr2_ids.intersection(_tf_ids))>0:
            has_polr2[cnt] = 1

    info_mat = np.concatenate((PE, num_tfs, has_ctcf, has_rad21, has_p300, has_polr2), axis=1)

    return M, info_mat


def run_PCA(M):

    pca_object = PCA()

    pca_object = pca_object.fit(M)
    explained_variances = pca_object.explained_variance_ratio_
    pca_matrix = pca_object.transform(M)

    return pca_matrix, explained_variances


def load_PCA():

    M, info_mat = get_loci_matrices()
    pca_matrix, explained_variances = run_PCA(M)

    return pca_matrix, info_mat, explained_variances


def run_LogReg(X, Y):

    clf = LogisticRegression(random_state=0).fit(X, Y)
    pred = clf.predict_proba(X)
    auc = roc_auc_score(Y, pred[:, 0])
    auc = np.round(auc, decimals=2)
    c = clf.intercept_[0]
    w1, w2 = clf.coef_.T
    c, m = -c / w2, -w1 / w2

    auc = 1 - auc if auc < 0.5 else auc

    return auc, c, m


def plot_pca_PE(data=None, x_ind=0, y_ind=1, save_file=None):

    if not save_file:
        return

    if not data:
        pca_matrix, info_mat, _ = load_PCA()
        return [pca_matrix, info_mat]

    [pca_matrix, info_mat] = data

    X = pca_matrix[:, [x_ind, y_ind]]
    Y = info_mat[:, 0]

    merged_data = np.concatenate((X, info_mat), axis=1)
    df = pandas.DataFrame(data=merged_data, columns=["x", "y", "PE", "tfs", "ctcf", "rad21", "p300", "polr2"])

    plt.figure(figsize=(3, 3))

    g = seaborn.scatterplot(data=df, x="x", y="y", hue="PE", hue_order=[0, 1], s=10)
    plt.legend([], [], frameon=False)
    l = g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.19), ncol=2, frameon=False, handletextpad=0.001, columnspacing=0.001)
    l.get_texts()[0].set_text('promoters')
    l.get_texts()[1].set_text('enhancers')

    auc, c, m = run_LogReg(X, Y)

    g.set_xlabel("PC%d (AUC=%.2f)" % (x_ind + 1, auc))
    g.set_ylabel("PC%d" % (y_ind + 1))
    plt.tight_layout()

    g.yaxis.set_label_coords(-0.15, 0.5)
    g.xaxis.set_label_coords(0.46, -0.11)

    x_vals = np.array(g.get_xlim())
    if x_ind == 1 and y_ind == 2:
        x_vals = [-3, 2]
    if x_ind == 1 and y_ind == 3:
        x_vals = [-3, 3]
    y_vals = c + m * x_vals
    g.plot(x_vals, y_vals, "--", color="red")

    plt.savefig(save_file, bbox_inches='tight')
    plt.close("all")


def plot_pca_1vs2_tfs(data=None):

    if not data:
        pca_matrix, info_mat, _ = load_PCA()
        return [pca_matrix, info_mat]

    pca_matrix, info_mat = data

    full_ix = np.arange(pca_matrix.shape[0])
    sub_ix = full_ix

    _data = pca_matrix[:, [0, 1]]
    _data = _data[sub_ix, :]
    _info_mat = info_mat[sub_ix, :]

    merged_data = np.concatenate((_data, _info_mat), axis=1)
    df = pandas.DataFrame(data=merged_data, columns=columns)

    plt.figure(figsize=(3, 3))

    g = seaborn.scatterplot(data=df, x="x", y="y", hue="tfs", s=8, legend="brief")
    g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=6, handletextpad=0.001, columnspacing=0.0001, fontsize=7, frameon=False)

    g.set_xlabel("PC1")
    g.set_ylabel("PC2")
    plt.tight_layout()

    g.yaxis.set_label_coords(-0.15, 0.5)
    g.xaxis.set_label_coords(0.46, -0.11)

    save_file = join(PLOTS_DIR, "PCA_1vs2_tfs_full.pdf")
    plt.savefig(save_file)
    plt.close("all")


def plot_pca_param(data=None, tf="polr2", x_ind=0, y_ind=1, save_file=None):

    if not save_file:
        return

    if tf == "polr2":
        tf_ind = 5
    if tf == "p300":
        tf_ind = 4
    if tf == "ctcf":
        tf_ind = 2
    if tf == "rad21":
        tf_ind = 2

    if not data:
        pca_matrix, info_mat, _ = load_PCA()
        return [pca_matrix, info_mat]

    [pca_matrix, info_mat] = data

    X = pca_matrix[:, [x_ind, y_ind]]
    Y = info_mat[:, tf_ind]

    merged_data = np.concatenate((X, info_mat), axis=1)
    df = pandas.DataFrame(data=merged_data, columns=["x", "y", "PE", "tfs", "ctcf", "rad21", "p300", "polr2"])

    plt.figure(figsize=(3, 3))

    g = seaborn.scatterplot(data=df, x="x", y="y", hue=tf, hue_order=[1.0, 0.0], s=10, palette="Set1")
    l = g.legend(loc="upper center", bbox_to_anchor=(0.5, 1.19), ncol=2,
                 frameon=False, handletextpad=0.001, columnspacing=0.001)

    if tf == "rad21":
        l.get_texts()[0].set_text("Cohesin")
        l.get_texts()[1].set_text("non-Cohesin")
    else:
        l.get_texts()[0].set_text(tf.upper())
        l.get_texts()[1].set_text("non-" + tf.upper())

    auc, c, m = run_LogReg(X, Y)

    g.set_xlabel("PC%d (AUC=%.2f)" % (x_ind + 1, auc))
    g.set_ylabel("PC%d" % (y_ind + 1))
    plt.tight_layout()

    g.yaxis.set_label_coords(-0.15, 0.5)
    g.xaxis.set_label_coords(0.46, -0.12)

    x_vals = np.array(g.get_xlim())

    if tf == "polr2":
        if x_ind == 0 and y_ind == 1:
            x_vals = [-5.63659173, 7]
        if x_ind == 1 and y_ind == 3:
            x_vals = [-3, 6]
        if x_ind == 1 and y_ind == 4:
            x_vals = [-5, 6]
        if x_ind == 4 and y_ind == 1:
            x_vals = [-2.5, 3.5]

    if tf == "p300":
        if x_ind == 0 and y_ind == 1:
            x_vals = [-4, 11]
        if x_ind == 0 and y_ind == 2:
            x_vals = [-1, 6]
        if x_ind == 1 and y_ind == 2:
            x_vals = [-3, 4.5]

    if tf == "ctcf":
        if x_ind == 3 and y_ind == 0:
            x_vals = [-1.8, 2.3]

    y_vals = c + m * x_vals
    g.plot(x_vals, y_vals, "--", color="black")

    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
    plt.close("all")


