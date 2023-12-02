#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")
import os.path
from os.path import join
import h5py
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.metrics import roc_auc_score, average_precision_score
import pickle
import argparse
from pathlib import Path

DATA_PATH = Path("../data")
DATASETS_DIR = DATA_PATH/"classification/datasets/features_datasets"
MODELS_DIR = DATA_PATH/"classification/models/feature_models"
MODELS_DIR.mkdir(exist_ok=True, parents=True)


def train_model_LogReg(cl="HepG2", ctr="proms", max_iter=100):

    data_file = DATASETS_DIR/f"{cl}_binary_{ctr}_400_features.hdf5"
    save_dir = MODELS_DIR/f"{cl}_{ctr}"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    model_file = join(save_dir, "LogReg.model.pk")
    auc_file = join(save_dir, "LogReg.AUC.txt")

    print(f"Model file: {model_file}")
    print(f"AUC file: {auc_file}")

    with h5py.File(data_file, "r") as inf:

        X_train = inf["train_data"][()]
        Y_train = inf["train_labels"][()]

        logreg = LogisticRegression(solver="saga",
                                    max_iter=max_iter,
                                    verbose=0,
                                    n_jobs=10)

        model = logreg.fit(X_train, Y_train)
        pickle.dump(model, open(model_file, 'wb'))

        X_test = inf["test_data"][()]
        Y_test = inf["test_labels"][()]
        Y_pred = model.predict_proba(X_test)

        auROC = roc_auc_score(Y_test, Y_pred[:, 1])
        auPRC = average_precision_score(Y_test, Y_pred[:, 1])

        with open(auc_file, "w") as _of:
            _of.write("auROC: %f\n" % auROC)
            _of.write("auPRC: %f\n" % auPRC)


def train_model_LogReg_partial(cl="HepG2", ctr="proms", max_iter=10**6, feature=0):

    data_file = DATASETS_DIR / f"{cl}_binary_{ctr}_400_features.hdf5"
    save_dir = MODELS_DIR / f"{cl}_{ctr}"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    model_file = join(save_dir, f"LogReg.model.feature_{feature}.pk")
    auc_file = join(save_dir, f"LogReg.AUC.feature.{feature}.txt")

    print(f"Model file: {model_file}")
    print(f"AUC file: {auc_file}")

    with h5py.File(data_file, "r") as inf:

        X_train = inf["train_data"][()][:, feature][..., np.newaxis]
        Y_train = inf["train_labels"][()]

        logreg = LogisticRegression(solver="saga",
                                    max_iter=max_iter,
                                    verbose=0,
                                    n_jobs=10)

        model = logreg.fit(X_train, Y_train)
        pickle.dump(model, open(model_file, 'wb'))

        X_test = inf["test_data"][()][:, feature][..., np.newaxis]
        Y_test = inf["test_labels"][()]
        Y_pred = model.predict_proba(X_test)

        auROC = roc_auc_score(Y_test, Y_pred[:, 1])
        auPRC = average_precision_score(Y_test, Y_pred[:, 1])

        with open(auc_file, "w") as _of:
            _of.write("auROC: %f\n" % auROC)
            _of.write("auPRC: %f\n" % auPRC)

        print(f"\t\tauROC: {auROC}, auPRC: {auPRC}")


def train_model_SVM_partial(cl="HepG2", ctr="proms", kernel="linear", feature=0):

    max_iter = 20*10**6 if ctr == "dhs" else -1

    data_file = DATASETS_DIR / f"{cl}_binary_{ctr}_400_features.hdf5"
    save_dir = MODELS_DIR / f"{cl}_{ctr}"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    model_file = join(save_dir, f"SVM.{kernel}.model.feature.{feature}.pk")
    auc_file = join(save_dir, f"SVM.{kernel}.AUC.feature.{feature}.txt")

    print(f"Model file: {model_file}")
    print(f"AUC file: {auc_file}")

    with h5py.File(data_file, "r") as inf:

        X_train = inf["train_data"][()][:, feature][..., np.newaxis]
        if feature in [1, 2]:
            X_train /= 100
        Y_train = inf["train_labels"][()]

        model = svm.SVC(kernel=kernel,
                        max_iter=max_iter,
                        probability=True,
                        verbose=True,
                        cache_size=1000,
                        random_state=1024)

        model.fit(X_train, Y_train)
        pickle.dump(model, open(model_file, 'wb'))

        X_test = inf["test_data"][()][:, feature][..., np.newaxis]
        if feature in [1, 2]:
            X_test /= 100

        Y_test = inf["test_labels"][()]
        Y_pred = model.predict_proba(X_test)

        auROC = roc_auc_score(Y_test, Y_pred[:, 1])
        auPRC = average_precision_score(Y_test, Y_pred[:, 1])

        with open(auc_file, "w") as _of:
            _of.write("auROC: %f\n" % auROC)
            _of.write("auPRC: %f\n" % auPRC)

        print(f"\tauROC: {auROC}, auPRC: {auPRC}")


def train_model_SVM(cl="HepG2", ctr="proms", kernel="linear"):

    max_iter = 10**6 if ctr == "dhs" else -1

    data_file = DATASETS_DIR / f"{cl}_binary_{ctr}_400_features.hdf5"
    save_dir = MODELS_DIR / f"{cl}_{ctr}"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    model_file = join(save_dir, f"SVM.{kernel}.model.pk")
    auc_file = join(save_dir, f"SVM.{kernel}.AUC.txt")

    print(f"Model file: {model_file}")
    print(f"AUC file: {auc_file}")

    with h5py.File(data_file, "r") as inf:

        X_train = inf["train_data"][()]
        X_train[:, [1, 2]] /= 100
        Y_train = inf["train_labels"][()]

        model = svm.SVC(kernel=kernel, max_iter=max_iter, probability=True, verbose=True)
        model.fit(X_train, Y_train)
        pickle.dump(model, open(model_file, 'wb'))

        X_test = inf["test_data"][()]
        X_test[:, [1, 2]] /= 100
        Y_test = inf["test_labels"][()]
        Y_pred = model.predict_proba(X_test)

        auROC = roc_auc_score(Y_test, Y_pred[:, 1])
        auPRC = average_precision_score(Y_test, Y_pred[:, 1])

        with open(auc_file, "w") as _of:
            _of.write("auROC: %f\n" % auROC)
            _of.write("auPRC: %f\n" % auPRC)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-cl',
                        choices=["HepG2", "K562"],
                        type=str,
                        default="HepG2")
    parser.add_argument('-ctr',
                        choices=["proms", "dhs", "re"],
                        type=str,
                        default="dhs")
    parser.add_argument('-model',
                        choices=["logreg", "svm"],
                        type=str,
                        default="logreg")
    parser.add_argument('-feature',
                        choices=[-1, 0, 1, 2, 3],
                        help="-1 for using all features.",
                        type=int,
                        default=-1)
    parser.add_argument('-kernel',
                        help="SVM kernel to use",
                        choices=["linear", "poly", "rbf", "sigmoid"],
                        type=str,
                        default="linear")

    args = parser.parse_args()

    print("argparse arguments:")
    for arg in vars(args):
        print(f"{arg}:\t{getattr(args, arg)}")
    print("\n")

    if args.model == "logreg":
        if args.feature == -1:
            train_model_LogReg(cl=args.cl, ctr=args.ctr)
        else:
            train_model_LogReg_partial(cl=args.cl, ctr=args.ctr, feature=args.feature)
    else:

        if args.feature == -1:
            train_model_SVM(cl=args.cl, ctr=args.ctr, kernel=args.kernel)
        else:
            train_model_SVM_partial(cl=args.cl, ctr=args.ctr, kernel=args.kernel, feature=args.feature)
