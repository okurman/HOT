import os.path
import sys
from os.path import join

import h5py
from sklearn import metrics
import numpy as np

PROJECT_DIR = ROOT_DIR


def extract_auPRC(cl="HepG2", w=1000):

    pred_file = os.path.join(PROJECT_DIR, f"DL_analysis/dl_runs/{cl}_{w}_14/predictions.hdf5")
    save_file = os.path.join(PROJECT_DIR, f"DL_analysis/dl_runs/{cl}_{w}_14/auPRC.txt")

    with h5py.File(pred_file, "r") as inf:

        Y = inf["Y"][()]
        Y_pred = inf["Y_pred"][()]

        aucs = [metrics.average_precision_score(Y[:, i], Y_pred[:, i]) for i in range(Y_pred.shape[1])]

        np.savetxt(save_file, aucs)

    print(save_file)
