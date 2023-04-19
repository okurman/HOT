import os.path
import sys
from os.path import join

import h5py
from keras.models import Sequential
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.constraints import max_norm
from keras.regularizers import L1L2
from keras.optimizers import Adam, Adadelta
from keras.callbacks import ModelCheckpoint, EarlyStopping
from sklearn import metrics
from sklearn.metrics import roc_curve
import numpy as np

MAX_NORM = 0.9
L1 = 5e-07
L2 = 1e-08
EPOCH = 200
BATCH_SIZE = 256

PROJECT_DIR = ROOT_DIR


def build_model(input_length=400, tasks=3):

    model = Sequential()
    model.add(Conv1D(input_shape=(input_length, 4),
                     filters=480,
                     kernel_size=9,
                     strides=1,
                     activation="relu",
                     kernel_constraint=max_norm(MAX_NORM),
                     kernel_regularizer=L1L2(L1, L2)))
    model.add(MaxPooling1D(pool_size=9, strides=3))
    model.add(Dropout(0.2))

    model.add(Conv1D(input_shape=(None, 480),
                     filters=480,
                     kernel_size=4,
                     strides=1,
                     activation="relu",
                     kernel_constraint=max_norm(MAX_NORM),
                     kernel_regularizer=L1L2(L1, L2)))
    model.add(MaxPooling1D(pool_size=4, strides=2))
    model.add(Dropout(0.2))

    model.add(Conv1D(input_shape=(None, 480),
                     filters=240,
                     kernel_size=4,
                     strides=1,
                     activation="relu",
                     kernel_constraint=max_norm(MAX_NORM),
                     kernel_regularizer=L1L2(L1, L2)))
    model.add(MaxPooling1D(pool_size=4, strides=2))
    model.add(Dropout(0.2))

    model.add(Conv1D(input_shape=(None, 240),
                     filters=320,
                     kernel_size=4,
                     strides=1,
                     activation="relu",
                     kernel_constraint=max_norm(MAX_NORM),
                     kernel_regularizer=L1L2(L1, L2)))
    model.add(MaxPooling1D(pool_size=4, strides=2))

    model.add(Flatten())
    model.add(Dense(units=180, activation="relu"))
    model.add(Dense(units=tasks))

    model.add(Activation("sigmoid"))

    return model


def run_model(data_file, model, save_dir):

    weights_file = join(save_dir, "weights.hdf5")

    # Adadelta is recommended to be used with default values
    opt = Adadelta()

    compiled_model = model
    compiled_model.compile(loss='binary_crossentropy', optimizer=opt, metrics=["accuracy"])

    with h5py.File(data_file, "r") as inf:

        X_train = inf["train_data"][()]
        Y_train = inf["train_labels"][()]
        X_validation = inf["validation_data"][()]
        Y_validation = inf["validation_labels"][()]
        X_test = inf["test_data"][()]
        Y_test = inf["test_labels"][()]

    _callbacks = []
    checkpoint = ModelCheckpoint(filepath=weights_file, save_best_only=True, verbose=1)
    _callbacks.append(checkpoint)
    earlystopping = EarlyStopping(monitor="val_loss", patience=15, verbose=1)
    _callbacks.append(earlystopping)

    compiled_model.fit(X_train,
                       Y_train,
                       batch_size=BATCH_SIZE,
                       epochs=EPOCH,
                       validation_data=(X_validation, Y_validation),
                       shuffle=True,
                       callbacks=_callbacks)

    Y_pred = compiled_model.predict(X_test)

    aucs = [metrics.roc_auc_score(Y_test[:, i], Y_pred[:, i]) for i in range(Y_pred.shape[1])]

    np.savetxt(join(save_dir, "auc.txt"), aucs)
    
    with h5py.File(join(save_dir, "predictions.hdf5"), "w") as of:
        of.create_dataset(name="Y", data=Y_test)
        of.create_dataset(name="Y_pred", data=Y_pred)



def launch_train(cl, input_length, tasks):

    if tasks == 3:

        data_file = join(PROJECT_DIR, "data/DL_analysis/datasets/%s_3class_%d.hdf5" % (cl, input_length))
        dl_dir = join(PROJECT_DIR, "data/DL_analysis/dl_runs/%s_%d" % (cl, input_length))
        if not os.path.exists(dl_dir):
            os.mkdir(dl_dir)
        model = build_model(input_length=input_length)

    else:

        data_file = join(PROJECT_DIR, "data/DL_analysis/datasets/%s_14class_%d.hdf5" % (cl, input_length))
        dl_dir = join(PROJECT_DIR, "data/DL_analysis/dl_runs/%s_%d_14" % (cl, input_length))
        if not os.path.exists(dl_dir):
            os.mkdir(dl_dir)
        model = build_model(input_length=input_length, tasks=14)

    # print "Launching training:"
    # print data_file
    # print dl_dir
    run_model(data_file, model, dl_dir)



if __name__ == "__main__":

    launch_train(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))