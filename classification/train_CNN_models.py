import os.path
import sys
from os.path import join

import h5py
from keras.models import Sequential, load_model
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.constraints import max_norm
from keras.regularizers import L1L2
from keras.optimizers import Adadelta
from keras.callbacks import ModelCheckpoint, EarlyStopping, CSVLogger
from sklearn import metrics
import numpy as np
import pandas as pd

MAX_NORM = 0.9
L1 = 5e-07
L2 = 1e-08
EPOCH = 200
BATCH_SIZE = 256
import matplotlib.pyplot as plt


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
    model_file = join(save_dir, "model.hdf5")
    model.save(model_file)

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
    history_log_file = os.path.join(save_dir, "training_history_log.tab")
    history_logger = CSVLogger(filename=history_log_file, separator="\t", append=True)
    _callbacks.append(history_logger)

    history = compiled_model.fit(X_train,
                       Y_train,
                       batch_size=BATCH_SIZE,
                       epochs=EPOCH,
                       validation_data=(X_validation, Y_validation),
                       shuffle=True,
                       callbacks=_callbacks)

    Y_pred = compiled_model.predict(X_test)

    aucs = []
    aucs.append(metrics.roc_auc_score(Y_test, Y_pred))
    aucs.append(metrics.average_precision_score(Y_test, Y_pred))

    with open(join(save_dir, "AUCs.txt"), "w") as of:
        of.write("auROC: %f\n" % aucs[0])
        of.write("auPRC: %f\n" % aucs[1])

    # with h5py.File(join(save_dir, "predictions.hdf5"), "w") as of:
    #     of.create_dataset(name="Y", data=Y_test, compression="gzip")
    #     of.create_dataset(name="Y_pred", data=Y_pred, compression="gzip")

    plot_file = join(save_dir, "training.pdf")
    generate_plot_history(history, plot_file)


def evaluate_model(data_file, run_dir):

    weights_file = join(run_dir, "weights.hdf5")
    model_file = join(run_dir, "model.hdf5")
    model = load_model(model_file)
    model.load_weights(weights_file)

    with h5py.File(data_file, "r") as inf:
        X_test = inf["test_data"][()]
        Y_test = inf["test_labels"][()]

    Y_pred = model.predict(X_test)

    aucs = []
    aucs.append(metrics.roc_auc_score(Y_test, Y_pred))
    aucs.append(metrics.average_precision_score(Y_test, Y_pred))

    with open(join(run_dir, "auc.txt"), "w") as of:
        of.write("auROC: %f\n" % aucs[0])
        of.write("auPRC: %f\n" % aucs[1])

        print("auROC: %f\n" % aucs[0])
        print("auPRC: %f\n" % aucs[1])

    with h5py.File(join(run_dir, "predictions.hdf5"), "w") as of:
        of.create_dataset(name="Y", data=Y_test, compression="gzip")
        of.create_dataset(name="Y_pred", data=Y_pred, compression="gzip")

    history_log_file = os.path.join(run_dir, "training_history_log.tab")
    plot_file = join(run_dir, "training.pdf")
    generate_plot_history(history_log_file, plot_file)


def generate_plot_history(history, save_file, title=None):

    if type(history) == str:
        pass
        history = pd.read_table(history)
    else:
        history = history.history

    loss = history["loss"]
    val_loss = history["val_loss"]
    acc = history["acc"]
    val_acc = history["val_acc"]

    x_range = np.arange(len(loss))

    fig, axs = plt.subplots(1, 2)

    axs[0].plot(x_range, loss, label="train")
    axs[0].plot(x_range, val_loss, label="val")
    axs[0].set_xlabel("Epochs")
    axs[0].set_ylabel("Loss")
    axs[0].set_title("Loss")
    axs[0].legend(loc="upper right")

    axs[1].plot(x_range, acc, label="train")
    axs[1].plot(x_range, val_acc, label="val")
    axs[1].set_xlabel("Epochs")
    axs[1].set_ylabel("Accuracy")
    axs[1].set_title("Accuracy")
    axs[1].legend(loc="upper left")

    if title:
        plt.suptitle(title)

    plt.tight_layout()
    plt.savefig(save_file)
    plt.close()


if __name__ == "__main__":

    data_file = sys.argv[1]
    save_dir = sys.argv[2]
    input_length = int(sys.argv[3])

    model = build_model(input_length=input_length, tasks=1)
    run_model(data_file, model, save_dir)
