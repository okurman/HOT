import os
from os.path import join

CELL_LINES = ["HepG2", "K562"]
LENS = [1000, 400]

# noinspection SmkAvoidTabWhitespace

SAVE_DIR = os.environ["HOT_DATA"] + "/classification/datasets/"
FASTA_DATASET_FILE_NAMES = ["controls_test.fa", "controls_train.fa", "hots_test.fa", "hots_train.fa"]


rule all:
	input:
		# fasta datasets
		expand(SAVE_DIR + "fasta_datasets/{cell_line}_binary_{region}_{_len}_lsgkm/{fname}",
		   cell_line=CELL_LINES,
		   region=["dhs", "proms", "re"],
		   _len=LENS,
		   fname=FASTA_DATASET_FILE_NAMES),
		# one-hot datasets
		expand(SAVE_DIR + "one_hot_datasets/{cell_line}_binary_{region}_{_len}.hdf5",
		   cell_line=CELL_LINES,
		   region=["dhs", "proms", "re"],
		   _len=LENS),
		# feature datasets
		expand(SAVE_DIR + "features_datasets/{cell_line}_binary_{region}_400_features.hdf5",
		   cell_line=CELL_LINES,
		   region=["dhs", "proms", "re"])


rule create_datasets_hepg2:
	output:
		expand(SAVE_DIR + "fasta_datasets/HepG2_binary_{region}_{_len}_lsgkm/{fname}",
			region=["dhs", "proms", "re"],
			_len=LENS,
			fname=FASTA_DATASET_FILE_NAMES),
		# one-hot datasets
		expand(SAVE_DIR + "one_hot_datasets/HepG2_binary_{region}_{_len}.hdf5",
			region=["dhs", "proms", "re"],
			_len=LENS),
		# feature datasets
		expand(SAVE_DIR + "features_datasets/HepG2_binary_{region}_400_features.hdf5",
			region=["dhs", "proms", "re"])
	threads: 4
	params:
		script = join(os.environ["HOT_CODE"], "data_prep/classification/create_datasets.py")
	run:
		shell("python {params.script} HepG2")



rule create_datasets_k562:
	output:
		expand(SAVE_DIR + "fasta_datasets/K562_binary_{region}_{_len}_lsgkm/{fname}",
			region=["dhs", "proms", "re"],
			_len=LENS,
			fname=FASTA_DATASET_FILE_NAMES),
		# one-hot datasets
		expand(SAVE_DIR + "one_hot_datasets/K562_binary_{region}_{_len}.hdf5",
			region=["dhs", "proms", "re"],
			_len=LENS),
		# feature datasets
		expand(SAVE_DIR + "features_datasets/K562_binary_{region}_400_features.hdf5",
			region=["dhs", "proms", "re"])
	threads: 4
	params:
		script = join(os.environ["HOT_CODE"], "data_prep/classification/create_datasets.py")
	run:
		shell("python {params.script} K562")

