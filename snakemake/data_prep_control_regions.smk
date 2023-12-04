import os
from os.path import join

CELL_LINES = ["HepG2", "K562"]
LENS = [1000, 400]

# noinspection SmkAvoidTabWhitespace

SAVE_DIR = os.environ["HOT_DATA"] + "/classification/datasets/"
FASTA_FILE = join(os.environ["HOT_DATA"], "src_files/hg19_files/hg19.fa")
FASTA_DATASET_FILE_NAMES = ["controls_test.fa", "controls_train.fa", "hots_test.fa", "hots_train.fa"]


rule all:
	input:
		FASTA_FILE,

		# control regions
		expand(SAVE_DIR + "control_regions/{cell_line}_{region}_{_len}.bed.gz",
				cell_line=CELL_LINES,
				region=["hots", "proms_controls", "re_controls", "dhs_controls"],
				_len=LENS),


rule download_fasta:
	output:
		join(os.environ["HOT_DATA"], "src_files/hg19_files/hg19.fa.gz")
	script:
		join(os.environ["HOT_CODE"], "data_prep/classification/create_datasets.py")


rule generate_control_regions_hepg2:
	input:
		FASTA_FILE
	output:
		expand(SAVE_DIR + "control_regions/HepG2_{region}_{_len}.bed.gz",
			region=["hots", "proms_controls", "re_controls", "dhs_controls"],
			_len=LENS)
	params:
		script=join(os.environ["HOT_CODE"], "data_prep/classification/generate_control_regions.py")
	run:
		shell("python {params.script} HepG2")


rule generate_control_regions_k562:
	input:
		FASTA_FILE
	output:
		expand(SAVE_DIR + "control_regions/K562_{region}_{_len}.bed.gz",
			region=["hots", "proms_controls", "re_controls", "dhs_controls"],
			_len=LENS)
	params:
		script=join(os.environ["HOT_CODE"], "data_prep/classification/generate_control_regions.py")
	run:
		shell("python {params.script} K562")

