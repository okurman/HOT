import os

CELL_LINES = ["HepG2", "K562", "H1"]

# noinspection SmkAvoidTabWhitespace

rule all:
	input:
		os.environ["HOT_DATA"] + "/phastCons/vertebrate",
		expand(os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.bed.gz", cell_line=CELL_LINES),
		expand(os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.proms.bed.gz", cell_line=CELL_LINES),
		expand(os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.noproms.bed.gz", cell_line=CELL_LINES)


rule define_hot_loci:
	params:
		loci=os.environ["HOT_CODE"] + "/data_prep/extract_loci.py",
		phastcons=os.environ["HOT_CODE"] + "/data_prep/phastCons.py"
	output:
		os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.proms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/{cell_line}_HOTs.noproms.bed.gz"
	run:
		shell("python {params.loci} {wildcards.cell_line}")
		shell("python {params.phastcons} {wildcards.cell_line} vertebrate")


rule download_phastcons:
	output:
		directory(os.environ["HOT_DATA"] + "/phastCons/vertebrate")
	params:
		script=os.environ["HOT_CODE"] + "/data_prep/phastCons_download.py"
	shell:
		"python {params.script} vertebrate"

