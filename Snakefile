import os

rule define_hot_loci:
	output:
		os.environ["HOT_DATA"] + "/HOTs/H1_HOTs.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/H1_HOTs.noproms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/H1_HOTs.proms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/HepG2_HOTs.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/HepG2_HOTs.noproms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/HepG2_HOTs.proms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/K562_HOTs.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/K562_HOTs.noproms.bed.gz",
		os.environ["HOT_DATA"] + "/HOTs/K562_HOTs.proms.bed.gz"

	script:
		"data_prep/extract_loci.py"


rule extract_phastcons_scores:
	output:
		os.environ["HOT_DATA"] + "/HOTs/HepG2_HOTs.bed.vertebrate.phastcons.gz",
		os.environ["HOT_DATA"] + "/HOTs/K562_HOTs.bed.vertebrate.phastcons.gz",
		os.environ["HOT_DATA"] + "/HOTs/H1_HOTs.bed.vertebrate.phastcons.gz"

	script:
		"data_prep/phastCons.py"


rule all:
	input:
		os.environ["HOT_DATA"]+"/plots/figure_1",
		os.environ["HOT_DATA"]+"/plots/figure_2",
		os.environ["HOT_DATA"]+"/plots/figure_3",
		os.environ["HOT_DATA"]+"/plots/figure_4",
		os.environ["HOT_DATA"]+"/plots/figure_5",
		os.environ["HOT_DATA"]+"/plots/figure_6",
		os.environ["HOT_DATA"]+"/plots/figure_7",
		os.environ["HOT_DATA"]+"/plots/figure_8"


rule create_figure_1:
	output:
		directory(os.environ["HOT_DATA"]+"/plots/figure_1")
	script:
		"plots/figure_1.py"


rule create_figure_2:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_2")
	script:
		"plots/figure_2.py"


rule create_figure_3:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_3")
	script:
		"plots/figure_3.py"


rule create_figure_4:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_4")
	script:
		"plots/figure_4.py"


rule create_figure_5:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_5")
	script:
		"plots/figure_5.py"


rule create_figure_6:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_6/")
	script:
		"plots/figure_6.py"


rule create_figure_7:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_7")
	script:
		"plots/figure_7.py"


rule create_figure_8:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_8")
	script:
		"plots/figure_8.py"


