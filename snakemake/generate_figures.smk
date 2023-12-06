import os
from os.path import join
from pathlib import Path
from collections import defaultdict

plot2files = defaultdict(list)
fname = Path(os.environ["HOT_CODE"])/"snakemake/figure_files.txt"
for line in open(fname).readlines():
	figure = line.split("/")[0]
	_path = Path(os.environ["HOT_DATA"])/"plots"/line.strip()
	plot2files[figure].append(_path)

DATA_PATH = Path(os.environ["HOT_DATA"])

rule all:
	input:
		[str(f) for files in plot2files.values() for f in files]


rule create_figure_1:
	output:
		expand(plot2files["figure_1"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_1.py")


rule create_figure_2:
	output:
		expand(plot2files["figure_2"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_2.py")


rule create_figure_3:
	output:
		expand(plot2files["figure_3"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_3.py")


rule create_figure_4:
	output:
		expand(plot2files["figure_4"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_4.py")


rule create_figure_5:
	output:
		expand(plot2files["figure_5"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_5.py")


rule create_figure_6:
	output:
		expand(plot2files["figure_6"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_6.py")


rule create_figure_7:
	output:
		expand(plot2files["figure_7"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_7.py")


rule create_figure_8:
	output:
		expand(plot2files["figure_8"])
	script:
		join(os.environ["HOT_CODE"], "plots/figure_8.py")

