import os
from os.path import join

FIGURE_DIRS = [os.environ["HOT_DATA"]+f"/plots/figure_{i}" for i in range(1, 9)]

rule all:
	input:
		expand(FIGURE_DIRS)

rule create_figure_1:
	output:
		directory(os.environ["HOT_DATA"]+"/plots/figure_1")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_1.py")


rule create_figure_2:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_2")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_2.py")


rule create_figure_3:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_3")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_3.py")


rule create_figure_4:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_4")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_4.py")


rule create_figure_5:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_5")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_5.py")


rule create_figure_6:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_6/")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_6.py")


rule create_figure_7:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_7")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_7.py")


rule create_figure_8:
	output:
		directory(os.environ["HOT_DATA"] + "/plots/figure_8")
	script:
		join(os.environ["HOT_CODE"], "plots/figure_8.py")

