Codebase of the manuscript:

## "Sequence characteristics and an accurate model of high-occupancy target loci in the human genome 
Sanjarbek Hudaiberdiev, Ivan Ovcharenko 

NCBI, NLM, NIH
Bethesda, MD 20894

https://www.biorxiv.org/content/10.1101/2023.02.05.527203v1
*********

### Getting started

Processed data and source files are available in the Zenodo repository (https://zenodo.org/records/10251023).  
Download the file `supplemental_files.tar.gz` and extract its content to the `data/src_files` directory

Create a python virtual environment and set up the running environment: \
```
conda env create -f environment.yml
conda activate hot_env
export HOT_DATA=$(readlink --canonicalize data)
export HOT_CODE=$(pwd)
```

### Define logarithmically binned DAP-bound loci and HOT loci 

`snakemake -c -R define_hot_loci`

This will populate with files the following the directories:
- `data/loci/`
- `data/log_bins/`
- `data/HOTs/`

`snakemake -c4 -R extract_phastcons_scores`
This will download phastCons scores and extract the conservation scores of the compared loci.

The HOT loci in three cell lines will be extracted to `data/HOTs/` with the names `HepG2_HOTs.bed.gz`, `K562_HOTs.bed.gz` and `H1_HOTs.bed.gz`.

### Generate figures in batch

`snakemake -c8 -R all`
This will create a directory for each of 8 figures under `data/HOTs/plots/` and generate the subplots depicted on the main text.
For questions about the figures in the supplemental figures please create an issue on this repo or reach out to the authors.


### Generate figures individually

***
#### Figure 1: Definition of HOT loci and compositional analyses

<img src="./src_figures/Figure1.png" width="480">

```
cd plots
python figure_1.py
```


***
#### Figure 2: PCA analysis of HOT loci.

<img src="./src_figures/Figure2.png" width="400">

```
cd plots
python figure_2.py
```


***
#### Figure 3: Hi-C analysis of the HOT loci.

<img src="./src_figures/Figure3.png" width="450">

```
cd plots
python figure_3.py
```


***
#### Figure 4: ChIP-seq signal strength analysis.

<img src="./src_figures/Figure4.png" width="450">

```
cd plots
python figure_4.py
```


***
#### Figure 5: Evolutionary conservation scores (phastCons, phyloP) and classification analysis.

<img src="./src_figures/Figure5.png" width="450">

```
cd data_pred
python phastCons.py
```

Note that this will first create a folder `phastCons` and download the phastCons files from UCSC Genome Browser Database and proceed to extract respective phastCons scores for each bin in `log_bins` directory. This will run approximately for 30 minutes with enough RAM available.

```
cd plots
python figure_5.py
```

By default, the data will be downloaded and processed for vertebrates. To re-create the conservation score analyses for placental mammals and primates, re-run the scripts `data_prep/phastCons.py`, `plots/figure_5.py` with parameter values `placentalMammals` and `primates` 

For Figure5-D, please refer to `Classification analysis of HOT loci`


***
#### Figure 6: Functional analysis (GO enrichment, housekeeping, tissue-specificity).

<img src="./src_figures/Figure6.png" width="450">

```
cd plots
python figure_6.py`
```

For Figure6-C refer to the main text and supplemental methods.

***
#### Figure 7: 

<img src="./src_figures/Figure7.png" width="450">

```
cd plots
python figure_7.py
```

***
#### Figure 8: 

<img src="./src_figures/Figure8.png" width="450">

```
cd plots
python figure_8.py
```