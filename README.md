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

`snakemake --cores 4 -s snakemake/data_prep.smk`

This will populate with files the following the directories:
- `data/loci/`
- `data/log_bins/`
- `data/HOTs/`

`snakemake -c4 -R extract_phastcons_scores`
This will download phastCons scores and extract the conservation scores of the compared loci.

The HOT loci in three cell lines will be extracted to `data/HOTs/` with the names `HepG2_HOTs.bed.gz`, `K562_HOTs.bed.gz` and `H1_HOTs.bed.gz`.

### Generate figures in batch

`snakemake --cores 5 -s snakemake/generate_figures.smk`

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

### Classification analysis

First, generate the control regions used in the study for classification analysis with:

`snakemake --cores 2 -s snakemake/data_prep_control_regions.smk`

This will create and populate the directory: 
 - `data/classification/datasets/control_regions`

Then, create the datasets with:
`snakemake --cores 2 -s snakemake/data_prep_classification.smk`

This will create and populate the directories: 
 - `data/classification/datasets/fasta_datasets`
 - `data/classification/datasets/features_datasets`
 - `data/classification/datasets/one_hot_datasets`

#### CNN classification

CNNs were trained using old versions of Tensorflow/Keras so we recommend using the singularity container provided in the Zenodo repository.
To train models, run:

```
cl=HepG2
ctr=dhs
seq_len=400
name=$cl"_binary_"$ctr"_"$seq_len
data_file=/classification/datasets/one_hot_datasets/$name.hdf5
save_dir=classification/models/CNN_models/$name

singularity exec --nv data/src_files/src_files/singularity_container.sif python \ 
classification/train_CNN_models.py $data_file $save_dir $seq_len
````

Training this model will require GPUs on the computer.
The file structure under `$save_dir` will contain files: `AUCs.txt  model.hdf5  training_history_log.tab  weights.hdf5`

Repeat the training procedure for all the values as follows:
```bash
for cl in HepG2 K562:
do
    for ctr in dhs re proms
    do
        for seq_len in 400 1000
        do
            name=$cl"_binary_"$ctr"_"$seq_len
            data_file=/classification/datasets/one_hot_datasets/$name.hdf5
            save_dir=classification/models/CNN_models/$name
            
            singularity exec --nv data/src_files/src_files/singularity_container.sif python \ 
            classification/train_CNN_models.py $data_file $save_dir $seq_len
        done
    done
done
```
