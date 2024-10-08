Codebase of the manuscript:

## "Functional characteristics and an accurate model of high-occupancy target loci in the human genome 
Sanjarbek Hudaiberdiev, Ivan Ovcharenko 

NCBI, NLM, NIH
Bethesda, MD 20894


https://www.biorxiv.org/content/10.1101/2023.02.05.527203v4
*********


# Table of contents
1. [Getting started](#gts)
2. [Define loci](#dfnlc)
3. [Extract phastCons scores](#phst)
4. [Generate figures in batch](#gnrtfgrsbtch)
5. [Generate individual figures](#gnrtfgrsind)
6. [Classification analysis](#cls)
   1. [Generate control regions](#cls-ctr)
   2. [Generate datasets](#cls-dts)
   3. [Train Convolutional Neural Networks](#cls-cnn)
   4. [Train LS-GKM models](#cls-lsgkm)
   5. [Train sequence feature (GC, CG, CpG, GCI) models](#cls-features)

### Getting started <a name="gts"></a>
___

Clone this codebase to a directory with at least 50 GB space and download the source files from Zenodo:

```bash

git clone https://github.com/okurman/HOT.git
mv HOT code
mkdir data

export HOT_DATA=$(readlink --canonicalize data)
export HOT_CODE=$(readlink --canonicalize code)

wget https://zenodo.org/records/13271790/files/supplemental_files.tar.gz
tar -xzf supplemental_files.tar.gz
mv src_files data/

cd code
conda env create -f environment.yml
conda activate hot_env
```

### Define logarithmically binned DAP-bound loci and HOT loci <a name="dfnlc"></a>
___

`snakemake --cores 4 -s snakemake/data_prep.smk`

Or individually by running the python scripts:
```bash
for cl in HepG2 K562 H1
do
  python data_prep/extract_loci.py $cl
done
```

This will populate with files the following the directories:
- `data/loci/`
- `data/log_bins/`
- `data/HOTs/`


### Generate figures in batch <a name="gnrtfgrsbtch"></a>
___

`snakemake --cores 5 -s snakemake/generate_figures.smk`

This will create a directory for each of 8 figures under `data/plots/` and generate the subplots depicted on the main text.
For questions about the figures in the supplemental figures please create an issue on this repo or reach out to the authors.

### Generate figures individually <a name="gnrtfgrsind"></a>
___

#### Figure 1: Definition of HOT loci and compositional analyses

<img src="./src_figures/Figure1.png" width="480">

```bash
python plots/figure_1.py
```

#### Figure 2: PCA analysis of HOT loci.

<img src="./src_figures/Figure2.png" width="400">

```bash
python plots/figure_2.py
```

#### Figure 3: Hi-C analysis of the HOT loci.

<img src="./src_figures/Figure3.png" width="450">

```bash
python plots/figure_3.py
```

#### Figure 4: ChIP-seq signal strength analysis.

<img src="./src_figures/Figure4.png" width="450">

```bash
python plots/figure_4.py
```

#### Figure 5: Evolutionary conservation scores (phastCons, phyloP) and classification analysis.

<img src="./src_figures/Figure5.png" width="450">

```bash
python data_prep/phastCons.py
```

Note that this will first create a folder `phastCons` and download the phastCons files from UCSC Genome Browser Database and proceed to extract respective phastCons scores for each bin in `log_bins` directory. This will run approximately for 30 minutes with enough RAM available.

```bash
python plots/figure_5.py
```

By default, the data will be downloaded and processed for vertebrates. To re-create the conservation score analyses for placental mammals and primates, re-run the scripts `data_prep/phastCons.py`, `plots/figure_5.py` with parameter values `placentalMammals` and `primates` 

For Figure5-D, please refer to `Classification analysis of HOT loci`


#### Figure 6: Functional analysis (GO enrichment, housekeeping, tissue-specificity).

<img src="./src_figures/Figure6.png" width="450">

```bash
python plots/figure_6.py
```

For Figure6-C refer to the main text and supplemental methods.

#### Figure 7: 

<img src="./src_figures/Figure7.png" width="450">

```bash
python plots/figure_7.py
```

#### Figure 8: 

<img src="./src_figures/Figure8.png" width="450">

```bash
python plots/figure_8.py
```

### Classification analysis <a name="cls"></a>
___
All the pre-trained models used in the analysis are available in the archived file in Zenodo repository:
-`models_trained.tar.gz`

The classification evaluation results (auROC and auPRC values) of each of the classification experiments are available in the file:
`data/src_files/classification_results.txt`

Below are the instructions for re-generating the whole classification analysis.

#### Generate control regions <a name="cls-ctr"></a>

First, generate the control regions used in the study for classification analysis with:
`snakemake --cores 2 -s snakemake/data_prep_control_regions.smk`

This will create and populate the directory: 
 - `data/classification/datasets/control_regions`

#### Prepare the datasets <a name="cls-dts"></a>

Then, create the datasets with:
`snakemake --cores 2 -s snakemake/data_prep_classification.smk`

This will create and populate the directories: 
- data/classification/datasets/fasta_datasets
- data/classification/datasets/features_datasets
- data/classification/datasets/one_hot_datasets

#### Train Convolutional Neural Networks <a name="cls-cnn"></a>

CNNs were trained using old versions of Tensorflow/Keras so we recommend using the singularity container provided in the Zenodo repository.
To train models, run:

```bash
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
            data_file=classification/datasets/one_hot_datasets/$name.hdf5
            save_dir=classification/models/CNN_models/$name
            
            singularity exec --nv data/src_files/src_files/singularity_container.sif python \ 
            classification/train_CNN_models.py $data_file $save_dir $seq_len
        done
    done
done
```

#### SVMs trained on DNA sequences <a name="cls-lsgkm"></a>

This will require the presence of the LS-GKM tool in the system path. Download and install it from:
https://github.com/Dongwon-Lee/lsgkm

Then run the following command to iteratively train all combinations of datasets:
```bash
for cl in HepG2 K562:
do
    for ctr in dhs re proms
    do
        for seq_len in 400 1000
        do
            python classification/train_lsgkm_models.py -cl $cl -ctr $ctr -len $seq_len
        done
    done
done
```

The results will be saved under the path:
`classification/models/lsgkm_models`

To train the models with different kernel versions supply `-kernel` parameter with the values from the range `1 - 7`. Default version will run with `kernel=4`, i.e. "gapped-kmers"

#### Classification using the sequence features <a name="cls-features"></a>

This part uses sequence features of GC, CG, CpG and CGI island coverage, instead of the direct sequences.

Run the following command to iteratively train all combinations of datasets:
```bash
for cl in HepG2 K562:
do
    for ctr in dhs re proms
    do    
      python classification/train_feature_models.py -cl $cl -ctr $ctr -model logreg -feature -1  
    done
done
```

`-feature` parameter will indicate on which feature to train the models. Value `-1` will ensure to train on all of the features at once. 

`-model` parameter can have self-explanatory `logreg` and `svm` values. In case of training the SVM models, the `-kernel` parameter can also be supplied to use different kernels with values `linear, rbf (radial basis function), poly (polynomial), sigmoid` 
