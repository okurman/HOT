#!/usr/bin/env bash

c=../data/src_files/singularity_container.sif

datasets_dir=../data/classification/datasets/one_hot_datasets/
models_dir=../data/classification/models/CNN_models/
[ -d $models_dir ] || mkdir -p $models_dir

for cl in HepG2 K562
do
  for ctr in dhs proms re
  do
    for len in 400 1000
    do

      data_file=$datasets_dir/$cl"_binary_"$ctr"_"$len".hdf5"
      [ -f $data_file ] || continue
      save_dir=$models_dir/$cl"_binary_"$ctr"_"$len
      [ -f $save_dir ] || mkdir $save_dir

      singularity exec --nv $c python train_CNN_models.py $data_file $save_dir $len

    done
  done
done

