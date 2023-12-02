#!/usr/bin/env bash

logs_dir=../data/classification/models/feature_models/logs/

[ -d $logs_dir ] || mkdir -p $logs_dir

for cl in HepG2 K562
do
  for ctr in dhs proms re
  do
    for feature in -1 0 1 2 3
    do
      for kernel in poly rbf sigmoid
      do

        name=$cl"_"$ctr"_"$feature"_"$kernel

        #########################################################
        ###### To submit the training jobs to compute farm ######
        #########################################################
        qsub  -v SGE_FACILITIES -v SGE_SUMMARY="stdout" -m n \
            -l h_vmem=5G,mem_free=5G,reserve_mem=5G,m_mem_free=5G,h_rt=10:0:0 \
            -N $name \
            -o $logs_dir/$name.out \
            -e $logs_dir/$name.err \
            -b y \
            run_python_script.sh train_feature_models.py -cl $cl -ctr $ctr -model svm -feature $feature -kernel $kernel

        # train_feature_models.py -cl $cl -ctr $ctr -model svm -feature $feature -kernel $kernel

      done
    done
  done
done

