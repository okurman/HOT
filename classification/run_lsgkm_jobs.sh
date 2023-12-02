#!/usr/bin/env bash

logs_dir=../data/classification/models/lsgkm_models/logs/
[ -d $logs_dir ] || mkdir -p $logs_dir

for cl in HepG2 K562
do
  for ctr in dhs proms re
  do
    for len in 400 1000
    do

      for kernel in {0..5}
      do

        #########################################################
        ###### To submit the training jobs to compute farm ######
        #########################################################
        #name=$cl"_"$ctr"_"$ctr"_"$len"_"$kernel
        #qsub  -v SGE_FACILITIES -v SGE_SUMMARY="stdout" -m n \
        #    -l h_vmem=5G,mem_free=5G,reserve_mem=5G,m_mem_free=5G,h_rt=10:0:0 \
        #    -N $name \
        #    -o $logs_dir/$name.out \
        #    -e $logs_dir/$name.err \
        #    -b y \
        #    run_python_script.sh train_lsgkm_models.py -cl $cl -ctr $ctr -len $len -kernel $kernel

        train_lsgkm_models.py -cl $cl -ctr $ctr -len $len -kernel $kernel

      done
    done
  done
done

