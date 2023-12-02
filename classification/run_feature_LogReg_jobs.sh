#!/usr/bin/env bash

for cl in HepG2 K562
do
  for ctr in dhs proms re
  do
    for feature in -1 0 1 2 3
    do
      python train_feature_models.py -cl $cl -ctr $ctr -model logreg -feature $feature
    done
  done
done


