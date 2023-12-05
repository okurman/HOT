#!/bin/bash

# provide a path to python when submitting to a cluster.
#p=[..]/python

p=`which python`

$p "$@"