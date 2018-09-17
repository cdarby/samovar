#!/bin/sh

INSTALLDIR="."

case "$1" in
  "generateVarfile")
    python pipeline/0.setup/generateVarfile.py "${@:2}"
    ;;
   "simulate")
    python pipeline/1.simulate/simulate.py "${@:2}"
    ;;
  "train")
    python pipeline/2.train/train.py "${@:2}"
    ;;
  "preFilter")
    python pipeline/3.classifyAndFilter/filter.py "${@:2}"
    ;;
  "classify")
    python pipeline/3.classifyAndFilter/classify.py "${@:2}"
    ;;
  "postFilter")
    python pipeline/3.classifyAndFilter/linkageFilter.py "${@:2}"
    ;;
  *)
    echo "Options: generateVarfile, simulate, train, preFilter, classify, postFilter (run without options to see arguments)"
    exit 1
    ;;
esac
