#!/bin/bash

JD=${1}
INPATH="/lustre/aoc/projects/hera/jaguirre/PolarizedMosaics/" #${2}
OUTPATH="/lustre/aoc/projects/hera/jaguirre/PolarizedMosaics/" #${3}

echo ${JD}

export JD
export INPATH
export OUTPATH

casa -c convert_to_ms_and_dirty_image.py