#!/bin/bash
PREFIX=$1
echo 'Using prefix '$PREFIX
python merge_pols_uv.py $PREFIX
echo 'Converting to uvfits'
# Can we change default output name?
miriad_to_uvfits.py $PREFIX'.HH.uv'
