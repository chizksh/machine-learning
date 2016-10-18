#!/bin/bash
INPUT=$1
OUTPUT=$2

ipython casoffinder_gen_inputfile.py $INPUT $OUTPUT
sed -i 1i"/data/GENOME_SEQ/chromosome/human_hg19\nNNNNNNNNNNNNNNNNNNNNNRG 2 1" $OUTPUT
