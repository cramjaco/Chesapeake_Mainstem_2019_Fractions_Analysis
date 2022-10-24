#!/bin/bash

# I'm just keeping notes here. This isn't designed to run
# place_seqs.py -s ./ASVs_NS20.fa -o out.tre -p 1 \
#               --intermediate intermediate/place_seqs

# following step
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1