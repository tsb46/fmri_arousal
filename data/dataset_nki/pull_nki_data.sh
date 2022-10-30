#!/bin/bash

mkdir -p raw
python download_rockland_raw_bids_ver2.py -a aws_links_sample.csv -o raw
python reorganize_nki.py
