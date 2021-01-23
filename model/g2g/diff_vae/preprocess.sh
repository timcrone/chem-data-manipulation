#!/bin/bash

python ../scripts/preprocess.py --train ../data/logp06/train_pairs.txt
python ../scripts/preprocess.py --train ../data/logp04/train_pairs.txt
python ../scripts/preprocess.py --train ../data/qed/train_pairs.txt
python ../scripts/preprocess.py --train ../data/drd2/train_pairs.txt
