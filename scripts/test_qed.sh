#!/bin/bash
#--enum_root --hidden_size 270 --embed_size 200
DIR=$1
modi=$2
f=$DIR/model.$modi
if [ -e $f ]; then
    echo $f
    python decode.py --test data/qed/covid_test.txt --vocab data/qed/covid_vocab.txt --model $f --num_decode 20  | python scripts/qed_score.py > covid_results_v3/results.$modi &
fi

