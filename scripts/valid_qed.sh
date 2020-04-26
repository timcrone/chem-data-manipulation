#!/bin/bash
#--enum_root --hidden_size 270 --embed_size 200
DIR=$1
ST=$2
ED=$3

for ((i=ST; i<=ED; i++)); do
    f=$DIR/model.$i
    if [ -e $f ]; then
        echo $f
        python decode.py --test data/qed/covid_valid.txt --vocab data/qed/covid_vocab.txt --model $f --num_decode 20  | python scripts/qed_score.py > covid_results_v3/results.$i &
    fi
done
