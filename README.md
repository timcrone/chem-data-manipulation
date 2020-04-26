# cse6250-team5-data
Code for project team 5

To see selections of data and trained models from repository, please clone https://github.gatech.edu/tcrone3/cse6250-team5-data.git and pull the master branch.

## Introduction
This repository contains code for preprocessing and handling the 827 million compounds in the raw ZINC data.  To generate the test, training, and chemical data for our project we used the tools in this folder in addition to standard Unix/Linux command-line stream processing tools.

## Hardware requirements
SSDs are strongly recommended.  The system will need at least 32GB of RAM for several of these steps, a CUDA-enabled GPU for others.  Most of the chemical calculations are CPU bound due to limitations in RDKit; however all cores will be used where possible.  It is strongly recommended that an additional SSD drive is configured as a swap disk for situations where memory is exceeded.

## Software requirements
All steps beyond the initial ZINC data download assume a Linux BASH-based environment.

Python modules variously required:
* PySpark
* RDKit
* PyTorch

Because RDKit requires a conda environment, Anaconda or Miniconda are a hard requirement.  Additionally, pyspark requires that 'python' and 'python3' be the same version.

In addition to this repository, you will need to load the Docker image of Schrodinger, Inc.'s gpusimilarity tool:
https://github.com/schrodinger/gpusimilarity/tree/master/python

Because this Docker image requires CUDA integration it will need to be loaded on a Linux system with a version of Docker above 19.03, as well as the Nvidia Docker tools.  Additionally the GPU will need more than 2GB of onboard RAM.

Model training was performed using NVIDIA Ti1070, 32GB of system RAM, and took about 2-3 hours per iteration.

## Installation
Use anaconda to install.
* `conda env create -f environment.yml`

## Steps
__Download the ZINC data from http://files.docking.org/3D/__
  * Use the wget command from tools/zinc-download to save the ZINC data. This should work in Windows and Linux.
  * The final download size is close to 250 GB.
  * Because of the way docking.org is set up you will only be able to download at a very low bitrate. It took about three weeks to download all the data.
  * If the download is interrupted you will need to delete \*.html in the tree to refresh the wget DB.
  * If you can preemptively narrow your target data subset you should instead use the Tranche browser http://zinc15.docking.org/tranches/home/ as it is reputed to be much faster.
  * Any output data in the repository is from our snapshot loaded during February / March 2020.

__Collate the tree data into several larger files.__
  * `cd <source>/files.docking.org/3D/`
  * `ls -1d ?? | xargs -n 1 -I {} <git-repo>/tools/concat-files.sh {}`
  * At this point the raw source data will no longer be used; it is recommended to compress and store them away for safekeeping.
  * Verify that there are 121 files, labeled in alphabetical sequence AA.tsv through KK.tsv.

__Use ChemStats to load and organize the molecule set, and to generate QED scores for all molecules.__
  * `cd <repo>/etl`
  * `find <source>/files.docking.org/3D/ -name \*.tsv -exec python3 chemstats.py {} \<out-path> \;`
  * This will output CSV files containing QED data in \<out-path>; these files are tab separated even though they are labeled CSV for convenience.
  * Further manipulation for the report was done by organizing the molecules by QED score (see etl/chemfilter.py).  However this is not necessary for the overall work flow.

__Identify the base compounds to be used for selection.__
  * We used prospective treatments for COVID-19 as the foundation for our search.
  * We were also able to do elementary statistical analysis of the overall data available to us.
  * These compounds need to be stored in a file in the format: `name SMILES`
  * Our compounds are in tools/mols-to-check

__Collate all compounds into manageable .smi-format files for use in gpusimilarity.__
  * `cd <out-path>`
  * `find . -name ??\*.csv -exec cut -f 1,2 \; >> mols.smi`
  * `split mols.smi mols -l 20000000  # 20 million rows per`
  
__Gzip all generated .smi files.__
 * `find . -name \*.smi -exec gzip {} \;`

__Start the gpusimilarity Docker image.__
  * On the host, copy \<our git repo>/tools/* to \<out-path>
  * `docker run --gpus all --net=host -v <out-path>:/mnt/fsim -it klorton/gpusimilarity:latest /bin/bash`

__Use the gpusimilarity Docker image to calculate fingerprints for the known molecules.__
  * It is strongly recommended to run the below in a screen session, as it will take a long while to run.
  * `wget https://bootstrap.pypa.io/get-pip.py ; python3 get-pip.py ; python3 -m pip install ipyparallel`
  * `cd /gpusimilarity/bld/python`
  * `ipcontroller --ip=* &`
  * `ipcluster start -n 8 &` then wait for "Engines appear to have started successfully"
  * `find /mnt/fsim/ -name \*.smi.gz -exec python3 gpusim_createdb.py {} {}.fsim \;`
  * Because this uses rdkit it is completely CPU bound.

__Use the gpusimilarity Docker image to calculate similarity values for the desired molecules.__
  * `cd /gpusimilarity/bld/python`
  * For a variety of reasons the base Docker image does not work for our purposes.  Copy our modified gpusim_search.py to continue this exercise.
  * `cp /mnt/fsim/gpusim_search.py .`
  * `cd /gpusimilarity/bld`
  * `cp /mnt/fsim/run-*.sh .`
  * `sh ./run-all-sims.sh`
  * You can now close the Docker image.

__The .sim files in \<out-path> are now named based on the desired chemicals.  If desired, these can be organized with `sort -k` and the like.__

__Collect an overall file by collating all the .sim files into an organized .list file.__
  * `find \<out-path> -name \*.sim -exec cat {} \; >> covid.list`

__Run etl/chemfinal.py to generate the training and testing data for the graph2graph models.__
  * `python3 chemfinal.py > modeldata.txt`
  * Since each line and pair is randomly selected, simply use `tail -n X modeldata.txt` and `head -n Y modeldata.txt` to collect training and testing data from the output.

The output of these commands, as well as some intermediate output, is stored in data/*.gz

__Run etl/make_sorted_qed_pairs.py to sort training data and generate train/valid/test/vocab partitions.__

  * `python etl/make_sorted_qed_pairs.py < data/pairs`
  * Out: `data/covid_train_pairs.txt`
  * Out: `data/covid_valid.txt`
  * Out: `data/covid_test.txt`
  * Out: `data/covid_mols.txt`

__Make custom vocabulary file.__
  * `python models/hg2g/get_vocab.py < data/covid_mols.txt > data/covid_vocab.txt`

__Preprocess raw smiles pairs into pytorch tensors.__
  * Use G2G or HG2G scripts for respective architectures
  * `python models/g2g/preprocess.py --train data/qed/covid_train_pairs.txt --vocab data/qed/covid_vocab.txt --ncpu 4 < data/qed/covid_train_pairs.txt`
  * `python models/hg2g/preprocess.py --train data/qed/covid_train_pairs.txt --vocab data/qed/covid_vocab.txt --ncpu 4 < data/qed/covid_train_pairs.txt`

__Train models__.
  * HG2G: `python models/hg2g/gnn_train.py --train train_processed/ --vocab data/covid_vocab.txt --save_dir model/`
  * G2G: `python models/g2g/diff_vae/vae_train.py --train processed/ --vocab ../data/covid_vocab.txt --save_dir model/  | tee model/train_logs`
  * HG2G Preprocess and train:
  *  `python get_vocab.py && python preprocess.py --train data/covid_train_pairs.txt --vocab data/covid_vocab.txt --ncpu 4 < data/covid_train_pairs.txt && mv tensors*.pkl covid_train_processed && python gnn_train.py --train covid_train_processed/ --vocab data/covid_vocab.txt --save_dir model/ --batch_size 32`

__Generate model predictions.__
  * Use G2G or HG2G scripts for respective architectures
  * Make predictions `python models/hg2g/decode.py --test ../data/test_covid.txt --vocab ../data/covid_vocab.txt --model model/model.iter-4 --use_molatt | python ../scripts/qed_score.py > results.covid.qed`
  * Quantify QED results `python ../scripts/qed_analyze.py < results.covid.qed`
  * Predict + Analyze `python decode.py --test ../data/qed/train_pairs_sub.txt --vocab ../data/qed/vocab.txt --model model/model.iter-4 --use_molatt | python ../scripts/qed_score.py > results.train_sub.qed`
  * Validate 12 models saved in `models_hg2g`: `scripts/valid_qed.sh models_hg2g/ 0 12`

__Make validation figures.__
  * input: `model/train_logs/hg2g_v3.log`
  * input: `model/validation/covid_results_v3/`
  * `python make_val_figure.py`
  * Output figure in `figures`

__Make comparison test figures.__
* input target file: `data/covid.csv`
* input predictions output 1: `model/test/covid_model_10_target_results.txt`
* input predictions output 2: `model/test/newmodel_1_target_results.txt`
* `python scripts/make_test_figure.py`
* Output figures and table in `figures`

