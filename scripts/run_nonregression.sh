#!/bin/bash

OS=$(uname -s)

make env.conda

if [ "$OS" == "Darwin" ]
then
    source ~/miniconda3/etc/profile.d/conda.sh
elif [ "$OS" == "Linux" ]
then
    source /builds/miniconda/etc/profile.d/conda.sh
else
    echo Unsupported OS $OS
    exit
fi

conda activate $WORKSPACE/env

if [ "$OS" == "Darwin" ]
then
    source "$(brew --prefix)/opt/modules/init/bash"
else
    source /usr/local/Modules/init/profile.sh
fi

module load clinica.all

make install

cd test

if [ "$OS" == "Darwin" ]
then
    WORKING_DIRECTORY=/Volumes/data/working_dir_mac
    INPUT_DATA_DIRECTORY=/Volumes/data_ci
    TMP_DIRECTORY=/Volumes/data/tmp
else
    WORKING_DIRECTORY=/mnt/data/ci/working_dir_linux
    INPUT_DATA_DIRECTORY=/mnt/data_ci
    TMP_DIRECTORY=/mnt/data/ci/tmp
fi

poetry run pytest --verbose \
    --working_directory=$WORKING_DIRECTORY \
    --input_data_directory=$INPUT_DATA_DIRECTORY \
    --basetemp=$TMP_DIRECTORY \
    --junitxml=./test-reports/non_regression_$MODALITY_mac.xml \
    --disable-warnings \
    ./nonregression/pipelines/$MODALITY

