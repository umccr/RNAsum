#!/usr/bin/env bash

### Clone the repo
git clone https://github.com/umccr/RNAseq-Analysis-Report

### Install conda
unset PYTHONPATH
unset CONDA_PREFIX
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda && rm miniconda.sh
export PATH=$(pwd)/miniconda/bin:$PATH
conda update conda

### Install environments
ENV_NAME=rnaseq-report
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME} --file RNAseq-Analysis-Report/envs/environment.yml
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$PATH
pip install -e rnaseq-report

### Create the loader script
ENV_NAME=rnaseq-report
cat <<EOT > load_rnaseq-report.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$(pwd)/miniconda/bin:\$PATH
export CONDA_PREFIX=$(pwd)/miniconda/envs/${ENV_NAME}
EOT