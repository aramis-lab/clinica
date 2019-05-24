#!/bin/bash
# A shell script to launnch clinica in CI machines

# Name of the Conda environment according to the branch
CLINICA_ENV_BRANCH="clinica_env_$BRANCH_NAME"

set -e

# Verify that the conda enviroment correponding to the branch exists, otherwise
# create it.
ENVS=$(conda env list | awk '{print $1}' )

if  [[ $ENVS = *"$CLINICA_ENV_BRANCH"* ]]
then
  echo "Create Conda environment..."
  conda env create --force --file environment.yml -n $CLINICA_ENV_BRANCH
else
  echo "Continue..."
fi;

# Activate conda environment
source activate $CLINICA_ENV_BRANCH
# Install conda source using pip
pip install --ignore-installed .
eval "$(register-python-argcomplete clinica)"
# Show clinica help message
clinica --help
# Desactivate conda environment
source deactivate
