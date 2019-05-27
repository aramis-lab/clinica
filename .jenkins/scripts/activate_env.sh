#!/bin/bash
# A shell script to launch clinica in CI machines

# Name of the Conda environment according to the branch
CLINICA_ENV_BRANCH="clinica_env_$BRANCH_NAME"

set -e

# Verify that the conda enviroment correponding to the branch exists, otherwise
# create it.
ENVS=$(conda env list | awk '{print $1}' )

if  [[ $ENVS = *"$CLINICA_ENV_BRANCH"* ]]
then
  echo "Conda environment exits, continue..."
else
  echo "Create Conda environment $CLINICA_ENV_BRANCH because it does not exist."
  conda env create --force --file environment.yml -n $CLINICA_ENV_BRANCH
fi;

# Activate conda environment
echo "Activate conda environment $CLINICA_ENV_BRANCH..."
source activate $CLINICA_ENV_BRANCH
