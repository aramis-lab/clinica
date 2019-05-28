#!/bin/bash
# A shell script to launch clinica in CI machines

# Name of the Conda environment according to the branch
CLINICA_ENV_BRANCH="clinica_env_$BRANCH_NAME"

set -e
set +x


# Verify that the conda enviroment correponding to the branch exists, otherwise
# create it.
ENVS=$(conda env list | awk '{print $1}' )
echo $ENVS

for ENV in $ENVS
do
  if  ! [[ "$ENV " == *"$CLINICA_ENV_BRANCH "* ]]
  then
    #echo "Conda env named $CLINICA_ENV_BRANCH not found, try next"
    continue
  else
    echo "Find Conda environment named $ENV, continue."
    exit 0
  fi;
done
echo "Conda env $CLINICA_ENV_BRANCH not found... Creating"
conda env create --force --file environment.yml -n $CLINICA_ENV_BRANCH
echo "Conda env $CLINICA_ENV_BRANCH was created."
