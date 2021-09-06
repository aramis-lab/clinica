#!/bin/bash
# A shell script to launch clinica in CI machines

# Name of the Conda environment according to the branch
CLINICA_ENV_BRANCH="clinica_env_$BRANCH_NAME"

set -e
set +x

ENV_EXISTS=0
# Verify that the conda enviroment correponding to the branch exists, otherwise
# create it.
ENVS=$(conda env list | awk '{print $1}' )
echo $ENVS

for ENV in $ENVS
do
  if  [[ "$ENV " == *"$CLINICA_ENV_BRANCH "* ]]
  then
    echo "Find Conda environment named $ENV, continue."
    ENV_EXISTS=1
    break
  fi;
done
if [ "$ENV_EXISTS" = 0 ]; then
  echo "Conda env $CLINICA_ENV_BRANCH not found... Creating"
  conda env create --force --file environment.yml -n $CLINICA_ENV_BRANCH
  echo "Conda env $CLINICA_ENV_BRANCH was created."
  eval "$(conda shell.bash hook)"
  conda activate $CLINICA_ENV_BRANCH
  pip install -r requirements-dev.txt
  echo "Developement requirements has been installed in  $CLINICA_ENV_BRANCH."
  conda deactivate
fi
