#! /bin/sh

#--------------------------------------#
# Clinica package creations ( wheels)  
#--------------------------------------#
#
# WARNING: Activate a conda environment with the right pip version.
# Use at your own risk.


CURRENT_DIR=$(pwd)
echo "$CURRENT_DIR"

# ensure we are in the right dir
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR"
echo "Entering ${SCRIPT_DIR}/../../"
cd "${SCRIPT_DIR}/../../"
ls 

# clean pycache stuff
rm -rf dist build clinica.egg-info/
find . -name "*__pycache__*" -exec rm {} -rf \;
find . -name "*.pyc*" -exec rm {} -rf \;

set -o errexit
set -e
# generate wheel
LC_ALL=C.UTF-8 poetry build
# come back to directory of
cd "$CURRENT_DIR"
