#!/bin/bash
set -e

if [ "$1" = 'bash' ]; then
    exec bash
fi

source ~/.bashrc
exec bash -c "clinica $@"
