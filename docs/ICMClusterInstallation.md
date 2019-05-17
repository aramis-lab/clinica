# Installation of Clinica on the ICM cluster

## Configuration file

### Bash / Zsh

Clinica dependencies are present in `/export/data/applications/clinica/Cluster`. Meanwhile, the `/export/applications/clinica/Cluster/dot_path` file contains the different `export PATH` needed to make Clinica run.

Either you copy the content of your `~/.bashrc` / `~/.zshrc` file (don't forget to unlog then relog to the Cluster to reload your configuration file) or you can simply do a `source /export/applications/clinica/Cluster/dot_path`.

### Matlab ?

`startup.m` to modify ?


### Check that `clinica` is loaded
If it worked, the `which clinica` command should return the following path: `/export/applications/aramis/Cluster/Miniconda2/bin/clinica`.

## Submit a clinica pipeline
Please read the [Wiki section](https://wiki.icm-institute.org/display/SIKB/CLUSTER) dedicated to the ICM cluster.
