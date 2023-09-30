# Installation of Clinica on the ICM cluster

## Load Clinica

Clinica is available through Module system.
After an `ssh` on the Cluster (`ssh login01` or `ssh login02`), you will land in your home directory.
Once done, type `module avail` to see all packages.
Finally, type:

```shell
module load clinica/aramis
```

## Configuration file

### Bash / Zsh

Clinica dependencies are present in `/export/data/applications/clinica/Cluster`.
Meanwhile, the `/export/applications/clinica/Cluster/dot_path` file contains the different `export PATH` needed to make Clinica run.

Either you copy the content of your `~/.bashrc` / `~/.zshrc` file (don't forget to unlog then relog to the Cluster to reload your configuration file) or you can simply do a `source /export/applications/clinica/Cluster/dot_path`.

### Matlab ?

`matlab/startup.m` to modify ?

## Submit a clinica pipeline

Please read the [Wiki section](https://wiki.icm-institute.org/display/SIKB/CLUSTER) dedicated to the ICM cluster.
