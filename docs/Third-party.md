<!-- markdownlint-disable MD046 -->
# Third-party software

## Converters

If you want to run the `convert <dataset>-to-bids` commands (e.g. `adni-to-bids`), you may have to install the [**dcm2niix**](https://github.com/rordenlab/dcm2niix), [**dcm2nii**](https://www.nitrc.org/frs/?group_id=152) and/or [**FreeSurfer**](http://surfer.nmr.mgh.harvard.edu/) tools.

|                   | dcm2nii | dcm2niix | FreeSurfer |
|:------------------|:-------:|:--------:|:----------:|
| `adni-to-bids`    |    x    |    x     |            |
| `aibl-to-bids`    |    x    |    x     |     x      |
| `nifd-to-bids`    |         |    x     |            |
| `oasis-to-bids`   |         |          |            |

Please refer to each tool’s website for installation instructions:

- [**dcm2niix**](https://github.com/rordenlab/dcm2niix)
  - Download [here](https://github.com/rordenlab/dcm2niix) and follow the installation instructions on the same page.
  - For Mac users: use Homebrew to install `dcm2niix` with `brew install dcm2niix`.
- [**dcm2nii**](https://www.nitrc.org/frs/?group_id=152) `dcm2nii` is included in an older release of MRIcron which can be downloaded [here](https://www.nitrc.org/frs/download.php/1976/mricronmac.zip).
Unpack the downloaded zip archive and move the `dcm2nii` binary somewhere within your registered `$PATH`.
- [**FreeSurfer 6.0**](http://surfer.nmr.mgh.harvard.edu/)
  - For Linux users, download and install FreeSurfer following the instructions on the [wiki](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall).
  Please note that on Ubuntu you will need to install the packages `tcsh` and `libjpeg62` (a `sudo apt-get install tcsh libjpeg62` should do the job).
  - For Mac users, download [here](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and follow the instructions on the FreeSurfer [wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/MacOsInstall).

Do not forget to check the installations following each tool’s guidelines.

## Pipeline-specific interfaces

Not all the following dependencies are necessary to install and run Clinica.
You may want to only install the software packages used by certain pipelines of Clinica.
Pipelines' specific dependencies are listed below:

|                          | ANTs | Convert3D | FreeSurfer | FSL | ITK | Matlab | MRtrix3 | PETPVC | SPM |
|:-------------------------|:----:|:---------:|:----------:|:---:|:---:|:------:|:-------:|:------:|:---:|
| `t1-volume-*`            |      |           |            |     |     |    x   |         |        |  x  |
| `t1-freesurfer`          |      |           |     x      |     |     |        |         |        |     |
| `dwi-preprocessing-*`    |   x  |     x     |            |  x  |     |        |    x    |        |     |
| `dwi-dti`                |   x  |           |            |  x  |     |        |    x    |        |     |
| `dwi-connectome`         |   x  |           |            |  x  |     |        |    x    |        |     |
| `pet-surface`            |      |           |            |  x  |  x* |        |         |   x*   |  x  |
| `pet-volume`             |      |           |            |     |  x* |    x   |         |   x*   |  x  |
| `statistics-surface`     |      |           |            |     |     |    x   |         |        |     |
| `machine-learning-*`     |      |           |            |     |     |        |         |        |     |

!!! note "CAT12 toolbox"
    Starting from Clinica `v0.3.7`, the [**CAT12**](http://dbm.neuro.uni-jena.de/cat/) toolbox is no longer needed for the `t1-volume` and `pet-volume` pipelines.
    For previous versions of Clinica, you will need to download the latest version of the toolbox [here](http://dbm.neuro.uni-jena.de/cat/index.html#DOWNLOAD) and follow the instructions to ensure that your `cat12` folder is located in your `spm/toolbox` folder.

_*You only need to install ITK if you plan to perform partial volume correction using PETPVC._

Depending on the architecture and OS of your system, setup of third party libraries can change.
Please refer to each tool’s website for installation instructions:

- [**ANTs v2.3.1**](http://stnava.github.io/ANTs/) Download [here](https://github.com/stnava/ANTs/releases) and follow the instructions on the ANTs [wiki](https://github.com/stnava/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS).
- [**Convert3D**](http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D) Download [here](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D).
- [**FreeSurfer**](http://surfer.nmr.mgh.harvard.edu/)
  - For Linux users, download and install FreeSurfer following the instructions on the [wiki](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall).
  Please note that on Ubuntu you will need to install the packages `tcsh` and `libjpeg62` ( a `sudo apt-get install tcsh libjpeg62` should do the job).
  - For Mac users, download [here](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and follow the instructions on the FreeSurfer [wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/MacOsInstall).
- [**FSL 6.0**](https://fsl.fmrib.ox.ac.uk/) Download [here](https://fsl.fmrib.ox.ac.uk/fsldownloads) and follow the instructions on the FSL wiki (this [page](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/Linux) for Linux users and this [page](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/MacOsX) for Mac users).
- [**ITK**](https://itk.org/) Follow the instructions on the ITK blog (this [page](https://blog.kitware.com/itk-packages-in-linux-distributions/) for Linux users and this [page](https://blog.kitware.com/kitware-packages-on-os-x-with-homebrew/) for Mac users).
- [**MRtrix3**](http://www.mrtrix.org) Follow the instructions on the MRtrix [website](https://www.mrtrix.org/download/)
  - For Linux users: use the official package for [Anaconda](https://www.mrtrix.org/download/linux-anaconda/) or follow this set of [instructions](https://mrtrix.readthedocs.io/en/latest/installation/build_from_source.html#linux) to build it from source.
  - For Mac users: use the official package for [Anaconda](https://www.mrtrix.org/download/macos-anaconda/) or the official [installer](https://www.mrtrix.org/download/macos-application/). There is also this [Homebrew formula](https://github.com/MRtrix3/homebrew-mrtrix3), although large dependencies such as XCode and Qt5 are required.
- [**Matlab**](https://fr.mathworks.com/products/matlab/)
- [**PETPVC 1.2.4**](https://github.com/UCL/PETPVC) Follow the instructions [here](https://github.com/UCL/PETPVC).
Do not forget to compile in RELEASE mode, otherwise, partial volume correction will be very slow.
- [**SPM12**](http://www.fil.ion.ucl.ac.uk/spm/) Download the latest version [here](http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip) and follow the instructions on the SPM wiki (this [page](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Linux) for Linux users and this [page](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)) for Mac users).
  - For systems running on MacOS Big Sur, a development version of SPM12 (check [here](https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/) for updates) as well as a more recent release of the MCR (minimum 2019a) are required.

Do not forget to check the installations following each tool’s guidelines.

## Environment variables setup

Edit the configuration file associated to your shell.
If you are using `bash`, it can be `~/.bashrc` or `~/.bash_profile`, if you are using `zsh`, it must be `~/.zshrc`.
Your file must look like as below.
Please be careful if you copy/paste this, some adjustments are needed.
You must specify on the corresponding lines the path to your different installations.

```bash
# System language setting:
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Miniconda
source /path/to/your/Miniconda/etc/profile.d/conda.sh

# ANTs
export ANTSPATH="/path/to/your/ANTs"
export PATH=${ANTSPATH}:${PATH}

# FreeSurfer
export FREESURFER_HOME="/Applications/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh &> /dev/null

# FSL
# Uncomment the line below if you are on Mac:
#export FSLDIR="/usr/local/fsl"
# Uncomment the line below if you are on Linux:
#export FSLDIR="/usr/share/fsl/6.0"
export PATH="${FSLDIR}/bin":${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh

# Matlab
export MATLAB_HOME="/path/to/your/matlab/bin/"
export PATH=${MATLAB_HOME}:${PATH}
export MATLABCMD="${MATLAB_HOME}/matlab"

# SPM
export SPM_HOME="/path/to/your/spm12"

# Dcm2nii
export PATH="/path/to/your/dcm2nii:$PATH"
```

<!-- # Autocomplete system
eval "$(register-python-argcomplete clinica)" -->

!!! Note
    `zsh` shell users will have to add this right before the last line of their configuration file to enable autocompletion in Clinica:

    ```bash
    autoload bashcompinit
    bashcompinit
    source ~/.bash_completion.d/python-argcomplete.sh
    ```

You must also add SPM to the MATLAB path variable if you installed it as a toolbox.
To do so, add the following line to your `startup.m` file located in your *initial working folder*, by default `~/Documents/MATLAB` (see [here](https://fr.mathworks.com/help/matlab/ref/startup.html) for more details).
If the file does not exist, you can create it and type inside:

```matlab
addpath('/path/to/your/spm12');
```

You can also replace the previous line by the following, assuming the `SPM_HOME` environment variable is set in your `~/.bashrc` file.

```matlab
[~, spmhome] = system('source ~/.bashrc > /dev/null; echo $SPM_HOME;');
spmhome = strsplit(spmhome,'\n');
addpath(spmhome{end-1});
```

!!! Note
    `zsh` shell users will have to replace `~/.bashrc` by `~/.zshrc`.
