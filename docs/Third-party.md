<!-- markdownlint-disable MD046 -->
# Third-party software

## Converters

Some converters require a recent version of [dcm2niix](https://github.com/rordenlab/dcm2niix) to transform DICOM files into NIfTI:

|                   | dcm2niix |
|:------------------|:--------:|
| `adni-to-bids`    |    x     |
| `aibl-to-bids`    |    x     |
| `nifd-to-bids`    |    x     |
| `oasis-to-bids`   |          |

Please check the installation instructions for all platforms [here](https://github.com/rordenlab/dcm2niix#install).

Clinica requires dcm2niix version `1.0.20190902` or later.

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
| `dwi-connectome`         |   x  |           |     x      |  x  |     |        |    x    |        |     |
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

### Environment variables setup

When installing some of the third party software, environment variables might be needed in order to configure some installations.

If you simply run the `export` and `source` commands in your terminal, these environment variables will be defined only for the duration of your session, and will be lost if you re-launch your terminal.

In order to define these variables permanently, you need to manually edit the configuration file associated to your shell.

If you are using `bash`, it can be `~/.bashrc` or `~/.bash_profile`, if you are using `zsh`, it should be `~/.zshrc`.

!!! note "test before"
    Please do not copy/paste the provided commands without adapting them to your system and without testing them.
    Most paths provided here require to be adapted to your system, and also depends on how you installed the software.
    A good approach is to verify that the different paths you want to assign to a variable exist first.
    If they do, then try running the `export` and `source` commands in your terminal, and verify that the software run as expected.
    You can also verify that the pipeline you want to use is also running as expected.
    If this works, then consider modifying your shell configuration file to have these variables automatically defined on every session.


### ANTs

#### Installation

To install `ANTs`, download it from [here](https://github.com/stnava/ANTs/releases) and follow the instructions on the `ANTs` [wiki](https://github.com/stnava/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS).

#### Configuration

We strongly recommend installing `ANTs >= 2.5.0` from which **no environment variable are needed**.

Nonetheless, if you are using an older version of `ANTs`, make sure to have the following environment variables defined:

```bash
export ANTSPATH="/path/to/your/ANTs/"
export PATH=${ANTSPATH}:${PATH}
```

### Convert3D

You can find more details about `Convert3D` [here](http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D).

#### Installation

You have two options to install `Convert3D`:

- [Use pre-built binaries](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D).
- [Use the official conda package](https://anaconda.org/conda-forge/convert3d).

### Freesurfer

You can find more details about `Freesurfer` [here](http://surfer.nmr.mgh.harvard.edu/).

#### Installation

##### On Linux

Download and install `FreeSurfer` following the instructions on the [wiki](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall).

!!! note
    Please note that on Ubuntu you will need to install the packages `tcsh` and `libjpeg62` ( a `sudo apt-get install tcsh libjpeg62` should do the job).

##### On MacOS

Download it from [here](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and follow the instructions on the `FreeSurfer` [wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/MacOsInstall).

#### Configuration

Make sure to have the following environment variables defined:

```bash
export FREESURFER_HOME="/Applications/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh &> /dev/null
```

### FSL

We recommend installing [**FSL 6.0**](https://fsl.fmrib.ox.ac.uk/).

#### Installation

##### On Linux

Download it from [here](https://fsl.fmrib.ox.ac.uk/fsldownloads) and follow the instructions on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/Linux).

##### On MacOS

Download it from [here](https://fsl.fmrib.ox.ac.uk/fsldownloads) and follow the instructions on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/MacOsX).

#### Configuration

##### On Linux

Make sure to have the following environment variables defined:

```bash
export FSLDIR="/usr/share/fsl/6.0"
export PATH="${FSLDIR}/bin":${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh
```

##### On MacOS

Make sure to have the following environment variables defined:

```bash
export FSLDIR="/usr/local/fsl"
export PATH="${FSLDIR}/bin":${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh
```

### ITK

You can find more details about `ITK` [here](https://itk.org/).

#### Installation

##### On Linux

Follow the instructions on the [ITK blog](https://blog.kitware.com/itk-packages-in-linux-distributions/).

##### On MacOS

Follow the instructions on the [ITK blog](https://blog.kitware.com/kitware-packages-on-os-x-with-homebrew/).

### MRtrix3

You can find more details about `MRtrix3` [here](http://www.mrtrix.org).

#### Installation

You can find the official instructions on the `MRtrix` [website](https://www.mrtrix.org/download/).

##### On Linux

You have basically two options:

- [Use the official conda package](https://www.mrtrix.org/download/linux-anaconda/).
- [Build MRtrix3 from source](https://mrtrix.readthedocs.io/en/latest/installation/build_from_source.html#linux).

!!! note
    Note that using the conda package should be easier in most cases.

##### On MacOS

- [Use the official conda package](https://www.mrtrix.org/download/macos-anaconda/).
- [Use the MacOS pre-compiled application package installer](https://www.mrtrix.org/download/macos-application/).
- [Use the Homebrew formula](https://github.com/MRtrix3/homebrew-mrtrix3) (although large dependencies such as `XCode` and `Qt5` are required).

!!! note
    As for Linux, note that using the conda package should be easier in most cases.

### Matlab

You can find more details about `Matlab` [here](https://fr.mathworks.com/products/matlab/).

!!! note
    Note that using `Matlab` requires having a valid license which might be available through your university or institution.

#### Configuration

Make sure to have the following environment variables defined:

```bash
export MATLAB_HOME="/path/to/your/matlab/bin/"
export PATH=${MATLAB_HOME}:${PATH}
export MATLABCMD="${MATLAB_HOME}/matlab"
```

### PETPVC

You can find more details about `PETPVC` [here](https://github.com/UCL/PETPVC).

#### Installation

You can find the official instructions in the README of [this page](https://github.com/UCL/PETPVC).

You have basically three options:

- Use pre-built binaries, available [here](https://github.com/UCL/PETPVC/releases).
- [Use the official conda package](https://anaconda.org/conda-forge/petpvc).
- [Build from source](https://github.com/UCL/PETPVC?tab=readme-ov-file#installation-from-source-instructions).

!!! note
    If building from source, do not forget to compile in RELEASE mode, otherwise, partial volume correction will be very slow.

### SPM12

You can find more details about `SPM12` [here](http://www.fil.ion.ucl.ac.uk/spm/).

Note that `SPM12` works with [Matlab](#matlab) such that clinica pipelines which require `SPM12`, will also need a `Matlab` installation.

If you cannot install `Matlab`, you can install [SPM standalone](#spm12-standalone).

#### Installation

##### On Linux

Download the latest version [here](http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip) and follow the instructions on the [SPM wiki](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Linux).

##### On MacOS

Download the latest version [here](http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip) and follow the instructions on the [SPM wiki](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)).

!!! note
    For systems running on MacOS Big Sur, a [development version of SPM12](https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/) as well as a more recent release of the MCR (minimum 2019a) are required.

#### Configuration

Make sure to have the following environment variable defined:

```bash
export SPM_HOME="/path/to/your/spm12"
```

You must also add `SPM` to the `MATLAB` path variable if you installed it as a toolbox.

To do so, add the following line to your `startup.m` file located in your *initial working folder*, by default `~/Documents/MATLAB` (see [here](https://fr.mathworks.com/help/matlab/ref/startup.html) for more details).

If the file does not exist, you can create it and type inside:

```matlab
addpath('/path/to/your/spm12');
```

You can also replace the previous line by the following, assuming the `$SPM_HOME` environment variable is set in your `~/.bashrc` file.

```matlab
[~, spmhome] = system('source ~/.bashrc > /dev/null; echo $SPM_HOME;');
spmhome = strsplit(spmhome,'\n');
addpath(spmhome{end-1});
```

!!! Note
    `zsh` shell users will have to replace `~/.bashrc` by `~/.zshrc`.

### SPM12 standalone

If you want to install `SPM12` without installing [Matlab](#matlab), you will need to install two things:

- The Matlab runtime (often abbreviated into MCR), for which no license is required.
- The SPM standalone itself.

#### Installation

You can find details on how to install these on [this page](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/standalone/).

#### Configuration

!!! note
    If you followed the installation instructions, you should have set the environment variable `$LD_LIBRARY_PATH`.

In addition, you need to define the following environment variables:

```bash
export MCR_HOME="/path/to/your/MCR/"
export SPMSTANDALONE_HOME="/path/to/your/spmstandalone/home/"
```

## Autocompletion

<!-- # Autocomplete system
eval "$(register-python-argcomplete clinica)" -->

!!! Note
    `zsh` shell users will have to add this right before the last line of their configuration file to enable autocompletion in Clinica:

    ```bash
    autoload bashcompinit
    bashcompinit
    source ~/.bash_completion.d/python-argcomplete.sh
    ```

