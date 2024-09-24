<!-- markdownlint-disable MD046 -->
# Third-party software

## Environment variables set-up

As you will see on this page, environment variables are often needed in order to configure some third-party software installations.

1. Run the provided `export` and `source` commands in your terminal. Environment variables will be defined only for the duration of your session, 
and will be lost if you re-launch your terminal.

    ??? info "Learn how to set environment variables"
        - `export var_name = value` : defines an environment variable as value and makes it accessible to other processes started from the terminal.
        - `source file` : used to read/execute a file, for example one where configurations are set

2. In order to define these variables permanently, you need to manually edit the configuration file associated to your shell.

    If you are using `bash`, it can be `~/.bashrc` or `~/.bash_profile`, if you are using `zsh`, it should be `~/.zshrc`.

!!! danger "Provided command lines should be adapted"
    **Please** do not copy/paste the provided commands without adapting them to your system and without testing them.
    Most paths provided here require to be adapted to your system and depend on how you installed the software.
    
    A good approach is to verify that the different paths you want to assign to a variable exist first : 

    → If they do, then try running the `export` and `source` **(step 1. above)** commands in your terminal. Verify that the software or the pipeline you want runs as expected.
    
    → If this works, then consider modifying your shell configuration file **(step 2. above)** to have these variables automatically defined on every session.


## Converters

Some converters require **dcm2niix** to transform DICOM files into NIfTI :

<div class="grid">
  <a href="./../../Converters/ADNI2BIDS/index.html" class="card">adni-to-bids</a>
  <a href="./../Converters/AIBL2BIDS/index.html" class="card">aibl-to-bids</a>
  <a href="./../Converters/GENFItoBIDS/index.html" class="card">genfi-to-bids</a>
  <a href="./../Converters/NIFD2BIDS/index.html" class="card">nifd-to-bids</a>
  <a href="./../Converters/UKBtoBIDS/index.html" class="card">ukb-to-bids</a>
</div>

### DCM2NIX

!!! warning "Version required"
    Clinica requires dcm2niix version `1.0.20190902` or later.

Please check the installation instructions for all platforms on [dcm2niix Git repository](https://github.com/rordenlab/dcm2niix#install).



## Pipelines

<div class="annotate" markdown>
Some, but not all pipelines use specific third-party software. Depending on your usage of Clinica, you will need to install additional packages.
Specific dependencies are described in the table below (1) :
</div>

1. If not listed, the pipeline does not require any additional dependency outside Clinica.

??? info "Clinica available pipelines"
    --8<-- "snippets/inventory_pipelines.md"

<div markdown="1" class="third_party_table">

|                     | ANTs | Convert3D | FreeSurfer | FSL | ITK | Matlab | MRtrix3 | PETPVC | SPM |
|:-------------------:|:----:|:---------:|:----------:|:---:|:---:|:------:|:-------:|:------:|:---:|
|    Anat > Linear    |  ✓∘  |           |            |     |     |        |         |        |     |
|    Anat > Volume    |      |           |            |     |     |   ✓    |         |        |  ✓  |
|  Anat > FreeSurfer  |      |           |     ✓      |     |     |        |         |        |     |
| DWI > Preprocessing |  ✓   |     ✓     |            |  ✓  |     |        |    ✓    |        |     |
|      DWI > DTI      |  ✓   |           |            |  ✓  |     |        |    ✓    |        |     |
|  DWI > Connectome   |      |           |     ✓      |  ✓  |     |        |    ✓    |        |     |
|    PET > Linear     |  ✓   |           |            |     |     |        |         |        |     |
|    PET > Surface    |      |           |     ✓      |  ✓  | ✓⟡  |   ✓    |         |   ✓⟡   |  ✓  |
|    PET > Volume     |      |           |            |     | ✓⟡  |   ✓    |         |   ✓⟡   |  ✓  |
|   Stats > Surface   |      |           |     ✓      |     |     |        |         |        |     |
|   Stats > Volume    |      |           |            |     |     |   ✓    |         |        |  ✓  |

</div>

- *✓∘ : for anatomical linear pipelines there is also the possibility to use ANTsPy instead of ANTs since Clinica `v0.9.0`* 
- *✓⟡ : you only need to install ITK if you plan to perform partial volume correction using PETPVC.*

??? warning "CAT12 toolbox and Clinica < `v0.3.7`"
    Starting from Clinica `v0.3.7`, the [**CAT12**](http://dbm.neuro.uni-jena.de/cat/) toolbox is no longer needed for the `t1-volume` and `pet-volume` pipelines.
    For previous versions of Clinica, you will need to download the latest version of the toolbox [here](http://dbm.neuro.uni-jena.de/cat/index.html#DOWNLOAD) and follow the instructions to ensure that your `cat12` folder is located in your `spm/toolbox` folder.

Depending on the architecture and OS of your system, setup of third party libraries can change.
Please refer to each tool’s website for installation instructions :
___

### ANTs

To install `ANTs`, download it from [ANTs release list](https://github.com/stnava/ANTs/releases) and follow the instructions on [ANTs wiki](https://github.com/stnava/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS).

=== "`ANTs >= 2.5.0` 🔺"
    We **strongly** recommend installing `ANTs >= 2.5.0` from which **no environment variable are needed**.

=== "`ANTs < 2.5.0`"
    Nonetheless, if you are using an older version of `ANTs`, make sure to have the following environment variables defined:
    
    ```{ .bash .copy }
    export ANTSPATH="/path/to/your/ANTs/"
    export PATH=${ANTSPATH}:${PATH}
    ```

___

### Convert3D

You can find more details about `Convert3D` on their [website](http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D). There are
two options to install it :

- [Use pre-built binaries](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D).
- [Use the official conda package](https://anaconda.org/conda-forge/convert3d).

___

### FreeSurfer

You can find more details about `FreeSurfer` on their [website](http://surfer.nmr.mgh.harvard.edu/). To install it :

=== "Linux"
    Download and install `FreeSurfer` following the instructions on the [wiki](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall).

    !!! warning
        Please note that on Ubuntu you will need to install the packages `tcsh` and `libjpeg62` ( a `sudo apt-get install tcsh libjpeg62` should do the job).

=== "MacOS"
    Download it from [here](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and follow the instructions on the `FreeSurfer` [wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/MacOsInstall).


Make sure to have the following environment variables defined:

```{ .bash .copy }
export FREESURFER_HOME="/Applications/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh &> /dev/null
```

___

### FSL

We recommend installing [**FSL 6.0**](https://fsl.fmrib.ox.ac.uk/).

=== "Linux"
    Download it from [here](https://fsl.fmrib.ox.ac.uk/fsldownloads) and follow the instructions on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/Linux).

    Make sure to have the following environment variables defined:
    ```{ .bash .copy }
    export FSLDIR="/usr/share/fsl/6.0"
    export PATH="${FSLDIR}/bin":${PATH}
    source ${FSLDIR}/etc/fslconf/fsl.sh
    ```

=== "MacOS"
    Download it from [here](https://fsl.fmrib.ox.ac.uk/fsldownloads) and follow the instructions on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/MacOsX).

    Make sure to have the following environment variables defined:
    ```{ .bash .copy }
    export FSLDIR="/usr/local/fsl"
    export PATH="${FSLDIR}/bin":${PATH}
    source ${FSLDIR}/etc/fslconf/fsl.sh
    ```

___

### ITK

You can find more details about `ITK` on their [website](https://itk.org/). To install it :

=== "Linux"
    Follow the instructions on the [ITK blog](https://blog.kitware.com/itk-packages-in-linux-distributions/).

=== "MacOS"
    Follow the instructions on the [ITK blog](https://blog.kitware.com/kitware-packages-on-os-x-with-homebrew/).

___

### MRtrix3

You can find more details about `MRtrix3` on their [website](http://www.mrtrix.org), including
official instructions for [downloading](https://www.mrtrix.org/download/).

=== "Linux"
    You have basically two options:

    - [Use the official conda package](https://www.mrtrix.org/download/linux-anaconda/).
    - [Build MRtrix3 from source](https://mrtrix.readthedocs.io/en/latest/installation/build_from_source.html#linux).

    !!! tip
        Note that using the conda package should be easier in most cases.

=== "MacOS"

    - [Use the official conda package](https://www.mrtrix.org/download/macos-anaconda/).
    - [Use the MacOS pre-compiled application package installer](https://www.mrtrix.org/download/macos-application/).
    - [Use the Homebrew formula](https://github.com/MRtrix3/homebrew-mrtrix3) (although large dependencies such as `XCode` and `Qt5` are required).

___

### Matlab

You can find more details about `Matlab` on their [website](https://fr.mathworks.com/products/matlab/).

!!! warning
    Note that using `Matlab` requires having a **valid license** which might be available through your university or institution.

Make sure to have the following environment variables defined:

```{ .bash .copy }
export MATLAB_HOME="/path/to/your/matlab/bin/"
export PATH=${MATLAB_HOME}:${PATH}
export MATLABCMD="${MATLAB_HOME}/matlab"
```

___

### PETPVC

You can find more details about `PETPVC` on their [website](https://github.com/UCL/PETPVC),
including official instructions for downloading on their [Github](https://github.com/UCL/PETPVC).

You have basically three options:

- Use pre-built binaries, available [here](https://github.com/UCL/PETPVC/releases).
- [Use the official conda package](https://anaconda.org/conda-forge/petpvc).
- [Build from source](https://github.com/UCL/PETPVC?tab=readme-ov-file#installation-from-source-instructions).

!!! tip
    If building from source, do not forget to compile in **RELEASE** mode, otherwise, partial volume correction will be very slow.

___

### SPM12

!!! warning "SPM12 and Matlab"
    Note that `SPM12` works with [Matlab](#matlab) such that Clinica pipelines which require `SPM12`, will also need a `Matlab` installation.
    If you cannot install `Matlab`, you can install [SPM standalone](#spm12-standalone).

You can find more details about `SPM12` on their [website](http://www.fil.ion.ucl.ac.uk/spm/). To install it :

=== "Linux"
    Download the latest version [here](http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip) and follow the instructions on the [SPM wiki](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Linux).

=== "MacOS"
    Download the latest version [here](http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip) and follow the instructions on the [SPM wiki](https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)).

    ??? warning "MacOS Big Sur"
        For systems running on MacOS Big Sur, a [development version of SPM12](https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/) as well as a more recent release of the MCR (minimum 2019a) are required.

Then :

1. Make sure to have the following environment variable defined:
    ```{ .bash .copy }
    export SPM_HOME="/path/to/your/spm12"
    ```

   2. You must also add `SPM` to the `MATLAB` path variable if you installed it as a toolbox.
      - First, verify a `startup.m` file exists in your *initial working folder*, by default `~/Documents/MATLAB` (see [here](https://fr.mathworks.com/help/matlab/ref/startup.html) for more details). 
        If the file does not exist, you can create it.
      - Then, type inside:
       ```matlab
       addpath('/path/to/your/spm12');
       ```
       **OR**, assuming the `$SPM_HOME` environment variable is set in your `~/.bashrc` file :
       ```matlab
       [~, spmhome] = system('source ~/.bashrc > /dev/null; echo $SPM_HOME;');
       spmhome = strsplit(spmhome,'\n');
       addpath(spmhome{end-1});
       ```

    !!! warning
        `zsh` shell users will have to replace `~/.bashrc` by `~/.zshrc`.

___

### SPM12 standalone

If you want to install `SPM12` without installing [Matlab](#matlab), you will need to install two things:

- The Matlab runtime (often abbreviated into MCR), for which no license is required.
- The SPM standalone itself.

You can find details on how to install these on [this page](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/standalone/).

!!! tip "Verify your environment variables"
    If you followed the installation instructions, you should have set the environment variable `$LD_LIBRARY_PATH`.

In addition, you need to define the following environment variables:

```bash
export MCR_HOME="/path/to/your/MCR/"
export SPMSTANDALONE_HOME="/path/to/your/spmstandalone/home/"
```

___

### Autocompletion (optional)

```{ .bash .copy }
eval "$(register-python-argcomplete clinica)"
```

!!! warning "Autocompletion for `zsh` shell users"
    `zsh` shell users will have to add this right before the last line of their configuration file to enable autocompletion in Clinica:

    ```{ .bash .copy }
    autoload bashcompinit
    bashcompinit
    source ~/.bash_completion.d/python-argcomplete.sh
    ```


--8<-- "snippets/known_issues.md:matlab"
