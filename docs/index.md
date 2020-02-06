# Clinica Documentation

## [What is Clinica?](WhatIsClinica)

## Installation

Clinica can be installed on Mac OS X and Linux (CentOS or Debian/Ubuntu) machines, and possibly on Windows computers with a Linux Virtual Machine. We assume that users installing and using Clinica are comfortable with using the command line.

<!-- !!! info "Clinica on GitHub!"
    Clinica has moved to GitHub, you can find the repo [here](https://github.com/aramis-lab/clinica/).-->

!!! info "New release: Clinica 0.3.2!"
    We are very pleased to announce the release 0.3.2 of Clinica. The release notes are available here: [v0.3.2](http://bit.ly/2tfZjvh).

<!-- ### Installing Clinica from source -->
  - [Installation](./Installation)
  - [Third-party software](./Third-party)  

<!-- ### Installing Clinica using Docker
Another way to install Clinica is to use [Docker](https://www.docker.com/what-docker). The installation procedure of the Clinica Docker image, which contains everything required to launch any pipeline of Clinica, is explained [here](https://gitlab.inria.fr/aramis/clinica_docker).    -->

<!-- ### Using Clinica on the ICM cluster
ICM members are encouraged to use the version of Clinica available on the cluster. Installation instructions are available [here](./ICMClusterInstallation). -->


## User documentation

### Clinica environment
- [Interacting with Clinica](InteractingWithClinica)
- [BIDS: the input data structure](BIDS)
- [CAPS: the processed data structure](CAPS/Introduction)

### Pipelines (`clinica run`)
- Anatomical MRI (T1-weighted)
    - `t1-volume` - [Processing of T1w MR images using SPM](Pipelines/T1_Volume): tissue segmentation and spatial normalization
    - `t1-freesurfer` - [Processing of T1w MR images using FreeSurfer](Pipelines/T1_FreeSurfer): cortical surface, subcortical structures and volumetrics
    <!--- `t1-freesurfer-longitudinal` - [Longitudinal processing of T1w MR images using FreeSurfer](Pipelines/T1_FreeSurfer): cortical surface, subcortical structures and volumetrics-->
- Diffusion MRI (DWI)
    - `dwi-preprocessing-*` - [DWI pre-processing](Pipelines/DWI_Preprocessing): correction of head motion, magnetic susceptibility, eddy current and bias field induced distortions
    - `dwi-dti` - [DTI scalar maps (FA, MD, AD, RD) and spatial normalization](Pipelines/DWI_DTI): extraction of DTI-based measures (FA, MD, AD, RD)
    - `dwi-connectome` - [Construction of structural connectome](Pipelines/DWI_Connectome): computation of fiber orientation distributions, tractogram and connectome
- PET
    - `pet-volume` - [Volume-based processing of PET images](Pipelines/PET_Volume): registration to T1w MRI, intensity normalization, partial volume correction and spatial normalization
    - `pet-surface` - [Surface-based processing of PET images](Pipelines/PET_Surface): projection of the PET signal onto the subject’s cortical surface
- Statistics
    - `statistics-surface` - [Surface-based mass-univariate analysis with SurfStat](Pipelines/Stats_Surface)
- Machine Learning
    - `machinelearning-prepare-spatial-svm` - [Prepare input data for spatially regularized SVM](Pipelines/MachineLearning_PrepareSVM)
    - [Classification based on machine learning](Pipelines/MachineLearning_Classification)

### Dataset converters (`clinica convert`)
- [Online neuroimaging databases (ADNI, AIBL, NIFD, OASIS) to BIDS converters](DatabasesToBIDS)

### I/O tools (`clinica iotools`)
- [Data handling tools for BIDS and CAPS compliant datasets](IO)

### Visualize pipeline outputs (`clinica visualize`)
Clinica allows visualization of the main outputs of some pipelines. Currently only supported for the [`t1-freesurfer` pipeline](Pipelines/T1_FreeSurfer).

## Clinica at conferences
Find on [this page](ClinicaConferences) the presentations and demo materials used when we showcase Clinica.

## Support
- [Report an issue on GitHub](https://github.com/aramis-lab/clinica/issues)
- Use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!

## License
Clinica is distributed under the terms of the MIT license given [here](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).

## Citing Clinica
For publications or communications using Clinica, please cite [[Routier et al](https://hal.inria.fr/hal-02308126/)] as well as the references mentioned on the wiki page of the pipelines you used. Each page includes text to cite the software packages that are used by Clinica (for example, citing SPM when using the `t1-volume` pipeline).

!!! info "Disclaimer"
    Clinica is a software for research studies. It is not intended for use in medical routine

---

![Clinica_Partners_Banner](img/Clinica_Partners_Banner.png)
