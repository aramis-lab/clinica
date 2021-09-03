<!-- markdownlint-disable MD007 -->
# Clinica Documentation

## [What is Clinica?](WhatIsClinica)

## Installation

Clinica can be installed on Mac OS X and Linux (CentOS or Debian/Ubuntu) machines,
and possibly on Windows computers with a Linux Virtual Machine.
We assume that users installing and using Clinica are comfortable with using the command line.

<!--!!! info "New release: Clinica 0.3.6!"
    We are very pleased to announce the release 0.3.6 of Clinica. The release notes are available here: [v0.3.6](http://bit.ly/2tfZjvh).-->

<!-- ### Installing Clinica from source -->

- [Installation](./Installation)
- [Third-party software](./Third-party)  

<!--
### Installing Clinica using Docker
Another way to install Clinica is to use [Docker](https://www.docker.com/what-docker).
The installation procedure of the Clinica Docker image, which contains everything required to launch any pipeline of Clinica, is explained [here](https://gitlab.inria.fr/aramis/clinica_docker).
-->

<!--
### Using Clinica on the ICM cluster
ICM members are encouraged to use the version of Clinica available on the cluster.
Installation instructions are available [here](./ICMClusterInstallation).
-->

## User documentation

### Clinica environment

- [Interacting with Clinica](InteractingWithClinica)
- [BIDS: the input data structure](BIDS)
- [CAPS: the processed data structure](CAPS/Introduction)

### Pipelines (`clinica run`)

- Anatomical MRI (T1-weighted)
    - `t1-linear` - [Linear processing of T1w MR images](Pipelines/T1_Linear): affine registration to the MNI standard space
    - `t1-volume` - [Processing of T1w MR images using SPM](Pipelines/T1_Volume): tissue segmentation and spatial normalization
    - `t1-freesurfer` - [Processing of T1w MR images using FreeSurfer](Pipelines/T1_FreeSurfer): cortical surface, subcortical structures and volumetrics
    - `t1-freesurfer-longitudinal` - [Longitudinal processing of T1w MR images using FreeSurfer](Pipelines/T1_FreeSurfer_Longitudinal): cortical surface, subcortical structures and volumetrics
- Diffusion MRI (DWI)
    - `dwi-preprocessing-*` - [DWI pre-processing](Pipelines/DWI_Preprocessing): correction of head motion, magnetic susceptibility, eddy current and bias field induced distortions
    - `dwi-dti` - [DTI scalar maps (FA, MD, AD, RD) and spatial normalization](Pipelines/DWI_DTI): extraction of DTI-based measures (FA, MD, AD, RD)
    - `dwi-connectome` - [Construction of structural connectome](Pipelines/DWI_Connectome): computation of fiber orientation distributions, tractogram and connectome

- PET
    - [Introduction to concepts used in the PET pipelines](Pipelines/PET_Introduction): partial volume correction and standardized uptake value ratio (SUVR) map computation
    - `pet-linear` - [Linear processing of PET images](Pipelines/PET_Linear): affine registration to the MNI standard space and intensity normalization
    - `pet-volume` - [Volume-based processing of PET images](Pipelines/PET_Volume): registration to T1w MRI, intensity normalization, partial volume correction and spatial normalization
    - `pet-surface` - [Surface-based processing of PET images](Pipelines/PET_Surface): projection of the PET signal onto the subject’s cortical surface
    - `pet-surface-longitudinal` - [Surface-based longitudinal processing of PET images](Pipelines/PET_Surface_Longitudinal): projection of the PET signal onto the subject’s cortical surface

- Statistics
    - `statistics-surface` - [Surface-based mass-univariate analysis with SurfStat](Pipelines/Stats_Surface)
    - `statistics-volume` - [Volume-based mass-univariate analysis with SPM](Pipelines/Stats_Volume)

- Machine Learning
    - `machinelearning-prepare-spatial-svm` - [Prepare input data for spatially regularized SVM](Pipelines/MachineLearning_PrepareSVM)
    - `machinelearning-classification` - [Classification based on machine learning](Pipelines/MachineLearning_Classification)

- Deep Learning
    - `deeplearning-prepare-data` - [Prepare input data for deep learning with PyTorch](Pipelines/DeepLearning_PrepareData)

### Dataset converters (`clinica convert`)

Clinica provides tools to curate several publicly available neuroimaging datasets and convert them to BIDS namely:

- `adni-2-bids` - [ADNI: Alzheimer’s Disease Neuroimaging Initiative](Converters/ADNI2BIDS)
- `aibl-2-bids` - [AIBL: Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing](Converters/AIBL2BIDS)
- `nifd-2-bids` - [NIFD: Neuroimaging in Frontotemporal Dementia](Converters/NIFD2BIDS)
- `oasis-2-bids` - [OASIS: Open Access Series of Imaging Studies](Converters/OASIS2BIDS)
- `oasis3-2-bids` - [OASIS-3: Longitudinal Neuroimaging, Clinical, and Cognitive Dataset for Normal Aging and Alzheimer’s Disease](Converters/OASIS3TOBIDS)

!!! note
    We provide converters for the datasets used in the [Aramis Lab](http://www.aramislab.fr/).
    Feel free to contact us if you are interested in another dataset or to contribute!

### I/O tools (`clinica iotools`)

- [Data handling tools for BIDS and CAPS compliant datasets](IO)

### Visualize pipeline outputs (`clinica visualize`)

Clinica allows visualization of the main outputs of some pipelines.
Currently only supported for the [`t1-freesurfer` pipeline](Pipelines/T1_FreeSurfer).

## Clinica at conferences

Find on [this page](ClinicaConferences) the presentations and demo materials used when we showcase Clinica.

## Support

- Check for [past answers](https://groups.google.com/forum/#!forum/clinica-user) in the old Clinica Google Group
- Start a [discussion](https://github.com/aramis-lab/clinica/discussions) on Github
- Report an [issue](https://github.com/aramis-lab/clinica/issues) on GitHub

## License

Clinica is distributed under the terms of the MIT license given
[here](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).

## Citing Clinica

For publications or communications using Clinica, please cite
[[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)]
as well as the references mentioned on the wiki page of the pipelines you used.
Each page includes text to cite the software packages that are used by Clinica
(for example, citing SPM when using the `t1-volume` pipeline).

!!! info "Disclaimer"
    Clinica is a software for research studies.
    It is not intended for use in medical routine.

---

![Clinica_Partners_Banner](img/Clinica_Partners_Banner.png)
