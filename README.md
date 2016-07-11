# Clinica

Welcome to the **Clinica** Software!

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Clinica](#clinica)
    - [Introduction](#introduction)
	- [Installation](#installation)
	- [Tasks performed by Clinica](#tasks-performed-by-clinica)
		- [Anatomical MRI](#anatomical-mri)
		- [Diffusion MRI](#diffusion-mri)
	- [Documentation](#documentation)
	- [License](#license)

<!-- /TOC -->

## <a name="introduction"></a> Introducation

**Clinica** is:
-  a set of neuroimaging pipelines build on top of [*nipype* framework](http://nipy.org/nipype/)
-  developped at [ICM](http://icm-institute.org) by the [Aramis Lab](http://www.aramislab.fr/)
-  written using Python 2.7
-  dependent of a set of tools like: FreeSurfer, SPM, ... TODO:adding all the dependences...

## <a name="installation"></a> Installation

Clinica should run on any Unix-like system. 

For installation instructions, please refer to this [Getting Started](https://gitlab.icm-institute.org/aramis/clinica/wikis/GettingStarted) page.

## <a name="tasks-performed-by-clinica"></a> Tasks performed by Clinica

The list of every example can be found [here](https://gitlab.icm-institute.org/aramis/clinica/wikis/ListOfExamples).

### <a name="anatomical-mri"></a> Anatomical MRI
* [Perform recon-all pipeline with T1 images in FreeSurfer](https://gitlab.icm-institute.org/aramis/clinica/wikis/T1_FreesurferReconAll)
* [Perform voxel-based morphometry with SPM](https://gitlab.icm-institute.org/aramis/clinica/wikis/T1_VoxelBasedMorphometry)

### <a name="diffusion-mri"></a> Diffusion MRI
* [Preprocesing of raw DWI dataset](https://gitlab.icm-institute.org/aramis/clinica/wikis/DWI_Preprocessing)
* [Processing of corrected DWI dataset (DTI, tractography)](https://gitlab.icm-institute.org/aramis/clinica/wikis/DWI_Processing)
* Postprocessing of tractogrammes ([computation of connectome](https://gitlab.icm-institute.org/aramis/clinica/wikis/connectome_construction)) and DTI derived scalar map ([atlas-based tracks scalar analysis](https://gitlab.icm-institute.org/aramis/clinica/wikis/DWI_WM_scalar_analysis))
* Statistical analysis of connectomes ([connection wise analysis](https://gitlab.icm-institute.org/aramis/clinica/wikis/connection_wise_analysis), [network-based statistics](https://gitlab.icm-institute.org/aramis/clinica/wikis/Network-based_statistics))

## <a name="documentation"></a> Documentation

To access a fully detailed documentation of **Clinica**... [TODO]

## <a name="license"></a> License

Clinica is a free software but..

TODO:each developer should add in this section the related licence he included in his pipeline! 


> **Notes:**
>
> _This wiki is under active development, further documentation will be added soon!_
