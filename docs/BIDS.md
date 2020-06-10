# BIDS

## Introduction

[**BIDS**](http://bids.neuroimaging.io/) (Brain Imaging Data Structure) is the standard adopted for the organisation of the datasets used by Clinica pipelines through the command line. BIDS provides a unified structure for organising and describing neuroimaging and behavioural data. The use of a standard like BIDS makes easier developing and distributing code that uses neuroimaging datasets. For this reason, when using Clinica pipelines from the command line, the input format of the dataset is required to be BIDS-compliant.


## An overview of the BIDS structure

Here is a general overview of the BIDS structure. If you need more details, please check the [documentation](https://bids-specification.readthedocs.io/en/latest/) on the [website](http://bids.neuroimaging.io/).

```Text
BIDS_Dataset/
├── sub-CLNC01/
│   │   ├── ses-M00/
│   │   │   ├── anat/
│   │   │   │   ├── sub-CLNC01_ses-M00_T1w.json
│   │   │   │   └── sub-CLNC01_ses-M00_T1w.nii.gz
│   │   │   ├── dwi/
│   │   │   │   ├── sub-CLNC01_ses-M00_dwi.bval
│   │   │   │   ├── sub-CLNC01_ses-M00_dwi.bvec
│   │   │   │   ├── sub-CLNC01_ses-M00_dwi.json
│   │   │   │   ├── sub-CLNC01_ses-M00_dwi.nii.gz
│   │   │   │   └── ...
│   │   │   ├── fmap/
│   │   │   │   ├── sub-CLNC01_ses-M00_phasediff.json
│   │   │   │   ├── sub-CLNC01_ses-M00_phasediff.nii.gz
│   │   │   │   ├── sub-CLNC01_ses-M00_magnitude1.nii.gz
│   │   │   │   └── sub-CLNC01_ses-M00_magnitude2.nii.gz
│   │   │   ├── func/
│   │   │   │   ├── sub-CLNC01_ses-M00_task­-rest_bold.json
│   │   │   │   ├── sub-CLNC01_ses-M00_task­-rest_bold.nii.gz
│   │   │   │   └── ...
│   │   │   └── sub-CLNC01_ses-M00_scans.tsv
│   │   ├── ses-M18/
│   │   │   └── ...
│   │   └── sub-CLNC01_sessions.tsv
├── sub-CLNC02/
│   └── ...
├── ...
└── participants.tsv
```
