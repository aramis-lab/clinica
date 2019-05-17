# CAPS

Clinica has its own specifications for storing processed data, called [**CAPS** (ClinicA Processed Structure)](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8), inspired from the BIDS structure.


## Motivation
The BIDS specifications currently do not provide specific rules for the processed data, so the goal of CAPS (designed by the [Aramis Lab](http://www.aramislab.fr/)) is to define a hierarchy for the Clinica processed data. The idea is to include in a single folder all the results of the different pipelines and organize the data following the main patterns of the BIDS specification.

!!! note
    There is an ongoing initiative called BIDS-derivatives that aim to provide a BIDS standard for processed data. However, in its current state, it is not well adapted to Clinica outputs. This is why Clinica has its own specifications called CAPS. We try to contribute to BIDS-derivatives and the two specifications should ultimately converge.

## General overview (example)
### At the subject level
This CAPS folder contains the outputs of the [`dwi-preprocessing-*`](../Pipelines/DWI_Preprocessing) and [`dwi-dti`](../Pipelines/DWI_DTI) of a fictional participant `CLNC01` at session `M00`:

```
CAPS
└── subjects
    ├── sub-CNLC01
    │   ├── ses-M00
    │   │   └── dwi
    │   │       ├── dti_based_processing
    │   │       │   ├── atlas_statistics
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-AD_statistics.tsv
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-FA_statistics.tsv
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-MD_statistics.tsv
    │   │       │   │   └── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-RD_statistics.tsv
    │   │       │   ├── native_space
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_model-dti_diffmodel.nii.gz
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_AD.nii.gz
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_DECFA.nii.gz
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_FA.nii.gz
    │   │       │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_MD.nii.gz
    │   │       │   │   └── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_RD.nii.gz
    │   │       │   └── normalized_space
    │   │       │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_AD.nii.gz
    │   │       │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_affine.mat
    │   │       │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_deformation.nii.gz
    │   │       │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_FA.nii.gz
    │   │       │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_MD.nii.gz
    │   │       │       └── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_RD.nii.gz
    │   │       └── preprocessing
    │   │           ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_brainmask.nii.gz
    │   │           ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.bval
    │   │           ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.bvec
    │   │           └── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.nii.gz
    │   └── ses-M18
    │       └── ...
    └── sub-CLNC02
        └── ...
```


### At the group level
This CAPS folder contains the outputs of a group comparison of patients with Alzheimer’s disease (`AD`) and healthy subjects (`Control`) thanks to the [`statistics-surface`](../Pipelines/Stats_Surface) pipeline. Results are stored under the group ID `ADvsHC`:

```
CAPS
└── groups
    └── group-ADvsHC
        └── statistics
            ├── participant.tsv
            └── surfstat_group_comparison
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_FDR.jpg
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_FDR.mat
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_TStatistics.jpg
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_TStatistics.mat
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_correctedPValue.jpg
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_correctedPValue.mat
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_uncorrectedPValue.jpg
                ├── group-ADvsHC_Control-lt-AD_measure-ct_fwhm-20_uncorrectedPValue.mat
                ├── group-ADvsHC_glm.json
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_FDR.jpg
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_FDR.mat
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_TStatistics.jpg
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_TStatistics.mat
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_correctedPValue.jpg
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_correctedPValue.mat
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_uncorrectedPValue.jpg
                ├── group-ADvsHC_AD-lt-Control_measure-ct_fwhm-20_uncorrectedPValue.mat
                └── matlab_output.log
```

!!! note
    This is just a general overview of the CAPS structure. For detailed information, please see the document [The ClinicA Processed Structure (CAPS) Specification](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8/edit#).
