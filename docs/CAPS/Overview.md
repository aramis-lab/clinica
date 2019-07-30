# Overview
### At the subject level
This CAPS folder contains the outputs of the [`dwi-preprocessing-*` pipeline](../../Pipelines/DWI_Preprocessing) and [`dwi-dti` pipeline](../../Pipelines/DWI_DTI) of a fictional participant `CLNC01` at session `M00`:
```
subjects/
└── sub-CLNC01/
    └── ses-M00/
        └── dwi/
            ├── dti_based_processing/
            │   ├── atlas_statistics/
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-FA_statistics.tsv
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-MD_statistics.tsv
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-AD_statistics.tsv
            │   │   └── sub-CNLC01_ses-M00_acq-axial_dwi_space-JHUTracts0_res-1x1x1_map-RD_statistics.tsv
            │   ├── native_space/
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_model-DTI_diffmodel.nii.gz
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_FA.nii.gz
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_MD.nii.gz
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_AD.nii.gz
            │   │   ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_RD.nii.gz
            │   │   └── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_DECFA.nii.gz
            │   └── normalized_space/
            │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_AD.nii.gz
            │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_affine.mat
            │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_deformation.nii.gz
            │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_FA.nii.gz
            │       ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_MD.nii.gz
            │       └── sub-CNLC01_ses-M00_acq-axial_dwi_space-MNI152Lin_res-1x1x1_RD.nii.gz
            └── preprocessing/
                ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_brainmask.nii.gz
                ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.bval
                ├── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.bvec
                └── sub-CNLC01_ses-M00_acq-axial_dwi_space-b0_preproc.nii.gz
```

## At the group level
This CAPS folder contains the outputs of a group comparison of patients with Alzheimer’s disease (`AD`) and healthy subjects (`HC`) thanks to the [`statistics-surface` pipeline](../../Pipelines/Stats_Surface). Results are stored under the group ID `ADvsHC`:

```
groups/
└── group-ADvsHC/
    └── statistics/
        ├── participant.tsv
        └── surfstat_group_comparison/
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_FDR.jpg
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_FDR.mat
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_TStatistics.jpg
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_TStatistics.mat
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_correctedPValue.jpg
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_correctedPValue.mat
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_uncorrectedPValue.jpg
            ├── group-ADvsHC_HC-lt-AD_measure-ct_fwhm-20_uncorrectedPValue.mat
            ├── group-ADvsHC_glm.json
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_FDR.jpg
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_FDR.mat
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_TStatistics.jpg
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_TStatistics.mat
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_correctedPValue.jpg
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_correctedPValue.mat
            ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_uncorrectedPValue.jpg
            └── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_uncorrectedPValue.mat
```
