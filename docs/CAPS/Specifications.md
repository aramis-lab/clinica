# Detailed file descriptions

In the following, brackets `[`/`]` will denote optional key/value pairs in the filename while accolades `{`/`}` will indicate a list of compulsory values (e.g. `hemi-{left|right}` means that the key `hemi` only accepts `left` or `right` as values).

##  T1 MRI data

### `t1-volume` pipeline - Volume-based processing of T1-weighted MR images

#### Segmentation
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── spm/
                └── segmentation/
                    ├── normalized_space/
                    │   ├── <source_file>_target-Ixi549Space_transformation-{inverse|forward}_deformation.nii.gz
                    │   ├── <source_file>_segm-<segm>_space-Ixi549Space_modulated-{on|off}_probability.nii.gz
                    │   └── <source_file>_space-Ixi549Space_T1w.nii.gz
                    ├── native_space/
                    │   └── <source_file>_segm-<segm>_probability.nii.gz
                    └── dartel_input/
                        └── <source_file>_segm-<segm>_dartelinput.nii.gz
```

The `modulated-{on|off}` key indicates if modulation has been used in SPM to compensate for the effect of spatial normalization.

The possible values for the `segm-<segm>` key/value are: `graymatter`, `whitematter`, `csf`, `bone`, `softtissue`, and `background`.

The T1 image in `Ixi549Space` (reference space of the TPM) is obtained by applying the transformation obtained from the SPM Segmentation routine to the T1 image in native space.


#### DARTEL
```
groups/
└── group-<group_label>/
    ├── group-<group_label>_subjects_visits_list.tsv
    └── t1/
        ├── group-<group_label>_iteration-<index>_template.nii.gz
        └── group-<group_label>_template.nii.gz
```

The final group template is `group-<group_label>_template.nii.gz`.

The `group-<group_label>_iteration-<index>_template.nii.gz` obtained at each iteration will only be used when obtaining flow fields for registering a new image into an existing template (SPM DARTEL Existing templates procedure).

!!! Note "Note for SPM experts"
    The original name of `group-<group_label>_iteration-<index>_template.nii.gz` is `Template<index>.nii`.


#### DARTEL to MNI
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── spm/
                └── dartel/
                    └── group-<group_label>/
                        ├── <source_file>_target-<group_label>_transformation-forward_deformation.nii.gz
                        └── <source_file>_segm-<segm>_space-Ixi549Space_modulated-{on|off}[_fwhm-<X>mm]_probability.nii.gz
```

#### Atlas statistics
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── spm/
                └── dartel/
                    └── group-<group_label>/
                        └── atlas_statistics/
                            └── <source_file>_space-<space>_map-graymatter_statistics.tsv
```
Statistics files (with `_statistics.tsv` suffix) are detailed in [appendix](#appendix-content-of-a-statistic-file).


### `t1-freesurfer` - FreeSurfer-based processing of T1-weighted MR images
The outputs of the `t1-freesurfer` pipeline are split into two subfolders, the first one containing the FreeSurfer outputs and a second with additional outputs specific to Clinica.

FreeSurfer outputs:
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── freesurfer_cross_sectional/
                └── sub-<participant_label>_ses-<session_label>/
                    ├── label/
                    ├── mri/
                    ├── scripts/
                    ├── stats/
                    └── surf/
```

Clinica additional outputs:
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── freesurfer_cross_sectional/
                └── regional_measures/
                    ├── <source_file>_parcellation-wm_volume.tsv
                    ├── <source_file>_segmentationVolumes.tsv
                    └── <source_file>_hemi-{left|right}_parcellation-<parcellation>_thickness.tsv
```

For the file: `*_hemi-{left|right}_parcellation-<parcellation>_thickness.tsv`, `_thickness` is just an example for the cortical thickness, we also have other measurements defined in the table below:


| Name                    | Suffix       | Description |
|-------------------------|--------------|-------------|
| Cortical thickness      | `_thickness` | Cortical thickness between pial surface and white surface|
| Cortical volume         | `_volume`    | Volume of gray matter|
| Cortical surface area   | `_area`      | Cortical surface area|
| Cortical mean curvature | `_meancurv`  | Mean curvature of cortical surface|


The `hemi-{left|right}` key/value stands for `left` or `right` hemisphere.

The possible values for the `parcellation-<parcellation>` key/value are: `desikan` (Desikan-Killiany Atlas), `destrieux` (Destrieux Atlas) and `ba` (Brodmann Area Maps). The TSV files for Brodmann areas contain a selection of regions, see this link for the content of this selection: [http://ftp.nmr.mgh.harvard.edu/fswiki/BrodmannAreaMaps](http://ftp.nmr.mgh.harvard.edu/fswiki/BrodmannAreaMaps)).

The details of the white matter parcellation of FreeSurfer can be found here: [https://surfer.nmr.mgh.harvard.edu/pub/articles/salat_2008.pdf](https://surfer.nmr.mgh.harvard.edu/pub/articles/salat_2008.pdf).


!!! Example "Example - Content of the TSV files"
    Content of `sub-CLNC01_ses­-M00_T1w_segmentationVolumes.tsv`:
    ```
    Measure:volume	Left-Lateral-Ventricle	Left-Inf-Lat-Vent	...
    /path/to/freesurfer/segmentation/	12345.6	12.334	...
    ```

    This file contains the volume of the different subcortical structures after segmentation.

    Content of `sub-CLNC01_ses­-M00_T1w_parcellation-wm_volume.tsv`:
    ```
    Measure:volume	wm-lh-bankssts	wm-lh-caudalanteriorcingulate ...
    /path/to/freesurfer/wm parcellation/ 2474.6	1863.7 ...
    ```
    This file contains the volume of the different white matter regions after parcellation.


    Content of `sub-CLNC01_ses­-M00_hemi-left_parcellation-desikan_thickness.tsv`:
    ```
    lh.aparc.thickness	lh_bankssts_thickness	lh_caudalanteriorcingulate_thickness …
    /path/to/freesurfer/cortical thickness/parcellation 2.048 2.892 …
    ```
    This file contains the cortical thickness in different regions of the Desikan atlas.


### `t1-freesurfer-longitudinal` – FreeSurfer-based longitudinal processing of T1-weighted MR images
The outputs of the `t1-freesurfer-longitudinal` pipeline are split into three subfolders, the first one containing the FreeSurfer unbiased template, the second containing the FreeSurfer longitudinal outputs and a third with additional outputs specific to Clinica.


FreeSurfer unbiased template:
```
subjects/
└── sub-<participant_label>/
    └── long-<long_label>/
        └── freesurfer_unbiased_template/
            └── sub-<participant_label>/
                ├── label/
                ├── mri/
                ├── scripts/
                ├── stats/
                └── surf/
```

FreeSurfer longitudinal outputs:
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── long-<long_label>/
                └── freesurfer_longitudinal/
                    └── sub-<participant_label>_ses-<session_label>.long.sub-<participant_label>/
                        ├── label/
                        ├── mri/
                        ├── scripts/
                        ├── stats/
                        └── surf/
```

Clinica additional outputs:
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── t1/
            └── long-<long_label>/
                └── freesurfer_longitudinal/
                    └── regional_measures/
                        ├── <source_file>_parcellation-wm_volume.tsv
                        ├── <source_file>_segmentationVolumes.tsv
                        └── <source_file>_hemi-{left|right}_parcellation-<parcellation>_thickness.tsv
```
where each file is explained in the `t1-freesurfer` sub-section.

!!! Note
    The naming convention `<subject_name>.long.<template_name>` is imposed by FreeSurfer.


## Diffusion imaging data
### `dwi-preprocessing-*` - Preprocessing of raw diffusion weighted imaging (DWI) datasets
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── dwi/
            └── preprocessing/
                ├── <source_file>_space-<space>_preproc.bval
                ├── <source_file>_space-<space>_preproc.bvec
                ├── <source_file>_space-<space>_preproc.nii.gz
                └── <source_file>_space-<space>_brainmask.nii.gz
```

The resulting DWI file after preprocessing. According to the subtype of pipeline run, `<space>` can be `T1w` ([`dwi-preprocessing-using-t1` pipeline](../../Pipelines/DWI_Preprocessing)) or `b0` ([`dwi-preprocessing-using-fieldmap` pipeline](../../Pipelines/DWI_Preprocessing)). A brain mask of the preprocessed file is provided.


### `dwi-dti` - DTI-based processing of corrected DWI datasets
```
subjects/
└── sub-<participant_label>/
    └── ses­-<session_label>/
        └── dwi/
            └── dti_based_processing/
                ├── native_space/
                │   ├── <source_file>_space-<space>_model-DTI_diffmodel.nii.gz
                │   ├── <source_file>_space-<space>_{FA|MD|AD|RD}.nii.gz
                │   └── <source_file>_space-<space>_DECFA.nii.gz
                │── normalized_space/
                │   ├── <source_file>_space-MNI152Lin_res-1x1x1_affine.mat
                │   ├── <source_file>_space-MNI152Lin_res-1x1x1_deformation.nii.gz
                │   └── <source_file>_space-MNI152Lin_res-1x1x1_{FA|MD|AD|RD}.nii.gz
                └── atlas_statistics/
                    └── <source_file>_space-<space>_res-1x1x1_map-{FA|MD|AD|RD}_statistics.tsv
```

The DTI is saved under the `*_model-DTI_diffmodel.nii.gz` filename. The different maps based on the DTI are the fractional anisotropy (`FA`), mean diffusivity (`MD`), axial diffusivity (`AD`) and radial diffusivity (`RD`) parametric maps, as well as the directionally-encoded colour (DEC) FA (`DECFA`) map.

Current atlases used for statistics are the 1mm version of `JHUDTI81`, `JHUTract0` and `JHUTract25` (see [Atlases page](../../Atlases) for further details).

Statistics files (with `_statistics.tsv` suffix) are detailed in [appendix](#appendix-content-of-a-statistic-file).


!!! Note
    The naming convention for suffixes follows the BIDS derivative specifications except for the statistics files (specific files for our needs) and the `_deformation.nii.gz` file (it is equivalent to the `_warp.nii.gz` file in the [BEP014 specifications](https://docs.google.com/document/d/11gCzXOPUbYyuQx8fErtMO9tnOKC3kTWiL9axWkkILNE)).


### `dwi-connectome` - Computation of structural connectome from corrected DWI datasets
```
subjects/
└── sub-<participant_label>/
    └── ses­-<session_label>/
        └── dwi/
            └── connectome_based_processing/
                ├── <source_file>_space-{b0|T1w}_model-CSD_diffmodel.nii.gz
                │── <source_file>_space-{b0|T1w}_model-CSD_tractography.tck
                └── <source_file>_space-{b0|T1w}_model-CSD_parcellation-{desikan|destrieux}_connectivity.tsv
```

The constrained spherical deconvolution (CSD) diffusion model is saved under the `*_model-CSD_diffmodel.nii.gz` filename. The whole-brain tractography is saved under the `*_tractography.tck` filename. The connectivity matrices are saved under `*_connectivity.tsv` filenames.

Current parcellations used for the computation of connectivity matrices are `desikan` and `desikan` (see [Atlases page](../../Atlases) for further details).

## Functional MRI data
### `fmri-preprocessing` - Preprocessing of raw functional MRI
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── fmri/
            └── preprocessing/
                ├── <source_T1w>_brainmask.nii.gz
                ├── <source_bold>_motion.tsv
                ├── <source_bold>_space-Ixi549Space_fwhm-8x8x8_preproc.nii.gz
                ├── <source_bold>_space-Ixi549Space_preproc.nii.gz
                ├── <source_bold>_space-meanBOLD_preproc.nii.gz
                └── <source_bold>_space-T1w_preproc.nii.gz
```

The `<source_bold>_motion.tsv` corresponds to the translations (`TransX`, `TransY`, and `TransZ`) in mm and the rotations (`RotX`, `RotY`, and `RotZ`) in each volume as compared to the first one, generated by the Realign tool of **SPM**. Thus, the size of the array should be `Mx6`, `M` being the number of volumes in the original fMRI file.

!!! note
    The naming convention for suffixes tried to follow the [BIDS Extension Proposal 12 (BEP012): BOLD processing derivatives](https://docs.google.com/document/d/16CvBwVMAs0IMhdoKmlmcm3W8254dQmNARo-7HhE-lJU/edit#) as much as possible.


## PET imaging data

### `pet-volume` - Volume-based processing of PET images
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── pet/
            └── preprocessing/
                └── group-<group_label>/
                    ├── <source_file>_space-T1w_pet.nii.gz
                    ├── <source_file>_space-T1w[_pvc-rbv]_pet.nii.gz
                    ├── <source_file>_space-Ixi549Space[_pvc-rbv]_pet.nii.gz
                    ├── <source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<suvr>_pet.nii.gz
                    ├── <source_file>_space-Ixi549Space_brainmask.nii.gz
                    └── <source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<suvr>_mask-brain_pet.nii.gz
```
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── pet/
            └── preprocessing/
                └── atlas_statistics/
                    └── <source_file>_space-<space>[_pvc-rbv]_suvr-<suvr>_statistics.tsv
```

The `_acq-<label>` key/value describes the radiotracer used for the PET acquisition (currently supported: `fdg` and `av45`).

The `[_pvc-rbv]` label is optional, depending on whether your image has undergone partial volume correction (region-based voxel-wise (RBV) method) or not.

The possible values for the `suvr-<suvr>` key/value are: `pons` for FDG-PET and `cerebellumPons` for different types of amyloid PET.

Statistics files (with `_statistics.tsv` suffix) are detailed in [appendix](#appendix-content-of-a-statistic-file).


### `pet-surface` - Surface-based processing of PET images
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── pet/
            └── surface/
                ├── atlas_statistics/
                │   └── <source_file>_task-<label>_acq-<label>_pet_space-<space>_pvc-iy_suvr-<suvr>_statistics.tsv
                ├── <source_file>_hemi-{left|right}_midcorticalsurface
                └── <source_file>_hemi-{left|right}_task-rest_acq-<label>_pet_space-<space>_suvr-<suvr>_pvc-iy_hemi-{left|right}_fwhm-<label>_projection.mgh
```

The `_acq-<label>` key/value describes the radiotracer used for the PET acquisition (currently supported: `fdg` and `av45`).

The `[_pvc-iy]` label describes the partial volume correction used in the algorithm for the projection (Iterative Yang).

The possible values for the `suvr-<suvr>` key/value are: `pons` for FDG-PET and `cerebellumPons` for different types of amyloid PET

The `fwhm` represents the FWHM (in mm) of the Gaussian filter applied to the data mapped onto the FsAverage surface. The different values are 0 (no smoothing), 5, 10, 15, 20, and 25.

Files with the `_midcorticalsurface` suffix represent the surface at equal distance between the white matter/gray matter interface and the pial surface (one per hemisphere).

Files with the `projection` suffix are PET data that can be mapped onto meshes. If `_space_fsaverage_` is in the name, it can be mapped either onto the white or pial surface of FsAverage. If `_space_native_` is in the name, it can be mapped onto the white or pial surface of the subject’s surface (`{l|r}h.white`, `{l|r}h.pial` files from `t1-freesurfer` pipeline).

Files with the `statistics` suffix are text files that display average PET values on either `_space-desikan` or `_space-destrieux` atlases. Example of this content can be found in [appendix](#appendix-content-of-a-statistic-file).


## Statistics

### `statistics-surface` - Surface-based mass-univariate analysis with SurfStat

#### Group comparison
```
groups/
└── group-<group_label>/
    └── statistics/
        ├── participants.tsv
        └── surfstat_group_comparison/
            ├── group-<group_label>_<group_1>-lt-<group_2>_measure-<measure>_fwhm-<label>_correctedPValue.jpg
            ├── group-<group_label>_<group_2>-lt-<group_1>_measure-<measure>_fwhm-<label>_correctedPValue.jpg
            ├── group-<group_label>_<group_1>-lt-<group_2>_measure-<measure>_fwhm-<label>_correctedPValue.mat
            ├── group-<group_label>_<group_2>-lt-<group_1>_measure-<measure>_fwhm-<label>_correctedPValue.mat
            ├── group-<group_label>_participants.tsv
            ├── group-<group_label>_output.log
            └── group-<group_label>_glm.json
```


In the case above, `_correctedPValue` indicates that these are maps of corrected p-values. Other types of maps include, but are not limited to:


| Name                | Suffix               | Description |
|---------------------|----------------------|-------------|
| Corrected p-value   | `_correctedPValue`   | Corrected P-values for vertices and clusters level based on random field theory|
| Uncorrected p-value | `_uncorrectedPValue` | Uncorrected P-value for a  generalized linear model|
| T-statistics        | `_TStatistics`       | T statistics for a generalized linear model|
| FDR                 | `_FDR`               | Q-values for False Discovery Rate of resels|



The `<group_1>-lt-<group_2>` means that the tested hypothesis is: the measurement of `<group_1>` is lower than (`lt`) that `<group_2>`.

The value for `measure` is `ct`, which corresponds to cortical thickness (currently, this is the only measure proposed by Clinica), and the value for `fwhm` corresponds to the size of the surface-based smoothing in mm and can be `5`, `10`, `15` or `20`.

The JPEG files are simple snapshots. The `*.mat` files can be read later by tools like PySurfer and Surfstat.


!!! Example
    ```
    groups/
    └── group-ADvsHC/
        └── statistics/
            ├── participants.tsv
            └── surfstat_group_comparison/
                ├── group-ADvsHC_AD-lt-HC_measure-ct_fwhm-20_correctedPValue.jpg
                ├── group-ADvsHC_participants.tsv
                ├── group-ADvsHC_output.log
                └── group-ADvsHC_glm.json
    ```

    Group comparison between patients with Alzheimer’s Disease (`group_1` = `AD`) and healthy subjects (`group_2` = `HC`). `ADvsHC` defines the `group_label`.

    The `group-ADvsHC_glm.json` contains the information for your generalized linear model, for example:

    ```javascript
    {
        "DesignMatrix": "1 + age + sex + group",
        "StringFormatTSV": "%s %f %f",
        "Contrast": "group",
        "ClusterThreshold": 0.001
    }
    ```

    This file describes the model that you want to create, you should include the factor and covariates in your generalized linear model as a column name in this TSV file. For example, the linear model formula is: `CorticalThickness = 1 + age + sex + group`, the contrasts (factors) `group`, `age` and `sex` are the covariates. Additional information is included in the log file.

    The content of `group-ADvsHC_participants.tsv` is:
    ```
    participant_id   session_id   sex      group   age
    sub-CLNC0001     ses-M00      Female   CN      71.1
    sub-CLNC0002     ses-M00      Male     CN      81.3
    sub-CLNC0003     ses-M00      Male     CN      75.4
    sub-CLNC0004     ses-M00      Female   CN      73.9
    sub-CLNC0005     ses-M00      Female   AD      64.1
    sub-CLNC0006     ses-M00      Male     AD      80.1
    sub-CLNC0007     ses-M00      Male     AD      78.3
    sub-CLNC0008     ses-M00      Female   AD      73.2
    ```

    (Note that to make the display clearer, the rows contain successive tabs, which should not happen in an actual TSV file.)

    The `group-<group_label>` key/value stands for the `group_label` for your analysis. It can be used to run different analyses for different subjects or different analyses for the same subjects.

    The example image here maps statistically significant differences in cortical thickness between a group of patients with Alzheimer’s disease and a group of healthy controls (yellow: correction at the vertex level; blue: correction at the cluster level).

    ![](../img/StatsSurfStat_images/ContrastNegative-CorrectedPValue.jpg)


#### Correlation analysis
```
groups/
└── group-<group_label>/
    └── statistics/
        ├── participants.tsv
        └── surfstat_correlation_analysis/
            ├── group-<group_label>_correlation-<label>_contrast-{negative|positive}_measure-<measure>_fwhm-<label>_correctedPValue.jpg
            ├── group-<group_label>_correlation-<label>_contrast-{negative|positive}_measure-<measure>_fwhm-<label>_correctedPValue.mat
            ├── group-<group_label>_output.log
            ├── group-<group_label>_participants.tsv
            └── group-<group_label>_glm.json
```

The `correlation-<label>` here describes the factor of the model which can be, for example, `age`. The `contrast-{negative|positive}` is the sign of the correlation you want to study, which can be `negative` or `positive`.

All other key/value pairs are defined in the same way as in the previous section.


<!--### Generalised Linear Model (GLM)
```
groups/
└── group-<group_label>/
    └── statistics/
        └── participants.tsv
            ├── surfstat_glm/
            └── group-<group_label>_measure-<measure>_fwhm-<label>_correctedPValue.jpeg
                ├── group-<group_label>_measure-<measure>_fwhm-<label>_correctedPValue.mat
                ├── group-<group_label>_output.log
                ├── group-<group_label>_participants.tsv
                └── group-<group_label>_glm.json
```

This section includes all the other situation for the generalized linear model. It follows the same spirit as the sections above.
-->

## Machine Learning
### `machinelearning-prepare-spatial-svm` - Prepare input data for spatially regularized SVM
```
subjects/
└── sub-<participant_label>/
    └── ses-<session_label>/
        └── machine_learning/
            └── input_spatial_svm/
                └── group-<group_label>/
                    ├── <source_file_t1w>_segm-{graymatter|whitematter|csf}_space-Ixi549Space_modulated-on_spatialregularization.nii.gz
                    └── <source_file_pet>_space-Ixi549Space[_pvc-rbv]_suvr-<suvr>_spatialregularization.nii.gz
```
```
groups/
└── group-<group_label>/
    └── machine_learning/
        └── input_spatial_svm/
            ├── group-<group_label>_space-Ixi549Space_gram.npy
            └── group-<group_label>_space-Ixi549Space_parameters.json
```


At the subject level, it contains SVM regularization of gray matter/white matter/CSF maps or PET data that accounts for the spatial and anatomical structure of neuroimaging data.

At the group level, it contains the Gram matrix with respect to gray matter/white matter/CSF maps needed for the SVM regularization and the information regarding the regularization. An example of JSON file is:


```javascript
{
    "MaxDeltaT": "0.0025",
    "Alpha": "0.0025", // Alpha such that: delta_t = MaxDeltaT * Alpha
    "Epsilon": "10E-6",
    "BoundaryConditions": "TimeInvariant",
    "SigmaLoc": "10",
    "TimeStepMax": "0.07760115580830161",
    "SpatialPrior": "Tissues (GM,WM,CSF)",
    "RegularizationType": "Fisher",
    "FWHM": "4",
}
```


## Appendix - Content of a statistic file
```
<source_file>_space-<space>_map-<map>_statistics.tsv
```

Statistic file for a given [atlas](../../Atlases). The TSV file summarizes regional volumes or averages for a given parametric map. With the help of pandas (Python library), it can be easily parsed for machine learning purposes.

Possible values for `_map-<map>` key/value are:

- For T1: `graymatter` (gray matter), `whitematter` (white matter),  `csf` (CSF) and  `ct` (cortical thickness)

- For DWI: `FA` (fractional anisotropy), `MD` (mean diffusivity, also called apparent diffusion coefficient), `AD` (axial diffusivity), `RD` (radial diffusivity), `NDI` (neurite density index), `ODI` (orientation dispersion index) and `FWF` (free water fraction).

- For PET: `fdg` (<sup>18</sup>F-Fluorodeoxyglucose), `av45` (<sup>18</sup>F-Florbetapir).

!!! Example
    Content of `sub-CLNC01_ses-M00_T1w_space-Hammers_map-graymatter_statistics.tsv`:
    ```
    index   label_name          mean_scalar
    0.0     Background          0.0011357992189
    1.0     Left Hippocampus    0.576250553131
    2.0     Right Hippocampus   0.681780695915
    3.0     Left Amygdala       0.577247679234
    ...
    ```
    (Note that to make the display clearer, the rows contain successive tabs, which should not happen in an actual TSV file.)
