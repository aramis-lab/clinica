<!-- markdownlint-disable MD007 MD046 -->
# `adni-2-bids` – Conversion of the Alzheimer’s Disease Neuroimaging Initiative (ADNI) to BIDS

!!! quote "Description adapted from the [ADNI website](http://adni.loni.usc.edu)"
    ADNI is a global research effort that actively supports the investigation and development of treatments that slow or stop the progression of Alzheimer's disease (AD).
    This multisite, longitudinal study assesses clinical, imaging, genetic and biospecimen biomarkers through the process of normal aging to mild cognitive impairment (MCI) and AD dementia.
    With established, standardized methods for imaging and biomarker collection and analysis, ADNI facilitates a way for scientists to conduct cohesive research and share compatible data with other researchers around the world.

    The ADNI study has four phases: ADNI1, ADNI GO, ADNI2 and ADNI 3.
    New participants were recruited across North America during each phase of the study and agreed to complete a variety of imaging and clinical assessments.
    Participants are followed and reassessed over time to track the pathology of the disease as it progresses.

    | Study characteristics   | ADNI 1 | ADNI GO | ADNI 2 | ADNI 3 |
    | :----------------------:|:------:|:-------:|:------:|:------:|
    | Primary goal            | Develop biomarkers as outcome measures for clinical trials | Examine biomarkers in earlier stages of disease | Develop biomarkers as predictors of cognitive decline, and as outcome measures | Study the use of tau PET and functional imaging techniques in clinical trials |
    | Duration / Start date   | 5 years / October 2004 | 2 years / September 2009 | 5 years / September 2011 | 5 years / September 2016 |
    | Cohort                  | 200 elderly controls + 400 MCI + 200 AD | Existing ADNI-1 + 200 early MCI | Existing ADNI-1 and ADNI-GO + 150 elderly controls + 100 early MCI + 150 late MCI + 150 AD | Existing ADNI-1, ADNI-GO, ADNI-2 + 133 elderly controls + 151 MCI + 87 AD |

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of the **dcm2niix** DICOM to NIfTI converter and **FreeSurfer**.
You can find how to install these software packages on the [installation](../../#installing-clinica-from-source) page.

## Downloading ADNI

To download the ADNI dataset you first need to register to the
[LONI Image & Data Archive (IDA)](https://ida.loni.usc.edu/login.jsp),
a secure research data repository, and then request access to the ADNI dataset through the submission of an
[online application form](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=ADNI).

In order to use the converter, you will need to download both the images and the clinical data.
To do so, from the [main page](https://ida.loni.usc.edu/login.jsp?returnPage=UserManagement.jsp&project=)
click on `PROJECTS` and `ADNI`. To download the imaging data, click on `Download` and choose `Image collections`.

In the `Advanced search` tab:
    1. pick the images you wish to download (for example tick `MRI` to download all the MR images),
    2. tick all the boxes in the `IMAGE TYPES` section (`Original`, `Pre-processed`, `Post-processed`),
    3. click on `SEARCH`.
In the `Advanced search results` tab, click `Select All` and `Add To Collection`.
Finally, in the `Data Collection` tab, select the collection you just created, tick `All` and
click on `Advanced download`.
We advise you to group files as 10 zip files.

To download the clinical data, click on `Download` and choose `Study Data`.
Select all the CSV files which are present in ALL by ticking `Select ALL tabular data` and click `Download`.

!!! note
    You do not have to modify the original folder name or rename the clinical data files before using the converter.

## Modalities supported

Currently, the modalities supported by our converter are:

|                                           | ADNI 1 | ADNI GO/2 | ADNI 3 |
| :----------------------------------------:|:------:|:---------:|:------:|
| T1-weighted MRI                           | ✓      | ✓         | ✓      |
| Diffusion weighted imaging (DWI)          | -      | ✓         | ✓      |
| FLAIR MRI                                 | -      | ✓         | ✓      |
| Functional MRI                            | -      | ✓         | ✓      |
| Fluorodeoxyglucose (FDG) PET              | ✓      | ✓         | ✓      |
| Pittsburgh compound B (PiB) PET (amyloid) | ✓      | -         | -      |
| Florbetapir (AV45) PET (amyloid)          | -      | ✓         | ✓      |
| Florbetaben (FBB) PET (amyloid)           | -      | -         | ✓      |
| Flortaucipir (AV1451) PET (tau)           | -      | -         | ✓      |
| Clinical data                             | ✓      | ✓         | ✓      |

To convert the imaging data to BIDS, a list of subjects with their sessions is first obtained from the `ADNIMERGE` spreadsheet.
This list is compared for each modality of interest to the list of scans available, as provided by modality-specific csv files (e.g. `MPRAGEMETA.csv`).
If the modality was acquired for a specific pair of subject-session, and
several scans and/or preprocessed images are available, only one is converted:

- **T1-weighted MRI** When several scans are available for a single session, the preferred scan (as identified in `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`) is chosen.
If a preferred scan is not specified then the higher quality scan (as defined in `MRIQUALITY.csv`) is selected.
If no quality control is found, then we choose the first scan for the visit.
Gradwarp, B1-inhomogeneity corrected and N3 bias field corrected images are selected.
1.5 T images are preferred for ADNI 1 since they are available for a larger number of patients.
- **DWI** We select images containing 'DTI' in the sequence name, that are not multiband, processed nor enhanced.
- **FLAIR** We select images containing 'FLAIR' in the sequence name, without multiplanar reconstruction (MPR).
- **fMRI** We select images containing 'MRI' in the sequence name, that are not multiband.
- **FDG, Amyloid and Tau PET** The images co-registered and averaged across time frames are selected.

For all imaging modalities, the scans failing quality control (if it was performed) are discarded.

As a final step in the conversion, images from some modalities are centered (currently, T1, FLAIR and the different PET tracers).
For each image, the coordinates of the origin are set to the center of the box containing the image data.
This allows other image processing pipelines in Clinica (mainly SPM based) to run without needing further image preprocessing.

Data that do not change over time, such as the subject's sex, education level or diagnosis at baseline, are obtained from the ADNIMERGE spreadsheet and gathered in the `participants.tsv` file, located at the top of the BIDS folder hierarchy.
The session-dependent data, such as the clinical scores, are obtained from specific CSV files (e.g. `MMSE.csv`) and gathered in `<participant_id>_session.tsv` files in each participant subfolder.
The clinical data being converted are defined in a spreadsheet (`clinical_specifications_adni.xlsx`) that is available with the code of the converter.
The user can easily modify this file if they want to convert additional clinical data.

For further details regarding clinica data, we recommend to look at the [ADNI Data Package
](https://adni.bitbucket.io/index.html) developped by the Alzheimer's Disease Neuroimaging Initiative.

## Using the converter

The converter can be run with the following command line:

```Text
clinica convert adni-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY
```

where:

- `DATASET_DIRECTORY` is the path to the downloaded ADNI images' directory.
Its content looks like:

```text
DATASET_DIRECTORY
├── 027_S_0074
│   ├── 3-plane_localizer
│   │   ├── ...
│   │   └── 2015-02-13_09_52_18.0
│   │       └── S249015
│   ├── ADNI_Brain_PET__Raw
│   │   ├── ...
│   │   └── 2019-01-23_15_54_06.0
│   │       └── I1119527
│   ├── ADNI_Brain_PET__Raw_AV45
│   │   ├── ...
│   │   └── 2015-04-01_16_18_44.0
│   │       └── I481838
│   ├── Axial_DTI
│   │   ├── ...
│   │   └── 2019-01-24_10_35_14.0
│   │       └── S788290
│   ├── ...
├── 041_S_1260
│   ├── ...
├── ...
```

- `CLINICAL_DATA_DIRECTORY` is the path to the directory where the CSV files with the clinical data are located.
Its content looks like:

```text
CLINICAL_DATA_DIRECTORY
├── ADAS_ADNI1.csv
├── ADAS_ADNIGO23.csv
├── ...
├── VISITS.csv
└── VITALS.csv
```

- `BIDS_DIRECTORY` is the path to the output directory, where the BIDS-converted version of ADNI will be stored.

### Optional parameters

The converter offers the possibility of converting only the clinical data (once the images are in the BIDS format) using the optional parameter `-c`.

Due to the high computational time required for converting all the modalities of the whole ADNI dataset, it is possible to convert a single modality at a time using the parameter `-m` with one of the following values:

- `T1` for T1-weighted MRI
- `FLAIR` for FLAIR MRI
- `DWI` for diffusion weighted imaging
- `fMRI` for resting-state functional MRI
- `PET_FDG` for Fluorodeoxyglucose (FDG) PET
- `PET_AMYLOID` for Pittsburgh compound B (PIB), Florbetapir (AV45) and Florbetaben (FBB) PET
- `PET_TAU` for Flortaucipir (AV1451) PET

It is also possible to provide the path to a .txt file with the list of subjects to convert using the optional parameter `--subjects_list`.
This file must contain one subject identifier per line.
The used format for the identifier is the one corresponding to column `PTID` in ADNIMERGE.
For example, we can provide a `subjects.txt` file with the following content:

```text
006_S_4713
006_S_4485
```

For more information about the optional parameters, you can type:

```Text
clinica convert adni-to-bids -h
```

??? failure "Known errors"
    **The following images are not converted since they generate errors.
    They will not be present in the BIDS output.**

    - **T1**
        - _Interslice distance varies in the volume (incompatible with NIfTI format):_
            - Subject sub-ADNI031S0830 for session ses-M48
            - Subject sub-ADNI100S0995 for session ses-M18
            - Subject sub-ADNI031S0867 for session ses-M48
            - Subject sub-ADNI100S0892 for session ses-M18
            - Subject sub-ADNI003S6264 for session ses-M12
        - _Image conversion does not generate an output file:_
            - Subject sub-ADNI029S0845 for session ses-M24
            - Subject sub-ADNI094S1267 for session ses-M24
            - Subject sub-ADNI029S0843 for session ses-M24
            - Subject sub-ADNI027S0307 for session ses-M48
            - Subject sub-ADNI057S1269 for session ses-M24
            - Subject sub-ADNI036S4899 for session ses-M03

    - **DWI**
        - _Interslice distance varies in the volume (incompatible with NIfTI format):_
            - Subject sub-ADNI006S6252 for session ses-M12
            - Subject sub-ADNI007S4611 for session ses-M03
            - Subject sub-ADNI010S0419 for session ses-M156
            - Subject sub-ADNI016S4638 for session ses-M00
            - Subject sub-ADNI021S5237 for session ses-M84
            - Subject sub-ADNI027S2245 for session ses-M120
            - Subject sub-ADNI027S5118 for session ses-M00
            - Subject sub-ADNI094S2238 for session ses-M48
            - Subject sub-ADNI099S6396 for session ses-M24
            - Subject sub-ADNI126S4507 for session ses-M96
            - Subject sub-ADNI126S4891 for session ses-M96
            - Subject sub-ADNI126S4896 for session ses-M96
            - Subject sub-ADNI126S6559 for session ses-M24
            - Subject sub-ADNI126S6724 for session ses-M12
            - Subject sub-ADNI127S0259 for session ses-M156
            - Subject sub-ADNI127S6203 for session ses-M24
            - Subject sub-ADNI127S6330 for session ses-M24
            - Subject sub-ADNI127S6512 for session ses-M24
            - Subject sub-ADNI127S6549 for session ses-M24
            - Subject sub-ADNI129S4287 for session ses-M00
            - Subject sub-ADNI129S6459 for session ses-M24
            - Subject sub-ADNI129S6763 for session ses-M12
            - Subject sub-ADNI129S6784 for session ses-M12
            - Subject sub-ADNI129S6830 for session ses-M00, ses-M12
            - Subject sub-ADNI135S6104 for session ses-M24
            - Subject sub-ADNI135S6446 for session ses-M24
            - Subject sub-ADNI301S6224 for session ses-M36
        - _Two images are generated and we can not choose the correct one:_
            - Subject sub-ADNI098S4018 for session ses-M00
            - Subject sub-ADNI098S4003 for session ses-M12
        - _Wrong b-val/b-vec files:_
            - Subject sub-ADNI029S4585 for session ses-M48
            - Subject sub-ADNI029S2395 for session ses-M60
            - Subject sub-ADNI029S0824 for session ses-M108
            - Subject sub-ADNI029S0914 for session ses-M108
            - Subject sub-ADNI029S4384 for session ses-M48
            - Subject sub-ADNI029S4385 for session ses-M48
            - Subject sub-ADNI094S4630 for session ses-M06
            - Subject sub-ADNI094S4649 for session ses-M06
            - Subject sub-ADNI029S5219 for session ses-M24
            - Subject sub-ADNI020S5203 for session ses-M72
            - Subject sub-ADNI020S6185 for session ses-M24
            - Subject sub-ADNI020S6227 for session ses-M24
            - Subject sub-ADNI020S6449 for session ses-M24
            - Subject sub-ADNI020S6513 for session ses-M12, ses-M24
            - Subject sub-ADNI021S0178 for session ses-M156
            - Subject sub-ADNI153S6336 for session ses-M12
            - Subject sub-ADNI153S6755 for session ses-M00
        - _Volume mismatch between .nii and .bvec / .bval files:_
            - Subject sub-ADNI006S6610 for session ses-M00
            - Subject sub-ADNI006S6682 for session ses-M00
            - Subject sub-ADNI006S6696 for session ses-M00
            - Subject sub-ADNI006S6770 for session ses-M00
            - Subject sub-ADNI027S6183 for session ses-M24
            - Subject sub-ADNI029S6289 for session ses-M00
            - Subject sub-ADNI123S6891 for session ses-M00
            - Subject sub-ADNI130S6043 for session ses-M00
            - Subject sub-ADNI130S6329 for session ses-M00
        - _Wrong output dimensions:_
            - Subject sub-ADNI027S2219 for session ses-M36 (256 x 256 x 2013)
            - Subject sub-ADNI129S2332 for session ses-M12 (256 x 256 x 1549)

    - **FLAIR**
        - _Interslice distance varies in the volume (incompatible with NIfTI format):_
            - Subject sub-ADNI141S0767 for session ses-M84
            - Subject sub-ADNI067S5205 for session ses-M00
            - Subject sub-ADNI127S4928 for session ses-M24
            - Subject sub-ADNI024S4674 for session ses-M06
            - Subject sub-ADNI123S2363 for session ses-M24
            - Subject sub-ADNI053S4578 for session ses-M48
            - Subject sub-ADNI128S4586 for session ses-M48
            - Subject sub-ADNI053S4813 for session ses-M48
            - Subject sub-ADNI053S5272 for session ses-M24
            - Subject sub-ADNI135S6284 for session ses-M12
            - Subject sub-ADNI027S5170 for session ses-M72
            - Subject sub-ADNI068S0127 for session ses-M180
            - Subject sub-ADNI068S2187 for session ses-M120


    - **fMRI**
        - _Image conversion does not generate an output file:_
            - Subject sub-ADNI006S4485 for session ses-M84
            - Subject sub-ADNI123S4127 for session ses-M96
        - _Interslice distance varies in the volume (incompatible with NIfTI format):_
            - Subject sub-ADNI016S6802 for session ses-M00
            - Subject sub-ADNI016S6816 for session ses-M00
            - Subject sub-ADNI126S4891 for session ses-M84

    - **FDG PET**
        - _Image conversion does not generate an output file:_
            - Subject sub-ADNI941S1195 for session ses-M48 folder is empty
            - Subject sub-ADNI005S0223 for session ses-M12 folder is empty
        - _Inconsistent output filename (`NONAME.nii`):_
            - Subject sub-ADNI037S1421 for session ses-M36
            - Subject sub-ADNI037S1078 for session ses-M36

    - **AV45 PET**
        - _Interslice distance varies in the volume (incompatible with NIfTI format):_
            - Subject sub-ADNI128S2220 for session ses-M48

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The ADNI data have been curated and converted to the Brain Imaging Data Structure (BIDS) format
    [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica
    [[Routier et al.](https://hal.inria.fr/hal-02308126/);
    [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

If you used the converter for DWI data, please also cite [[Wen et al., 2020](https://doi.org/10.1007/s12021-020-09469-5)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).

## (Advanced) Appendix - How ADNI is converting into BIDS?

For all the imaging modalities, the `adni-2-bids` converter performs three steps:

1. Image selection from subject and imaging metadata
2. Path extraction
3. Paths to BIDS conversion

??? abstract "List of the CSV files used for the conversion"
    - For all modalities: `ADNIMERGE.csv`
    - T1-weighted MRI: `MPRAGEMETA.csv`, `MRIQUALITY.csv`, `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`
    - Diffusion weighted imaging (DWI): `MRILIST.csv`, `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`
    - FLAIR MRI: `MRILIST.csv`, `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`
    - Functional MRI: `MRILIST.csv`, `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`, `MAYOADIRL_MRI_QUALITY_ADNI3.csv`
    - Fluorodeoxyglucose (FDG) PET: `PET_META_LIST.csv`, `PETQC.csv`, `PETQC3.csv`
    - Pittsburgh compound B (PIB) PET: `PET_META_LIST.csv`, `PIBQC.csv`
    - Florbetapir (AV45) and Florbetaben (FBB) PET: `PET_META_LIST.csv`, `AV45QC.csv`, `AMYQC.csv`
    - Flortaucipir (AV1451) TAU PET: `PET_META_LIST.csv`, `TAUQC.csv`, `TAUQC3.csv`

!!! note "BIDS participant ID: `sub-ADNI<PTID>`"

    As the [BIDS specification](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#participant-names-and-other-labels) only allows letters and numbers in the participant names, underscores were deleted from the original ADNI ID (`PTID` in `ADNIMERGE.csv`).

    Example: if the PTID of a subject is `011_S_0002` the BIDS ID will be `sub-ADNI011S0002`.

### Step 1: Image selection

The following steps are followed for each modality:

- First, we need a list of subjects.
It can be provided by the user or obtained from the `ADNIMERGE` file.
- We want to obtain a list of images containing for each of these subjects one image per visit, if available.
- We need a file containing a list of scans and their metadata for the current modality (we will call it `metadata_file` from now onwards) and a file containing quality control of scans (`qc_file`).

The general principle is, for each subject and each visit, to filter the `metadata_file` according to a certain image selection criterion to obtain the existing scans for the subject at the visit (usually filtering the name of the sequence).
Then the quality of the scans is checked using `qc_file`.
If there are several scans for a visit, we choose the "preferred scan" if indicated, or the scan with better QC score, or the last performed scan if no QC is available, in that order.

!!! Note
    Since we are constantly filtering pandas Dataframes, we have to check regularly if the resulting dataframe is empty (e.g. `if filtered_df.empty:` ...).

!!! Note
    Since there are many missing values in the CSV files we use, we have to check regularly if the values are nan (e.g. `import pandas as pd; if pd.isna(field_value):` ...).

Once a scan has been selected, subject, session, image identifiers and metadata are saved and added to the image list.
Known conversion exceptions are removed from the list.

??? abstract "Criteria for T1-weighted MRI"
    - First, since the `ADNIMERGE` and `MPRAGEMETA` files have different notations for the visits, a correspondence must be established.
    - For each subject, we pair the closest dates from the two files as the same visit (`visits_to_timepoints_t1`)
    - For each visit, the image is processed differently according to which cohort it belongs.
        - If the visit occurs in ADNI 1, GO or 2 (`adni1go2_image`):
            - Filter out images which do not pass QC check based on `MRIQUALITY.csv`.
            - Since there might be several acquisitions (each acquisition corresponds to a different series ID, which is the same for all the processed images associated to the corresponding original acquisition), we have to choose a series.
            - The first option is to look at which one has been further processed (preferred series).
            The last processing step is scaling.
            So we look for a ‘Processed’ image with its sequence name ending by ‘Scaled’ (preferred_processed_scan
            - If not found, we look for images with the previous processing step, N3 bias field correction.
            So we look for a `Processed` image with its sequence name ending by `N3m` (`preferred_processed_scan`)
            - If there is no N3 processed image, we can only use the original image (original_image).
                - We keep images that are MPRAGE (`mprage|mp-rage|mp rage` in the sequence, but no `2`, to keep the first acquisition) or SPGR (`spgr` in sequence name, but no `acc`, since for ADNI 2 we have a single standard acquisition and sometimes a second accelerated acquisition).
            - If there are images with different magnetic field strengths it means the image was acquired during ADNI1, so we filter and keep preferred 1.5 T.
            - If we have several processed images, we keep the one with the lower series ID (acquired first).
            - We save as sequence name of the image we want to use, the current sequence up to the point after where the processed scan contains N3m or N3. We do this because we do not want to keep Scaled images
            - !!! Warning
                The imageID saved in the t1_paths.csv file corresponds to that of the scaled image
        - If the visit occurs in ADNI 3 (`adni3_image`), the selection criteria are (for the current visit and patient):
            - Original images (All the previous preprocessing steps are now done directly on the scanners, so we only need the original image.)
            - Containing `accel` in the sequence name (In the ADNI 3 T1 protocol all the sequences are `accelerated`.)
            - Not containing `_nd` in the sequence name.
              > [The `_ND` suffix is an automated output from some Siemens scanners when distortion correction is applied. The ND stands for No Distortion correction, so they system provides the distortion corrected images and the non corrected images.](http://adni.loni.usc.edu/support/experts-knowledge-base/question/?QID=1191)
            - QC check is done with `select_scan_from_qc`.
        - To check QC in the case of an original image for ADNI1, GO and 2, or always for ADNI 3, the function `select_scan_from_qc` is used:
            - We separate the scans according to the preferred magnetic field strength
                - If there are 3T scans, we check for QC
                    - For the provided scans we try to see if there is a “selected series"
                    - If not, we get the series with the best QC (lower), if QC passed
                    - Otherwise we take the first scan (`select_scan_no_qc`)
                - If there are 1.5T scans, there is no available QC.Then we take the first scan (`select_scan_no_qc`)

??? abstract "Criteria for diffusion weighted imaging (DWI)""
    - The image sequence names in the `MRILIST.csv` file are filtered to keep only DTI that are not multiband, processed or enhanced images (sequence name contains `dti` and does not contain `MB`, `ADC`, `FA`, `TRACEW`, `Enhanced` or `Reg`).
    - For each subject, since the `ADNIMERGE` and `MRILIST` files have different notations for the visits, a correspondence must be established.
    For each subject, we pair the closest dates from the two files as the same visit (`visits_to_timepoints_t1`).
    - For each visit, the images are filtered to keep only images for the current visit.
    - This list of images for the current visit is passed to the function `dwi_image`.
        - The best image is selected (`select_image_qc` in `adni_utils.py`).
          - For the provided scans we try to see if there is a "selected series".
          - If not, we get the series with the higher QC (if QC passed),
          - Otherwise we take the last scan (max(no_qc_ids)).
        - The image metadata is extracted from `MRILIST`.
    - The image data is added to the list of images.

??? abstract "Criteria for FLAIR MRI"
    - The image sequence names in the `MRILIST.csv` file are filtered to keep only FLAIR that are not MPR images (sequence name contains `flair` and does not contain `MPR`).
    - For each subject, since the `ADNIMERGE` and `MRILIST` files have different notations for the visits, a correspondence must be established.
    For each subject, we pair the closest dates from the two files as the same visit (`visits_to_timepoints_t1`).
    - For each visit, the images are filtered to keep only images for the current visit.
    - This list of images for the current visit is passed to the function `flair_image`.
        - The best image is selected (`select_image_qc` in `adni_utils.py`).
            - For the provided scans we try to see if there is a "selected series".
            - If not, we get the series with the higher QC (if QC passed).
            - Otherwise we take the last scan (max(no_qc_ids)).
        - The image metadata is extracted from `MRILIST`.
    - The image data is added to the list of images.

??? abstract "Criteria for functional MRI"
    - QC files are filtered to keep only entries corresponding to fMRI scans. We keep:
        - `MAYOADIRL_MRI_IMAGEQC_12_08_15` rows containing `fMRI` as `series_type`,
        - `MAYOADIRL_MRI_QUALITY_ADNI3` rows containing `EPB` as `SERIES_TYPE`.
        - Resulting entries from both QC files are concatenated in one dataframe.
    - The image sequence names in the `MRILIST.csv` file are filtered to keep only fMRI that are not multiband images (sequence name contains `MRI` and does not contain `MB`).
    - For each subject, since the `ADNIMERGE` and `MRILIST` files have different notations for the visits, a correspondence must be established.
    For each subject, we pair the closest dates from the two files as the same visit (`visits_to_timepoints_t1`).
    - For each visit, the images are filtered to keep only images for the current visit.
    - This list of images for the current visit is passed to the function `fmri_image`.
        - The best image is selected (`select_image_qc` in `adni_utils.py`).
            - For the provided scans we try to see if there is a "selected series".
            - If not, we get the series with the higher QC (if QC passed).
            - Otherwise we take the last scan (max(no_qc_ids)).
        - The image metadata is extracted from MRILIST.
    - The image data is added to the list of images.

??? abstract "Criteria for Fluorodeoxyglucose (FDG) PET"
    - For each subject, the QC files are filtered to keep only entries corresponding to the current subject that pass the QC (PASS == 1) .
        - Resulting entries from both QC files (`PETQC` for ADNI1, GO, 2 and `PETC3` for ADNI3) are concatenated in one dataframe.
    - Then we apply the function `get_images_pet` in `adni_utils.py`:
        - For each unique visit in the list of images passing QC, the images are filtered to keep only images for the current visit.
        - We are looking to get the metadata of the last original image that passed QC.
        If there are several scans for a timepoint we order them so we start by the image_id of the scan acquired last (higher LONIUID).
            - We iterate over QC entries ordered by LONIUID (higher first).
            For each we are looking for Original images, that passed QC (image id corresponding to current QC entry LONIUID), acquired at the same date as the current entry, not containing `early` in the sequence name.
            - We stop the loop at the moment we find the first matching image.
        - With the original image metadata, we can have the series id (it is the same for all the images coming from the same original image and independently of the preprocessing).
        - We see if there are images with the same series id as the original image and the desired preprocessing sequence ("Co-registered, Averaged").
        - If an explicit "Co-registered, Averaged" image does not exist, it means that the original image is already in that preprocessing stage, so we keep the original image.
        - The selected image data is added to the list of images for the subject.
    - The list of images for each subject is added to the list of images to convert.

??? abstract "Criteria for Pittsburgh compound B (PIB) PET"
    - For each subject, the QC file `PIBQC.csv` is filtered to keep only entries corresponding to the current subject that pass the QC (PASS == 1) .
    - Then we apply function the `get_images_pet` in `adni_utils.py`:
        - Functioning is the same as described above for FDG PET but we look for a different sequence ("PIB Co-registered, Averaged").
    - The list of images for each subject is added to the list of images to convert.

??? abstract "Criteria for Florbetapir (AV45) PET and Florbetaben (FBB) PET"
    - For each subject, the QC files are filtered to keep only entries corresponding to the current subject that pass the QC (PASS == 1) .
        - Resulting entries from both QC files (`AV45QC.csv` for AV45 for ADNI1, GO, 2 and `AMYQC.csv` for ADNI3) are concatenated in one dataframe.
    - Then we apply the function `get_images_pet` in `adni_utils.py`:
        - Functioning is the same as described above for FDG PET but we look for different sequences ("AV45 Co-registered, Averaged" and "FBB Co-registered, Averaged").
        - Inside the function the tracer is determined based on the sequence name.
    - The list of images for each subject is added to the list of images to convert.

??? abstract "Criteria for Flortaucipir (AV1451) TAU PET"
    - For each subject, the QC files are filtered to keep only entries corresponding to the current subject that pass the QC (PASS == 1) .
        - Resulting entries from both QC files (`TAUQC.csv` for ADNI1, GO, 2 and `TAUQC3.csv` for ADNI3) are concatenated in one dataframe.
    - Then we apply the function `get_images_pet` in `adni_utils.py`:
        - Functioning is the same as described above for FDG PET but we look for a different sequence ("AV1451 Co-registered, Averaged").
    - The list of images for each subject is added to the list of images to convert.

### Step 2: Path extraction

In this step the input is a pandas dataframe of images containing metadata.
Usually: `Subject_ID`, `VISCODE`, `Visit`, `Sequence`, `Scan_Date`, `Study_ID`, `Series_ID`, `Image_ID`.
This dataframe is passed to the function `find_image_path(images, source_dir, modality, prefix, id_field)`
in the `adni_utils.py` file.
The source directory containing the downloaded ADNI images is used as the base path for the paths we will create next.

For each image in the dataframe, we create the path to the folder of the subject and the corresponding image sequence (after escaping special characters).
Inside these folders, there are folders named as timestamps corresponding to the dates of the different scans with the same sequence name for the current subject.
Inside each of these "timestamp" folders there is another folder with the corresponding series ID or image ID (depending on the modality).
Since the exact timestamp of a scan is hard to obtain, we need to iterate through the "timestamp" folders to see which one contains the folder named as the identifier we are looking for.
Once we have this folder we can see if there is an image inside it and if it is a NIfTI file or a series of DICOM files.
If the image is in NIfTI format, the path to the file is saved.
If it is a DICOM image, the path to the folder is saved.
If at some point we do not find a corresponding folder for the subject, the sequence, or the image or series identifier, we save an empty path.

The result of this step is a `MODALITY_paths.tsv` file (e.g. `t1_paths.tsv`) containing the list of image metadata, paths and whether the images are DICOM or NIfTI.
It will be located in the BIDS output folder in a directory called `conversion_info`.

!!! note "Use of symlinks"
    If a symlink is used for the folder of the imaging data, the `MODALITY_paths.tsv` files will include
    the paths using the symlink (and not the real path).

### Step 3: Paths to BIDS conversion

In this step the images in the dataframe are going to be converted into NIfTI images stored in a BIDS tree hierarchy.
Most of the modalities use the function `paths_to_bids(images, bids_dir, modality)` that creates multiple parallel processes to accelerate the conversion of the images.
The main steps are, for each image (inside the function `create_file`):

- Creation of the appropriate filename and output path for the image.

- Creation of the corresponding folder in the BIDS hierarchy to store the image with the correct BIDS name.

- If the image is DICOM:
    - Check if there are several DICOM images in the path and creation of a new temporary folder containing only the image we want to convert;
    - Conversion using `dcm2niix` (more recent converter);
    - Check if all the needed files are generated, usually a NIfTI file (there might be JSON files depending on the image modality);
    - If conversion failed try `dcm2nii` (older version of converter);
    - If failed (no NIfTI image) we move on to the next image.

- Now that we have a NIfTI image (either in the ADNI folder or a newly converted one), we will "center" the image if the modality requires it (currently done for T1, FLAIR and PETs).
We will use NiBabel to load the image, find the closest canonical representation, and then, we will set the coordinates of the origin, to the center of the box containing the image data (`center_nifti_origin(nifti_file, output_image)`).
