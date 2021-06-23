<!-- markdownlint-disable MD046 -->
# `oasis-2-bids` – Conversion of the Open Access Series of Imaging Studies (OASIS) to BIDS

!!! quote "Description reproduced from the [OASIS' webpage](http://oasis-brains.org/)"
    The Open Access Series of Imaging Studies (OASIS) is a project aimed at making MRI data sets of the brain freely available to the scientific community.
    By compiling and freely distributing MRI data sets, we hope to facilitate future discoveries in basic and clinical neuroscience.
    OASIS is made available by the Washington University Alzheimer’s Disease Research Center, Dr. Randy Buckner at the Howard Hughes Medical Institute (HHMI) at Harvard University, the Neuroinformatics Research Group (NRG) at Washington University School of Medicine, and the Biomedical Informatics Research Network (BIRN).

    The "Cross-sectional MRI Data in Young, Middle Aged, Nondemented and Demented Older Adults" set consists of a cross-sectional collection of 416 subjects aged 18 to 96.
    For each subject, 3 or 4 individual T1-weighted MRI scans obtained in single scan sessions are included.
    The subjects are all right-handed and include both men and women.
    100 of the included subjects over the age of 60 have been clinically diagnosed with very mild to moderate Alzheimer’s disease (AD).
    Additionally, a reliability data set is included containing 20 nondemented subjects imaged on a subsequent visit within 90 days of their initial session.

  For more information about the images and the dataset you can read the [OASIS Fact Sheet](http://www.oasis-brains.org/pdf/oasis_cross-sectional_facts.pdf).

## Dependencies

If you installed the core of Clinica, this converter needs no further dependencies.

## Downloading OASIS

The OASIS to BIDS converter requires the user to have downloaded the OASIS-1 (also called *Cross-sectional MRI Data in Young, Middle Aged, Nondemented and Demented Older Adults*) imaging and clinical data. To do so, visit the [OASIS website](http://www.oasis-brains.org/), click on `DATASETS` then `OASIS-1`. You can download the imaging data via FTP or XNAT (if you click on `Browse Data`). Download the clinical data by clicking on `Download CSV`.

!!! note
    You do not have to modify the original folder name or rename the clinical data file before using the converter.

!!! warning
    We do not currently support the conversion of OASIS-2 and OASIS-3.

## Modalities supported

Currently, the converter converts only the T1-weighted MRI images and clinical data.

To convert the imaging data to BIDS, the list of subjects is obtained from the downloaded folders.
For each subject, among the multiple T1w MR images available, we select the average of the motion-corrected co-registered individual images resampled to 1 mm isotropic voxels, located in the `SUBJ_111` subfolder.
The clinical data being converted are defined in a spreadsheet (`clinical_specifications.xlsx`) available with the code of the converter, which the user can modify.

## Using the converter

The converter can be run with the following command line:

```Text
clinica convert oasis-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original OASIS images' directory.
Its content looks like:

```text
DATASET_DIRECTORY
├── OAS1_0001_MR1
│   ├── FSL_SEG
│   ├── PROCESSED
│   │   └── MPRAGE
│   │       ├── SUBJ_111
│   │       └── T88_111
│   │           └── t4_files
│   └── RAW
├── ...
```

- `CLINICAL_DATA_DIRECTORY` is the path to the directory containing the CSV file called `oasis_cross-sectional.csv`.

- `BIDS_DIRECTORY` is the path to the output directory, where the BIDS-converted version of OASIS will be stored.

!!! note
    In order to improve the readability, the BIDS subject ID is different from the original OASIS ID and is defined as follows:

    ```Text
    sub-OASIS1+ original numerical ID of the subject
    ```

    !!! example
        If the original subject ID is `OAS1_0001_MR1`, since the numerical ID is `0001`, the final BIDS ID will be `sub-OASIS10001`.

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The OASIS data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica [[Routier et al.](https://hal.inria.fr/hal-02308126/); [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
