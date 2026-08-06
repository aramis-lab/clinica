<!-- markdownlint-disable MD046 -->
# `oasis3-to-bids` – Conversion of the Cross-sectional MRLongitudinal MRI Data (OASIS-2) to BIDS

??? quote "Description reproduced from the [OASIS' webpage](https://sites.wustl.edu/oasisbrains/)"
    The Open Access Series of Imaging Studies (OASIS) is a project aimed at making MRI data sets of the brain freely available to the scientific community.
    By compiling and freely distributing MRI data sets, we hope to facilitate future discoveries in basic and clinical neuroscience.
    OASIS is made available by the Washington University Alzheimer’s Disease Research Center, Dr. Randy Buckner at the Howard Hughes Medical Institute (HHMI) at Harvard University, the Neuroinformatics Research Group (NRG) at Washington University School of Medicine, and the Biomedical Informatics Research Network (BIRN).

    This set consists of a longitudinal collection of MPRAGE T1w for 150 subjects aged 60 to 96. Each subject was scanned on two or more visits,
    separated by at least one year for a total of 373 imaging sessions. For each subject, 3 or 4 individual T1-weighted MRI scans
    obtained in single scan sessions are included. The subjects are all right-handed and include both men and women.
    72 of the subjects were characterized as nondemented throughout the study. 64 of the included subjects were characterized as demented
    at the time of their initial visits and remained so for subsequent scans, including 51 individuals with mild to moderate Alzheimer’s disease.
    Another 14 subjects were characterized as nondemented at the time of their initial visit and were subsequently characterized as demented at a later visit.

## Dependencies

If you installed the core of Clinica, this converter needs no further dependencies.

## Downloading OASIS-2

The OASIS2 to BIDS converter requires the user to have downloaded the OASIS-2 (also called *Cross-sectional MRLongitudinal MRI Data in
Nondemented and Demented Older AdultsI Data in Young, Middle Aged, Nondemented and Demented Older Adults*) imaging and clinical data.

To do so, visit the OASIS2 website page [specific to OASIS-2](https://sites.wustl.edu/oasisbrains/home/oasis-2/), and go to `Download Instructions`.

- Download demographic data by clicking on `OASIS-2 Longitudinal Subject Data`, then `Demographic Data`.
- Download image data by clicking on `OASIS-2 Longitudinal Scan Data` then any archive link.

!!! warning
    All subjects to be converted should be placed in the same folder, which means the archives downloaded from the website should be unzipped in the same folder.

   
## Using the converter

The converter can be run with the following command line:

```Text
clinica convert oasis2-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original OASIS-2 imaging directory, which content should look like:

    ```text
    DATASET_DIRECTORY
    .
    ├── OAS2_0100_MR1
    │   └── RAW
    │       ├── ...
    │       ├── mpr-4.nifti.hdr
    │       └── mpr-4.nifti.img
    ├── OAS2_0100_MR2
    │   ├── OLD
    │   └── RAW
    │       ├── ...
    │       ├── mpr-3.nifti.hdr
    │       └── mpr-3.nifti.img
    ├── OAS2_0100_MR3
    │   ├── OLD
    │   └── RAW
    │       ├── ...
    │       ├── mpr-3.nifti.hdr
    │       └── mpr-3.nifti.img
    ├── OAS2_0101_MR1
    │   └── RAW
    │       ├── ...
    │       ├── mpr-4.nifti.hdr
    │       └── mpr-4.nifti.img
    ├── OAS2_0101_MR2
    │   ├── OLD
    │   └── RAW
    │       ├── ...
    │       ├── mpr-4.nifti.hdr
    │       └── mpr-4.nifti.img

    ```

- `CLINICAL_DATA_DIRECTORY` is the path to the directory containing the clinical excel file.

- `BIDS_DIRECTORY` is the path to the output directory where the BIDS-converted version of OASIS-2 will be stored.


--8<-- "snippets/converters_options.md"

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The OASIS-2 data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica [[Routier et al.](https://hal.inria.fr/hal-02308126/); [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
