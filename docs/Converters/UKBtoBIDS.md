<!-- markdownlint-disable MD046 -->
# `ukb-to-bids` – Conversion of the UK Biobank (UKB) to BIDS

!!! quote "Description reproduced from the [UK Biobank webpage](https://www.ukbiobank.ac.uk/)"
    UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and health information from half a million UK participants. The database is regularly augmented with additional data and is globally accessible to approved researchers undertaking vital research into the most common and life-threatening diseases. It is a major contributor to the advancement of modern medicine and treatment and has enabled several scientific discoveries that improve human health

## Dependencies

If you installed the core of Clinica, this converter needs no further dependencies.

## Downloading UK Biobank

The UK Biobank to BIDS converter requires the user to have downloaded the UK Biobank imaging and clinical data.

Once you have access to it, the imaging data can be downloaded and extracted using this tutorial: https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf, and the clinical data can be extracted using the information from this repository https://github.com/kenhanscombe/ukbtools.


Be careful about two things: 
    - the converter only takes in the nifti provided by UK Biobank, not the dicoms.
    - the converter will not work if you do not have among the clinical data the year of birth, the age at recruitment and at sessions and the sex of the subject.

## Supported modalities

Please note that this converter processes the following modalities : 
- T1W
- T2 Flair
- DWI
- SWI
- tfMRI
- rsfMRI

!!! Chosen files
    Whenever possible, we use the rawest files in order not to impose opinions onto the user. When available we get the associated json.
    - T1W -> T1.nii.gz: cropped so there is no neck and defaced. The rawest image would be the simply defaced one, but since we have no interest in the neck we don't use it. We do not want any other corrections.
    - T2 Flair -> T2_FLAIR: cropped so there is no neck and defaced. The reasons for choosing this image are the same as for T1.
    - DWI -> AP.nii.gz/PA.nii.gz + bvals and bvecs. It is the rawest images we can get.
    - rsfMRI -> rfMRI.nii.gz
    - tfMRI -> tfMRI.nii.gz or dicom equivalent to get the json associated
    - SWI -> SWI.nii.gz which is the combined coil version. We would get a rawer version, but dcm2niix seems to have trouble handling this format, so it is simpler to go with this version. Plus, swi is not fully integrated to bids so changes may be coming, hence a simple version for now.
    


## Using the converter

The converter can be run with the following command line:

```Text
clinica convert ukb-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original UK BIobank imaging directory, which content should look like:

    ```text
    DATASET_DIRECTORY
    ├── 1000223_20227_2_0.zip
    ├── 1000223_20249_2_0.zip
    ├── 1000223_20250_2_0.zip
    ├── 1000223_20251_2_0.zip
    ├── 1000223_20252_2_0.zip
    ├── 1000223_20253_2_0.zip
    ├── 3338337_20251_2_0.zip
    └── 5566112_20253_2_0.zip
    ```

- `CLINICAL_DATA_DIRECTORY` is the path to the directory containing the clinical CSV file.

- `BIDS_DIRECTORY` is the path to the output directory where the BIDS-converted version of UK Biobank will be stored.

!!! note
    In order to improve the readability, the BIDS subject ID is different from the original UK Biobank ID and is defined as follows:

    ```Text
    sub-UKB+ original numerical ID of the subject
    ```

    !!! example
        If the original subject ID is `0001`, the final BIDS ID will be `sub-UKB0001`.

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The UK Biobank data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica [[Routier et al.](https://hal.inria.fr/hal-02308126/); [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
