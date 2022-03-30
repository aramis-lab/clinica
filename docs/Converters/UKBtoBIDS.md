<!-- markdownlint-disable MD046 -->
# `ukb-to-bids` – Conversion of the UK Biobank (UKB) to BIDS

!!! quote "Description reproduced from the [UK Biobank webpage](https://www.ukbiobank.ac.uk/)"
    UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and health information from half a million UK participants. The database is regularly augmented with additional data and is globally accessible to approved researchers undertaking vital research into the most common and life-threatening diseases. It is a major contributor to the advancement of modern medicine and treatment and has enabled several scientific discoveries that improve human health.

## Dependencies

If you installed the core of Clinica, this converter needs no further dependencies.

## Downloading UK Biobank

The UK Biobank to BIDS converter requires the user to have applied to the data and to have downloaded the UK Biobank imaging and clinical data.

Once you have access to it, the imaging data can be downloaded and extracted using this tutorial: https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf, and the clinical data can be extracted using the information from this repository https://github.com/kenhanscombe/ukbtools.

Be careful, the converter only takes in the nifti provided by UK Biobank, not the dicoms.
Also, the clinical_data.tsv needs to contain these informations about the subject: sex, year of birth, age at recruitment, age at sessions.

## Supported modalities

Please note that this converter processes the following modalities : 
- T1W
- T2 Flair
- DWI
- SWI
- tfMRI
- rsfMRI

!!! Chosen files
    Whenever possible, we use the rawest files so as to leave it to the user to choose the processing needed. When available we get the associated json.
| Modality    | Chosen image(s) | Justification |
| :----------:|:---------------:|:-------:|
| T1W                     | T1.nii.gz       | Defaced and cropped so there is no neck. The rawest image would be the simply defaced one, but since we have no interest in the neck we don't use it. We do not want any other corrections. |
| T2 Flair    | T2_FLAIR.nii.gz | Defaced and cropped so there is no neck. The reasons for choosing this image are the same as for T1. |
| DWI         | AP/PAnii.gz     | It is the rawest images we can get. bvals and bvec are also available. |
| rsfMRI      | rsfMRI.nii.gz   |  |
| tfMRI       | tfMRI.nii.gz    |  |
| SWI         | SWI.nii.gz      | Combined coil version. We would get a rawer version, but dcm2niix has trouble handling the slices direction, so it is simpler to go with this version. Plus, swi is not fully integrated to bids so changes may be coming, hence a simple version for now. | 
    

## Using the converter

The converter can be run with the following command line:

```Text
clinica convert ukb-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original UK BIobank imaging directory, whose content should look like:

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
