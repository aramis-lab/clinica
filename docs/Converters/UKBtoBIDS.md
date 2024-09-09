<!-- markdownlint-disable MD046 -->
# `ukb-to-bids` – Conversion of the UK Biobank (UKB) to BIDS

??? quote "Description reproduced from the [UK Biobank webpage](https://www.ukbiobank.ac.uk/)"
    UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and health information from half a million UK participants. The database is regularly augmented with additional data and is globally accessible to approved researchers undertaking vital research into the most common and life-threatening diseases. It is a major contributor to the advancement of modern medicine it and has led to the discovery of several scientific advances and numerous treatments to improve human health.

## Downloading the Data

The UKB to BIDS converter assumes that the user has already got access to the data and has downloaded both imaging and clinical data locally.

!!! danger "Organising data with the aim of using the converter"
    Be careful, the file {==clinical_data.csv==}  needs to contain the following information about the subject :
    *sex, year of birth, age at recruitment, age at sessions*.
    There are two columns as there are two imaging sessions.
    The names of the columns should be kept as they are when downloaded using `ukbconv` using options `csv`.



## Using the converter
### Dependencies

If you [installed the core of Clinica](../Installation.md#install-clinica), this converter needs the [dcm2niix]((../Third-party.md#converters)) package.

### Supported modalities

<div class="annotate" markdown>
Please note that this converter processes the following modalities (1) :
</div>

1. Whenever possible, the converter uses the rawest files available. This decision allows the user to choose the processing they need. If possible the converter also gets the associated json.

| Modality    | Chosen image | Format |                                                                                                                                                              Justification                                                                                                                                                               |
| :----------:|:---------------:|:-------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| T1W         | T1.nii.gz       | nifti |                                                                  Defaced and cropped so there is no neck. The rawest image would be the simply defaced one, but brain studies usually are not interested in the region of the neck. No other corrections are included.                                                                   |
| T2 Flair    | T2_FLAIR.nii.gz | nifti |                                                                                                                                           Defaced and cropped so there is no neck. Same as T1.                                                                                                                                           |
| DWI         | AP/PA.nii.gz     | nifti |                                                                                                                                          Rawest existing images with available bvals and bvec.                                                                                                                                           |
| rsfMRI      | rsfMRI.dcm | dicom |                                                                                                               The nifti doesn't always include a json and in consideration for its usage, the dicom is converted instead.                                                                                                                |
| tfMRI       | tfMRI.dcm   | dicom |                                                                                                                                                             Same as rsfMRI.                                                                                                                                                              |
| SWI         | SWI.nii.gz      | nifti | Combined coil version. We would get a rawer version, but `dcm2niix` has trouble handling the slices direction, so it is simpler to go with this version. In addition, SWI modality is not fully integrated to BIDS specification and some changes may be coming. A simple version is available in the current version of this converter. |


### Understanding the command line

The converter can be run with the following command line:

```{ .bash .copy }
clinica convert ukb-to-bids DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY [OPTIONS]
```

where:

<div class="grid" markdown>

=== "Imaging data :"

    - `DATASET_DIRECTORY` is the path to the
    original UK Biobank imaging directory.

    - `BIDS_DIRECTORY` is the path to the
    output directory where the BIDS-converted
    version of UK Biobank will be stored.


```title="DATASET_DIRECTORY Organisation"
    DATASET_DIRECTORY
    ├── 1000223_20227_2_0.zip
    ├── 1000223_20249_2_0.zip
    ├── 1000223_20250_2_0.zip
    ├── 3338337_20251_2_0.zip
    └── 5566112_20253_2_0.zip
```

=== "Clinical data :"

    - `CLINICAL_DATA_DIRECTORY` is the path to the directory
    containing the clinical CSV file.

```title="CLINICAL_DATA_DIRECTORY Organisation"
    CLINICAL_DATA_DIRECTORY
    ├── clinical_data.csv
    ├── ...
```

</div>



!!! note
    In order to improve the readability, the BIDS subject ID is different from the original UK Biobank ID and is defined as follows:

    ```Text
    sub-UKB+ original numerical ID of the subject
    ```

    !!! example
        If the original subject ID is `0001`, the final BIDS ID will be `sub-UKB0001`.


--8<-- "snippets/converters_options.md"


## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The UK Biobank data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica [[Routier et al.](https://hal.inria.fr/hal-02308126/); [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
