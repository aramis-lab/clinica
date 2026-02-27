<!-- markdownlint-disable MD046 -->
# `genfi-to-bids` – Conversion of the GENFI to BIDS

!!! quote "Description reproduced from the [GENFI webpage](https://www.genfi.org)"
    The Genetic Frontotemporal dementia Initiative (GENFI) is a group of research centres across Europe and Canada with expertise in familial FTD, and is co-ordinated by Professor Jonathan Rohrer at University College London. GENFI is the largest genetic FTD consortium to date and currently consists of sites across the UK, Netherlands, Belgium, France, Spain, Portugal, Italy, Germany, Sweden, Denmark, Finland and Canada. The aim of the study is to understand more about genetic FTD, particularly in those who have mutations in the progranulin (GRN), microtubule-associated protein tau (MAPT) and chromosome 9 open reading frame 72 (C9orf72) genes. GENFI investigates both people who have developed symptoms and also people who have a risk of developing symptoms in the future because they carry an abnormal genetic mutation. By studying these individuals who are destined to develop the disease later in life we can understand the development from the very earliest changes. The key objectives of GENFI are therefore to develop markers which help identify the disease at its earliest stage as well as markers that allow the progression of the disease to be tracked. We are now collaborating closely with other similar studies around the world through the FTD Prevention Initiative. Through this worldwide initiative we are working with pharmaceutical companies to help design clinical trials for genetic FTD.

## Dependencies

If you only [installed the core of Clinica](../Software/Installation.md), this pipeline needs the installation of the [**dcm2niix**](../Software/Third-party.md#dcm2nix) DICOM to NIfTI converter.

## Supported modalities

Please note that this converter currently processes the following modalities : 

- T1W
- T2W
- DWI
- Fieldmaps
- rsfMRI

## Downloading GENFI

To download GENFI in a way that you can convert it, you need to make sure of two things: 
- Only select the "simplify downloaded archive structure" in download data options.
- Select only "DICOM" in "Select Image Data" section "Scans Format".

![](../img/GENFI_download/selection.png)

## Using the converter

The converter can be run with the following command line:

```Text
clinica convert genfi-to-bids [OPTIONS] DATASET_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original GENFI imaging directory, whose content should look like:

    ```text
    DATASET_DIRECTORY
    ├── C9ORF001-01-MR00
    │   ├── 1
    │   ├── 11
    │   ├── 12
    │   └── 13
    ├── C9ORF001-11
    │   ├── 1
    │   ├── 2
    │   ├── 4
    │   ├── 5
    │   ├── 6
    │   ├── 7
    │   ├── 8
    │   └── 9
    └── GRN001-01-MR00
        ├── 1
        ├── 10
        ├── 8
        └── 9
    ```

- `BIDS_DIRECTORY` is the path to the output directory where the BIDS-converted version of GENFI will be stored, whose content should look like:

    ```text
    BIDS_DIRECTORY
    ├── dataset_description.json
    ├── participants.tsv
    ├── README
    ├── sub-C9ORF001
    │   ├── ses-01
    │   │   ├── anat
    │   │   │   ├── sub-C9ORF001_ses-01_run-01_T1w.json
    │   │   │   ├── sub-C9ORF001_ses-01_run-01_T1w.nii.gz
    │   │   │   ├── ...
    │   │   ├── dwi
    │   │   │   ├── sub-C9ORF001_ses-01_run-01_dwi.bval
    │   │   │   ├── sub-C9ORF001_ses-01_run-01_dwi.bvec
    │   │   │   ├── ...
    │   │   ├── func
    │   │   │   ├── sub-C9ORF001_ses-01_task-rest_run-01_bold.json
    │   │   │   └── sub-C9ORF001_ses-01_task-rest_run-01_bold.nii.gz
    │   │   └── sub-C9ORF001_ses-01_scans.tsv
    │   ├── ses-11
    │   │   ├── ...
    │   └── sub-C9ORF001_sessions.tsv
    ├── sub-GRN001
    │   ├── ses-01
    │   │   ├── ...
    │   └── sub-GRN001_sessions.tsv
    └── ...
    ```

    Sessions respect the naming convention used in GENFI clinical data. In particular, for `ses-XY`, `X` represents the GENFI phase, and `Y` represents the visit number.

- `OPTIONS`:
    - `--clinical-data-dir/-cdd` is the path to the clinical data directory. Allows the user to add mandatory clinical data to `participants.tsv` and `sessions.tsv`. Mandatory clinical data are the following :
        - For `participants.tsv` : `blinded_code`, `blinded_family`, `blinded_site`, `gender`.
        - For `sessions.tsv` : `age_at_visit`, `date_of_assessment`, `diagnosis`, `education`, `ftld-cdr-nm-global`, `genetic_group`, `genetic_status_1`, `genetic_status_2`, `visit`. 
        <br> All the remaining clinical data is optional and added through the use of the `-full` flag.
    - `--clinical-data-txt/-cdt` is a txt file containing the additional fields the user wants. The available data can be retrieved from the specification file `full_specs.csv`, located in the Clinica installation directory at : `your_path_to_clinica/clinica/converters/genfi_to_bids/specifications/full_specs.csv`. The txt file should be written one field per line such as in the example below :
        ```text
        diagnosis_1
        diagnosis_1.1
        diagnosis_2
        diagnosis_3
        diagnosis_4
        diagnosis_5
        diagnosis_child1
        diagnosis_child2
        digit_symbol
        disinhibition
        dob
        drc_qc
        drug_history
        ...
        ```
    If the `-full` flag is used, this option is considered redundant and will be ignored.
    - `-gif` allows the user to add all the clinical data related to the imaging volumes (GIF, Geodesic Information Flow) to `session.tsv`. The added clinical data also contain the mandatory ones.
    - `-full` allows the user to add all clinical data (mandatory, GIF, and the remaining optional ones) to `sessions.tsv`.



!!! note
    In order to improve the readability, the BIDS subject ID is the genetic group concatenated with the original GENFI ID and is defined as follows:

    ```Text
    sub-GRN/C9ORF/MAPT+ original numerical ID of the subject
    ```

    !!! example
        If the original subject ID is `0001`, the final BIDS ID will be `sub-GRN0001`.

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The GENFI data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica [[Routier et al.](https://hal.inria.fr/hal-02308126/); [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
