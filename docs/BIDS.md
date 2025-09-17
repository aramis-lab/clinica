<!-- markdownlint-disable MD007 -->
# BIDS

## Introduction

[**BIDS**](http://bids.neuroimaging.io/) (Brain Imaging Data Structure) is the standard adopted for the organisation of the datasets used by Clinica pipelines through the command line.
BIDS provides a unified structure for organising and describing neuroimaging and behavioural data.
The use of a standard like BIDS makes easier developing and distributing code that uses neuroimaging datasets.
For this reason, when using Clinica pipelines from the command line, the input format of the dataset is required to be BIDS-compliant.

## An overview of the BIDS structure

Here is a general overview of the BIDS structure.
If you need more details, please check the [documentation](https://bids-specification.readthedocs.io/en/latest/) on the [website](http://bids.neuroimaging.io/).

```Text
BIDS_Dataset/
├── dataset_description.json
├── sub-CLNC01/
│   │   ├── ses-M000/
│   │   │   ├── anat/
│   │   │   │   ├── sub-CLNC01_ses-M000_T1w.json
│   │   │   │   └── sub-CLNC01_ses-M000_T1w.nii.gz
│   │   │   ├── dwi/
│   │   │   │   ├── sub-CLNC01_ses-M000_dwi.bval
│   │   │   │   ├── sub-CLNC01_ses-M000_dwi.bvec
│   │   │   │   ├── sub-CLNC01_ses-M000_dwi.json
│   │   │   │   ├── sub-CLNC01_ses-M000_dwi.nii.gz
│   │   │   │   └── ...
│   │   │   ├── fmap/
│   │   │   │   ├── sub-CLNC01_ses-M000_phasediff.json
│   │   │   │   ├── sub-CLNC01_ses-M000_phasediff.nii.gz
│   │   │   │   ├── sub-CLNC01_ses-M000_magnitude1.nii.gz
│   │   │   │   └── sub-CLNC01_ses-M000_magnitude2.nii.gz
│   │   │   ├── func/
│   │   │   │   ├── sub-CLNC01_ses-M000_task­-rest_bold.json
│   │   │   │   ├── sub-CLNC01_ses-M000_task­-rest_bold.nii.gz
│   │   │   │   └── ...
│   │   │   ├── pet/
│   │   │   │   ├── sub-CLNC01_ses-M000_trc-11CPIB_pet.json
│   │   │   │   ├── sub-CLNC01_ses-M000_trc-11CPIB_pet.nii.gz
│   │   │   │   ├── sub-CLNC01_ses-M000_trc-18FFDG_pet.json
│   │   │   │   ├── sub-CLNC01_ses-M000_trc-18FFDG_pet.nii.gz
│   │   │   └── sub-CLNC01_ses-M000_scans.tsv
│   │   ├── ses-M018/
│   │   │   └── ...
│   │   └── sub-CLNC01_sessions.tsv
├── sub-CLNC02/
│   └── ...
├── ...
└── participants.tsv
```

!!! note "PET modality"
    Clinica early adopted [BIDS Extension Proposal (BEP) 9](https://docs.google.com/document/d/1mqMLnxVdLwZjDd4ZiWFqjEAmOmfcModA_R535v3eQs0/edit) regarding PET modality before its official adoption in BIDS version 1.6.0.

    Since version 0.6, Clinica is now compliant with the official specifications for PET modality.
    This involves three main changes:

    - `task-rest` key/entity is now optional and ignored by Clinica
    - PET acquisitions with different tracers now use the new tracer entity `trc` instead of
    BIDS acquisition entity `acq`.
    - PET files must use tracer entity `trc` (optional in BIDS but necessary for Clinica pipelines) and
    label names follow the proposed convention in BIDS
    (pib -> 11CPIB, av45 -> 18FAV45, fbb -> 18FFBB, fdg -> 18FFDG, flute -> 18FFMM, tau -> 18FAV1451).


## The dataset_description.json file

### Specifications

This file MUST be present at the root of a BIDS or CAPS dataset and MUST contain the following minimal information:

```json
{
	"Name": "name identifier for the dataset",
	"BIDSVersion": "1.7.0",
	"CAPSVersion": "1.0.0",
	"DatasetType": "derivative"
}
```

- `Name`: String identifier of the dataset. It can be the name of your study for example. By default Clinica generates a random UUID for this field. When running a pipeline which will create a new CAPS dataset, you can use the `--caps-name` option to provide a name (see [common options](/Software/InteractingWithClinica.md#common-options-for-pipelines)). If the CAPS dataset already exist, the existing name will be kept.
- `BIDSVersion`: The version number of the BIDS specifications that the BIDS input dataset is using when this CAPS dataset was generated.
- `CAPSVersion`: The version number of the CAPS specifications used for this dataset.
- `DatasetType`: Either "raw" or "derivative". For a CAPS dataset this should always be "derivative" as it contains processed data.

In addition, the `dataset_description.json` file MAY contain a `Processing` key which is a list of objects describing the different processing pipelines that were run on this CAPS.
Here is an example for a CAPS dataset containing the outputs of two pipelines: `t1-linear` and `pet-linear`:

```json
{
    "Name": "e6719ef6-2411-4ad2-8abd-da1fd8fbdf32",
    "BIDSVersion": "1.7.0",
    "CAPSVersion": "1.0.0",
    "DatasetType": "derivative",
    "Processing": [
        {
            "Name": "t1-linear",
            "Date": "2024-08-06T10:28:21.848950",
            "Author": "ci",
            "Machine": "ubuntu",
            "InputPath": "/mnt/data_ci/T1Linear/in/bids"
        },
        {
            "Name": "pet-linear",
            "Date": "2024-08-06T10:36:27.403373",
            "Author": "ci",
            "Machine": "ubuntu",
            "InputPath": "/mnt/data_ci/PETLinear/in/bids"
        }
    ]
}
```

A `Processing` is described with the following fields:

- `Name`: The name of the processing. For Clinica pipelines, this is the name of the pipeline.
- `Date`: This date is in iso-format and indicates when the processing was run.
- `Author`: This indicates the user name which triggered the processing.
- `Machine`: This indicates the name of the machine on which the processing was run.
- `InputPath`: This is the full path (on the machine on which the processing was run) to the input dataset of the processing.

### Potential problems

The `dataset_description.json` file for BIDS and CAPS datasets was introduced in Clinica `0.9.0`.

This means that results obtained with prior versions of Clinica do not have this file automatically generated.
Clinica will interpret this as a `<1.0.0` dataset and should error with a suggestion of a minimal `dataset_description.json` that you should add to your dataset.
In this situation, create this new file with the suggested content and re-start the pipeline.

You might also see the following error message:

```
Impossible to write the 'dataset_description.json' file in <FOLDER> because it already exists and it contains incompatible metadata.
```

This means that you have version mismatch for the BIDS and/or CAPS specifications.
That is, the versions indicated in the input (or output) dataset(s) does not match the versions currently used by Clinica.
If this happens, it is recommended to re-run the conversion or the pipeline which initially generated the dataset with the current version of Clinica.



## Validation of BIDS datasets

[bids-validator](https://github.com/bids-standard/bids-validator) can be run to ensure that a dataset is [BIDS](glossary.md#bids)-compliant.

Clinica provides tools to curate several publicly available neuroimaging datasets and convert them to [BIDS](http://bids.neuroimaging.io/).

Datasets currently supported can be found [here](../#dataset-converters-clinica-convert).

We decided to ignore several warnings and errors detected by the validator.
These are listed in the `.bids-validator-config.json` and `.bidsignore` files at the root of each [BIDS](http://bids.neuroimaging.io/) folder.
These files are automatically generated by Clinica converters to ignore the following issues:

- Won't fix errors:
    - All participants do not have the same sessions or modalities (`MISSING_SESSION` / `INCONSISTENT_SUBJECTS`).
    - [JSON](glossary.md#json) files with column description are not generated as [TSV](glossary.md#tsv) files already explicit the purpose of the columns in metadata files (`CUSTOM_COLUMN_WITHOUT_DESCRIPTION`).
    - The optional field `authors` is not filled in `dataset_description.json` (`NO_AUTHORS`).
    - The folder `conversion_info` is not BIDS-compliant (`conversion_info/` in `.bidsignore`).
- fMRI-specific: these errors are not fixed as there are no preprocessing pipelines for this modality in Clinica yet
(`SLICE_TIMING_NOT_DEFINED` / `NIFTI_PIXDIM4` / `BOLD_NOT_4D` / `REPETITION_TIME_MUST_DEFINE` / `TASK_NAME_MUST_DEFINE`).
- Probable DICOM errors (`NIFTI_UNIT` / `INCONSISTENT_PARAMETERS`).
