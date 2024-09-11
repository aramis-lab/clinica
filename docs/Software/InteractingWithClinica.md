<!-- markdownlint-disable MD046 -->
# Interacting with Clinica

## Preparing your data

Clinica pipelines require to have your data organized in the [BIDS format](../BIDS.md).
BIDS is currently becoming the standard for data organization in the brain imaging community and we strongly recommend to use it.

If your dataset does not follow this standard, you will need to **convert it** :

- If your data are in [DICOM format](https://www.dicomstandard.org), you can use one of the converters referenced on the [BIDS website](https://bids.neuroimaging.io/benefits.html#converters).
- Otherwise, Clinica includes converters for public datasets :

??? info "Clinica available converters"
    --8<-- "snippets/inventory_converters.md"

!!! warning "Clinica and cross-sectional BIDS datasets"
    If you run Clinica with a dataset containing **no timepoints** e.g.:

    ```Text
    BIDS
    └── sub-CLNC0001
        ├── anat
        │   └── sub-CLNC0001_T1w.nii.gz
        └── pet
            ├── sub-CLNC0001_trc-18FFDG_pet.json
            └── sub-CLNC0001_trc-18FFDG_pet.nii.gz
    ```
    Clinica will propose you to create a new BIDS dataset with a fake timepoint.
    This will result in a new dataset resembling to:
    ```Text
    BIDS
    └── sub-CLNC0001
        └── ses-M000
            ├── anat
            │   └── sub-CLNC0001_ses-M000_T1w.nii.gz
            └── pet
                ├── sub-CLNC0001_ses-M000_trc-18FFDG_pet.json
                └── sub-CLNC0001_ses-M000_trc-18FFDG_pet.nii.gz
    ```
    You should **accept** the creation of this dataset, else Clinica will fail.

!!! tip
     If you need to create [BIDS](http://bids.neuroimaging.io/) compliant datasets or need tutorials on [BIDS](http://bids.neuroimaging.io/), you can look at this [BIDS Starter Kit](https://github.com/INCF/bids-starter-kit/).

## Clinica command-line interface

Clinica's main usage is through command-line.

!!! tip "See the list of available commands on your terminal"
    - Do not hesitate to use `-h` or `--help` after any component of the command line,
    ex : `clinica -h`, `clinica run -h` ... to see what options are available.
    - Clinica supports autocompletion for zsh shell users : to see the list of commands, simply type `clinica` followed by ++tab++.

In general, a Clinica command-line has the following syntax:

```bash
clinica category_of_command command argument options
```

where the arguments are usually your input/output folders, and where the options look like `--flag_1 option_1 --flag_2 option_2`.

!!! note "Components order"
    Please note that the ordering of options on the command-line is not important,
    whereas arguments must be given in the exact order specified in the documentation (or in the command line helper).

## Categories of command line

The command-line `clinica` has been divided into four main categories :

### `clinica run`

This category allows the user to run the different image processing and analysis pipelines using the following syntax:

```bash
clinica run modality-pipeline bids_directory caps_directory options
```

"modality" is a prefix that corresponds to the data modality (e.g. T1, DWI, fMRI, PET) or to the category of processing (machine learning, statistics...).
If you execute `clinica run --help`, you can see the list of `modality-pipeline` available :

??? info "Clinica available pipelines"
    --8<-- "snippets/inventory_pipelines.md"

!!! tip "Clinica run logs"
    Clinica run logs are written in the current working directory by default. A different directory may be specified by setting the `CLINICA_LOGGING_DIR` environment variable.

### `clinica convert`

These tools allow you to convert unorganized datasets from publicly available neuroimaging studies into a [BIDS](http://bids.neuroimaging.io/) hierarchy.

Clinica currently includes some converters for public datasets :

??? info "Clinica available converters"
    --8<-- "snippets/inventory_converters.md"

### `clinica iotools`

`iotools` is a set of tools that allows the user to handle [BIDS](http://bids.neuroimaging.io) and [CAPS](../CAPS/Introduction.md) datasets.
It allows generating lists of subjects or merging all tabular data into a single TSV file for analysis with external statistical software packages.
See [here](../IO.md) for more details.

### `clinica generate` (for developers)

This category allows developers to generate the skeleton for a new pipeline.
The syntax is:

```bash
clinica generate template "Modality My Pipeline" -d output_folder
```

## The main arguments

### `BIDS_DIRECTORY` and/or `CAPS_DIRECTORY`

Running a pipeline involves most of the time these two parameters:

- `BIDS_DIRECTORY`, which is the input folder containing the dataset in a [BIDS](../BIDS.md) hierarchy;
- `CAPS_DIRECTORY`, which is the output folder containing the expected results in a [CAPS](../CAPS/Introduction.md) hierarchy.
It can be also the input folder containing the dataset in a [CAPS](../CAPS/Introduction.md) hierarchy.

### `GROUP_LABEL`

You will see the `GROUP_LABEL` argument when working on any group-wise analysis (e.g. template creation from a list of subjects, statistical analysis).
This is simply a label name that will define the  group of subjects used for this analysis.
It will be written in your output [CAPS](../CAPS/Introduction.md) folder, for possible future reuses.
For example, an `AD` group ID label could be used when creating a template for a group of Alzheimer’s disease patients.
Any time you would like to use this `AD` template you will need to provide the group ID used to identify the pipeline output obtained from this group.
You might also use `CNvsAD`, for instance, as group ID for a statistical group comparison between patients with Alzheimer's disease (`AD`) and cognitively normal (`CN`) subjects.

--8<-- "snippets/pipelines_options.md"

--8<-- "snippets/converters_options.md"

--8<-- "snippets/known_issues.md:matlab"


## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
