<!-- markdownlint-disable MD046 -->
# Interacting with Clinica

## Preparing your data

The easiest way to use Clinica is to have your data organized using the [BIDS standard](http://bids.neuroimaging.io/).
[BIDS](http://bids.neuroimaging.io/) is currently becoming the standard for data organization in the brain imaging community and we strongly recommend to use it.

If your dataset does not follow this standard, you will need to convert it:

- If your data are in DICOM format, you can use one of the converters referenced on the [BIDS website](https://bids.neuroimaging.io/benefits.html#converters).
- Otherwise, Clinica includes converters for public datasets such as [ADNI](http://adni.loni.usc.edu/), [AIBL](https://aibl.csiro.au), [NIFD](http://4rtni-ftldni.ini.usc.edu/) and [OASIS](http://www.oasis-brains.org).
See [here](../DatabasesToBIDS) for more details.

!!! note "Regarding cross-sectional BIDS datasets"
    If you run Clinica with a dataset containing no timepoints e.g.:

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
    This will result to a new dataset ressembling to:
    ```Text
    BIDS
    └── sub-CLNC0001
        └── ses-M00
            ├── anat
            │   └── sub-CLNC0001_ses-M00_T1w.nii.gz
            └── pet
                ├── sub-CLNC0001_ses-M00_trc-18FFDG_pet.json
                └── sub-CLNC0001_ses-M00_trc-18FFDG_pet.nii.gz
    ```

!!! tip
     If you need to create BIDS compliant datasets or need tutorials on BIDS, you can look at this [BIDS Starter Kit](https://github.com/INCF/bids-starter-kit/).

## Clinica command-line interface

Clinica's main usage is through command-line.
Clinica supports autocompletion: to see the list of commands, simply type `clinica` followed by ++tab++.

In general, a Clinica command-line has the following syntax:

```bash
clinica category_of_command command argument options
```

where the arguments are usually your input/output folders, and where the options look like `--flag_1 option_1 --flag_2 option_2`.

Please note that the ordering of options on the command-line is not important, whereas arguments must be given in the exact order specified in the documentation (or in the command line helper).

## Categories of command line

The command-line `clinica` has been divided into four main categories.

### `clinica run`

This category allows the user to run the different image processing and analysis pipelines using the following syntax:

```bash
clinica run modality-pipeline bids_directory caps_directory -tsv my_participants.tsv
```

"modality" is a prefix that corresponds to the data modality (e.g. T1, DWI, fMRI, PET) or to the category of processing (machine learning, statistics...).
If you execute `clinica run --help`, you can see the list of `modality-pipeline` available: they correspond to the different pipelines displayed on the [main page of the documentation](..).

<!-- ### clinica visualize

!!! note
    We are currently rewriting this section. We will update this section ASAP. -->

### `clinica convert`

These tools allow you to convert unorganized datasets from publicly available neuroimaging studies into a BIDS hierarchy.
Clinica currently includes converters for [ADNI](http://adni.loni.usc.edu/), [AIBL](https://aibl.csiro.au) and [OASIS](http://www.oasis-brains.org).
See [here](../DatabasesToBIDS) for more details.

### `clinica iotools`

`iotools` is a set of tools that allows the user to handle [BIDS](http://bids.neuroimaging.io) and [CAPS](../CAPS/Introduction) datasets.
It allows generating lists of subjects or merging all tabular data into a single TSV file for analysis with external statistical software packages.
See [here](../IO) for more details.

### `clinica generate` (for developers)

This category allows developers to generate the skeleton for a new pipeline.
The syntax is:

```bash
clinica generate template "Modality My Pipeline" -d output_folder
```

## The main arguments

### `BIDS_DIRECTORY` and/or `CAPS_DIRECTORY`

Running a pipeline involves most of the time these two parameters:

- `BIDS_DIRECTORY`, which is the input folder containing the dataset in a [BIDS](../BIDS) hierarchy;
- `CAPS_DIRECTORY`, which is the output folder containing the expected results in a [CAPS](../CAPS/Introduction) hierarchy.
It can be also the input folder containing the dataset in a [CAPS](../CAPS/Introduction) hierarchy.

### `GROUP_LABEL`

You will see the `GROUP_LABEL` argument when working on any group-wise analysis (e.g. template creation from a list of subjects, statistical analysis).
This is simply a label name that will define the  group of subjects used for this analysis.
It will be written in your output CAPS folder, for possible future reuses.
For example, an `AD` group ID label could be used when creating a template for a group of Alzheimer’s disease patients.
Any time you would like to use this `AD` template you will need to provide the group ID used to identify the pipeline output obtained from this group.
You might also use `CNvsAD`, for instance, as group ID for a statistical group comparison between patients with Alzheimer's disease (`AD`) and cognitively normal (`CN`) subjects.

## Common options

### `-tsv` / `--subjects_sessions_tsv`

The `-tsv` flag allows you to specify in a TSV file the participants belonging to your subset.
For instance, running the [FreeSurfer pipeline](../Pipelines/T1_FreeSurfer) on T1w MRI can be done using :

```shell
clinica run t1-freesurfer path/to/my/bids/dataset path/where/results/will/be/stored -tsv my_list_of_subjects.tsv
```

where your TSV file looks as follows:

```text
participant_id  session_id
sub-CLNC0001    ses-M00
sub-CLNC0001    ses-M18
sub-CLNC0001    ses-M36
sub-CLNC0002    ses-M00
sub-CLNC0002    ses-M18
sub-CLNC0002    ses-M36
sub-CLNC0003    ses-M00
```
<!-- Note that to make the display clearer, the rows contain successive tabs, which should not happen in an actual TSV file. -->

### `-wd` / `--working_directory`

In every pipeline, a working directory can be specified.
This directory gathers all the inputs and outputs of the different steps of the pipeline.
It is then very useful for the debugging process.
It is specially useful in the case where your pipeline execution crashes and you relaunch it with the exact same parameters, allowing you to continue from the last successfully executed node. <!--If you do not specify any working directory, a temporary one will be created, then deleted at the end if everything went well.--> For the pipelines that generate many files, such as `dwi-preprocessing` (especially if you run it on multiple subjects), a specific drive/partition with enough space can be used to store the working directory.

### `-np` / `--n_procs`

The `--n_procs` flag allows you to exploit several cores of your machine to run pipelines in parallel, which is very useful when dealing with numerous subjects and multiple sessions.
Thanks to Nipype, even for a single subject, a pipeline can be run in parallel by exploiting the cores available to process simultaneously independent sub-parts.

If you do not specify `-np` / `--n_procs` flag, Clinica will detect the number of threads to run in parallel and propose the adequate number of threads to the user.

## :warning: Known issues

Matlab and SPM12 (whose implementation is based on Matlab) can sometimes randomly crash, causing a rather unreadable error in the console.
Those events are unpredictable.
In case it occurs to you, please do the following:

- Check that you have a valid Matlab license.
- Before relaunching the command line, be sure to remove the content of the working directory (if you specified one).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
