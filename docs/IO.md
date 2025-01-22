<!-- markdownlint-disable MD046 -->
# Data handling tools

This page describes data handling tools provided by Clinica for [BIDS](http://bids.neuroimaging.io) and [CAPS](../CAPS/Introduction) compliant datasets.
These tools provide easy interaction mechanisms with datasets, including generating subject lists or merging all tabular data into a single TSV for analysis with external statistical software tools.

## `describe` - Describe a BIDS or CAPS dataset

This tool describes a BIDS or CAPS dataset, basically parsing the `dataset_description.json` file and displaying the information in the console in a nice way.

!!! note
    This tool has been added in Clinica `0.9.0` and is not available in older versions.

```shell
clinica iotools describe DATASET_PATH
```

where:

- `DATASET_PATH`: path to the BIDS or CAPS dataset to be described.

**Examples:**

With a BIDS dataset:

```shell
clinica iotools describe ./bids
                         Dataset information
╭───────────────────────────────┬──────┬─────────────────────────────╮
│                          Name │ Type │ BIDS Specifications Version │
├───────────────────────────────┼──────┼─────────────────────────────┤
│ The mother of all experiments │  raw │ 1.6.0                       │
╰───────────────────────────────┴──────┴─────────────────────────────╯
```

With a CAPS dataset:

```shell
clinica iotools describe ./caps
                                               Dataset information
╭──────────────────────────────────────┬────────────┬─────────────────────────────┬─────────────────────────────╮
│                                 Name │       Type │ BIDS Specifications Version │ CAPS Specifications Version │
├──────────────────────────────────────┼────────────┼─────────────────────────────┼─────────────────────────────┤
│ 05c5794b-2d20-4217-b727-c215d079ab35 │ derivative │ 1.9.0                       │                       1.0.0 │
╰──────────────────────────────────────┴────────────┴─────────────────────────────┴─────────────────────────────╯
                                                                                           Processing information
╭────────────────────────────┬────────────────────────────┬───────────────────┬─────────────────┬────────────────────────────────────────────┬────────────────┬────────────────────────────────────────────╮
│                       Name │                       Date │            Author │         Machine │                                  InputPath │ ClinicaVersion │                               Dependencies │
├────────────────────────────┼────────────────────────────┼───────────────────┼─────────────────┼────────────────────────────────────────────┼────────────────┼────────────────────────────────────────────┤
│ dwi-preprocessing-using-t1 │ 2024-09-01 13:41:05.529799 │ nicolas.gensollen │ UMR-COLLI-MP050 │ /Users/nicolas.gensollen/dwi_datasets/bids │          0.9.0 │                Dependencies                │
│                            │                            │                   │                 │                                            │                │ ╭───────────┬──────────────┬─────────────╮ │
│                            │                            │                   │                 │                                            │                │ │      Name │ VersionCons… │ InstalledV… │ │
│                            │                            │                   │                 │                                            │                │ ├───────────┼──────────────┼─────────────┤ │
│                            │                            │                   │                 │                                            │                │ │      ants │      >=2.3.0 │       2.5.0 │ │
│                            │                            │                   │                 │                                            │                │ │       fsl │      >=6.0.1 │       6.0.2 │ │
│                            │                            │                   │                 │                                            │                │ │    mrtrix │      >=0.0.0 │       3.0.3 │ │
│                            │                            │                   │                 │                                            │                │ │ convert3d │      >=0.0.0 │       1.0.0 │ │
│                            │                            │                   │                 │                                            │                │ ╰───────────┴──────────────┴─────────────╯ │
╰────────────────────────────┴────────────────────────────┴───────────────────┴─────────────────┴────────────────────────────────────────────┴────────────────┴────────────────────────────────────────────╯
```

## `create-subjects-visits` - Generate the list all subjects and visits of a given dataset

A TSV file with two columns (`participant_id` and `session_id`) containing the list of visits for each subject can be created as follows:

```shell
clinica iotools create-subjects-visits BIDS_DIRECTORY OUTPUT_TSV
```

where:

- `BIDS_DIRECTORY`: input folder of a BIDS compliant dataset,
- `OUTPUT_TSV`: output TSV file containing the subjects with their sessions.

Here is an example of the file generated by this tool:

```text
participant_id   session_id
sub-01           ses-M000
sub-02           ses-M024
sub-03           ses-M024
...
```

!!! note
    The format of the participant ID and the session ID follows the [BIDS standard](http://bids.neuroimaging.io/bids_spec1.0.0.pdf).

!!! example

    ```shell
    clinica iotools create-subjects-visits /home/ADNI_BIDS/ adni_participants.tsv
    ```

## `check-missing-modalities` - Check missing modalities for each subject

Starting from a BIDS compliant dataset, this command creates:

1. `<prefix>_ses-<session_label>.tsv`: TSV files for each session available with the list 
   of the modalities found for each subject.
2. `<prefix>_summary.txt`: a text file containing the number and the percentage of modalities missing for each session.
3. `analysis.txt`: a text file in which a table is written per session. This table contains the number of
images per modality per diagnosis when the column `diagnosis` is available in the session-
level files of the BIDS directory.

If no value for `<prefix>` is specified by the user, the default will be `missing_mods`.

```shell
clinica iotools check-missing-modalities [OPTIONS] BIDS_DIRECTORY OUTPUT_DIRECTORY
```

where:

- `BIDS_DIRECTORY`: input folder of a BIDS compliant dataset
- `OUTPUT_DIRECTORY`: output folder
- `-op` / `--output_prefix` (Optional):  prefix used for the name of the output files.
If not specified the default value will be `missing_mods`

If, for example, only the session M00 is available and the parameter `-op` is not specified, the command will create the files:

- `missing_mods_ses-M000.tsv`
- `missing_mods_summary.txt`.

The content of `missing_mods_ses-M000.tsv` will look like:

```Text
participant_id   T1w   DWI
sub-01           1       1
sub-02           1       0
sub-03           1       0
```

Where the column `participant_id` contains all the subjects found and the following columns correspond to the list of all the modalities available for the given dataset.
The availability is expressed by a boolean value.

The nomenclature of the modalities tries to follow, as much as possible, the one proposed by the BIDS standard.

!!! example

    ```shell
    clinica iotools check-missing-modalities /Home/ADNI_BIDS/ /Home/
    clinica iotools check-missing-modalities /Home/ADNI_BIDS/ /Home/ -op new_name
    ```

## `check-missing-processing` - Check missing processing in a CAPS directory

Starting from a CAPS compliant dataset, this command creates a TSV file with columns
`participant_id`, `session_id` and names corresponding to steps of
`t1-volume`, `t1-freesurfer`, `t1-linear`, `pet-volume` and `pet-surface`.
For PET pipelines one column is created per tracer and the PVC option is considered for `pet-volume`.

```shell
clinica iotools check-missing-processing BIDS_DIRECTORY CAPS_DIRECTORY OUTPUT_FILE
```

where:

- `BIDS_DIRECTORY`: input folder of a BIDS compliant dataset
- `CAPS_DIRECTORY`: input folder of a [CAPS](../CAPS/Introduction) compliant dataset
- `OUTPUT_FILE`: output file path (filename included).

The content of `output_file` will look like:

```Text
participant_id   session_id     t1-linear   ...     pet-volume_trc-<tracer>_group-<group_label>_pvc-{True|False}
sub-01           ses-M000        1                   1
sub-01           ses-M012        1                   0
sub-02           ses-M000        0                   0
```

- columns associated with `pet-volume` outputs will specify the PET tracer,
the group label and if a PVC correction was performed.
- columns associated with `t1-volume` outputs will specify the group label and which steps of `t1-volume` were performed.
- columns associated with `pet-surface` outputs will specify the PET tracer used.


## `merge-tsv` - Gather BIDS and CAPS data into a single TSV file

[BIDS](http://bids.neuroimaging.io) and [CAPS](../CAPS/Introduction) datasets are composed of multiple TSV files for the different subjects and sessions.
While this has some advantages, it may not be convenient when performing statistical analyses (with external statistical software tools for instance).
This command merges all the TSV files into a single larger TSV file and can be run with the following command line:

```shell
clinica iotools merge-tsv [OPTIONS] BIDS_DIRECTORY OUTPUT_TSV
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](http://bids.neuroimaging.io) hierarchy.
- `OUTPUT_TSV` is the path of the output TSV file.
If a directory is specified instead of a file name, the default name for the file created will be `merge-tsv.tsv`.

The optional arguments allow the user to also merge data from a CAPS directory, which will be concatenated to the BIDS summary.
The main optional arguments are the following:

- `-caps`: input folder of a [CAPS](../CAPS/Introduction) compliant dataset

If a CAPS folder is given, data generated by the pipelines of Clinica (regional measures) will be merged to the output file, and a summary file containing the names of the atlases merged will be generated in the same folder.

- `-tsv`: input list of subjects and sessions

If an input list of subjects and sessions is given, the merged file will only gather information from the pairs of subjects and sessions specified.

!!! example

    ```shell
    clinica iotools merge-tsv /Home/ADNI_BIDS /Home/merge-tsv.tsv -caps /Home/ADNI_CAPS -tsv /Home/list_subjects.tsv
    ```

    The output file will contain one row for each visit:

    ```Text
    participant_id   session_id   date_of_birth   ...   ..._ROI-0   ..._ROI-1  ...
    sub-01           ses-M000      25/04/41        ...   9.824750    0.023562
    sub-01           ses-M018      25/04/41        ...   8.865353    0.012349
    sub-02           ses-M000      09/01/91        ...   9.586342    0.027254
    ...
    ```


!!! Note for t1-volume and pet-volume pipelines
    The suffix "_intensity" is added systematically to the atlas statistics of t1-volume and pet-volume pipelines.

A complete list of optional arguments can be obtained with the command line `clinica merge-tsv --help`

## `center-nifti` - Center NIfTI files of a BIDS directory

Your [BIDS](http://bids.neuroimaging.io) dataset may contain NIfTI files where the origin does not correspond to the center of the image (i.e. the anterior commissure).
SPM is especially sensitive to this case and segmentation procedures may result in blank images or even fail.
To mitigate this issue we offer a simple tool that converts generates from your BIDS a new dataset with centered NIfTI files for the selected modalities.

!!! warning "By default :"

    - This tool will only center **T1w** images.
    - Only NIfTI volumes whose center is at **more than 50 mm** from the origin of the world coordinate system are centered. This threshold has been chosen empirically after a set of experiments to determine at which distance from the origin SPM segmentation and coregistration procedures stop working properly.

```shell
clinica iotools center-nifti [OPTIONS] BIDS_DIRECTORY OUTPUT_BIDS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](http://bids.neuroimaging.io) hierarchy.
- `OUTPUT_BIDS_DIRECTORY` is the output path to the new version of your BIDS dataset, with faulty NIfTI centered.
This folder can be empty or nonexistent.

Optional arguments:

- `--modality` is a case-insensitive parameter that defines which modalities are converted.

    !!! tip "How to use :"
        - If you want to convert T1w images only, do not use the option.
        - If you want to convert all types of pet, use `--modality pet` or `-m pet`
        - If you want to convert 18FFDG_PET, use `-m 18ffdg_pet`
        - If you want to convert both pet and T1, use `-m T1 -m pet`

        Basically, the software searches for the modality key inside the filename. Understanding this, you can now center any modality you want!

- `--center_all_files` is an option that forces Clinica to center all the files of the modalities selected with the `--modality` flag.

!!! note
    The images contained in the input `bids_directory` folder that do not need to be centered will also be copied to the output folder `new_bids_directory`.

     The list of the converted files will appear in a text file in `new_bids_directory/centered_nifti_list_TIMESTAMP.txt`.
