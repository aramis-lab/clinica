# `check-missing-modalities` - Check missing modalities for each subject

Starting from a BIDS compliant dataset, this command creates:

1. `<prefix>_ses-<session_label>.tsv`: TSV files for each session available with the list of the modalities found for each subject.
2. `<prefix>_summary.txt`: a text file containing the number and the percentage of modalities missing for each session.
3. `analysis.txt`: a text file in which a table is written per session. This table contains the number of images per modality per diagnosis when the column `diagnosis` is available in the session-level files of the BIDS directory.

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
