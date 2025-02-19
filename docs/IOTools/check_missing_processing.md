# `check-missing-processing` - Check missing processing in a CAPS directory

Starting from a CAPS compliant dataset, this command creates a TSV file with columns `participant_id`, `session_id` and names corresponding to steps of `t1-volume`, `t1-freesurfer`, `t1-linear`, `pet-volume` and `pet-surface`.

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

- columns associated with `pet-volume` outputs will specify the PET tracer, the group label and if a PVC correction was performed.
- columns associated with `t1-volume` outputs will specify the group label and which steps of `t1-volume` were performed.
- columns associated with `pet-surface` outputs will specify the PET tracer used.
