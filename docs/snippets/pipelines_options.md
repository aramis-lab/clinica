## Common options for pipelines

### `-tsv` / `--subjects_sessions_tsv`

This flag allows you to specify in a TSV file the participants belonging to your subset.
For instance, running the [FreeSurfer pipeline](/Pipelines/T1_FreeSurfer.md) on T1w MRI can be done using :

```shell
clinica run t1-freesurfer BIDS_PATH OUTPUT_PATH -tsv my_subjects.tsv
```

where your TSV file looks as follows:

```text
participant_id  session_id
sub-CLNC0001    ses-M000
sub-CLNC0001    ses-M018
sub-CLNC0002    ses-M000
...
```

!!! warning "Writing the TSV"
    To make the display clearer the rows here contain successive tabs but that should not happen in an actual TSV.

### `-wd` / `--working_directory`

In every pipeline, a working directory can be specified.
This directory gathers all the inputs and outputs of the different steps of the pipeline.
It is then very useful for the debugging process.
In the case where your pipeline execution crashes, you can relaunch it with the exact same parameters which will allow you to continue from the last successfully executed node.

!!! info "No working directory"
    If you do not specify any working directory, a temporary one will be created, then deleted at the end if everything went well.

For the pipelines that generate many files, such as `dwi-preprocessing` (especially if you run it on multiple subjects), a specific drive/partition with enough space can be used to store the working directory.

### `-np` / `--n_procs`

This flag allows you to exploit several cores of your machine to run pipelines in parallel, which is very useful when dealing with numerous subjects and multiple sessions.
Thanks to [Nipype](https://nipype.readthedocs.io/en/latest/), even for a single subject, a pipeline can be run in parallel by exploiting the cores available to process simultaneously independent sub-parts.

If you do not specify `-np` / `--n_procs` flag, Clinica will detect the number of threads to run in parallel and propose the adequate number of threads to the user.

### `-cn` / `--caps-name`

Use this option if you want to specify the name of the CAPS dataset that will be used inside the `dataset_description.json` file, at the root of the CAPS folder (see [CAPS Specifications](../CAPS/Specifications.md#the-dataset_descriptionjson-file) for more details). This works if this CAPS dataset does not exist yet, otherwise the existing name will be kept. 
