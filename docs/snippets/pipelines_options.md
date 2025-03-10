--8<-- [start:all]
- `-tsv` / `--subjects_sessions_tsv`

This flag allows you to specify in a TSV file the participants belonging to your subset.
For instance, running the [FreeSurfer pipeline](/Pipelines/T1_FreeSurfer.md) on T1w MRI can be done using :

```shell
clinica run t1-freesurfer BIDS_PATH OUTPUT_PATH -tsv my_subjects.tsv
```

<div class="grid" markdown>

=== "TSV Example :"
    ```{ .text .copy }
    participant_id  session_id
    sub-CLNC0001    ses-M000
    sub-CLNC0001    ses-M018
    sub-CLNC0002    ses-M000
    ```
    
!!! warning "Creating the TSV"
    To make the display clearer the rows here contain successive tabs but that should not happen in an actual TSV.
</div>

- `-wd` / `--working_directory`

By default when running a pipeline, a temporary working directory is created. This directory stores all the intermediary inputs and outputs of the different steps of the pipeline. If everything goes well, the output directory is eventually created and the working directory is deleted. 

With this option, a working directory of your choice can be specified. It is very useful for the debugging process or if your pipeline crashes. Then, you can relaunch it with the exact same parameters which will allow you to continue from the last successfully executed node.
For the pipelines that generate many files, such as `dwi-preprocessing` (especially if you run it on multiple subjects), a specific drive/partition with enough space can be used to store the working directory.

- `-np` / `--n_procs`

This flag allows you to exploit several cores of your machine to run pipelines in parallel, which is very useful when dealing with numerous subjects and multiple sessions.
Thanks to [Nipype](https://nipype.readthedocs.io/en/latest/), even for a single subject, a pipeline can be run in parallel by exploiting the cores available to process simultaneously independent sub-parts.

If you do not specify `-np` / `--n_procs` flag, Clinica will detect the number of threads to run in parallel and propose the adequate number of threads to the user.

- `-cn` / `--caps-name`

Use this option if you want to specify the name of the CAPS dataset that will be used inside the `dataset_description.json` file, at the root of the CAPS folder (see [CAPS Specifications](../CAPS/Specifications.md#the-dataset_descriptionjson-file) for more details). This works if this CAPS dataset does not exist yet, otherwise the existing name will be kept. 
--8<-- [end:all]


--8<-- [start:freesurfer]
- `-ap`/`--atlas_path` : In case you wish to use another atlas, specify its folder path with this option. Your atlas will need to be in FreeSurfer `gcs` format (e.g `hemisphere.atlasname_6p0.gcs`). The results will be stored in the same folder as the original results with additional files in `labels`, `stats` and `regional measures`.
- `-overwrite`/`--overwrite_outputs` : Force the overwrite of output files in the CAPS folder with this option.
--8<-- [end:freesurfer]

--8<-- [start:linear]
- `-ui`/`--uncropped_image` : By default, output images are cropped to a **fixed** matrix size of 169×208×179, 1 mm isotropic voxels. This allows reducing the computing power required when training deep learning models afterwards.
Use this option if you do not want to get cropped images.
- `--random_seed` : By default, results are not deterministic. Use this option if you want to obtain a deterministic output. The value you set corresponds to the random seed used by ANTs. This option requires [ANTs](../Software/Third-party.md#ants) version `2.3.0` onwards and is also compatible with [ANTsPy](https://antspyx.readthedocs.io/en/latest/index.html).
--8<-- [end:linear]

--8<-- [start:antspy]
- `--use-antspy` : By default, the pipeline is running with [ANTs](../Software/Third-party.md#ants). Use this flag option if you want to use [ANTsPy](https://antspyx.readthedocs.io/en/latest/index.html) instead.
--8<-- [end:antspy]

--8<-- [start:pet_recon]
- `-rec`/`--reconstruction_method`: Select only images based on a [specific reconstruction method](/Pipelines/PET_Introduction.md#reconstruction-methods).
--8<-- [end:pet_recon]
