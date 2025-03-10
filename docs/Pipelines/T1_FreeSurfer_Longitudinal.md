<!-- markdownlint-disable MD046-->
# `t1-freesurfer-longitudinal` – FreeSurfer-based longitudinal processing of T1-weighted MR images

This pipeline processes a series of images acquired at different time points for the same subject with the longitudinal FreeSurfer stream
[[Reuter et al., 2012](http://dx.doi.org/10.1016/j.neuroimage.2012.02.084)]
to increase the accuracy of volume and thickness estimates.
It does so in a single command consisting of two consecutive steps:

- **Template creation** This corresponds to the
[FreeSurfer `recon-all -base`](https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing)
command.
It is executed on the subject to produce an unbiased template image using robust and
inverse consistent registration
[[Reuter et al., 2010](http://dx.doi.org/10.1016/j.neuroimage.2010.07.020)].

- **Longitudinal correction** This template is then used as an initialization to a sequence of processing steps triggered by the command
[FreeSurfer `recon-all -long`](https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing)
to perform segmentation, extract surfaces and derive measurements at each time point.

## Prerequisites

The pipeline requires a prior run of the cross-sectional [`t1-freesurfer`](../T1_FreeSurfer) pipeline.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of [FreeSurfer 6.0](../Software/Third-party.md#freesurfer) on your computer.

## Running the pipeline

The pipeline can be run with the following command line:

```Text
clinica run t1-freesurfer-longitudinal [OPTIONS] CAPS_DIRECTORY 
```

where:

- `CAPS_DIRECTORY` is the input/output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.
If you want to run the pipeline on a subset of your dataset, you can use the `-tsv` flag to specify in a TSV file the participants and the corresponding sessions of interest.

with specific options :
- -ap/--atlas_path : In case you wish to use another atlas, specify its folder path with this option. Your atlas will need to be in FreeSurfer `gcs` format (e.g `hemisphere.atlasname_6p0.gcs`). The results will be stored in the same folder as the original results with additional files in `labels`, `stats` and `regional measures`.
- -overwrite/--overwrite_outputs : Force the overwrite of output files in the CAPS folder with this option.

??? info "Optional parameters common to all pipelines"
    --8<-- "snippets/pipelines_options.md"

!!! note "Computational time"
    The computational time for one subject is around 6-8 hours (creation of the unbiased template) + 2-5 hours per corresponding session.
    The code execution speed depends on your CPU and the quality of your input T1 volumes.
    Please be aware that even though the pipeline runs in parallel, processing many subjects and sessions (e.g. ADNI dataset) is time consuming.


??? warning "Case when longitudinal correction is performed on macOS"
    If your run the `t1-freesurfer-longitudinal` pipeline on macOS, you will see warning messages when longitudinal correction is performed e.g.:

    ```text
    [19:29:11] Needs to create a $SUBJECTS_DIR folder in /var/folders/m_/j76n37kn4vs6zsj8fq0qcgy0000dn9/T/tmpbrbg9w1n for sub-01 | ses-2011 | long-20112015 (macOS case).
    ```

    When longitudinal correction is performed on macOS, FreeSurfer (`recon-all -long`) may crash if it has to handle a very long path.
    The workaround we are currently using is that FreeSurfer will be run in a temporary folder (e.g. `/tmp/tmp<hash>`) instead of `<path_to_wd>/t1-freesurfer-longitudinal-correction/ReconAll`.
    Then, the results will be copied to the working directory before the temporary folder is deleted.

??? warning "Case where one session is used for a participant"
    If your CAPS directory contains a participant with one session e.g.:

    ```text
    CAPS_DIRECTORY
    └── subjects
     ├── sub-CLNC01
     │   ├── ses-M000
     │   │   └── t1
     │   │       └── freesurfer_cross_sectional
     │   └── ses-M018
     │       └── t1
     │           └── freesurfer_cross_sectional
     └── sub-CLNC02
         └── ses-M000
             └── t1
                 └── freesurfer_cross_sectional
    ```

    You will see this type of message when running Clinica:

    ```console
    $ clinica run t1-freesurfer-template CAPS -np 2 -wd <path_to_wd>
    The pipeline will be run on the following 2 participant(s):
        sub-CLNC01 | ses-M018, ses-M000 | long-M000M018
        sub-CLNC02 | ses-M000 | long-M000
    List available in <path_to_wd>/t1-freesurfer-template/participants.tsv
    The pipeline will last approximately 10 hours per participant.
    [13:33:43] sub-CLNC02 | long-M000 has only one time point. Needs to create a $SUBJECTS_DIR folder in /tmp/tmpe7ztq9hq
    [13:33:43] Running pipeline for sub-CLNC01 | long-M000M018
    [13:33:43] Running pipeline for sub-CLNC02 | long-M000
    [19:51:18] sub-CLNC01 | long-M000M018 has completed
    [20:15:04] Segmentation of sub-CLNC02 | long-M000 has moved to working directory and $SUBJECTS_DIR folder (/tmp/tmpe7ztq9hq) was deleted
    [20:15:05] sub-CLNC02 | long-M000 has completed
    [20:15:09] The t1-freesurfer-template pipeline has completed. You can now delete the working directory (<path_to_wd>/t1-freesurfer-template).
    ```

    When one session is used for template creation, FreeSurfer (`recon-all -base`) may crash if it has to handle a very long path.
    The workaround we are currently using is that when one time point is detected for a given participant, FreeSurfer will be run in a temporary folder (e.g. `/tmp/tmp<hash>`) instead of `<path_to_wd>/t1-freesurfer-template/ReconAll`.
    Then, the results will be copied to the working directory before the temporary folder is deleted.

## Outputs

### Template creation

Results stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#t1-freesurfer-longitudinal-freesurfer-based-longitudinal-processing-of-t1-weighted-mr-images):
`subjects/<participant_id>/<long_id>/freesurfer_unbiased_template/<participant_id>_<long_id>`.

`<long_label>` is an identifier defined by concatenating all the sessions associated with the current `<participant_id>` (e.g. if the template for participant `sub-CLNC01` is built from sessions `M00`, `M01`, `M05`, then `<long_label>` will be `M00M01M05`).
See [CAPS specifications](../../CAPS/Introduction/#subject-and-group-naming) for full definition and example of `<long_id>`.

This folder contains the standard output structure of the `recon-all` command (`label/`, `mri/`, `surf/`, etc.) already explained in the [`t1-freesurfer`](../T1_FreeSurfer) pipeline.

### Longitudinal correction

Results are stored in the following folder of
[CAPS hierarchy](../../CAPS/Specifications/#t1-freesurfer-longitudinal-freesurfer-based-longitudinal-processing-of-t1-weighted-mr-images):
`subjects/<participant_id>/<session_id>/t1/<long_id>/freesurfer_longitudinal/<participant_id>_<session_id>.long.<participant_id>_<long_id>`.

Similar to the template creation folder, the longitudinal folder contains the standard output structure of the `recon-all` command (`label/`, `mri/`, `surf/`, etc.).
Among the files generated by FreeSurfer, you may be interested in the following outputs:

- `*/mri/aseg.mgz`: subcortical segmentation volume after correction with the unbiased template
- `*/mri/wm.mgz`: white matter mask after correction with the unbiased template
- `*/mri/brainmask.mgz`: skull-stripped volume after correction with the unbiased template
- `*/surf/{l|r}.white`: white surface between white matter and gray matter after correction with the unbiased template
- `*/surf/{l|r}.pial`: pial surface between gray matter and CSF after correction with the unbiased template

(where `*` stands for `<participant_id>_<session_id>.long.<participant_id>_<long_id>`)

More details regarding the `recon-all` output files can be found on the [FreeSurfer website](https://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllOutputFiles).

<!-- TODO: Add note regarding TSV files generated in this sub-section -->

!!! note
    The full list of features extracted from the FreeSurfer pipeline can be found in the
    [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#t1-freesurfer-longitudinal-freesurfer-based-longitudinal-processing-of-t1-weighted-mr-images).

<!-- ## Visualization of the results

!!! note
    The visualization command is not available for the moment. Please come back later, this section will be updated ASAP. -->

## Describing this pipeline in your paper

!!! cite "Example of paragraph (short version):"
    These results have been obtained using the `t1-freesurfer-longitudinal` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    This pipeline is a wrapper of different tools of the FreeSurfer software
    (<http://surfer.nmr.mgh.harvard.edu/>)
    [[Fischl et al., 2012](http://dx.doi.org/10.1016/j.neuroimage.2012.01.021)].
    This processing builds a subject-dependent template space and
    extracts volume and thickness estimates in this space at different points in time.

??? cite "Example of paragraph (long version):"
    These results have been obtained using the `t1-freesurfer-longitudinal` pipeline of Clinica
    [Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675).
    This pipeline is a wrapper of different tools of the FreeSurfer software
    [[Fischl et al., 2012](http://dx.doi.org/10.1016/j.neuroimage.2012.01.021)],
    which is documented and freely available for download online (<http://surfer.nmr.mgh.harvard.edu/>).
    The technical details of the procedures concerned with longitudinal analysis are described in prior publications
    [[Reuter et al., 2010](https://doi.org/10.1016/j.neuroimage.2010.07.020);
    [Reuter et al., M00](http://dx.doi.org/10.1016/j.neuroimage.M00.02.076);
    [Reuter et al., 2012](http://dx.doi.org/10.1016/j.neuroimage.2012.02.084)].
    The pipeline processes a series of images acquired at different time points for the same subject.
    It first produces an unbiased (with respect to any time point) template volume, and then, for each time point, uses the template as an initialisation (tailored to the subject) for the FreeSurfer cortical reconstruction process.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/GHAXT4R5).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).

## Advanced usage

The two main processing steps of the `t1-freesurfer-longitudinal` pipeline can be performed individually:

- **Template creation**

    Command line:

    ```Text
    clinica run t1-freesurfer-template [OPTIONS] CAPS_DIRECTORY
    ```

- **Longitudinal correction**

    Command line:

    ```Text
    clinica run t1-freesurfer-longitudinal-correction CAPS_DIRECTORY
    ```
