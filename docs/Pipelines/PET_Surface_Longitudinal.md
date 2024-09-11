<!-- markdownlint-disable MD046 -->
# `pet-surface-longitudinal` - Surface-based longitudinal processing of PET images

This pipeline performs several processing steps for the longitudinal analysis of PET data on the cortical surface
[[Marcoux et al., 2018](https://doi.org/10.3389/fninf.2018.00094)]:

- co-registration of PET and T1-weighted MR images;
- intensity normalization;
- partial volume correction;
- robust projection of the PET signal onto the subject’s cortical surface;
- spatial normalization to a template;
- atlas statistics.

This pipeline relies mainly on tools from **[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)** and
**[PETPVC](https://github.com/UCL/PETPVC)** [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)].

The only difference with the [`pet-surface`](../PET_Surface) pipeline is that the subject’s cortical surface is obtained with the [`t1-freesurfer-longitudinal`](../T1_FreeSurfer_Longitudinal) pipeline and
not the [`t1-freesurfer`](../T1_FreeSurfer) pipeline.

!!! note "Clinica & BIDS specifications for PET modality"
    Since Clinica `v0.6`, PET data following the official specifications in BIDS version 1.6.0 are now compatible with Clinica.
    See [BIDS](../../BIDS) page for more information.

## Prerequisite

You need to have performed the [`t1-freesurfer-longitudinal`](../T1_FreeSurfer_Longitudinal) pipeline on your T1-weighted MR images.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of [FreeSurfer 6.0](../Software/Third-party.md#freesurfer), [FSL 6.0](../Software/Third-party.md#fsl), and [PETPVC 1.2.4](../Software/Third-party.md#petpvc) (which depends on [ITK 4](../Software/Third-party.md#itk)) on your computer.

In addition, you also need to either install [SPM12](../Software/Third-party.md#spm12) and [Matlab](../Software/Third-party.md#matlab), or [spm standalone](../Software/Third-party.md#spm12-standalone).


## Running the pipeline

The pipeline can be run with the following command line:

```shell
clinica run pet-surface-longitudinal [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY ACQ_LABEL
                                     {pons|cerebellumPons|pons2|cerebellumPons2} PVC_PSF_TSV
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](../../BIDS) hierarchy.
- `CAPS_DIRECTORY` acts both as an input folder (where the results of the `t1-freesurfer-longitudinal` pipeline are stored) and
as the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.
- `ACQ_LABEL` is the label given to the PET acquisition, specifying the tracer used (`trc-<acq_label>`).
- The reference region is used to perform intensity normalization (i.e. dividing each voxel of the image by the average uptake in this region) resulting in a standardized uptake value ratio (SUVR) map.
It can be `cerebellumPons` or `cerebellumPons2 (used for amyloid tracers) or `pons` or `pons2` (used for FDG).
- `PVC_PSF_TSV` is the TSV file containing the `psf_x`, `psf_y` and `psf_z` of the PSF for each PET image.
More explanation is given in [PET Introduction](../PET_Introduction) page.

!!! info
    Since the release of Clinica v0.3.8, the handling of PSF information has changed.
    In previous versions of Clinica, each BIDS-PET image had to contain a JSON file with the `EffectiveResolutionInPlane` and `EffectiveResolutionAxial` fields corresponding to the PSF in mm.
    `EffectiveResolutionInPlane` is replaced by both `psf_x` and `psf_y`, `EffectiveResolutionAxial` is replaced by `psf_z` and the `acq_label` column has been added.
    Additionally, the SUVR reference region is now a compulsory argument: it will be easier for you to modify Clinica if you want to add a custom reference region ([PET Introduction](../PET_Introduction) page).
    Choose `cerebellumPons` for amyloid tracers or `pons` for FDG to have the previous behavior.

Pipeline options:

- `-np`: This parameter specifies the number of threads to run in parallel.
We recommend using `your_number_of_cpu - 1`.
Please note that PETPVC is extremely demanding in terms of resources and
may cause the pipeline to crash if many subjects happen to be partial volume corrected at the same time
(Error : `Failed to allocate memory for image`).
To mitigate this issue, you can do the following:

    **1)** Use a working directory when you launch Clinica.

    **2)** If the pipeline crashed, just launch again the command (while giving the same working directory).
    The whole processing will continue where it left (you can reduce the number of threads to run in parallel the second time).

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with Clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `clinica run pet-surface-longitudinal --help` to see the full list of parameters.

!!! failure "Known error on macOS"
    If you are running `pet-surface-longitudinal` on macOS, we noticed that if the path to the CAPS is too long, the pipeline fails when the `gtmseg` command from FreeSurfer is executed.
    This generates crash files with `gtmseg` in the filename, for instance:

    ```console
    $ nipypecli crash crash-20210404-115414-sheldon.cooper-gtmseg-278e3a57-294f-4121-8a46-9975801f24aa.pklz
    [...]
    Abort
    ERROR: mri_gtmseg --s sub-ADNI011S4105_ses-M000 --usf 2 --o gtmseg.mgz --apas apas+head.mgz --no-subseg-wm --no-keep-cc --no-keep-hypo
    gtmseg exited with errors
    Standard error:
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M000/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53505/tmp.mgh' (type = MGH )                       [ ok ]
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53727/tmp.mgh' (type = MGH )                       [ ok ]
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M000/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53946/tmp.mgh' (type = MGH )                       [ ok ]
    Return code: 1
    ```

    This is under investigation (see [Issue #119](https://github.com/aramis-lab/clinica/issues/119) for details) and will be solved as soon as possible.

!!! warning "Case where several longitudinal IDs are present"
    If a subject has more than two longitudinal IDs (e.g. `long-M000M018` and `long-M000M018M036`),
    Clinica is not currently able to deal with this particular case.
    This will be fixed in the future.

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#pet-surface-longitudinal-surface-based-longitudinal-processing-of-pet-images):
`subjects/<participant_id>/<session_id>/pet/<long_id>/surface`

The main output files are (where `*` stands for `<participant_id>_<session_id>_<long_id>`):

- `*_trc-<label>_pet_space-<label>_suvr-<label>_pvc-iy_hemi-<label>_fwhm-<value>_projection.mgh`:
PET data that can be mapped onto meshes.
If the `space` is `fsaverage`, it can be mapped either onto the white or pial surface of FsAverage.
If the `space` is `native`, it can be mapped onto the white or pial surface of the subject’s surface (i.e. `{l|r}h.white`, `{l|r}h.pial` files from the `t1-freesurfer-longitudinal` pipeline).
- `*_hemi-{left|right}_midcorticalsurface`:
surface at equal distance between the white matter/gray matter interface and the pial surface (one per hemisphere).
- `atlas_statistics/*_trc-<label>_pet_space-<label>_pvc-iy_suvr-<label>_statistics.tsv`:
TSV files summarizing the regional statistics on the labelled atlases (Desikan and Destrieux).

!!! note
    The full list of output files from the `pet-surface-longitudinal` pipeline can be found in [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#pet-surface-longitudinal-surface-based-longitudinal-processing-of-pet-images).

## Describing this pipeline in your paper

!!! cite "Example of paragraph:"
    These results have been obtained using the `pet-surface-longitudinal` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675);
    [Marcoux et al., 2018](https://doi.org/10.3389/fninf.2018.00094)].
    The subject’s PET image was registered to the T1-weighted MRI using spmregister
    ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)) and intensity normalized using
    the [pons | pons and cerebellum] from the Pick atlas in MNI space as reference region
    (registration to MNI space was performed using
    [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)).
    Partial volume correction was then performed using the iterative Yang algorithm implemented in
    [PETPVC](https://github.com/UCL/PETPVC)
    [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)]
    with regions obtained from gtmseg ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)).
    Based on the subject’s white surface and cortical thickness estimated with the longitudinal pipeline of FreeSurfer
    [[Fischl et al., 2012](http://dx.doi.org/10.1016/j.neuroimage.2012.01.021)] (`t1-freesurfer-longitudinal`),
    seven surfaces for each hemisphere were computed,
    ranging from 35% to 65% of the gray matter thickness.
    The partial volume corrected data were projected onto these meshes and
    the seven values were averaged, giving more weight to the vertices near the center of the cortex.
    Finally, the projected PET signal in the subject’s native space was
    spatially normalized to the standard space of FsAverage
    ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)).

!!! tip
    Easily access the papers cited on this page on
    [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/RGVVHC5W).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
