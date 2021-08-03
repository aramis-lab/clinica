<!-- markdownlint-disable MD046-->
# `dwi-preprocessing-*` – Preprocessing of raw diffusion weighted imaging (DWI) datasets

This pipeline corrects DWI datasets for motion, eddy current, magnetic susceptibility and bias field distortions.

In particular, head-motion and eddy current corrections are performed with the `eddy` tool [[Andersson et al., 2016a](https://dx.doi.org/10.1016%2Fj.neuroimage.2015.10.019); [Andersson et al., 2016b](https://doi.org/10.1016/j.neuroimage.2016.06.058)] from **FSL** [[Jenkinson et al., 2011](https://doi.org/10.1016/j.neuroimage.2011.09.015)].
Depending on the data available (see below), the magnetic susceptibility correction can also be performed with the `eddy` tool or as a separate step.
Finally, bias field correction is performed using the N4 algorithm [[Tustison et al., 2010](https://dx.doi.org/10.1109/TMI.2010.2046908)] from **ANTs** [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)] by
computing a single multiplicative bias field from the corrected b0 image(s) implemented in **MRtrix3** [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)].

!!! note "Notes regarding the fieldmaps for the magnetic susceptibility correction"
    Depending on the type of extra acquisition available (see details in the 'Fieldmap data' section of the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#fieldmap-data)), different magnetic susceptibility corrections are performed.
    Two approaches are currently implemented in Clinica:

    - Phase difference image and at least one magnitude image (case 1 in the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#case-1-phase-difference-map-and-at-least-one-magnitude-image)).

    - No extra data: this is the case with the public [Alzheimer’s Disease Neuroimaging Initiative (ADNI)](http://adni.loni.usc.edu/) dataset for instance.
    The phase unwrapping is simulated thanks to a non-linear registration towards the T1-weighted image, which does not suffer from these artifacts.

    The cases 2, 3 and 4 of the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#fieldmap-data) aka the "two phase maps and two magnitude images",
    the "direct field mapping" (showing the field inhomogeneity in each voxel) and the "multiple phase encoded directions" (`topup`) approaches are not currently implemented in Clinica.

## Dependencies

<!--If you installed the docker image of Clinica, nothing is required.-->

If you only installed the core of Clinica, the `dwi-preprocessing-*` pipeline needs the installation of **ANTs**, **FSL**, and **MRtrix3** on your computer.
Extra installation of **Convert3D** will be needed for `dwi-preprocessing-using-t1` pipeline.
You can find how to install these software packages on the [third-party](../../Third-party) page.

## Running the pipeline

The pipeline can be run with the following command lines depending on the data you have:

```Text
clinica run dwi-preprocessing-using-t1 [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY
```

```Text
clinica run dwi-preprocessing-using-phasediff-fmap [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](../../BIDS) hierarchy.
- `CAPS_DIRECTORY` is the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.

Please note that you will need the `PhaseEncodingDirection` and `TotalReadoutTime` metadata fields in the JSON file associated to your DWI images
(see [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#diffusion-imaging-data) for more details).
For the `dwi-preprocessing-using-phasediff-fmap`pipeline, you will also need the  `EchoTime1` and `EchoTime2` metadata fields in the JSON file associated to your fieldmap images
(see [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#phase-difference-image-and-at-least-one-magnitude-image)).
Without these metadata fields, the pipelines will not run.

If you want to run the pipeline on a subset of your BIDS dataset, you can use the `-tsv` flag to specify in a TSV file the participants belonging to your subset.

!!! tip "Decreasing computation time"
    By default, the `eddy` tool of FSL uses OpenMP for parallel computing.
    CUDA can be used instead to speed up processing.
    To do so, see instructions in [FSL `eddy` wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#The_eddy_executables) for installation and configuration of `eddy_cuda`.
    To check that CUDA has been installed properly, type the `eddy_cuda` command.
    If this type of message appears:

    ```Text
    $ eddy_cuda
    eddy_cuda9.1: error while loading shared libraries: libcudart.so.9.1: cannot open shared object file: No such file or directory
    ```
    it means that something went wrong during the installation or configuration of CUDA.
    Otherwise, you should see the help of the command and you can now add the `--use_cuda` flag in Clinica!

!!! tip
    If your b0 images are not identical to 0 (e.g. 5 or 10, you can check this information by opening your `*.bval` file), you can use the `--low_bval` parameter to consider these images as b0s.

!!! note
    Please note that the `dwi-preprocessing-using-t1` pipeline will generate many temporary files.
    It can generate up to 20 GB of data for a dataset with 45 volumes.
    The size of the working directory is proportional to the number of directions available in the DWI dataset.

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#dwi-preprocessing-preprocessing-of-raw-diffusion-weighted-imaging-dwi-datasets):
`subjects/<participant_id>/<session_id>/dwi/preprocessing`.

The output files are:

- `<source_file>_space-{T1w|b0}_preproc.{bval|bvec|nii.gz}`:
corrected DWI dataset where the first volume of the dataset is the reference b0.
- `<source_file>_space-{T1w|b0}_brainmask.nii.gz`:
brain extracted image based on the reference b0.

!!! note
    If you ran the preprocessing pipeline using the T1-weighted image, your DWI dataset will be registered with the T1w image.
    Otherwise, your DWI dataset will be registered with the reference b0 image.

## Describing this pipeline in your paper

??? cite "Example of paragraph for the `dwi-preprocessing-using-t1` pipeline"
    These results have been obtained using the `dwi-preprocessing-using-t1` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675);
    [Wen et al., 2020](https://doi.org/10.1007/s12021-020-09469-5)].
    For each subject, all the b0 images were rigidly registered to the first b0 image and averaged to create the b0 reference.
    Then, the raw DWIs were corrected for eddy current induced distortions and subject movements by simultaneously modelling the effects of diffusion eddy currents and movements on the image.
    This step was performed using the FSL `eddy` tool
    [[Andersson et al., 2016a](https://dx.doi.org/10.1016%2Fj.neuroimage.2015.10.019)] from FSL [[Jenkinson et al., 2012](https://doi.org/10.1016/j.neuroimage.2011.09.015)] with the replace outliers (`--repol`) option [[Andersson et al., 2016b](https://doi.org/10.1016/j.neuroimage.2016.06.058)].
    To correct for susceptibility induced distortions, the T1w MRI was used as fieldmap data.
    The reference b0 image was first skull-stripped.
    To obtain the transformation flow from the native diffusion space to the T1w MRI space, the skull-stripped b0 image was registered to the T1w MRI in two steps:
    first a rigid registration using the FSL `flirt` tool and then a non-linear registration using the SyN algorithm implemented in ANTs [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)].
    SyN is an inverse-consistent registration algorithm allowing EPI induced susceptibility artifact correction [[Leow et al., 2007](https://doi.org/10.1109/TMI.2007.892646)].
    The resulting transformation was applied to the DWIs to correct for the susceptibility induced distortions and the diffusion weighting directions were appropriately updated [[Leemans & Jones, 2009](https://doi.org/10.1002/mrm.21890)].
    Lastly, the DWI volumes were corrected for nonuniform intensity using the ANTs N4 bias field correction algorithm [[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)] and
    implemented in MRtrix3 [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)] software.
    A single multiplicative bias field from the reference b0 image was estimated,
    as suggested in [[Jeurissen et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.07.061)].
    Average b0 image was finally computed on corrected DWI volume in order to extract brain mask with FSL `bet` [[Smith, 2002](https://doi.org/10.1002/hbm.10062)].

??? cite "Example of paragraph for the `dwi-preprocessing-using-phasediff-fmap` pipeline"
    These results have been obtained using the `dwi-preprocessing-using-fmap` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    For each subject, a brain mask was estimated from the DWI volumes
    [[Dhollander et al., 2016](https://www.researchgate.net/publication/307863133_Unsupervised_3-tissue_response_function_estimation_from_single-shell_or_multi-shell_diffusion_MR_data_without_a_co-registered_T1_image)]
    and a first run of the `eddy` tool [[Andersson et al., 2016a](https://dx.doi.org/10.1016%2Fj.neuroimage.2015.10.019)] from FSL [[Jenkinson et al., 2012](https://doi.org/10.1016/j.neuroimage.2011.09.015)] with the replace outliers (`--repol`) option [[Andersson et al., 2016b](https://doi.org/10.1016/j.neuroimage.2016.06.058)] was used to align all the b0 images, average them, and create the reference b0 image.
    The fieldmap image was calibrated with the FSL `fugue`/`prelude` tools [[Jenkinson, 2002](https://doi.org/10.1002/mrm.10354)] and registered to the reference b0 image.
    Head motion, eddy current and magnetic susceptibility distortions were simultaneously estimated and corrected using a second run of `eddy` with the calibrated fieldmap.
    Lastly, the DWI volumes were corrected for nonuniform intensity using the ANTs N4 bias field correction algorithm [[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)].
    A single multiplicative bias field from the corrected reference b0 image was estimated,
    as suggested in [[Jeurissen et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.07.061)] and
    implemented in MRtrix3 [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)] software.
    Average b0 image was finally computed on corrected DWI volume in order to extract brain mask with FSL `bet` [[Smith, 2002](https://doi.org/10.1002/hbm.10062)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/BJV73LU7).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
