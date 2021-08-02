<!-- markdownlint-disable MD046-->
# `dwi-preprocessing-*` – Preprocessing of raw diffusion weighted imaging (DWI) datasets

This pipeline corrects DWI datasets for motion, eddy current, magnetic susceptibility and bias field distortions, assuming that the data have been acquired using an EPI sequence.

To that aim, a reference b0 image is computed: if only one b0 is present, it will act as the b0 reference, otherwise, the reference will correspond to the average of the different b0s after registration to the first b0.
Head-motion correction is performed by rigidly registering each DWI volume to the reference b0 and updating the b-vectors accordingly.
These rigid transformations initialize the affine registrations between each DWI volume and the reference b0, for eddy current correction.
The magnetic susceptibility correction is then performed using different strategies depending on the data available (see below), and all these corrections are combined into a single deformation field to minimize the number of interpolations.
Finally, bias field correction is performed using the ANTs N4 bias correction algorithm.

!!! note "Notes concerning the fieldmaps for the magnetic susceptibility corrections"
    Depending on the type of extra acquisition available (see details e.g. in the fieldmap section of the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#fieldmap-data)), different magnetic susceptibility corrections are performed.
    Two approaches are currently implemented in Clinica:

    - No extra data: this is the case with the public [Alzheimer’s Disease Neuroimaging Initiative (ADNI)](http://adni.loni.usc.edu/) dataset for instance.
    The phase unwrapping is simulated thanks to a non-linear registration towards the T1-weighted image, which does not suffer from these artifacts.

    - Phase difference image and at least one magnitude image (case 1 in the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#fieldmap-data))

    The cases 2, 3 and 4 of the [BIDS specifications](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#fieldmap-data) aka the “two phase images and two magnitude images”, the “single, real fieldmap image” (showing the field inhomogeneity in each voxel) and the “multiple phase encoded directions” (e.g. as implemented in FSL `topup`) approaches are currently not implemented in Clinica.

## Dependencies

<!--If you installed the docker image of Clinica, nothing is required.-->

If you only installed the core of Clinica, this pipeline needs the installation of
**ANTs v2.3.1** and **FSL 6.0** on your computer.
You can find how to install these software packages on the [third-party](../../Third-party) page.

The pipeline can be run with the following command lines depending on the type of  magnetic susceptibility correction performed:

```Text
clinica run dwi-preprocessing-using-t1 [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY
```

```Text
clinica run dwi-preprocessing-using-phasediff-fieldmap [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](../../BIDS) hierarchy.
- `CAPS_DIRECTORY` is the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.

If you want to run the pipeline on a subset of your BIDS dataset, you can use the `-tsv` flag to specify in a TSV file the participants belonging to your subset.

!!! tip
    If your b0 images are not identical to 0 (e.g. 5 or 10, you can check this information by opening your `*.bval` file), you can use the `--low_bval` parameter to consider these images as b0s.

!!! warning
    Please note that this pipeline will generate many temporary files.
    It can generate 3-4 GB of data for a dataset with 70 volumes.
    The more directions you have in the DWI dataset, the bigger the working directory will be.

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#dwi-preprocessing-preprocessing-of-raw-diffusion-weighted-imaging-dwi-datasets):
`subjects/<participant_id>/<session_id>/dwi/preprocessing`.

The output files are:

- `<source_file>_space-{T1w|b0}_preproc.{bval|bvec|nii.gz}`:
corrected DWI dataset where the first volume of the dataset is the reference b0.
- `<source_file>_brainmask.nii.gz`:
brain extracted image based on the reference b0.

!!! note
    If you ran the preprocessing pipeline using the T1-weighted image, your DWI dataset will be registered with the T1w image.
    Otherwise, your DWI dataset will be registered with the reference b0 image.

## Describing this pipeline in your paper

??? cite "Example of paragraph for the `dwi-preprocessing-using-t1` pipeline"
    These results have been obtained using the `dwi-preprocessing` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    For each subject, all raw DWI volumes were rigidly registered (6 degrees of freedom (dof)) to the reference b0 image (DWI volume with no diffusion sensitization) to correct for head motion.
    The diffusion weighting directions were appropriately updated [[Leemans & Jones, 2009](https://doi.org/10.1002/mrm.21890)].
    An affine registration (12 dof) was then performed between each DWI volume and
    the reference b0 to correct for eddy current distortions.
    These registrations were done using the FSL flirt tool (www.fmrib.ox.ac.uk/fsl).
    To correct for echo-planar imaging (EPI) induced susceptibility artifacts, the skull-stripped b0 images were registered to the T1-weighted image using the ANTs SyN registration algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004); [Leow et al., 2007](https://doi.org/10.1109/TMI.2007.892646)].
    The resulting deformation fields were then applied to the DWI volumes to align them with the T1 image.
    Finally, the DWI volumes were corrected for nonuniform intensity using the ANTs N4 bias correction algorithm [[Tustison et al.,, 2010](https://doi.org/10.1109/TMI.2010.2046908)].
    A single multiplicative bias field from the reference b0 image was estimated, as suggested in [[Jeurissen et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.07.061)].

??? cite "Example of paragraph for the `dwi-preprocessing-using-phasediff-fieldmap` pipeline"
    These results have been obtained using the `dwi-preprocessing` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    For each subject, all raw DWI volumes were rigidly registered (6 degrees of freedom (dof)) to the reference b0 image (DWI volume with no diffusion sensitization) to correct for head motion.
    The diffusion weighting directions were appropriately updated
    [[Leemans & Jones, 2009](https://doi.org/10.1002/mrm.21890)].
    An affine registration (12 dof) was then performed between each DWI volume and
    the reference b0 to correct for eddy current distortions.
    These registrations were done using the FSL flirt tool (www.fmrib.ox.ac.uk/fsl).
    To correct for echo-planar imaging (EPI) induced susceptibility artifacts, the fieldmap image was used as proposed by [[Jezzard and Balaban, 1995](https://doi.org/10.1002/mrm.1910340111)] with the FSL prelude/fugue tools.
    Finally, the DWI volumes were corrected for nonuniform intensity using the ANTs N4 bias correction algorithm [[Tustison et al.,, 2010](https://doi.org/10.1109/TMI.2010.2046908)].
    A single multiplicative bias field from the reference b0 image was estimated, as suggested in [[Jeurissen et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.07.061)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/BJV73LU7).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
