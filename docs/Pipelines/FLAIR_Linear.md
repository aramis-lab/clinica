# `flair-linear` - Affine registration of FLAIR images to the MNI standard space

This pipeline performs a set of steps in order to affinely align FLAIR
images to the MNI space using the [ANTs](http://stnava.github.io/ANTs/)
software package [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)].
These steps include: bias field correction using N4ITK
[[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)];
affine registration to the
[MNI152NLin2009cSym](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#template-based-coordinate-systems)
template [Fonov et al.,
[2011](https://doi.org/10.1016/j.neuroimage.2010.07.033),
[2009](https://doi.org/10.1016/S1053-8119(09)70884-5)] in MNI space with the
SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)];
cropping of the registered images to remove the background.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of
**ANTs** on your computer.
You can find how to install this software package on the [third-party](../../Third-party) page.

## Running the pipeline

The pipeline can be run with the following command line:

```Text
clinica run flair-linear [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a
[BIDS](../../BIDS) hierarchy.
- `CAPS_DIRECTORY` is the output folder containing the results in a
[CAPS](../../CAPS/Introduction) hierarchy.

On default, cropped images (matrix size 169×208×179, 1 mm isotropic voxels) are
generated to reduce the computing power required when training deep learning models.
Use the option `--uncropped_image` if you do not want to crop the image.

!!! note
    The arguments common to all Clinica pipelines are described in
    [Interacting with clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `clinica run flair-linear --help` to see the full
    list of parameters.

## Outputs

Results are stored in the following folder of the [CAPS
hierarchy](../../CAPS/Specifications/#flair-linear-affine-registration-of-flair-images-to-the-mni-standard-space):
`subjects/<participant_id>/<session_id>/flair_linear` with the following outputs:

- `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_FLAIR.nii.gz`:
FLAIR image affinely registered to the
[`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html).
- (optional)
  `<source_file>_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_FLAIR.nii.gz`: FLAIR
  image registered to the [`MNI152NLin2009cSym`
  template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html)
  and cropped.
- `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_affine.mat`:
affine transformation estimated with [ANTs](https://stnava.github.io/ANTs/).


## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    These results have been obtained using the `flair-linear` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675);
    [Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)].
    More precisely, bias field correction was applied using the N4ITK method
    [[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)].
    Next, an affine registration was performed using the SyN algorithm
    [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)]
    from ANTs [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)]
    to align each image to the MNI space with the ICBM 2009c nonlinear symmetric template
    [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033),
    [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)].
    (Optional) The registered images were further cropped to remove the background
    resulting in images of size 169×208×179, with 1 mm isotropic voxels.

!!! tip
    Easily access the papers cited on this page on
    [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/collections/8B2R2826).
