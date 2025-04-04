<!-- markdownlint-disable MD033 MD046-->
# `pet-linear` - Linear processing of PET images

This pipeline performs spatial normalization to the MNI space and intensity normalization of [PET](../glossary.md#pet) images.

Its steps include:

- clipping of the background to remove noisy values which can deteriorate further processing.
- affine registration to the [MNI152NLin2009cSym](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#template-based-coordinate-systems) template [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033), [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)] in MNI space with the SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)] from the [ANTs](http://stnava.github.io/ANTs/) software package [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)]
- intensity normalization using the average [PET](../glossary.md#pet) uptake in reference regions resulting in a standardized uptake value ratio ([SUVR](../glossary.md#suvr)) map
- cropping of the registered images to remove the background

!!! note "Clinica & BIDS specifications for PET modality"
    Since Clinica `v0.6`, PET data following the official specifications in BIDS version 1.6.0 are now compatible with Clinica.
    See [BIDS](../../BIDS) page for more information.

## Prerequisites

You need to have performed the [`t1-linear`](./T1_Linear.md) pipeline on your T1-weighted MR images.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of [ANTs](../Software/Third-party.md#ants) on your computer.

## Running the pipeline

The pipeline can be run with the following command line:

```shell
clinica run pet-linear [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY ACQ_LABEL
                       {pons|cerebellumPons|pons2|cerebellumPons2}
```

where:

--8<-- "snippets/cmd_inputs.md:bids_caps"
--8<-- "snippets/cmd_inputs.md:acq"
--8<-- "snippets/cmd_inputs.md:region"

with specific options :

--8<-- "snippets/pipelines_options.md:linear"
--8<-- "snippets/pipelines_options.md:pet_recon"
- `--save_pet_in_t1w_space`: Use this option to save the [PET](../glossary.md#pet) image in the T1w space after rigid transformation.

??? info "Optional parameters common to all pipelines"
    --8<-- "snippets/pipelines_options.md:all"

--8<-- "snippets/known_issues.md:several_pet"

!!! tip
    Do not hesitate to type `clinica run pet-linear --help` to see the full list of parameters.

## Outputs

Results are stored in the following folder of the [CAPS hierarchy](../CAPS/Specifications.md#pet-imaging-data): `subjects/<participant_id>/<session_id>/pet_linear`.

The main output files are:

- `<source_file>_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-<label>_pet.nii.gz`: [PET](../glossary.md#pet) [SUVR](../glossary.md#suvr) image registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html) and cropped.
- `<source_file>_space-T1w_rigid.mat`: rigid transformation between the [PET](../glossary.md#pet) and T1w images estimated with [ANTs](https://stnava.github.io/ANTs/).
- (optional) `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_pet.nii.gz`: [PET](../glossary.md#pet) [SUVR](../glossary.md#suvr) image affinely registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html) (i.e. not cropped).
- (optional) `<source_file>_space-T1w_pet.nii.gz`: [PET](../glossary.md#pet) image affinely registered to the associated T1w image.

## Going further

You can now use the [ClinicaDL framework](https://clinicadl.readthedocs.io/) presented in [[Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)] for classification based on deep learning methods.

## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    These results have been obtained using the `pet-linear` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    This pipeline first performs intra-subject rigid registration of the PET
    image into the space of the subject’s T1-weighted (T1w) MR image using the
    SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)]
    from ANTs [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)].
    The PET to T1w transformation is then composed with T1w to ICBM 2009c
    nonlinear symmetric template affine transformation, obtained with
    `t1-linear`, to transport the PET image to the MNI space [Fonov et al.,
    [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033),
    [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)].
    The PET image is further intensity normalized using the average PET uptake
    in a reference region ([pons | pons + cerebellum]), resulting in a
    standardized uptake value ratio ([SUVR](../glossary.md#suvr)) map.
    The [PET](../glossary.md#pet) [SUVR](../glossary.md#suvr) image in MNI space is finally cropped to remove the background,
    resulting in images of size 169×208×179 with 1 mm isotropic voxels.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/collections/8AEDUMZB).
