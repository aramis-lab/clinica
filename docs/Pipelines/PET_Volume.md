<!-- markdownlint-disable MD046 -->
# `pet-volume` – Volume-based processing of PET images

This pipeline performs several processing steps on [PET](../glossary.md#pet) data in voxel space, which include:

- intra-subject registration of the PET image into the space of the subject’s T1-weighted MR image using [SPM](http://www.fil.ion.ucl.ac.uk/spm/);
- (optional) partial volume correction ([PVC](../glossary.md#pvc)) using the [PETPVC toolbox](https://github.com/UCL/PETPVC) [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)];
- inter-subject spatial normalization of the PET image into MNI space based on the DARTEL deformation model of SPM [[Ashburner, 2007](http://dx.doi.org/10.1016/j.neuroimage.2007.07.007)];
- intensity normalization using the average PET uptake in reference regions resulting in a standardized uptake value ratio ([SUVR](../glossary.md#suvr)) map;
- parcellation into anatomical regions based on an atlas and computation of average values within each region.
The list of available atlases can be found [here](../Atlases.md).

!!! note "Clinica & BIDS specifications for PET modality"
    Since Clinica `v0.6`, PET data following the official specifications in BIDS version `1.6.0` are now compatible with Clinica.
    See [BIDS](../BIDS.md) page for more information.

## Prerequisite

You need to have performed the [`t1-volume`](./T1_Volume.md) pipeline on your T1-weighted MR images.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of either [SPM12](../Software/Third-party.md#spm12) and [Matlab](../Software/Third-party.md#matlab), or [spm standalone](../Software/Third-party.md#spm12-standalone).

In addition, if you want to apply partial volume correction (PVC) on your PET data, you will need to install [PETPVC 1.2.4](../Software/Third-party.md#petpvc), which depends on [ITK 4](../Software/Third-party.md#itk).

## Running the pipeline

The pipeline can be run with the following command line:

```shell
clinica run pet-volume [OPTIONS] BIDS_DIRECTORY CAPS_DIRECTORY GROUP_LABEL ACQ_LABEL
                       {pons|cerebellumPons|pons2|cerebellumPons2} 
```

where:

--8<-- "snippets/cmd_inputs.md:bids_caps"
- `GROUP_LABEL` is the label of the group that is associated to the DARTEL template that you had created when running the [`t1-volume`](./T1_Volume.md) pipeline.
- `ACQ_LABEL` is the label given to the PET acquisition, specifying the tracer used (`trc-<acq_label>`).
- The reference region is used to perform intensity normalization (i.e. dividing each voxel of the image by the average uptake in this region) resulting in a standardized uptake value ratio ([SUVR](../glossary.md#suvr)) map.
It can be `cerebellumPons` or `cerebellumPons2` (used for amyloid tracers) or `pons` or `pons2` (used for FDG).

with specific options : 

- `--pvc_psf_tsv`: TSV file containing the `psf_x`, `psf_y` and `psf_z` of the PSF for each PET image. More explanation is given in [PET Introduction](./PET_Introduction.md#partial-volume-correction-pvc) page.
    
    ??? info "Clinica v0.3.8+"
        Since the release of Clinica v0.3.8, the handling of PSF information in the TSV file has changed: `fwhm_x`, `fwhm_y`, `fwhm_z` columns have been replaced by `psf_x`, `psf_y`, `psf_z` and the `acq_label` column has been added.
        Additionally, the [SUVR](../glossary.md#suvr) reference region is now a compulsory argument: it will be easier for you to modify Clinica if you want to add a custom reference region ([PET Introduction](../PET_Introduction) page).
        Choose `cerebellumPons` for amyloid tracers or `pons` for FDG to have the previous behavior.

- `--smooth`: a list of integers specifying the different isotropic full width at half maximum ([FWHM](../glossary.md#fwhm)) in millimeters to smooth the image. Default value is: 0, 8 (both without smoothing and with an isotropic smoothing of 8 mm)
- `--reconstruction_method`: Select only images based on a [specific reconstruction method](./PET_Introduction.md#reconstruction-methods).

??? info "Optional parameters common to all pipelines"
    --8<-- "snippets/pipelines_options.md:all"

--8<-- "snippets/known_issues.md:several_pet"

!!! tip
    Do not hesitate to type `clinica run pet-volume --help` to see the full list of parameters.

## Outputs

Results are stored in the following folder of the [CAPS hierarchy](../CAPS/Specifications.md#pet-volume---volume-based-processing-of-pet-images): `subjects/<participant_id>/<session_id>/pet/preprocessing`.

The main output files are:

- `<source_file>_space-T1w[_pvc-rbv]_pet.nii.gz`: PET image registered into the T1-weighted MRI native space.

- `<source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<label>_mask-brain[_fwhm-<X>mm]_pet.nii.gz`: standard uptake value ratio ([SUVR](../glossary.md#suvr)) [PET](../glossary.md#pet) image in MNI space, masked to keep only the brain, and optionally smoothed.

- `atlas_statistics/<source_file>_space-<space>[_pvc-rbv]_suvr-<label>_statistics.tsv`: [TSV](../glossary.md#tsv) files summarizing the regional statistics on the labelled atlas `<space>`.

!!! note
    The `[_pvc-rbv]` label indicates whether the [PET](../glossary.md#pet) image has undergone partial value correction (region-based voxel-wise method) or not.

    The full list of output files from the pet-volume pipeline can be found in the [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#pet-volume-volume-based-processing-of-pet-images).

## Going further

- You can use standardized uptake value ratio ([SUVR](../glossary.md#suvr)) maps to perform group comparison with the [`statistics-volume` pipeline](../Stats_Volume).
- You can perform classification based on [machine learning](../MachineLearning_Classification), as showcased in the [AD-ML framework](https://github.com/aramis-lab/AD-ML).

## Describing this pipeline in your paper

!!! cite "Example of paragraph:"
    These results have been obtained using the `pet-volume` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675);
    [Samper et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].
    This pipeline first performs intra-subject registration of the PET image into
    the space of the subject’s T1-weighted MR image using
    [SPM](http://www.fil.ion.ucl.ac.uk/spm/).
    [The PET image is corrected for partial volume effects using the
    [PETPVC toolbox](https://github.com/UCL/PETPVC)
    [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)].]
    The PET image is then spatially normalized into MNI space using the DARTEL deformation model of SPM, and intensity normalized using the average PET uptake in a reference region ([pons | pons + cerebellum]).
    Finally, the average PET uptake is computed for a set of regions obtained from
    different atlases in MNI space [Tzourio-Mazoyer et al.,
    [2002](http://dx.doi.org/10.1006/nimg.2001.0978),
    [2015](http://dx.doi.org/10.1016/j.neuroimage.2015.07.075);
    [Joliot et al., 2015](http://dx.doi.org/10.1016/j.jneumeth.2015.07.013);
    [Hammers et al., 2003](http://dx.doi.org/10.1002/hbm.10123);
    [Gousias et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.11.034);
    [Shattuck et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.09.031);
    [CAT12](http://dbm.neuro.uni-jena.de/cat/)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/INDXD9QQ).
