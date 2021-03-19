# `pet-linear` - Affine registration of PET images to the MNI standard space

This pipeline performs a set of steps in order to affinely align PET images to the MNI space using the [ANTs](http://stnava.github.io/ANTs/) software package [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)]. These steps include: 

* affine registration to the [MNI152NLin2009cSym](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#template-based-coordinate-systems) template [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033), [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)] in MNI space with the SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)]; 
* intensity normalization using the average PET uptake in reference regions resulting in a standardized uptake value ratio (SUVR) map; 
* cropping of the registered images to remove the background.

## Prerequisites
You need to have performed the [`t1-linear`](../T1_Linear) pipeline on your T1-weighted MR images.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of **ANTs** on your computer.
You can find how to install this software package on the [third-party](../../Third-party) page.

## Running the pipeline

The pipeline can be run with the following command line:

```text
clinica run pet-linear <bids_directory> <caps_directory> <acq_label> <suvr_reference_region>
```

where:

- `bids_directory` is the input folder containing the dataset in a [BIDS](docs/BIDS) hierarchy;
- `caps_directory` is the output folder containing the results in a [CAPS](docs/CAPS) hierarchy;
- `acq_label` is the label given to the PET acquisition, specifying the tracer used (`acq-<acq_label>`). It can be for instance 'fdg' for fluorodésoxyglucose or 'av45' for florbetapir;
- `suvr_reference_region` is the reference region used to perform intensity normalization (i.e. dividing each voxel of the image by the average uptake in this region) resulting in a standardized uptake value ratio (SUVR) map. It can be `cerebellumPons` (used for amyloid tracers) or `pons` (used for FDG).

On default, cropped images (matrix size 169×208×179, 1 mm isotropic voxels) are generated to reduce the computing power required when training deep learning models. Use the `--uncropped_image` flag if you do not want to crop the image.

The pipelines also offers the possibility to save the PET image registered in the T1w space through a rigid transformation. This an intermediate step that is executed in the pipeline. Use the `--save_pet_in_t1w_space` flag to save this image.

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](docs/CAPS):
`subjects/<participant_id>/<session_id>/pet_linear`.

The main output files are:

* `<source_file>_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-<label>_pet.nii.gz`: PET image registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html) and cropped.
- `<source_file>_space-T1w_rigid.mat`: rigid transformation estimated with [ANTs](https://stnava.github.io/ANTs/).
- (optional) `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_pet.nii.gz`: PET image affinely registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html) and uncropped.
- (optional) `<source_file>_space-T1w_pet.nii.gz`: PET image affinely registered to the assiciated T1w image.

## Going further

You can now run the [`deeplearning-prepare-data` pipeline](../DeepLearning_PrepareData) to prepare images to be used with the PyTorch library
[[Paszke et al., 2019]](https://papers.nips.cc/paper/9015-pytorch-an-imperative-style-high-performance-deep-learning-library)
for classification based on deep learning using the [AD-DL framework](https://github.com/aramis-lab/AD-DL) presented in [[Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)].

## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    These results have been obtained using the `pet-linear` pipeline of Clinica
    [[Routier et al](https://hal.inria.fr/hal-02308126/);
    [Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)].
    More precisely, a rigid registration PET image to T1w image was performed using the SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)] from ANTs [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)]. Then, this registration was combined with the affine registration from `t1-linear` to align the PET images in MNI space with the ICBM 2009c nonlinear symmetric template  [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033), [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)]. Next image intensity has been normalized using the average PET uptake in reference regions resulting in a standardized uptake value ratio (SUVR) map. The registered and normalized images were further cropped to remove the background resulting in images of size 169×208×179, with 1 mm isotropic voxels.
!!! tip
   Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/RGVVHC5W).


