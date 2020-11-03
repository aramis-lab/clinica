# `t1-linear` - Affine registration of T1w images to the MNI standard space

This pipeline performs a set of steps in order to affinely align T1-weighted MR images to the MNI space using the [ANTs](http://stnava.github.io/ANTs/) software package [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)]. These steps include: bias field correction using N4ITK [[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)]; affine registration to the [MNI152NLin2009cSym](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#template-based-coordinate-systems) template [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033), [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)] in MNI space with the SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)]; cropping of the registered images to remove the background.

This pipeline was designed as a prerequisite for the [`deeplearning-prepare-data` pipeline](../DeepLearning_PrepareData) and deep learning classification algorithms presented in [[Wen et al., 2020](https://arxiv.org/abs/1904.07773)].

## Dependencies
If you only installed the core of Clinica, this pipeline needs the installation of **ANTs** on your computer. You can find how to install this software package on the [third-party](../../Third-party) page.

## Running the pipeline
The pipeline can be run with the following command line:
```Text
clinica run t1-linear <bids_directory> <caps_directory>
```
where:

- `bids_directory` is the input folder containing the dataset in a [BIDS](http://www.clinica.run/doc/BIDS/) hierarchy.
- `caps_directory` is the output folder containing the results in a [CAPS](http://www.clinica.run/doc/CAPS/Introduction) hierarchy.

On default, cropped images (matrix size 169×208×179, 1 mm isotropic voxels) are generated to reduce the computing power required when training deep learning models. Use `--uncropped_image` flag if you do not want to crop the image.

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).


## Outputs
Results are stored in the following folder of the [CAPS hierarchy](../../CAPS/Specifications/#t1-linear-affine-registration-of-t1w-images-to-the-mni-standard-space): `subjects/<participant_id>/<session_id>/t1_linear` with the following outputs:

- `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz`: T1w image affinely registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html).
- (optional) `<source_file>_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_T1w.nii.gz`: T1w image registered to the [`MNI152NLin2009cSym` template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html) and cropped.
- `<source_file>_space-MNI152NLin2009cSym_res-1x1x1_affine.mat`: affine transformation estimated with [ANTs](https://stnava.github.io/ANTs/).

## Going further

- You can now run the [`deeplearning-prepare-data` pipeline](../DeepLearning_PrepareData) to prepare images to be used with the PyTorch library [[Paszke et al., 2019]](https://papers.nips.cc/paper/9015-pytorch-an-imperative-style-high-performance-deep-learning-library) for classification based on deep learning using the [AD-DL framework](https://github.com/aramis-lab/AD-DL) presented in [[Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)].

## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    These results have been obtained using the `t1-linear` pipeline of Clinica [[Routier et al](https://hal.inria.fr/hal-02308126/); [Wen et al., 2020](https://doi.org/10.1016/j.media.2020.101694)]. More precisely, bias field correction was applied using the N4ITK method [[Tustison et al., 2010](https://doi.org/10.1109/TMI.2010.2046908)]. Next, an affine registration was performed using the SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)] from ANTs [[Avants et al., 2014](https://doi.org/10.3389/fninf.2014.00044)] to align each image to the MNI space with the ICBM 2009c nonlinear symmetric template  [Fonov et al., [2011](https://doi.org/10.1016/j.neuroimage.2010.07.033), [2009](https://doi.org/10.1016/S1053-8119(09)70884-5)]. (Optional) The registered images were further cropped to remove the background resulting in images of size 169×208×179, with 1 mm isotropic voxels.


!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/collections/8B2R2826).
