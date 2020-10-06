# `pet-volume` – Volume-based processing of PET images


This pipeline performs several processing steps on PET data in voxel space, which include:

- intra-subject registration of the PET image into the space of the subject’s T1-weighted MRI image using [SPM](http://www.fil.ion.ucl.ac.uk/spm/);
- (optional) partial volume correction (PVC) using the [PETPVC toolbox](https://github.com/UCL/PETPVC) [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)];
- inter-subject spatial normalization of the PET image into MNI space based on the DARTEL deformation model of SPM [[Ashburner, 2007](http://dx.doi.org/10.1016/j.neuroimage.2007.07.007)];
- intensity normalization using the average PET uptake in reference regions resulting in a standardized uptake value ratio (SUVR) map;
- parcellation into anatomical regions based on an atlas and computation of average values within each region. The list of available atlases can be found [here](../../Atlases).


## Prerequisite
You need to have performed the [`t1-volume`](../T1_Volume) pipeline on your T1-weighted MR images.


## Dependencies
<!--- If you installed the docker image of Clinica, nothing is required.-->

- If you only installed the core of Clinica, this pipeline needs the installation of **SPM12**. You can find how to install these software packages on the [third-party](../../Third-party) page.

- If you want to apply partial volume correction (PVC) on your PET data, you will need to install **PETPVC 1.2.4**, which depends on **ITK 4**. More information on the [third-party](../../Third-party) page.



## Running the pipeline
The pipeline can be run with the following command line:

```Text
clinica run pet-volume <bids_directory> <caps_directory> <group_id>
```
where:

- `bids_directory` is the input folder containing the dataset in a [BIDS](../../BIDS) hierarchy.
- `caps_directory` acts both as an input folder (where the results of the `t1-volume-*` pipeline are stored) and as the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.
- `group_id` is the ID of the group that is associated to the DARTEL template that you had created when running the `t1-volume-*` pipeline.

Pipeline options:

- `--pet_tracer`: type of PET image to process. Possible values are `fdg` and `av45`. Default value is `fdg`.
- `--smooth`: a list of integers specifying the different isomorphic full width at half maximum (FWHM) in millimeters to smooth the image. Default value is: 0, 8 (both without smoothing and with an isomorphic smoothing of 8 mm)
- `--pvc_fwhm`: TSV file containing the `fwhm_x`, `fwhm_y` and `fwhm_z` of the PSF for each PET image. More explanation below.



!!! note "Partial volume correction"
    To correct for [partial volume effects](http://www.turkupetcentre.net/petanalysis/image_pve.html), the pipeline uses the [region-based voxel-wise (RBV) correction](http://doc.pmod.com/pneuro/8893.htm) implemented in the [PETPVC toolbox](https://github.com/UCL/PETPVC).
    You need to specify in a TSV file the full width at half maximum (FWHM), in millimeters, of the [point spread function (PSF)](https://en.wikipedia.org/wiki/Point_spread_function) associated with your data, in the x, y and z directions. For instance, if the FWHM of the PSF associated with your first image is 8 mm along the x axis, 9 mm along the y axis, and 10 mm along z axis, the first row of your TSV file will look like this:
    ```
    participant_id    session_id     fwhm_x    fwhm_y    fwhm_z
    sub-CLNC0001      ses-M00        8         9         10
    sub-CLNC0002      ses-M00        7         6         5
    sub-CLNC0003      ses-M00        6         6         6
    ```

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `clinica run pet-volume --help` to see the full list of parameters.


## Outputs
Results are stored in the following folder of the [CAPS hierarchy](../../CAPS/Specifications/#pet-volume-volume-based-processing-of-pet-images): `subjects/sub-<participant_label>/ses-<session_label>/pet/preprocessing`.

The main output files are:

- `<source_file>_space-T1w[_pvc-rbv]_pet.nii.gz`: PET image registered into the T1w native space.

- `<source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<label>_mask-brain[_fwhm-<X>mm]_pet.nii.gz`: standard uptake value ratio (SUVR) PET image in MNI space, masked to keep only the brain, and optionally smoothed.

- `atlas_statistics/<source_file>_space-<space>[_pvc-rbv]_suvr-<label>_statistics.tsv`: TSV files summarizing the regional statistics on the labelled atlas <space\>.

!!! note
    The [pvc-rbv] label indicates whether the PET image has undergone partial value correction (region-based voxel-wise method) or not.

    The full list of output files from the pet-volume pipeline can be found in the [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#pet-volume-volume-based-processing-of-pet-images).


## Going further

- You can use standardized uptake value ratio (SUVR) maps to perform group comparison with the [`statistics-volume` pipeline](../Stats_Volume).
- You can perform classification based on [machine learning](../MachineLearning_Classification), as showcased in the [AD-ML framework](https://github.com/aramis-lab/AD-ML).


## Describing this pipeline in your paper

!!! cite "Example of paragraph:"
    Theses results have been obtained using the `pet-volume` pipeline of Clinica [[Routier et al](https://hal.inria.fr/hal-02308126/); [Samper et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)]. This pipeline first performs intra-subject registration of the PET image into the space of the subject’s T1-weighted MRI image using [SPM](http://www.fil.ion.ucl.ac.uk/spm/). [The PET image is corrected for partial volume effects using the [PETPVC toolbox](https://github.com/UCL/PETPVC) [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)]. The PET image is then spatially normalized into MNI space using the DARTEL deformation model of SPM, and intensity normalized using the average PET uptake in a reference region ([pons | pons + cerebellum]). Finally, the average PET uptake is computed for a set of regions obtained from different atlases in MNI space [Tzourio-Mazoyer et al., [2002](http://dx.doi.org/10.1006/nimg.2001.0978), [2015](http://dx.doi.org/10.1016/j.neuroimage.2015.07.075); [Joliot et al., 2015](http://dx.doi.org/10.1016/j.jneumeth.2015.07.013); [Hammers et al., 2003](http://dx.doi.org/10.1002/hbm.10123); [Gousias et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.11.034); [Shattuck et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.09.031); [CAT12](http://dbm.neuro.uni-jena.de/cat/)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/INDXD9QQ).

## Support

-   You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
-   Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
