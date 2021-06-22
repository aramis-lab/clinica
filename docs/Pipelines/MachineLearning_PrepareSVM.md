# `machinelearning-prepare-spatial-svm` - Prepare input data for spatially regularized SVM

This pipeline allows the preparation of T1-weighted MRI and PET data to perform classification with a support vector machine (SVM) with spatial and anatomical regularization [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)].
In this approach, the standard regularization of the SVM is replaced with a regularization that accounts for the spatial and anatomical structure of neuroimaging data.
More specifically, it is regularized with respect to the tissue maps (gray matter, white matter, cerebrospinal fluid [CSF]).
As a result, the decision function learned by the algorithm will be more regular and anatomically interpretable.

## Prerequisites

You need to execute the [`t1-volume`](../T1_Volume) pipeline to run the pipeline on T1-weighted MRI data, and the [`t1-volume`](../T1_Volume) + [`pet-volume`](../PET_Volume) pipelines to apply it to PET data.

## Dependencies

If you installed the core of Clinica, this pipeline needs no further dependencies.

## Running the pipeline

The pipeline can be run with the following command line:

```shell
clinica run machinelearning-prepare-spatial-svm [OPTIONS] CAPS_DIRECTORY GROUP_LABEL {t1-volume|pet-surface}
```

where:

- `CAPS_DIRECTORY` is the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy
- `GROUP_LABEL` is the user-defined identifier for the provided group of subjects
- The third positional argument can be `t1-volume` to use tissue maps or `pet-volume` to use standardized uptake value ratio (SUVR) maps.

Pipeline options if you use inputs from the `pet-volume` pipeline:

- `--acq_label`: name of the label given to the PET acquisition, specifying the tracer used (`acq-<acq_label>`).
- `--suvr_reference_region`: reference region used to perform intensity normalization
(i.e. dividing each voxel of the image by the average uptake in this region) resulting in a SUVR map.
It can be `cerebellumPons` (used for amyloid tracers) or `pons` (used for FDG).
- `--use_pvc_data`: use PET data with partial value correction (by default, PET data with no PVC are used)

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `machinelearning-prepare-spatial-svm --help` to see the full list of parameters.

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#machinelearning-prepare-spatial-svm-prepare-input-data-for-spatially-regularized-svm):
`subjects/<participant_id>/<session_id>/machine_learning/input_spatial_svm/<group_id>/`
and `groups/<group_id>/machine_learning/input_spatial_svm/`.

The main output files in the `subjects` subfolder are:

- `<source_file>_segm-{graymatter|whitematter|csf}_space-Ixi549Space_modulated-on_spatialregularization.nii.gz`:
SVM regularization that accounts for the spatial and anatomical structure of neuroimaging data for gray matter, white matter or CSF maps.
- `<source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<label>_spatialregularization.nii.gz`:
SVM regularization of PET data that accounts for the spatial and anatomical structure of neuroimaging data.

!!! note
    The full list of output files can be found in the [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#machinelearning-prepare-spatial-svm-prepare-input-data-for-spatially-regularized-svm).

## Going further

- You can now perform classification based on machine learning using the [AD-ML framework](https://github.com/aramis-lab/AD-ML).

## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    The classification was performed using a spatially regularized support vector machine (SVM), as proposed in [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)] and implemented in Clinica
    [[Routier et al](https://hal.inria.fr/hal-02308126/)].
    In this approach, the standard regularization of the SVM is replaced with a regularization that accounts for the spatial and anatomical structure of neuroimaging data.
    More specifically, we used the Fisher regularization and tissue maps (gray matter, white matter and cerebrospinal fluid) as spatial priors.
    The decision function of the SVM is made regular with respect to these tissues and is thus easier to interpret in terms of anatomical regions.
    To that purpose, feature maps were preprocessed using the `machinelearning-prepare-spatial-svm` pipeline of Clinica.
    The classification was then performed using an SVM on the preprocessed feature maps.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/1517933/aramis_clinica/items/collectionKey/78RQYITS).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).

## Advanced usage

The approach is general and can make use of different types of spatial and anatomical regularizations, introduce different types of spatial priors and varying amounts of regularization.
These different aspects are described in details in [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)].
Currently, this pipeline implements only one type of regularization (Fisher regularization), which is the most general one and should fit the vast majority of purposes.
As for the type of spatial prior, the pipeline currently only uses tissue maps (gray matter, white matter and CSF).
The decision function of the SVM is made regular with respect to these tissues.
Other types of priors (such as atlases of anatomical regions) are currently not available and might be implemented in future releases.
Finally, the amount of regularization can be changed using the `fwhm` option.
The default value is 4 mm.
In practice, we found this value to be optimal.
We therefore do not recommend to change it unless you have a specific reason to do so.
