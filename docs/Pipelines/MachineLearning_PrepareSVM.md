# `machinelearning-prepare-spatial-svm` - Prepare input data for spatially regularized SVM

This pipeline allows the preparation of T1 MRI and PET data to perform classification with an SVM (support vector machine) with spatial and anatomical regularization [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)]. In this approach, the standard regularization of the SVM is replaced with a regularization that accounts for the spatial and anatomical structure of neuroimaging data. More specifically, it is regularized with respect to the tissue maps (GM, WM, CSF). As a result, the decision function learned by the algorithm will be more regular and anatomically interpretable.

## Prerequisites
You need to execute the [`t1-volume`](../T1_Volume) pipeline to run the pipeline on T1 MRI data, and the [`t1-volume`](../T1_Volume) + [`pet-volume`](../PET_Volume) pipelines to apply it to PET data.

## Dependencies
If you installed the core of Clinica, this pipeline needs no further dependencies.

## Running the pipeline
The pipeline can be run with the following command line:
```
clinica run machinelearning-prepare-spatial-svm caps_directory group_id
```
where:

- `caps_directory` is the output folder containing the results in a [CAPS](../../CAPS) hierarchy
- `group_id` is the user-defined identifier for the provided group of subjects

Pipeline options:

- `image_type`: possible options are `t1` and `pet` depending on the imaging modality considered. Default value: `t1`. When this flag is set to `pet`, the `pet_tracer` must be specified too (default value: `fdg`, other possible value: `av45`).

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `machinelearning-prepare-spatial-svm  --help` to see the full list of parameters.

## Outputs

Results are stored in the following folder of the [CAPS hierarchy](../../CAPS): `subjects/sub-<participant_label>/ses-<session_label>machine_learning/input_spatial_svm/group-<group_label>/` and `groups/group-<group_label>machine_learning/input_spatial_svm/`.

The main output files in `subjects` subfolder are:

  - `<source_file>_segm-[graymatter|whitematter|csf]_space-Ixi549Space_modulated-on_spatialregularization.nii.gz`: SVM regularization that accounts for the spatial and anatomical structure of neuroimaging data for gray matter, white matter or CSF maps.
  - `<source_file>_space-Ixi549Space[_pvc-rbv]_suvr-<label>_spatialregularization.nii.gz`: SVM regularization of PET data that accounts for the spatial and anatomical structure of neuroimaging data.


!!! note
    The full list of output files can be found in the [The ClinicA Processed Structure (CAPS) Specification](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8/edit#).

## Describing this pipeline in your paper

!!! cite "Example of paragraph"
    The classification was performed using a spatially regularized support vector machine (SVM), as proposed in [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)]. In this approach, the standard regularization of the SVM is replaced with a regularization that accounts for the spatial and anatomical structure of neuroimaging data. More specifically, we used the Fisher regularization and tissue maps (gray matter, white matter and cerebrospinal fluid) as spatial priors. The decision function of the SVM is made regular with respect to these tissues and is thus easier to interpret in terms of anatomical regions. To that purpose, feature maps were preprocessed using the `machinelearning-prepare-spatial-svm` pipeline of Clinica. The classification was then performed using an SVM on the preprocessed feature maps.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/1517933/aramis_clinica/items/collectionKey/78RQYITS).

## Advanced usage
The approach is general and can make use of different types of spatial and anatomical regularizations, introduce different types of spatial priors and varying amounts of regularization. These different aspects are described in details in [[Cuingnet et al, 2013](https://doi.org/10.1109/TPAMI.2012.142)]. Currently, this pipeline implements only one type of regularization (Fisher regularization), which is the most general one and should fit the vast majority of purposes.
As for the type of spatial prior, the pipeline currently only uses tissue maps (gray matter, white matter and cerebrospinal fluid). The decision function of the SVM is made regular with respect to these tissues. Other types of priors (such as atlases of anatomical regions) are currently not available and might be implemented in future releases. Finally, the amount of regularization can be changed using the  `fwhm` option. The default value is 4 mm. In practice, we found this value to be optimal. We therefore do not recommend to change it unless you have a specific reason to do so.
