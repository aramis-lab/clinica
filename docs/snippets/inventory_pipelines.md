
- Anatomical MRI
    - `t1-linear` - [Linear processing of T1w MR images](./../Pipelines/T1_Linear): affine registration to the MNI standard space of T1 images
    - `flair-linear` - [Linear processing of FLAIR images](./../Pipelines/FLAIR_Linear.md): affine registration to the MNI standard space of FLAIR images
    - `t1-volume` - [Processing of T1w MR images using SPM](./../Pipelines/T1_Volume.md): tissue segmentation and spatial normalization
    - `t1-freesurfer` - [Processing of T1w MR images using FreeSurfer](./../Pipelines/T1_FreeSurfer.md): cortical surface, subcortical structures and volumetrics
    - `t1-freesurfer-longitudinal` - [Longitudinal processing of T1w MR images using FreeSurfer](./../Pipelines/T1_FreeSurfer_Longitudinal.md): cortical surface, subcortical structures and volumetrics
- Diffusion MRI (DWI)
    - `dwi-preprocessing-*` - [DWI pre-processing](./../Pipelines/DWI_Preprocessing.md): correction of head motion, magnetic susceptibility, eddy current and bias field induced distortions
    - `dwi-dti` - [DTI scalar maps (FA, MD, AD, RD) and spatial normalization](./../Pipelines/DWI_DTI.md): extraction of DTI-based measures (FA, MD, AD, RD)
    - `dwi-connectome` - [Construction of structural connectome](./../Pipelines/DWI_Connectome.md): computation of fiber orientation distributions, tractogram and connectome

- PET
    - [Introduction to concepts used in the PET pipelines](./../Pipelines/PET_Introduction.md): partial volume correction and standardized uptake value ratio (SUVR) map computation
    - `pet-linear` - [Linear processing of PET images](./../Pipelines/PET_Linear.md): affine registration to the MNI standard space and intensity normalization
    - `pet-volume` - [Volume-based processing of PET images](./../Pipelines/PET_Volume.md): registration to T1w MRI, intensity normalization, partial volume correction and spatial normalization
    - `pet-surface` - [Surface-based processing of PET images](./../Pipelines/PET_Surface.md): projection of the PET signal onto the subject’s cortical surface
    - `pet-surface-longitudinal` - [Surface-based longitudinal processing of PET images](./../Pipelines/PET_Surface_Longitudinal.md): projection of the PET signal onto the subject’s cortical surface

- Statistics
    - `statistics-surface` - [Surface-based mass-univariate analysis with SurfStat](./../Pipelines/Stats_Surface.md)
    - `statistics-volume` - [Volume-based mass-univariate analysis with SPM](./../Pipelines/Stats_Volume.md)

- Machine Learning
    - `machinelearning-prepare-spatial-svm` - [Prepare input data for spatially regularized SVM](./../Pipelines/MachineLearning_PrepareSVM.md)
    - `machinelearning-classification` - [Classification based on machine learning](./../Pipelines/MachineLearning_Classification.md)

- Deep learning
    - You can use the [ClinicaDL framework](https://clinicadl.readthedocs.io/) for the reproducible processing of neuroimaging data with deep learning methods.
