# Classification based on machine learning

Clinica provides a modular way to perform classification based on machine learning. To build its own classification pipeline, the user can combine three modules based on [scikit-learn](http://scikit-learn.org/stable/index.html):

  - Input (e.g. gray matter maps obtained from T1-weighted MRI images, FDG PET images)
  - Algorithm (e.g. support vector machine, logistic regression, random forest)
  - Validation (e.g. K-fold cross validation, repeated K-fold cross validation, repeated hold-out validation)

## Prerequisites
You need to have performed the [`t1-volume`](../T1_Volume) pipeline on your T1-weighted MRI images and/or the [`pet-volume`](../PET_Volume) pipeline on your PET images.


## Dependencies
If you installed the core of Clinica, this pipeline needs no further dependencies.


## Classification modules

### Input
Two classes corresponding to the voxel-based and the region-based approaches are implemented in `input.py`:

  - `CAPSVoxelBasedInput`: all the voxel of the image are used as features.
  - `CAPSRegionBasedInput`: a list of values stored in a TSV file is used as features. This list corresponds to PET or T1 image intensities averaged over a set of regions obtained from a brain parcellation when running the SPM pipeline.

!!! note
    The atlases that can be used for the region-based approaches are listed [here](../../Atlases).

### Algorithm
Three classes corresponding to the machine learning-based classification algorithms are implemented in `algorithm.py`:

  - `DualSVMAlgorithm`: support vector machine (SVM) algorithm (input: all the data available or a kernel that can be pre-computed)
  - `LogisticReg`: logistic regression algorithm (input: all the data available)
  - `RandomForest`: random forest algorithm (input: all the data available)

Each algorithm implements a grid search approach to choose the best parameters for the classification by looking at the value of the balanced accuracy. The area under curve (AUC) is also reported. The labels are automatically assigned based on the `diagnoses_tsv` file.


### Validation
Three classes corresponding to the validation strategies are implemented in `validation.py`:

  - `KFoldCV`: K-fold cross validation
  - `RepeatedKFoldCV`: repeated K-fold cross validation
  - `RepeatedHoldOut`: repeated hold-out validation

The input is the name of the classification algorithm used.


## Running your pipeline
No matter the combination of modules chosen, the inputs necessary are:

- `caps_directory`: the folder containing the results of the SPM pipeline (where TSV files are stored)
- `subjects_visits_tsv`: the TSV file containing the participant_id and the session_id
- `diagnoses_tsv`: a TSV file where the diagnosis for each participant (identified by a participant ID) is reported (i.e. AD, CN). It allows the algorithm to perform the dual classification (between the two labels reported).
Example of a diagnosis TSV file:
````
participant_id    diagnosis
sub-CLNC0001      AD
sub-CLNC0002      CN
sub-CLNC0003      AD
sub-CLNC0004      AD
sub-CLNC0005      CN
````
- `group_id`: the ID of the group of subjects studied
- `image_type`: a value to set the modality studied ("T1" or "fdg")
- `output_dir`: the directory where outputs are saved
- `atlas`: the name of the atlas used for the brain parcellation in case of a region-based approach
- `fwhm`: the FWHM value in mm used in the SPM pipeline
- `modulated`: a flag to indicate if when running the SPM pipeline the image has been modulated or not ("on", "off")
- `pvc`: type of PVC if used during the preprocessing of the PET images (e.g. "RBV")
- `precomputed_kernel`: to load the precomputed kernel if it exists
- `mask_zeros`: a flag to indicate if zero-valued voxels should be taken into account for the classification ("True", "False")
- `n_iterations`: number of times a task is repeated
- `grid_search_folds`: number of folds to use for the hyper-parameter grid search (e.g. 10)
- `c_range`: range used to select the best value for the C parameter, in the logspace
- `n_threads`: number of threads used if run in parallel
- `test_size`: percentage (between 0 and 1) representing the size of the test set for each shuffle split
- `balanced`:  option to balance the weights according to the number of samples
- `penalty`: type of penalty ("l2" or "l1")

!!! tip
    Usage examples are available in `ml_workflows.py`.


## Output
Results are saved in the output folder following this hierarchy:
```
└── <image-type>
    ├── region_based
    |    └── atlas-<atlas-id>
    |        └── <machine-learning-algorithm>
    |             └── <task1>_vs_<task2>
    |                 ├── classifier
    |                 |    └── iteration-<iteration-number>
    |                 |        ├── mean_results.tsv
    |                 |        ├── results.tsv
    |                 |        └── subjects.tsv
    |                 ├── best_parameters.json
    |                 ├── dual_coefficients.txt
    |                 ├── intersect.txt
    |                 ├── support_vector_indices.json
    |                 ├── weights.nii.gz
    |                 └── weights.txt
    └── voxel_based
        └── smoothing-<fwhm>
            └── <machine-learning-algorithm>
                └── <task1>_vs_<task2>
                    ├── classifier
                    |    └── iteration-<number-iteration>
                    |        ├── mean_results.tsv
                    |        ├── results.tsv
                    |        └── subjects.tsv
                    ├── best_parameters.json
                    ├── dual_coefficients.txt
                    ├── intersect.txt
                    ├── support_vector_indices.json
                    ├── weights.nii.gz
                    └── weights.txt
```

If `image_type` is `fdg`:
```
└── <image-type>
    └── region_based/voxel_base
        └── pvc-<pvc>
            └── ...
```


## Describing this pipeline in your paper

!!! cite "Example of paragraph:"
		These results have been obtained using the machine learning-based classification modules of Clinica. Clinica provides a modular way to perform classification based on machine learning by combining different inputs (e.g. gray matter maps obtained from T1-weighted MRI images, FDG PET images), algorithms (e.g. support vector machine, logistic regression, random forest) and validation strategies (e.g. K-fold cross validation, repeated K-fold cross validation, repeated hold-out validation). These modules rely on [scikit-learn](http://scikit-learn.org/stable/index.html).

## Support

-   You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
-   Report an issue on [GitLab](https://gitlab.icm-institute.org/aramislab/clinica/issues).
