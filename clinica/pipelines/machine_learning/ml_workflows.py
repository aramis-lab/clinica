import numpy as np

from clinica.pipelines.machine_learning import algorithm, base, input, validation


class VoxelBasedKFoldDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        output_dir,
        fwhm=0,
        modulated="on",
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        mask_zeros=True,
        n_threads=15,
        n_folds=10,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(VoxelBasedKFoldDualSVM, self).__init__(
            input.CAPSVoxelBasedInput,
            validation.KFoldCV,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class VoxelBasedRepKFoldDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        output_dir,
        fwhm=0,
        modulated="on",
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        mask_zeros=True,
        n_threads=15,
        n_iterations=100,
        n_folds=10,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(VoxelBasedRepKFoldDualSVM, self).__init__(
            input.CAPSVoxelBasedInput,
            validation.RepeatedKFoldCV,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class VoxelBasedRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        output_dir,
        fwhm=0,
        modulated="on",
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        mask_zeros=True,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super().__init__(
            input.CAPSVoxelBasedInput,
            validation.RepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class VertexBasedRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        output_dir,
        image_type="PET",
        acq_label=None,
        suvr_reference_region=None,
        fwhm=20,
        precomputed_kernel=None,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-10, 2, 1000),
        splits_indices=None,
    ):

        super(VertexBasedRepHoldOutDualSVM, self).__init__(
            input.CAPSVertexBasedInput,
            validation.RepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class RegionBasedRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(RegionBasedRepHoldOutDualSVM, self).__init__(
            input.CAPSRegionBasedInput,
            validation.RepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class RegionBasedRepHoldOutLogisticRegression(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(RegionBasedRepHoldOutLogisticRegression, self).__init__(
            input.CAPSRegionBasedInput,
            validation.RepeatedHoldOut,
            algorithm.LogisticReg,
            locals(),
            output_dir,
        )


class RegionBasedRepHoldOutRandomForest(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        n_estimators_range=(100, 200, 400),
        max_depth_range=[None],
        min_samples_split_range=[2],
        max_features_range=("auto", 0.25, 0.5),
        splits_indices=None,
    ):

        super(RegionBasedRepHoldOutRandomForest, self).__init__(
            input.CAPSRegionBasedInput,
            validation.RepeatedHoldOut,
            algorithm.RandomForest,
            locals(),
            output_dir,
        )


class RegionBasedLearningCurveRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        n_learning_points=10,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
    ):

        super(RegionBasedLearningCurveRepHoldOutDualSVM, self).__init__(
            input.CAPSRegionBasedInput,
            validation.LearningCurveRepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class VoxelBasedLearningCurveRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        output_dir,
        fwhm=0,
        modulated="on",
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        mask_zeros=True,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        n_learning_points=10,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
    ):

        super(VoxelBasedLearningCurveRepHoldOutDualSVM, self).__init__(
            input.CAPSVoxelBasedInput,
            validation.LearningCurveRepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class RegionBasedRepKFoldDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        n_folds=10,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(RegionBasedRepKFoldDualSVM, self).__init__(
            input.CAPSRegionBasedInput,
            validation.RepeatedKFoldCV,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class CAPSTsvRepHoldOutDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        dataset,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(CAPSTsvRepHoldOutDualSVM, self).__init__(
            input.CAPSTSVBasedInput,
            validation.RepeatedHoldOut,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


class CAPSTsvRepHoldOutRandomForest(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        atlas,
        dataset,
        output_dir,
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        n_threads=15,
        n_iterations=100,
        test_size=0.3,
        grid_search_folds=10,
        balanced=True,
        n_estimators_range=(100, 200, 400),
        max_depth_range=[None],
        min_samples_split_range=[2],
        max_features_range=("auto", 0.25, 0.5),
        splits_indices=None,
    ):

        super(CAPSTsvRepHoldOutRandomForest, self).__init__(
            input.CAPSTSVBasedInput,
            validation.RepeatedHoldOut,
            algorithm.RandomForest,
            locals(),
            output_dir,
        )


# SVM reg


class VoxelBasedREGRepKFoldDualSVM(base.MLWorkflow):
    def __init__(
        self,
        caps_directory,
        subjects_visits_tsv,
        diagnoses_tsv,
        group_label,
        image_type,
        output_dir,
        fwhm=0,
        modulated="on",
        acq_label=None,
        suvr_reference_region=None,
        use_pvc_data=False,
        precomputed_kernel=None,
        mask_zeros=True,
        n_threads=15,
        n_iterations=100,
        n_folds=10,
        test_size=0.1,
        grid_search_folds=10,
        balanced=True,
        c_range=np.logspace(-6, 2, 17),
        splits_indices=None,
    ):

        super(VoxelBasedREGRepKFoldDualSVM, self).__init__(
            input.CAPSTSVBasedInput,
            validation.RepeatedKFoldCV,
            algorithm.DualSVMAlgorithm,
            locals(),
            output_dir,
        )


# TSV


class TsvRepHoldOutRandomForest(base.MLWorkflow):
    def __init__(
        self,
        data_tsv,
        columns,
        output_dir,
        n_threads=20,
        n_iterations=250,
        test_size=0.2,
        grid_search_folds=10,
        balanced=True,
        n_estimators_range=(100, 200, 400),
        max_depth_range=[None],
        min_samples_split_range=[2],
        max_features_range=("auto", 0.25, 0.5),
        splits_indices=None,
        inner_cv=False,
    ):

        super(TsvRepHoldOutRandomForest, self).__init__(
            input.TsvInput,
            validation.RepeatedHoldOut,
            algorithm.RandomForest,
            locals(),
            output_dir,
        )
