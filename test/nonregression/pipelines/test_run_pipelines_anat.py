# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import pardir
from test.nonregression.testing_tools import *

# Determine location for working_directory
warnings.filterwarnings("ignore")


def test_run_T1FreeSurferCrossSectional(cmdopt):
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer.t1_freesurfer_pipeline import T1FreeSurfer

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1FreeSurfer")

    # Remove potential residual of previous tests
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1FreeSurfer"))

    parameters = {"recon_all_args": "-qcache", "skip_question": False}

    pipeline = T1FreeSurfer(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
        base_dir=join(working_dir, "T1FreeSurfer"),
    )
    pipeline.base_dir = join(working_dir, "T1FreeSurfer")
    pipeline.run(bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id, sess_id):
        import os

        output_folder = os.path.join(
            "caps", "subjects", part_id, sess_id, "t1", "freesurfer_cross_sectional"
        )
        return output_folder

    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(path_to_caps_fs("sub-01", "ses-2011"), "regional_measures"),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(path_to_caps_fs("sub-01", "ses-2011"), "sub-01_ses-2011", "label"),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(path_to_caps_fs("sub-01", "ses-2011"), "sub-01_ses-2011", "mri"),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(path_to_caps_fs("sub-01", "ses-2011"), "sub-01_ses-2011", "stats"),
    )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1FreeSurfer"), recreate=False)


def test_run_T1VolumeTissueSegmentation(cmdopt):
    import os
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1VolumeTissueSegmentation")
    clean_folder(join(working_dir, "T1VolumeTissueSegmentation"))
    clean_folder(join(root, "out", "caps"))

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1VolumeTissueSegmentation"),
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    out_file = join(
        root,
        "out/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/"
        + "sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz",
    )
    if not os.path.exists(out_file):
        raise IOError(
            "Pipeline did not produce file: "
            + out_file
            + ". Consider rerunning test_run_T1VolumeTissueSegmentation"
        )

    ref_file = join(
        root,
        "ref/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/"
        + "sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz",
    )

    assert likeliness_measure(out_file, ref_file, (1e-1, 0.02), (0.4, 0.01))

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeTissueSegmentation"), recreate=False)


def test_run_T1VolumeCreateDartel(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1VolumeCreateDartel")

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, "T1VolumeCreateDartel"))
    clean_folder(join(root, "out", "caps"), recreate=False)
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {"group_label": "UnitTest"}
    # Instantiate pipeline
    pipeline = T1VolumeCreateDartel(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1VolumeCreateDartel"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    out_data_template = join(
        root, "out/caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz"
    )
    ref_data_template = join(root, "ref/group-UnitTest_template.nii.gz")
    assert likeliness_measure(
        out_data_template, ref_data_template, (1e-3, 0.1), (1e-2, 0.1)
    )

    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_forward_def = [
        join(
            root,
            "out",
            "caps",
            "subjects",
            sub,
            "ses-M00",
            "t1",
            "spm",
            "dartel",
            "group-UnitTest",
            sub
            + "_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        )
        for sub in subjects
    ]
    ref_data_forward_def = [
        join(
            root,
            "ref",
            sub
            + "_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        )
        for sub in subjects
    ]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(
            out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1)
        )

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeCreateDartel"), recreate=False)


def test_run_T1VolumeDartel2MNI(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1VolumeDartel2MNI")

    # Remove potential residual of previous UT
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeDartel2MNI"))

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {"group_label": "UnitTest"}
    # Instantiate pipeline and run()
    pipeline = T1VolumeDartel2MNI(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1VolumeDartel2MNI"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_GM_MNI = [
        join(
            root,
            "out",
            "caps",
            "subjects",
            sub,
            "ses-M00",
            "t1",
            "spm",
            "dartel",
            "group-UnitTest",
            sub
            + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz",
        )
        for sub in subjects
    ]
    ref_data_GM_MNI = [
        join(
            root,
            "ref",
            sub
            + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz",
        )
        for sub in subjects
    ]
    for i in range(len(out_data_GM_MNI)):
        assert likeliness_measure(
            out_data_GM_MNI[i], ref_data_GM_MNI[i], (1e-4, 0.15), (1, 0.02)
        )

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeDartel2MNI"), recreate=False)


def test_run_T1VolumeRegisterDartel(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1VolumeExistingDartel")
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeExistingDartel"))

    # Copy necessary data to run pipeline
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    # Instantiate and run pipeline
    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeRegisterDartel(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1VolumeExistingDartel"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_forward_def = [
        join(
            root,
            "out",
            "caps",
            "subjects",
            sub,
            "ses-M00",
            "t1",
            "spm",
            "dartel",
            "group-UnitTest",
            sub
            + "_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        )
        for sub in subjects
    ]
    ref_data_forward_def = [
        join(
            root,
            "ref",
            sub
            + "_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        )
        for sub in subjects
    ]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(
            out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1)
        )

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeExistingDartel"), recreate=False)


def test_run_T1VolumeParcellation(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import numpy as np
    import pandas as pds

    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1VolumeParcellation")
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeParcellation"))

    # Copy data for use of pipeline
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    # Instantiate pipeline
    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeParcellation(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1VolumeParcellation"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_files = [
        join(
            root,
            "out/caps/subjects/sub-ADNI018S4696/ses-M00/t1/spm/dartel/group-UnitTest/atlas_statistics",
            "sub-ADNI018S4696_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-"
            + atlas
            + "_map-graymatter_statistics.tsv",
        )
        for atlas in pipeline.parameters["atlases"]
    ]
    ref_files = [
        join(
            root,
            "ref/sub-ADNI018S4696_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-"
            + atlas
            + "_map-graymatter_statistics.tsv",
        )
        for atlas in pipeline.parameters["atlases"]
    ]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep="\t")
        ref_csv = pds.read_csv(ref_files[i], sep="\t")
        assert np.allclose(
            np.array(out_csv.mean_scalar),
            np.array(ref_csv.mean_scalar),
            rtol=1e-8,
            equal_nan=True,
        )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1VolumeParcellation"), recreate=False)


def test_run_T1Linear(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_linear.t1_linear_pipeline import T1Linear

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1Linear")

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, "T1Linear"))
    clean_folder(join(root, "out", "caps"), recreate=False)

    parameters = {"uncropped_image": False}
    # Instantiate pipeline
    pipeline = T1Linear(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1Linear"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    out_folder = join(root, "out")
    ref_folder = join(root, "ref")

    compare_folders(out_folder, ref_folder, shared_folder_name="caps")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1Linear"), recreate=False)


def test_run_T1FreeSurferTemplate(cmdopt):
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    # sub-01 was duplicated into to sub-02 with one session in order to test the "one time point" case
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_pipeline import (
        T1FreeSurferTemplate,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1FreeSurferTemplate")

    # Remove potential residual of previous tests
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1FreeSurferTemplate"))

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    pipeline = T1FreeSurferTemplate(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1FreeSurferTemplate"),
    )
    pipeline.base_dir = join(working_dir, "T1FreeSurferTemplate")
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 2}, bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id, long_id):
        import os

        output_folder = os.path.join(
            "caps", "subjects", part_id, long_id, "freesurfer_unbiased_template"
        )
        return output_folder

    for (p_id, l_id) in zip(["sub-01", "sub-02"], ["long-20112015", "long-2011"]):
        compare_folders(
            join(root, "out"),
            join(root, "ref"),
            join(path_to_caps_fs(p_id, l_id), p_id + "_" + l_id, "label"),
        )
        compare_folders(
            join(root, "out"),
            join(root, "ref"),
            join(path_to_caps_fs(p_id, l_id), p_id + "_" + l_id, "mri"),
        )
        compare_folders(
            join(root, "out"),
            join(root, "ref"),
            join(path_to_caps_fs(p_id, l_id), p_id + "_" + l_id, "stats"),
        )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1FreeSurferTemplate"), recreate=False)


def test_run_T1FreeSurferLongitudinalCorrection(cmdopt):
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "T1FreeSurferLongitudinalCorrection")

    # Remove potential residual of previous tests
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "T1FreeSurferLongitudinalCorrection"))

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "T1FreeSurferLongitudinalCorrection"),
    )
    pipeline.base_dir = join(working_dir, "T1FreeSurferLongitudinalCorrection")
    pipeline.run(bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id, sess_id, long_id):
        import os

        output_folder = os.path.join(
            "caps",
            "subjects",
            part_id,
            sess_id,
            "t1",
            long_id,
            "freesurfer_longitudinal",
        )
        return output_folder

    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(
            path_to_caps_fs("sub-01", "ses-2011", "long-20112015"), "regional_measures"
        ),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(
            path_to_caps_fs("sub-01", "ses-2011", "long-20112015"),
            "sub-01_ses-2011.long.sub-01_long-20112015",
            "label",
        ),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(
            path_to_caps_fs("sub-01", "ses-2011", "long-20112015"),
            "sub-01_ses-2011.long.sub-01_long-20112015",
            "mri",
        ),
    )
    compare_folders(
        join(root, "out"),
        join(root, "ref"),
        join(
            path_to_caps_fs("sub-01", "ses-2011", "long-20112015"),
            "sub-01_ses-2011.long.sub-01_long-20112015",
            "stats",
        ),
    )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(
        join(working_dir, "T1FreeSurferLongitudinalCorrection"), recreate=False
    )
