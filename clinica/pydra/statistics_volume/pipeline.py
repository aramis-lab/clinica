import typing as ty
from os import PathLike
from pathlib import PurePath

import pydra
from pydra import Workflow

from clinica.pydra.engine import clinica_io


@clinica_io
def build_core_workflow(name: str = "core", parameters={}) -> Workflow:
    """doc"""
    from os.path import dirname, join

    from pydra.tasks.nipype1.utils import Nipype1Task

    import clinica.pydra.statistics_volume.task as utils

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            (
                "pet_volume",
                dict,
                {"acq_label": parameters["acq_label"]},
                {"group_label": parameters["group_label_dartel"]},
                {"suvr_reference_region": parameters["suvr_reference_region"]},
                {"use_brainmasked_image": True},
                {"use_pvc_data": parameters["use_pvc_data"]},
                {"fwhm": parameters["full_width_at_half_maximum"]},
                {"mandatory": True},
            ),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(name, input_spec=input_spec)

    wf.split(("pet_volume"))

    wf.add(
        utils.unzip_nii(
            name="unzip_nii",
            in_file=wf.lzin.pet_volume,
        )
    )
    wf.add(
        utils.get_group_1_and_2(
            name="get_groups",
            tsv=parameters["tsv_file"],
            contrast=parameters["contrast"],
            # idx_group1=wf.02_get_groups.lzout.first_group_idx,
            # idx_group2=wf.02_get_groups.lzout.second_group_idx,
            # class_names=wf.02_get_groups.lzout.class_names,
        )
    )
    # 1. Model creation
    # We use overwrite option to be sure this node is always run so that it can delete the output dir if it
    # already exists (this may cause error in output files otherwise)
    wf.add(
        utils.model_creation(
            name="model_creation",
            tsv=parameters["tsv_file"],
            contrast=parameters["contrast"],
            template_file=join(dirname(__file__), "template_model_creation.m"),
            file_list=wf.unzip_nii.lzout.unzipped_nii,
            idx_group1=wf.get_groups.lzout.first_group_idx,
            idx_group2=wf.get_groups.lzout.second_group_idx,
            # output
            # script_file=wf.model_creation.lzout.current_model,
            # covariate=wf.model_creation.lzout.covariates,
        )
    )
    wf.add(
        utils.run_m_script(
            name="run_spm_model_creation",
            m_file=wf.model_creation.lzout.current_model,
            # output
            # spm_mat=wf.run_spm_model_creation.lzout.output_mat_file,
        )
    )
    # 2. Model estimation
    wf.add(
        utils.estimate(
            name="model_estimation",
            # interface=utils.estimate
            mat_file=wf.run_spm_model_creation.lzout.output_mat_file,
            template_file=join(dirname(__file__), "template_model_estimation.m"),
            # output
            # script_file=wf.model_estimation.lzout.current_model_estimation,
        )
    )
    wf.add(
        utils.run_m_script(
            name="run_spm_model_estimation",
            m_file=wf.model_estimation.lzout.current_model_estimation,
            # output
            # spm_mat=wf.run_spm_model_estimation.lzout.output_mat_file,
        )
    )
    # 3. Contrast
    wf.add(
        utils.contrast(
            name="model_contrast",
            mat_file=wf.run_spm_model_estimation.lzout.output_mat_file,
            template_file=join(dirname(__file__), "template_model_contrast.m"),
            covariates=wf.model_creation.lzout.covariates,
            class_names=wf.get_groups.lzout.class_names,
            # output
            # script_file=wf.model_contrast.lzout.current_model_estimation,
        )
    )
    wf.add(
        utils.run_m_script(
            name="run_spm_model_contrast",
            m_file=wf.model_contrast.lzout.current_model_estimation,
            # output
            # spm_mat=wf.run_spm_model_contrast.lzout.output_mat_file,
        )
    )
    # 4. Results
    wf.add(
        utils.results(
            name="model_result_no_correction",
            mat_file=wf.run_spm_model_contrast.lzout.output_mat_file,
            template_file=join(dirname(__file__), "template_model_results.m"),
            method="none",
            threshold=parameters["cluster_threshold"],
            # output
            # script_file=wf.06_1_model_result_no_correction.lzout.current_model_result,
        )
    )
    wf.add(
        utils.run_m_script(
            name="run_spm_model_result_no_correction",
            m_file=wf.model_result_no_correction.lzout.current_model_result,
            # output
            # spm_mat=wf.run_spm_model_result_no_correction.lzout.output_mat_file,
        )
    )
    wf.add(
        utils.read_output(
            name="read_output_node",
            spm_mat=wf.run_spm_model_result_no_correction.lzout.output_mat_file,
            class_names=wf.get_groups.lzout.class_names,
            covariates=wf.model_creation.lzout.covariates,
            group_label=parameters["group_label"],
            fwhm=parameters["full_width_at_half_maximum"],
            measure=parameters["measure_label"],
        )
    )
    wf.set_output(
        [
            ("spmT_0001", wf.read_ouput_node.lzout.spmT_0001),
            ("spmT_0002", wf.read_ouput_node.lzout.spmT_0002),
            ("spm_figures", wf.read_ouput_node.lzout.spm_figures),
            ("variance_of_error", wf.read_ouput_node.lzout.variance_of_error),
            ("resels_per_voxels", wf.read_ouput_node.lzout.resels_per_voxels),
            ("mask", wf.read_ouput_node.lzout.mask),
            ("regression_coeff", wf.read_ouput_node.lzout.regression_coeff),
            ("contrasts", wf.read_ouput_node.lzout.contrasts),
        ]
    )

    return wf
