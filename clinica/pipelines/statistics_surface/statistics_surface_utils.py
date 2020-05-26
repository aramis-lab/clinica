# coding: utf8


def get_t1_freesurfer_custom_file():
    import os
    custom_file = os.path.join(
        '@subject',
        '@session',
        't1',
        'freesurfer_cross_sectional',
        '@subject_@session',
        'surf',
        '@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
    )
    return custom_file


def get_fdg_pet_surface_custom_file():
    import os
    custom_file = os.path.join(
        '@subject',
        '@session',
        'pet',
        'surface',
        '@subject_@session_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-@hemi_fwhm-@fwhm_projection.mgh'
    )
    return custom_file


def init_input_node(parameters, base_dir, subjects_visits_tsv):
    """Initialize the pipeline.

    This function will:
        - Create `surfstat_results_dir` in `base_dir`/<group_id> for SurfStat;
        - Save pipeline parameters in JSON file;
        - Copy TSV file with covariates;
        - Print begin execution message.
    """
    import os
    import errno
    import json
    import shutil
    from clinica.utils.ux import print_begin_image
    from clinica.pipelines.statistics_surface.statistics_surface_utils import create_glm_info_dictionary

    group_id = 'group-' + parameters['group_label']

    # Create surfstat_results_dir for SurfStat
    surfstat_results_dir = os.path.join(base_dir, group_id)
    try:
        os.makedirs(surfstat_results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:  # EEXIST: folder already exists
            raise e

    # Save pipeline parameters in JSON file
    glm_dict = create_glm_info_dictionary(parameters)
    json_filename = os.path.join(surfstat_results_dir, group_id + '_glm.json')
    with open(json_filename, 'w') as json_file:
        json.dump(glm_dict, json_file, indent=4)

    # Copy TSV file with covariates
    tsv_filename = os.path.join(surfstat_results_dir, group_id + '_covariates.tsv')
    shutil.copyfile(subjects_visits_tsv, tsv_filename)

    # Print begin message
    list_keys = ['AnalysisType', 'DesignMatrix', 'Contrast', 'FWHM', 'ClusterThreshold']
    list_values = [
        parameters['glm_type'],
        parameters['design_matrix'],
        parameters['contrast'],
        str(parameters['full_width_at_half_maximum']),
        str(parameters['cluster_threshold'])
    ]
    group_id = 'group-' + parameters['group_label']
    print_begin_image(group_id, list_keys, list_values)

    return parameters['group_label'], surfstat_results_dir


def run_matlab(caps_dir,
               output_dir,
               subjects_visits_tsv,
               pipeline_parameters):
    """
    Wrap the call of SurfStat using clinicasurfstat.m Matlab script.

    Args:
        caps_dir (str): CAPS directory containing surface-based features
        output_dir (str): Output directory that will contain outputs of clinicasurfstat.m
        subjects_visits_tsv (str): TSV file containing the GLM information
        pipeline_parameters (dict): parameters of StatisticsSurface pipeline
    """
    import os
    import sys
    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
    import clinica.pipelines as clinica_pipelines
    from clinica.utils.check_dependency import check_environment_variable

    path_to_matlab_script = os.path.join(os.path.dirname(clinica_pipelines.__path__[0]), 'lib', 'clinicasurfstat')
    freesurfer_home = check_environment_variable('FREESURFER_HOME', 'FreeSurfer')

    MatlabCommand.set_default_matlab_cmd(
        get_matlab_command()
    )
    matlab = MatlabCommand()

    # Add the dynamic traits
    # opengl_trait = traits.Bool(
    #     True, argstr='-nosoftwareopengl', usedefault=True, desc='Switch on hardware openGL', nohash=True
    # )
    # matlab.input_spec.add_trait(matlab.input_spec(), 'nosoftwareopengl', opengl_trait())
    if sys.platform.startswith('linux'):
        # Bug(JW): for my laptop, it does not work, but the command line does have the flag -nosoftwareopengl,
        # TODO: We should try on other computer's matlab to check if this flag works!
        matlab.inputs.args = '-nosoftwareopengl'
    matlab.inputs.paths = path_to_matlab_script

    matlab.inputs.script = """
    clinicasurfstat('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, '%s', %.3f, '%s', %.3f, '%s', %.3f);
    """ % (
        os.path.join(caps_dir, 'subjects'),
        output_dir,
        subjects_visits_tsv,
        pipeline_parameters['design_matrix'],
        pipeline_parameters['contrast'],
        pipeline_parameters['str_format'],
        pipeline_parameters['glm_type'],
        pipeline_parameters['group_label'],
        freesurfer_home,
        pipeline_parameters['custom_file'],
        pipeline_parameters['feature_label'],
        'sizeoffwhm', pipeline_parameters['full_width_at_half_maximum'],
        'thresholduncorrectedpvalue', 0.001,
        'thresholdcorrectedpvalue', 0.05,
        'clusterthreshold', pipeline_parameters['cluster_threshold']
        )
    # This will create a file: pyscript.m , the pyscript.m is the default name
    matlab.inputs.mfile = True
    # This will stop running with single thread
    matlab.inputs.single_comp_thread = False
    matlab.inputs.logfile = 'group-' + pipeline_parameters['group_label'] + '_matlab.log'

    # cprint("Matlab logfile is located at the following path: %s" % matlab.inputs.logfile)
    # cprint("Matlab script command = %s" % matlab.inputs.script)
    # cprint("MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread)
    # cprint("MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command())
    # if sys.platform.startswith('linux'):
    #     cprint("MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args)
    matlab.run()

    return output_dir


def create_glm_info_dictionary(pipeline_parameters):
    """Create dictionary containing the GLM information that will be stored in a JSON file."""
    out_dict = {
        'AnalysisType': pipeline_parameters['glm_type'],
        'DesignMatrix': pipeline_parameters['design_matrix'],
        'StringFormatTSV': pipeline_parameters['str_format'],
        'Contrast': pipeline_parameters['contrast'],
        'GroupLabel': pipeline_parameters['group_label'],
        'FWHM': pipeline_parameters['full_width_at_half_maximum'],
        'ThresholdUncorrectedPvalue': 0.001,
        'ThresholdCorrectedPvalue': 0.05,
        'ClusterThreshold': pipeline_parameters['cluster_threshold']
    }
    return out_dict


def save_to_caps(source_dir, caps_dir, overwrite_caps, pipeline_parameters):
    """Save `source_dir`/ to CAPS folder.

    This function copies outputs of `source_dir`/ >to
    `caps_dir`/groups/<group_id>/<statistics>/surfstat_<glm_type>/

    The `source_dir`/ folder should contain the following elements:
        - group-<group_label>_<group_1_or_2>-lt-<group_1_or_2>_measure-<measure>_fwhm-<label>_suffix.ext
    or
        - group-<group_label>_correlation-<label>_contrast-{-|+}_measure-<measure>_fwhm-<label>_suffix.ext
    and
        - group-<group_label>_covariates.tsv
        - group-<group_label>_glm.json

    Raise:
        NotImplementedError: If overwrite_caps=True.
    """
    import os
    import shutil
    from clinica.utils.ux import print_end_image

    group_id = 'group-' + pipeline_parameters['group_label']

    if pipeline_parameters['glm_type'] == "group_comparison":
        surfstat_folder = 'surfstat_' + pipeline_parameters['glm_type']
    elif pipeline_parameters['glm_type'] == "correlation":
        surfstat_folder = 'surfstat_' + pipeline_parameters['glm_type'] + '_analysis'
    else:
        raise NotImplementedError("The other GLM situations have not been implemented in this pipeline.")

    destination_dir = os.path.join(
        os.path.expanduser(caps_dir),
        'groups',
        group_id,
        'statistics',
        surfstat_folder
    )

    if overwrite_caps:
        raise NotImplementedError('save_to_caps(overwrite_caps=True) not implemented')
    shutil.copytree(source_dir, destination_dir, symlinks=True)
    print_end_image(group_id)
