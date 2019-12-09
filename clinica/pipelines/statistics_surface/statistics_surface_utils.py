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


def prepare_data(input_directory, subjects_visits_tsv, group_label, glm_type):
    """Fetch all the intermediate variables for this workflow.

    Args:
        input_directory: CAPS directory
        subjects_visits_tsv: TSV file defining the GLM model
        group_label: Current group name for this analysis
        glm_type: Based on the hypothesis, you should define one of the glm types, "group_comparison", "correlation"

    Returns:

    """
    import os
    from shutil import copy
    import clinica.pipelines as clp
    import pandas as pd

    path_to_matlab_script = os.path.join(os.path.dirname(clp.__path__[0]), 'lib', 'clinicasurfstat')

    # CAPS input and output vars
    input_directory = os.path.expanduser(input_directory)
    surfstat_input_dir = os.path.join(input_directory, 'subjects')

    group_id = 'group-' + group_label
    statistics_dir_tsv = os.path.join(input_directory, 'groups', group_id, 'statistics', 'participant.tsv')

    if glm_type == "group_comparison":
        output_directory = os.path.join(input_directory, 'groups', group_id, 'statistics', 'surfstat_group_comparison')
    elif glm_type == "correlation":
        output_directory = os.path.join(input_directory, 'groups', group_id, 'statistics', 'surfstat_correlation_analysis')
    else:
        raise NotImplementedError("The other GLM situations have not been implemented in this pipeline.")

    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except BaseException:
            raise OSError("SurfStat: can't create destination directory (%s)!" % output_directory)

    # Copy the subjects_visits_tsv to the result folder
    # First, check if the subjects_visits_tsv has the same info with the participant.tsv in the folder of statistics.
    # If the participant TSV does not exit, copy subjects_visits_tsv in the folder of statistics too,
    # if it is here, compare them.
    if not os.path.isfile(statistics_dir_tsv):
        copy(subjects_visits_tsv, statistics_dir_tsv)
    else:
        # Compare the two TSV files
        participant_df = pd.io.parsers.read_csv(statistics_dir_tsv, sep='\t')
        subjects_visits_tsv_df = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
        participant_list = list(participant_df.participant_id)
        subjects_visits_tsv_list = list(subjects_visits_tsv_df.participant_id)
        dif_list = list(set(participant_list) - set(subjects_visits_tsv_list))
        try:
            len(dif_list) == 0
        except BaseException:
            raise ValueError("It seems that this round of analysis does not contain the same subjects"
                             "where you want to put the results, please check it!")

        group_tsv = 'group-' + group_label + '_participants.tsv'
        copied_tsv = os.path.join(output_directory, group_tsv)
        copy(subjects_visits_tsv, copied_tsv)
    # Point to the path to the json file
    out_json = os.path.join(output_directory, 'group-' + group_label + '_glm.json')

    freesurfer_home = os.environ["FREESURFER_HOME"]

    return path_to_matlab_script, surfstat_input_dir, output_directory, freesurfer_home, out_json


def run_matlab(input_directory,
               output_directory,
               subjects_visits_tsv,
               pipeline_parameters,
               freesurfer_home,
               path_to_matlab_script,
               ):
    """
    Wrap the Matlab script of SurfStat.

    Args:
        input_directory (str): surfstat_input_dir where containing all the subjects' output in CAPS directory
        output_directory (str): output folder to contain the result in CAPS folder
        subjects_visits_tsv (str): TSV file containing the GLM information
        pipeline_parameters (dict): parameters of StatisticsSurface pipeline
        freesurfer_home (str): the environmental variable $FREESURFER_HOME
        path_to_matlab_script (str): path to find the matlab script

    Returns:

    """
    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
    from os.path import join
    import sys
    from clinica.utils.stream import cprint

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
        input_directory,
        output_directory,
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
        'thresholduncorrectedpvalue', pipeline_parameters['threshold_uncorrected_pvalue'],
        'thresholdcorrectedpvalue', pipeline_parameters['threshold_corrected_pvalue'],
        'clusterthreshold', pipeline_parameters['cluster_threshold']
        )
    # This will create a file: pyscript.m , the pyscript.m is the default name
    matlab.inputs.mfile = True
    # This will stop running with single thread
    matlab.inputs.single_comp_thread = False
    matlab.inputs.logfile = join(output_directory, "matlab_output.log")
    cprint("Matlab logfile is located at the following path: %s" % matlab.inputs.logfile)
    cprint("Matlab script command = %s" % matlab.inputs.script)
    cprint("MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread)
    cprint("MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command())
    if sys.platform.startswith('linux'):
        cprint("MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args)
    out = matlab.run()
    return out


def create_glm_info_dictionary(pipeline_parameters):
    """Create dictionary containing the GLM information that will be stored in a JSON file."""
    json_dict = {
        'AnalysisType': pipeline_parameters['glm_type'],
        'DesignMatrix': pipeline_parameters['design_matrix'],
        'StringFormatTSV': pipeline_parameters['str_format'],
        'Contrast': pipeline_parameters['contrast'],
        'GroupLabel': pipeline_parameters['group_label'],
        'FWHM': pipeline_parameters['full_width_at_half_maximum'],
        'ThresholdUncorrectedPvalue': pipeline_parameters['threshold_uncorrected_pvalue'],
        'ThresholdCorrectedPvalue': pipeline_parameters['threshold_corrected_pvalue'],
        'ClusterThreshold': pipeline_parameters['cluster_threshold']
    }

    return json_dict
