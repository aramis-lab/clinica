# coding: utf8

import clinica.engine as ce


class T1FreeSurferVisualizer(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-freesurfer'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Cross-sectional pre-processing of T1w images with FreeSurfer:\nhttp://clinica.run/doc/Pipelines/T1_FreeSurfer/'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("participant_id",
                                  help='Participant ID (e.g. sub-CLNC01).')
        clinica_comp.add_argument("session_id",
                                  help='Session ID (e.g. sub-M00).')

    def run_command(self, args):
        """
        """
        import os
        from clinica.utils.stream import cprint
        from clinica.utils.check_dependency import check_freesurfer

        check_freesurfer()

        participant_id = args.participant_id
        session_id = args.session_id
        subject_id = participant_id + '_' + session_id

        cprint('Visualizing outputs of t1-freesurfer pipeline for ' + participant_id + ' at session ' + session_id)
        caps_participant = os.path.join(
            args.caps_directory, 'subjects', participant_id, session_id,
            't1', 'freesurfer_cross_sectional', subject_id
        )
        command_line = \
            'freeview' \
            ' -v' \
            ' %s/mri/T1.mgz' \
            ' %s/mri/aseg.mgz:colormap=lut:opacity=0.2' \
            ' -f' \
            ' %s/surf/lh.white:edgecolor=blue' \
            ' %s/surf/lh.pial:edgecolor=green' \
            ' %s/surf/rh.white:edgecolor=blue' \
            ' %s/surf/rh.pial:edgecolor=green' % \
            (caps_participant, caps_participant,
             caps_participant, caps_participant, caps_participant, caps_participant)

        os.system(command_line)
