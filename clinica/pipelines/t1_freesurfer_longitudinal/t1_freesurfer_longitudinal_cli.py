# coding: utf8

import clinica.engine as ce


class T1FreeSurferLongitudinalCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "t1-freesurfer-longitudinal"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Longitudinal pre-processing of T1w images with FreeSurfer:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer_Longitudinal/"
        )

    def define_options(self):
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_overwrite_flag=True)

    def run_command(self, args):
        """Run the pipeline with defined args."""
        import datetime
        import os

        from clinica.utils.longitudinal import get_participants_long_id
        from clinica.utils.participant import get_subject_session_list
        from clinica.utils.stream import cprint

        from .longitudinal_utils import save_part_sess_long_ids_to_tsv
        from .t1_freesurfer_longitudinal_correction_cli import (
            T1FreeSurferLongitudinalCorrectionCLI,
        )
        from .t1_freesurfer_template_cli import T1FreeSurferTemplateCLI

        cprint(
            "The t1-freesurfer-longitudinal pipeline is divided into 2 parts:\n"
            "\tt1-freesurfer-unbiased-template pipeline: Creation of unbiased template\n"
            "\tt1-freesurfer-longitudinal-correction pipeline: Longitudinal correction\n"
        )

        if not self.absolute_path(args.subjects_sessions_tsv):
            l_sess, l_part = get_subject_session_list(
                self.absolute_path(args.caps_directory), None, False, False
            )
            l_long = get_participants_long_id(l_part, l_sess)
            now = datetime.datetime.now().strftime("%H%M%S")
            args.subjects_sessions_tsv = now + "_participants.tsv"
            save_part_sess_long_ids_to_tsv(
                l_part, l_sess, l_long, os.getcwd(), args.subjects_sessions_tsv
            )

        cprint("Part 1/2: Running t1-freesurfer-unbiased-template pipeline")
        unbiased_template_cli = T1FreeSurferTemplateCLI()
        unbiased_template_cli.run_command(args)

        cprint("Part 2/2 Running t1-freesurfer-longitudinal-correction pipeline")
        longitudinal_correction_cli = T1FreeSurferLongitudinalCorrectionCLI()
        longitudinal_correction_cli.run_command(args)
