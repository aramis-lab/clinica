# coding: utf8

import clinica.engine as ce


class T1FreeSurferCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 't1-freesurfer'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = 'Cross-sectional pre-processing of T1w images with FreeSurfer:\nhttp://clinica.run/doc/Pipelines/T1_FreeSurfer/'

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-raa", "--recon_all_args",
                              help='Additional flags for recon-all command line (default: --recon_all_args "-qcache")')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel')

    def run_command(self, args):
        """Run the pipeline with defined args."""
        import os
        import datetime
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.stream import cprint
        from .t1_freesurfer_pipeline import T1FreeSurfer
        from tempfile import mkdtemp

        pipeline = T1FreeSurfer(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
        )

        pipeline.parameters = {
            'recon_all_args': args.recon_all_args or '-qcache'
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        try:
            if args.n_procs:
                pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
            else:
                pipeline.run()
        except RuntimeError as e:
            # Check that it is a Nipype error
            if 'Workflow did not execute cleanly. Check log for details' in str(e):
                from clinica.iotools.grabcaps import CAPSLayout
                import re

                def read_tsv(tsv_file):
                    import pandas
                    df = pandas.io.parsers.read_csv(tsv_file, sep='\t')
                    return list(df.participant_id), list(df.session_id)

                # Extract subject IDs given to the pipeline
                [participant_ids, session_ids] = read_tsv(os.path.join(
                    self.absolute_path(args.working_directory),
                    pipeline.__class__.__name__,
                    'participants.tsv')
                )
                participant_labels = '|'.join(sub[4:] for sub in participant_ids)
                session_labels = '|'.join(ses[4:] for ses in session_ids)
                # Extract subject IDs generated by the pipeline on CAPS folder
                caps_layout = CAPSLayout(self.absolute_path(args.caps_directory))
                caps_files = caps_layout.get(freesurfer_file='recon-all.log', return_type='file',
                                             subject=participant_labels, session=session_labels)
                # Extract missing subject IDs
                subject_ids = [participant_ids[i] + '_' + session_ids[i] for i
                               in range(len(participant_ids))]
                id_caps_files = [
                    re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)', file).group()
                    for file in caps_files]
                id_missing_subjects = list(set(subject_ids) - set(id_caps_files))
                missing_caps = ', '.join(
                    id_missing_subjects[i].split('_')[0][4:] + '|' + id_missing_subjects[i].split('_')[1][4:]
                    for i in range(len(id_missing_subjects)))

                # Display summary of the pipeline
                now = datetime.datetime.now().strftime('%H:%M:%S')
                cprint('\n%s[%s] The %s pipeline finished with errors.%s\n' %
                       (Fore.RED, now, self._name, Fore.RESET))
                cprint('%sCAPS outputs were not found for %s subject(s): %s%s\n' %
                       (Fore.RED, len(id_missing_subjects), missing_caps, Fore.RESET))
                cprint('%sError details can be found either by opening the log file (%s) or '
                       'by opening the crash file(s) with the following command(s):%s' %
                       (Fore.YELLOW, args.logname, Fore.RESET))
                log_file = open(args.logname, "r")
                import re
                for line in log_file:
                    if re.match("(.*)crashfile:(.*)", line):
                        cprint(Fore.YELLOW + line.replace('\t crashfile:', '- nipypecli crash').replace('\n', '') + Fore.RESET)
                cprint('\n%sIf your pipeline crashed due to lack of space of network issues, re-run the pipeline with '
                       'the working directory (-wd %s).\nOtherwise, you can delete it.%s' %
                       (Fore.YELLOW, os.path.abspath(args.working_directory), Fore.RESET))
                # Force the display of "Documentation can be found..."
                raise ClinicaException('')
            else:
                raise e
        else:
            now = datetime.datetime.now().strftime('%H:%M:%S')
            cprint('%s[%s]%s The %s pipeline has completed. You can now delete the working directory (%s).' %
                   (Fore.GREEN, now, Fore.RESET, self._name, os.path.join(os.path.abspath(args.working_directory), pipeline.__class__.__name__)))