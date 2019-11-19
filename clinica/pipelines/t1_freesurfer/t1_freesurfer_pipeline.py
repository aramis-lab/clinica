# coding: utf8

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config

import clinica.pipelines.engine as cpe

cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class T1FreeSurfer(cpe.Pipeline):
    """FreeSurfer-based processing of T1-weighted MR images.

    Returns:
        A clinica pipeline object containing the T1FreeSurfer pipeline.
    """

    def get_subject_session_list(self, input_dir, tsv_file, is_bids_dir, base_dir):
        """Parse a BIDS or CAPS directory to get the subjects and sessions."""
        from colorama import Fore
        from clinica.utils.stream import cprint
        from clinica.utils.io import extract_subjects_sessions_from_filename
        from clinica.utils.inputs import clinica_file_reader
        import clinica.utils.input_files as input_files

        super(T1FreeSurfer, self).get_subject_session_list(input_dir, tsv_file, is_bids_dir, base_dir)

        if not tsv_file:
            # We only extract <participant_id>_<session_id> having T1w files
            t1w_files = clinica_file_reader(self.subjects, self.sessions,
                                            self.bids_directory, input_files.T1W_NII, False)
            self._subjects, self._sessions = extract_subjects_sessions_from_filename(t1w_files)

        # Display image(s) already present in CAPS folder
        t1_freesurfer_files = clinica_file_reader(self.subjects, self.sessions,
                                                  self.caps_directory, input_files.T1_FS_DESTRIEUX, False)

        processed_participants, processed_sessions = extract_subjects_sessions_from_filename(t1_freesurfer_files)
        if len(processed_participants) > 0:
            cprint("%sClinica found %s image(s) already processed in CAPS directory:%s" %
                   (Fore.YELLOW, len(processed_participants), Fore.RESET))
            for p_id, s_id in zip(processed_participants, processed_sessions):
                cprint("%s\t%s | %s%s" % (Fore.YELLOW, p_id, s_id, Fore.RESET))
            if self.overwrite_caps:
                output_folder = "<CAPS>/subjects/sub-<participant_id>/ses-<session_id>/t1/freesurfer_cross_sectional"
                cprint("%s\nOutput folders in %s will be recreated.\n%s" % (Fore.YELLOW, output_folder, Fore.RESET))
            else:
                cprint("%s\nImage(s) will be ignored by Clinica.\n%s" % (Fore.YELLOW, Fore.RESET))
                input_ids = [p_id + '_' + s_id for p_id, s_id in zip(self.subjects, self.sessions)]
                processed_ids = [p_id + '_' + s_id for p_id, s_id in zip(processed_participants, processed_sessions)]
                to_process_ids = list(set(input_ids) - set(processed_ids))
                self.subjects, self.sessions = extract_subjects_sessions_from_filename(to_process_ids)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the T1FreeSurfer pipeline is:
                * t1w (str): Path of the T1 weighted image in BIDS format

        Returns:
            A list of (string) input fields name.
        """
        return ['t1w']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Note:
            The list of outputs of the T1FreeSurfer pipeline is:
                * image_id (str): Image ID (e.g. sub-CLNC01_ses-M00)

        Returns:
            A list of (string) output fields name.
        """
        return ['image_id']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.

        Raise:
            ClinicaBIDSError: If there are duplicated files or missing files for any subject
        """
        import os
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (print_images_to_process, print_no_image_to_process)
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_NII
        from clinica.iotools.utils.data_handling import check_volume_location_in_world_coordinate_system
        from clinica.utils.io import save_participants_sessions

        # Inputs from anat/ folder
        # ========================
        # T1w file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.bids_directory,
                                            T1W_NII)
        except ClinicaException as e:
            err_msg = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n' + str(e)
            raise ClinicaBIDSError(err_msg)

        # Save subjects to process in <WD>/T1FreeSurfer/participants.tsv
        save_participants_sessions(self.subjects, self.sessions, os.path.join(self.base_dir, self.__class__.__name__))

        print_images_to_process(self.subjects, self.sessions)
        cprint('The pipeline will last approximately 10 hours per image.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )
        check_volume_location_in_world_coordinate_system(t1w_files, self.bids_directory)
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from .t1_freesurfer_utils import save_to_caps
        import os

        save_to_caps = npe.Node(interface=nutil.Function(
            input_names=['source_dir', 'image_id', 'caps_dir', 'overwrite_caps'],
            output_names=['image_id'],
            function=save_to_caps),
            name='SaveToCaps')
        save_to_caps.inputs.source_dir = os.path.join(self.base_dir, 'T1FreeSurfer', 'ReconAll')
        save_to_caps.inputs.caps_dir = self.caps_directory
        save_to_caps.inputs.overwrite_caps = self.overwrite_caps

        self.connect([
            (self.output_node, save_to_caps, [('image_id', 'image_id')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import os
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.freesurfer.preprocess import ReconAll
        from .t1_freesurfer_utils import (init_input_node, write_tsv_files)

        # Nodes declaration
        # =================
        # Initialize the pipeline
        #   - Extract <image_id> (e.g. sub-CLNC01_ses-M00) T1w filename;
        #   - Check FOV of T1w;
        #   - Create <subjects_dir> folder in <WD>/T1FreeSurfer/ReconAll/<image_id>/;
        #   - Print begin execution message.
        init_input = npe.Node(
            interface=nutil.Function(
                input_names=['t1w', 'recon_all_args', 'output_dir'],
                output_names=['image_id', 't1w', 'flags', 'subjects_dir'],
                function=init_input_node),
            name='0-InitPipeline')
        init_input.inputs.recon_all_args = self.parameters['recon_all_args']
        init_input.inputs.output_dir = os.path.join(self.base_dir, 'T1FreeSurfer', 'ReconAll')

        # Run recon-all command
        # FreeSurfer segmentation will be in <subjects_dir>/<image_id>/
        recon_all = npe.Node(interface=ReconAll(),
                             name='1-SegmentationReconAll')
        recon_all.inputs.directive = 'all'

        # Generate TSV files containing a summary of the regional statistics
        # in <subjects_dir>/regional_measures
        create_tsv = npe.Node(interface=nutil.Function(
            input_names=['subjects_dir', 'image_id'],
            output_names=['image_id'],
            function=write_tsv_files),
            name='2-CreateTsvFiles')

        # Connections
        # ===========
        self.connect([
            # Get <image_id> from input_node and print begin message
            (self.input_node, init_input, [('t1w', 't1w')]),
            # Run recon-all command
            (init_input, recon_all, [('subjects_dir', 'subjects_dir'),
                                     ('t1w', 'T1_files'),
                                     ('image_id', 'subject_id'),
                                     ('flags',  'flags')]),
            # Generate TSV files
            (init_input, create_tsv, [('subjects_dir', 'subjects_dir')]),
            (recon_all,  create_tsv, [('subject_id', 'image_id')]),
            # Output node
            (recon_all, self.output_node, [('subject_id', 'image_id')]),
        ])
